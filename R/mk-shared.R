# shared library loads and helper function for Mk model project

require(ape)
require(phangorn)
require(castor)
require(igraph)
require(ggplot2)
require(ggraph)
require(stringr)
require(ggpubr)
require(ggtree)

## Functions from the above packages called explicitly
## ape:
##     rphylo
##     stree
##     drop.tip

## phangorn:
##     getRoot

## castor:
##     fit_mk


# binary to decimal function
BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

# decimal to binary function
DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

# decimal to binary function, returning a numerical vector
DecToBinV <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(s)
}

# function to create a plot label given some statistics and labels
titlestr = function(expt, fit = c("rev", "irrev"), stats.df) {
  fit = match.arg(fit)
  if(fit == "rev") { lead.str = "Reversible" } else {lead.str = "Irreversible"}
  t.str = paste(c(lead.str, " fit, simplified AIC ~ ", round(stats.df$AIC.reduced, digits=2),
                  " (full ", round(stats.df$AIC, digits=2), ")"),
                collapse = "")
  return(t.str)
}

# this function converts a species name string from the Newick format which Common Taxonomy Tree gives us into a simpler lower-case, no quotes version comparable to Kostas' dataset
convname = function(str) {
  return(tolower(gsub("\'", "", gsub("_", " ", str))))
}

#' Make an indexing matrix containing required parameter structure
#'
#' Construct matrix assigning indices to elements describing possible transitions
#'
#' @param L (integer) number of features
#' @param reversible (Boolean, default TRUE) whether to create a reversible model or not
#' @param to.nullify (n x 2 matrix) from-to transitions to set to zero
#'
#' @return a matrix with the required parameter index structure
#' @export
#' @examples
#' mk_index_matrix(3)
#'
mk_index_matrix = function(L, reversible=TRUE, to.nullify=matrix(nrow=0,ncol=2)) {
  index_matrix = matrix(rep(0, (2**L)**2), nrow=2**L)
  index = 1
  # go through binary strings and put indices on transitions to Hamming +1 neighbours
  for(i in 0:(2**L-1)) {
    b = as.numeric(strsplit(DecToBin(i, L),"")[[1]])
    for(j in 1:L) {
      if(b[j] == 0) {
        c = b
        c[j] = 1
        partner = BinToDec(c)
        if(!any(to.nullify[,1]==i+1 & to.nullify[,2]==partner+1)) {
          index_matrix[i+1,partner+1] = index
          index = index+1
        }
        if(reversible == TRUE) {
          # add reverse transition if required
          if(!any(to.nullify[,1]==partner+1 & to.nullify[,2]==i+1)) {
            index_matrix[partner+1,i+1] = index
            index = index+1
          }
        }
      }
    }
  }

  return(index_matrix)
}


## FIXME: should be possible to add some minimal error checking
##        in mk_pull_transitions
## if reversible = FALSE but any non-zero entry in lower-triangular
##                       of transition_matrix: give warning
## if reversible = TRUE  but no non-zero entry in lower-triangular
##                       of transition_matrix: give warning


#' Extract transitions as data frame from fitted Mk model
#'
#' @param fitted.obj (hypermk model fit) the fitted model
#' @param reversible (Boolean, default TRUE) whether to use a reversible model or not
#'
#' @return a dataframe containing the transitions and their probabilities
#' @export
#' @examples
#' m = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = mk_infer_cross_sectional(m)
#' mk_pull_transitions(fit)
#'
mk_pull_transitions = function(fitted.obj, reversible=TRUE) {
  # deal either with the castor fit or the wrapped version from mk.inference being passed
  if("fitted_mk" %in% names(fitted.obj)) {
    m = fitted.obj$fitted_mk$transition_matrix
  } else {
    m = fitted.obj$transition_matrix
  }
  mk_df = data.frame()
  if(reversible == TRUE) {
    for(i in 1:(nrow(m))) {
      for(j in 1:(ncol(m))) {
        if(m[i,j] != 0) {
          mk_df = rbind(mk_df, data.frame(From=i-1, To=j-1, Rate=m[i,j]/sum(m[i,seq(1,ncol(m))[seq(1,ncol(m))!=i]])))
        }
      }
    }
  } else {
    # loop over transitions
    for(i in 1:(nrow(m)-1)) {
      ## FIXME: reversible case returns -1 for self -> self transitions
      ## but irreversible doesn't.
      ## If next loop went from i:ncol(m) it would also return -1 in self-transitions,
      ## like reversible case.
      ## Alternatively, for reversible, delete all rows where From == To.
      for(j in (i+1):ncol(m)) {
        if(m[i,j] != 0) {
          mk_df = rbind(mk_df, data.frame(From=i-1, To=j-1, Rate=m[i,j]/sum(m[i,(i+1):ncol(m)])))
        }
      }
    }
  }
  return(mk_df)
}

## FIXME: mk_pull_transitions has a single function, with argument reversible.

## In that function, it should be possible to add some minimal error checking
## - if reversible = FALSE and any non-zero lower-triangular: stop
## - if reversible = TRUE  and no  non-zero lower-triangular: give warning

#' Simulate evolutionary trajectories on a fitted Mk model
#'
#' @param fitted.obj (hypermk model fit) the fitted model
#' @param L (default read from model fit) number of features. The ability to specify this is required for a lower-level use of this function and shouldn't be necessary.
#' @param reversible (default read from model fit) reversible? The ability to specify this is required for a lower-level use of this function and shouldn't be necessary.
#' @param nwalker (integer, default 1e4) number of walkers to simulate
#'
#' @return a dataframe containing the simulated transitions
#' @export
#' @examples
#' m = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = mk_infer_cross_sectional(m)
#' mk_simulate_fluxes(fit)
#'
mk_simulate_fluxes = function(fitted.obj, L = 0, reversible = FALSE, nwalker = 10000) {
  # deal either with the castor fit or the wrapped version from mk.inference being passed
  if("fitted_mk" %in% names(fitted.obj)) {
    L = fitted.obj$L
    fit.mk = fitted.obj$fitted_mk
    reversible = fitted.obj$reversible
  } else {
    fit.mk = fitted.obj
  }
  # to get flux matrix we'll simulate random walkers on the transition matrix
  if(reversible == TRUE) {
    # set up data frame containing transitions and fluxes
    mk.rev.df = mk_pull_transitions(fitted.obj, reversible=TRUE)
    mk.rev.df = mk.rev.df[mk.rev.df$From != mk.rev.df$To,]
    mk.rev.df$Rate[mk.rev.df$Rate == Inf] = 10*max(mk.rev.df$Rate[mk.rev.df$Rate != Inf])
    mk.rev.df$Flux = 0
    # simulate walkers starting from 0^L
    for(walk in 1:nwalker) {
      state = 0
      for(t in 1:(2*L)) {
        # pull row of transition matrix
        trans = fit.mk$transition_matrix[state+1,]
        trans[state+1] = 0
        if(sum(trans)==0) {break}
        # sample the next transition
        this.out = sample(0:(2**L-1), size=1, prob=trans)

        ref = which(mk.rev.df$From == state & mk.rev.df$To == this.out)
        mk.rev.df$Flux[ref] = mk.rev.df$Flux[ref]+1
        state = this.out
      }
    }
    return(mk.rev.df)
  }
  else {

    mk.irrev.df = mk_pull_transitions(fitted.obj, reversible=FALSE)
    mk.irrev.df = mk.irrev.df[mk.irrev.df$From != mk.irrev.df$To,]
    mk.irrev.df$Rate[mk.irrev.df$Rate == Inf] = 10*max(mk.irrev.df$Rate[mk.irrev.df$Rate != Inf])
    mk.irrev.df$Flux = 0
    # simulate walkers starting from 0^L
    for(walk in 1:nwalker) {
      state = 0
      for(t in 1:L) {
        # get possible out transitions
        outs = which(mk.irrev.df$From == state)
        if(length(outs) == 0) {
          break
        } else if(length(outs) == 1) {
          this.out = outs[1]
        } else {
          # sample a possible out transition
          this.out = sample(outs, size=1, prob=mk.irrev.df$Rate[outs])
        }
        mk.irrev.df$Flux[this.out] = mk.irrev.df$Flux[this.out]+1
        state = mk.irrev.df$To[this.out]
      }
    }
    return(mk.irrev.df)
  }
}

#' Cast cross-sectional data into appropriate format for Mk model
#'
#' @param state.list (matrix) observation data. Can be binary or decimal (see below)
#' @param L (integer) number of features to consider
#' @param decimal.labels (Boolean, default FALSE) if FALSE, interpret state.list as a set of rows,
#' each containing a binary description of the state. if TRUE, interpret state.list as 1-indexed
#' decimal representations of states.
#'
#' @return a named list containing a "tree" structure and priors on the tips (as required for Mk model)
#' @export
#' @examples
#' m = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' dataset = mk_cross_sectional(m, 3)
#'
mk_cross_sectional = function(state.list, L, decimal.labels = FALSE) {
  if(decimal.labels == FALSE) {
    state.list = apply(state.list, 1, BinToDec)
    tip.states = state.list+1
  } else {
    tip.states = state.list
  }
  my.tree = tip.priors = vector("list", length(state.list))
  base.tree = ape::stree(2, type = "star")
  # initialise prior matrix for each tree with uniform prior over second tips
  zero.mat = matrix(0, nrow=2, ncol=2**L)
  zero.mat[2,] = 1/(2**L)
  for(i in 1:length(my.tree)) {
    my.tree[[i]] = base.tree
    tip.priors[[i]] = zero.mat
  }
  # see comment above. now we construct a list of 2-tip trees, one for each observation
  mk.tree = my.tree
  # example set of tip states


  # Minimal error checking
  if (any(tip.states > (2**L)))
    stop("At least one tip state has a value outside of possible range: ",
         "[0, (2**L) - 1]. \n",
         "(state.list argument: 0-indexed decimal representation)")

  # enforce deterministic prior for each cross-sectional observation
  for(i in 1:length(tip.states)) {
    tip.priors[[i]][1,tip.states[i]] = 1
  }
  cs.out = list(tree=mk.tree,
                tips=tip.priors)
  return(cs.out)
}

#' Cast phylogeny and set of tip labels data into appropriate format for Mk model
#'
#' @param state.list (matrix) observation data. Can be binary or decimal (see below)
#' @param tree (tree) describing the phylogenetic relationship of observations
#' @param L (integer) number of features to consider
#' @param decimal.labels (Boolean, default FALSE) if FALSE, interpret state.list as a set of rows,
#' each containing a binary description of the state. if TRUE, interpret state.list as 1-indexed
#' decimal representations of states.
#'
#' @return a named list containing a "tree" structure and priors on the tips (as required for Mk model)
#' @export
#' @examples
#' m = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' tree = ape::rphylo(3, 1, 1)
#' dataset = mk_phylogeny_precise(m, tree, 3)
mk_phylogeny_precise = function(state.list, tree, L, decimal.labels = FALSE) {
  if(decimal.labels == FALSE) {
    state.list = apply(state.list, 1, BinToDec)
    tip.states = state.list+1
  } else {
    tip.states = state.list
  }
  mk.tree = tree
  mk.tree$tip.label = tip.states

  # Minimal error checking
  if (any(tip.states > (2**L)))
    stop("At least one tip state has a value outside of possible range: ",
         "[0, (2**L) - 1]. \n",
         "(state.list argument: 0-indexed decimal representation)")

  cs.out = list(tree=mk.tree,
                tips=tip.states)
  return(cs.out)
}


## FIXME: it should be possible to write a wrapper for cross-sectional that
##  takes a data set (as 0-index decimal or matrix of subjects-by-alterations)
##  and a TRUE/FALSE for reversible and does all the rest
##  (calls mk_cross_sectional, mk_index_matrix, fit_mk, mk_pull_transitions,
##   mk_simulate_fluxes)
##  But maybe it is not worth it and even hides the key message
##  (cross-sectional can be thought of as a kind of phylog. data and can be
##   analyzed with Mk)


## FIXME: why use rates = TRUE if it seems emphasis is on fluxes,
##        (e.g., AIC.reduced) and rates can give impossible transitions
##        (transitions from states that cannot be reached)

#' Plot a collection of transitions as a graph. Mainly meant for lower-level functionality.
#'
#' @param trans.f (data frame) transitions, with From/To/Flux/Rate columns
#' @param L (integer) number of features to consider
#' @param rates (Boolean, default FALSE) use rates (as opposed to fluxes)?
#'
#' @return a ggraph plot
#' @export
#' @examples
#' m = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = mk_infer_cross_sectional(m)
#' plot_hypercube2(fit$mk_fluxes, fit$L)
#'
plot_hypercube2 = function(trans.f, L, rates=FALSE) {
  trans.f$Change = ""
  for(i in 1:nrow(trans.f)) {
    if(trans.f$From[i] < trans.f$To[i]) { trans.f$Change[i] = paste0("+", L-log2(trans.f$To[i]-trans.f$From[i]))
    } else { trans.f$Change[i] = paste0("-", L-log2(trans.f$From[i]-trans.f$To[i]))
    }
  }
  trans.g = igraph::graph_from_data_frame(trans.f[trans.f$Flux > 0,])
  bs = unlist(lapply(as.numeric(igraph::V(trans.g)$name), DecToBin, len=L))
  igraph::V(trans.g)$binname = bs #V(trans.g)$name
  layers = stringr::str_count(bs, "1")
  if(rates == FALSE) {
    g.flux = ggraph::ggraph(trans.g) +
      ggraph::geom_edge_arc(strength=0.1,ggplot2::aes(alpha=sqrt(Flux), label=Change),
                    arrow=ggplot2::arrow(length=ggplot2::unit(0.2, "inches"), type="closed")) +
      ggraph::geom_node_label(ggplot2::aes(label=binname), size=2.5, alpha=0.7) +
      ggraph::scale_edge_width(limits=c(0,NA))+ ggplot2::theme_void() +
      ggplot2::theme(legend.position="none", plot.title = ggplot2::element_text(size = 10))
    return(g.flux)
  } else {
    g.prob = ggraph::ggraph(trans.g) +
      ggraph::geom_edge_arc(strength=0.1,ggplot2::aes(alpha=sqrt(Rate), label=Change),
                    arrow=ggplot2::arrow(length=ggplot2::unit(0.2, "inches"), type="closed")) +
      ggraph::geom_node_label(ggplot2::aes(label=binname), size=2.5, alpha=0.7) +
      ggraph::scale_edge_width(limits=c(0,NA))+ ggplot2::theme_void() +
      ggplot2::theme(legend.position="none", plot.title = ggplot2::element_text(size = 10))
    return(g.prob)
  }
}

#' Plot the output of a fitted hypermk model (one mk.inference experiment)
#'
#' @param fitted.obj (hypermk model) the fitted model
#' @param flux.threshold (double, default 0) threshold to prune low fluxes
#' @param full.title (Boolean, default FALSE) give full title information
#'
#' @return a ggraph plot
#' @export
#' @examples
#' m = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = mk_infer_cross_sectional(m)
#' mk.inference.plot(fit, full.title=TRUE)
#'
mk.inference.plot = function(fitted.obj, flux.threshold = 0, full.title = FALSE) {
  L = log2(fitted.obj$fitted_mk$Nstates)
  AIC.val = round(fitted.obj$fitted_mk$AIC, digits=1) #AIC.rev - 2*length(which(combined.obj$mk.out.rev$mk_fluxes$Flux==0))
  title.val = paste0("AIC = ", AIC.val, collapse = "")
  if(full.title == TRUE) {
    k = (fitted.obj$fitted_mk$AIC + 2*fitted.obj$fitted_mk$loglikelihood)/2
    title.val = paste0("n = ", fitted.obj$n,
    ", log L = ", round(fitted.obj$fitted_mk$loglikelihood, digits=2),
    ", k = ", k, ",\n",
    "AIC = ", AIC.val,
    ", AICc = ", round(AIC.val + (k**2 + k)/(fitted.obj$n - k - 1), digits=2))
  }
  flux.threshold.scale = flux.threshold*max(fitted.obj$mk_fluxes$Flux)
  return(plot_hypercube2(fitted.obj$mk_fluxes[fitted.obj$mk_fluxes$Flux > flux.threshold.scale,], L) +
           ggplot2::ggtitle(title.val) )
}

#' General wrapper function for HyperMk inference
#'
#' takes tree, number of features, whether to use priors (vs specific states),
#' tip labels (priors or states), and whether to run reversible model or not
#'
#' because setting up the data is a bit awkward, this is mainly meant for lower-level functionality
#'
#' returns the fitted model with summary dataframes of transitions and fluxes
#'
#' @param mk.tree (tree) a phylogenetic tree (which can be used to reflect independence)
#' @param L (integer) number of features to consider
#' @param use.priors (Boolean) whether we are passing prior distributions on tips or specific observations
#' @param tips (matrix) XXX
#' @param reversible (Boolean) reversible model or irreversible?
#' @param optim_max_iterations (integer, default 200) number of iterations before termination
#' @param Ntrials (integer, default 1) number of repeated trials for the inference process
#' @param optim_algorithm ("optim" or "nlminb", default "optim") which algorithm to use
#' @param Nthreads (integer, default 1) number of threads to use
#' @param to.nullify (n x 2 matrix) from-to transitions to set to zero
#' @param flux.samples (integer, default 10000) number of samples to take in characterising flux
#'
#' @return a fitted hypermk model
#' @export
#' @examples
#' m = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' data = mk_cross_sectional(m, 3)
#' fit = mk.inference(data$tree, 3, use.priors=TRUE, data$tips, reversible = TRUE)
#'
mk.inference = function(mk.tree, L, use.priors, tips, reversible,
                        optim_max_iterations = 200, Ntrials = 1,
                        optim_algorithm = c("optim", "nlminb"), Nthreads = 1,
                        to.nullify = matrix(nrow=0, ncol=2),
                        flux.samples = 10000) {
  optim_algorithm = match.arg(optim_algorithm)
  # for cross-sectional data and uncertain data, tips = tip.priors, use.priors = TRUE
  # otherwise, tips = tip.states, use.priors = FALSE

  # construct matrix describing possible transitions
  index_matrix = mk_index_matrix(L, reversible=reversible, to.nullify=to.nullify)

  # do the Mk model fitting
  # remember the (deterministic) prior on the root state! this is important

  message("fitting model")
  if(use.priors == TRUE) {
    # specify priors, rather than precise states, on the tips of the tree
    fitted_mk = castor::fit_mk(mk.tree, 2**L,
                               tip_priors=tips,
                               optim_algorithm = optim_algorithm,
                               rate_model=index_matrix,
                               root_prior=c(1,rep(0, 2**L-1)),
                               optim_max_iterations = optim_max_iterations,
                               Ntrials = Ntrials,
                               Nthreads = Nthreads)
  } else {
    # specify precise states
    fitted_mk = castor::fit_mk(mk.tree, 2**L,
                               tip_states=tips,
                               optim_algorithm = optim_algorithm,
                               rate_model=index_matrix,
                               root_prior=c(1,rep(0, 2**L-1)),
                               optim_max_iterations = optim_max_iterations,
                               Ntrials = Ntrials,
                               Nthreads = Nthreads)
  }

  nsample = length(tips)
  if(fitted_mk$converged != TRUE) {
    message("WARNING: Mk model fit didn't converge!")
  }

  message("simulating fluxes")
  # convert inferred rate matrix into transition set
  mk_df = mk_pull_transitions(fitted_mk, reversible = reversible)
  # and simulate fluxes through this transition set
  mk_fluxes = mk_simulate_fluxes(fitted_mk, L, reversible = reversible, nwalker=flux.samples)

  # return a list of useful info
  l.return = list(fitted_mk = fitted_mk,
                  reversible = reversible,
                  mk_df = mk_df,
                  mk_fluxes = mk_fluxes,
                  tree = mk.tree, n = nsample, reversible = reversible,
                  L = L, tips = tips, use.priors=use.priors)

  return(l.return)
}


#' Simple wrapper to immediately do inference on a matrix of cross-sectional observations
#'
#' @param m (matrix) matrix of binary observations
#' @param reversible (Boolean, default TRUE) reversible model?
#' @param ... other parameters to mk.inference
#'
#' @return a fitted hypermk model
#' @export
#' @examples
#' m = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' fit = mk_infer_cross_sectional(m)
#'
mk_infer_cross_sectional = function(m, reversible = TRUE, ...) {
  L = ncol(m)
  # cast the matrix into a form that fit_mk will take: a collection of binary trees with root 0, one unspecified tip, and one tip corresponding to the observation
  cs.data = mk_cross_sectional(m, L)

  mk.tree = cs.data$tree
  tip.priors = cs.data$tips
  mk.fit =  mk.inference(mk.tree, L,
                         use.priors=TRUE, tip.priors,
                         reversible = reversible, ...)
  return(mk.fit)
}

#' Simple wrapper to immediately do inference on a matrix of observations linked by a given phylogeny
#'
#' @param m (matrix) matrix of binary observations
#' @param tree (tree) phylogeny linking observations
#' @param reversible (Boolean, default TRUE) reversible model?
#' @param ... other parameters to mk.inference
#'
#' @return a fitted hypermk model
#' @export
#' @examples
#'
#' m = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' tree = ape::rphylo(3, 1, 1)
#' fit = mk_infer_phylogenetic(m, tree)
#'
mk_infer_phylogenetic = function(m, tree, reversible = TRUE, ...) {
  L = ncol(m)

  # cast the matrix into a form that fit_mk will take: a phylogeny with tip labels
  phy.data = mk_phylogeny_precise(m, tree, L)

  mk.tree = phy.data$tree
  tip.states = phy.data$tips
  mk.fit =  mk.inference(mk.tree, L,
                         use.priors=FALSE, tip.states,
                         reversible = reversible, ...)
  return(mk.fit)
}

#' Prune a fitted hypermk model by removing zero-flux edges
#'
#' @param fitted.obj (hypermk model) the fitted model
#' @param flux.threshold (double, default 0) threshold to prune low fluxes
#'
#' @return a fitted, pruned model
#' @export
#' @examples
#' m = matrix(c(0,0,1, 0,1,1, 1,1,1), ncol=3, nrow=3)
#' tree = ape::rphylo(3, 1, 1)
#' fit = mk_infer_phylogenetic(m, tree)
#' fit.pruned = mk_prune_model(fit)
#'
mk_prune_model = function(fitted.obj, flux.threshold=0) {
  to.nullify = fitted.obj$mk_fluxes[which(fitted.obj$mk_fluxes$Flux<=flux.threshold),1:2]+1

  pruned.fit = mk.inference(fitted.obj$tree, fitted.obj$L,
                            fitted.obj$use.priors, fitted.obj$tips,
                            to.nullify = to.nullify,
                            reversible = fitted.obj$reversible)
  return(pruned.fit)
}

## FIXME:
## Questions:
##  - why are there sometimes estimates of rates for impossible transitions?
##  - would we want to remove those rates from the mk_pull_transitions
##    output? (i.e., remove anything involving a non-reachable state
##    --- I think this is easy in irreversible case,
##    and doable in the reversible---); though maybe looking at fluxes
##    is good enough? Would removing them make the simulate_fluxes code
##    faster? Would this speed increase even matter?


## FIXME: is there a -1 in the rate of self-transitions in reversible
## model and not irreversible? I understand why they are -1
## but why do something different for reversible and irreversible?
## See note in mk-shared, function mk_pull_transitions.
