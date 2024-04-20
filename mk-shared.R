# shared library loads and helper function for Mk model project

library(ape)
library(phangorn)
library(castor)
library(igraph)
library(ggplot2)
library(ggraph)
library(stringr)
library(ggpubr)
library(ggtree)

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

# construct matrix assigning indices to elements describing possible transitions
mk_index_matrix = function(L, reversible=TRUE) {
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
        index_matrix[i+1,partner+1] = index
        index = index+1
        if(reversible == TRUE) {
          # add reverse transition if required
          index_matrix[partner+1,i+1] = index
          index = index+1
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


# extract transitions as data frame from fitted Mk model
mk_pull_transitions = function(fit.mk, reversible=TRUE) {
  m = fit.mk$transition_matrix
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

mk_simulate_fluxes = function(fit.mk, L, reversible=TRUE, nwalker = 10000) {
  # to get flux matrix we'll simulate random walkers on the transition matrix
  if(reversible == TRUE) {
    # set up data frame containing transitions and fluxes
    mk.rev.df = mk_pull_transitions(fit.mk, reversible=TRUE)
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
    
    mk.irrev.df = mk_pull_transitions(fit.mk, reversible=FALSE)
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

# cast cross-sectional data into appropriate format for Mk model
# requires matrix of binary states or 1-indexed, decimal state representation; indexes from 1 for Mk fit
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

# cast phylogeny + set of tip labels data into appropriate format for Mk model
# requires matrix of binary states or 1-indexed, decimal state representation; indexes from 1 for Mk fit
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



# adapted HyperTraPS plot for HyperHMM outputs
plot.hypercube = function(trans.p, bigL, node.labels = TRUE, use.probability = FALSE) {
  ### produce hypercube subgraph
  trans.g = igraph::graph_from_data_frame(trans.p[!is.infinite(trans.p$Probability),])
  bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
  V(trans.g)$binname = V(trans.g)$name
  layers = stringr::str_count(bs, "1")
  this.plot =  ggraph(trans.g, layout="sugiyama", layers=layers) + 
    geom_edge_arc(aes(edge_width=Probability, edge_alpha=Probability), strength=0.1) +
    scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,0.5)) 
  if(node.labels == TRUE) {
    this.plot = this.plot + geom_node_label(aes(label=binname),angle=45,size=2) 
  }
  return(this.plot)
}


## FIXME: why use rates = TRUE if it seems emphasis is on fluxes,
##        (e.g., AIC.reduced) and rates can give impossible transitions
##        (transitions from states that cannot be reached)
plot.hypercube2 = function(trans.f, bigL, rates=FALSE) {
  trans.f$Change = ""
  for(i in 1:nrow(trans.f)) {
    if(trans.f$From[i] < trans.f$To[i]) { trans.f$Change[i] = paste0("+", bigL-log2(trans.f$To[i]-trans.f$From[i]))  
    } else { trans.f$Change[i] = paste0("-", bigL-log2(trans.f$From[i]-trans.f$To[i])) 
    }
  }
  trans.g = igraph::graph_from_data_frame(trans.f[trans.f$Flux > 0,])
  bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
  V(trans.g)$binname = bs #V(trans.g)$name
  layers = str_count(bs, "1")
  if(rates == FALSE) {
    g.flux = ggraph(trans.g) + 
      geom_edge_arc(strength=0.1,aes(alpha=sqrt(Flux), label=Change), 
                    arrow=arrow(length=unit(0.2, "inches"), type="closed")) +
      geom_node_label(aes(label=binname), size=2.5, alpha=0.7) + 
      scale_edge_width(limits=c(0,NA))+ theme_void() +
      theme(legend.position="none", plot.title = element_text(size = 10))
    return(g.flux)
  } else {
    g.prob = ggraph(trans.g) + 
      geom_edge_arc(strength=0.1,aes(alpha=sqrt(Rate), label=Change), 
                    arrow=arrow(length=unit(0.2, "inches"), type="closed")) +
      geom_node_label(aes(label=binname), size=2.5, alpha=0.7) + 
      scale_edge_width(limits=c(0,NA))+ theme_void() +
      theme(legend.position="none", plot.title = element_text(size = 10))
    return(g.prob)
  }
}

# general wrapper function for HyperMk inference
# takes tree, number of features, whether to use priors (vs specific states), tip labels (priors or states), and whether to run reversible model or not
# returns the fitted model with summary dataframes of transitions and fluxes
mk.inference = function(mk.tree, L, use.priors, tips, reversible, optim_max_iterations = 200, Ntrials = 1, optim_algorithm = c("optim", "nlminb"), Nthreads = 1)
{
  optim_algorithm = match.arg(optim_algorithm)
  # for cross-sectional data and uncertain data, tips = tip.priors, use.priors = TRUE
  # otherwise, tips = tip.states, use.priors = FALSE
  
  # construct matrix describing possible transitions
  index_matrix = mk_index_matrix(L, reversible=reversible)
  
  # do the Mk model fitting
  # remember the (deterministic) prior on the root state! this is important
  
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
  
  if(fitted_mk$converged != TRUE) {
    message("WARNING: Mk model fit didn't converge!")
  }
  
  # convert inferred rate matrix into transition set
  mk_df = mk_pull_transitions(fitted_mk, reversible = reversible)
  # and simulate fluxes through this transition set
  mk_fluxes = mk_simulate_fluxes(fitted_mk, L, reversible = reversible)
  
  # return a list of useful info
  l.return = list(fitted_mk = fitted_mk, 
                  reversible = reversible, 
                  mk_df = mk_df, 
                  mk_fluxes = mk_fluxes)  
  
  return(l.return)
}


# simple wrapper to immediately do inference on a matrix of cross-sectional observations
mk_infer_cross_sectional = function(m, reversible = TRUE) {
  L = ncol(m)
  # cast the matrix into a form that fit_mk will take: a collection of binary trees with root 0, one unspecified tip, and one tip corresponding to the observation
  cs.data = mk_cross_sectional(m, L)
  
  mk.tree = cs.data$tree
  tip.priors = cs.data$tips
  mk.fit =  mk.inference(mk.tree, L, 
                         use.priors=TRUE, tip.priors, 
                         reversible = reversible)
  return(mk.fit)
}


