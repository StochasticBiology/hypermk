# collection of HyperMk experiments

source("mk-shared.R")
require(parallel)

# populate a named list with data and visualisations corresponding to synthetic test and scientific cases
# argument specifies the particular case to produce
setup.data = function(expt) {
  set.seed(1)
  
  # cross-sectional cases -- use the mk_cross_sectional function to get a set of trees and tip priors for a given matrix
  if(expt == "cross.sectional.single" | 
     expt == "cross.sectional.many" | 
     expt == "cross.sectional.cross" |
     expt == "ovarian") {
    cross.sectional = TRUE
    uncertain = FALSE
    
    if(expt == "cross.sectional.single") {
      # A single cross-sectional observation with binary genotype 101
      # (0-indexed decimal: 5).
      # Not shown in figures in the paper.
      L = 3
      m = matrix(c(1,0,1), ncol=L, byrow=TRUE)
    }
    if(expt == "cross.sectional.many") {
      # Three cross-sectional observations, with binary genotypes
      # 001, 011, 111 (0-indexed decimal: 1, 3, 7, respectively)
      # Not shown in figures in the paper.
      L = 3
      m = matrix(c(0,0,1,
                   0,1,1,
                   1,1,1), ncol=L, byrow = TRUE)
    }
    if(expt == "cross.sectional.cross") {
      # cross-sectional data, supporting two competing pathways (set up for L=3)
      # fig-2.png; Figure 4 of current ms
      L = 4
      
      m = matrix(c(0,0,0,1,
                   0,0,1,1,
                   0,1,1,1,
                   1,0,0,0,
                   1,1,0,0,
                   1,1,1,0), ncol=L, byrow=TRUE)
    }
    
    if(expt == "ovarian") {
      L = 5
      m = read.csv("Data/ovarian.csv", header=FALSE)[,1:L]
    }
    
    # cast the matrix into a form that fit_mk will take: a collection of binary trees with root 0, one unspecified tip, and one tip corresponding to the observation
    cs.data = mk_cross_sectional(m, L)
    
    mk.tree = cs.data$tree
    tip.priors = cs.data$tips
    
    # just for plotting -- record feature sets
    # get barcodes from (converted) 0-indexed decimal states
    # (barcodes are not used for analysis per se, but for plots)
    tip.states = apply(m, 1, BinToDec)+1
    barcodes = unlist(lapply(tip.states-1, DecToBin, L))
    
    # construct tables of observed barcodes and their decimals in the dataset  
    b.stats = as.data.frame(table(barcodes))
    data.plot = ggplot(b.stats, aes(x=barcodes, y=Freq)) + geom_col() +
      theme_light() + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
      xlab("Observations") + ylab("Count") +
      scale_y_continuous(breaks = seq(0, max(b.stats$Freq), by = 1)) 
    data.plot.nb = ggplot() + geom_blank()
  }
  
  ####### random tree with random, single-pathway dynamics
  if(expt == "single" || expt == "single.rev" || expt == "single.uncertain") {
    # single and single.rev in fig-1.png 1, Figure 3 of current ms.
    # single.uncertain in fig-2.png, Figure 4 of current ms.
    L = 5
    cross.sectional = FALSE
    
    # parameterisation for tree construction
    tree.size = 64
    birth.rate = 1
    death.rate = 0.1
    # accumulation rate for features (and loss rate, for reversible setup)
    accumulation.rate = 1.2
    loss.rate = 1
    
    # create random phylogeny with tree.size nodes from birth-death process parameterised as above
    my.tree = ape::rphylo(tree.size, birth=birth.rate, death=death.rate)
    my.tree$node.label = as.character(1:my.tree$Nnode)
    tree.labels = c(my.tree$tip.label, my.tree$node.label)
    
    # generate state for all nodes traversing tree breadth-first
    # and setting the state of the child nodes according to
    # accumulation (and, if reversible, loss.rate) and branch length
    my.root = phangorn::getRoot(my.tree)
    to.do = c(my.root)
    # initialise state list
    x = list()
    x[[my.root]] = rep(0,L)
    # while we still have vertices to simulate
    while(length(to.do) > 0) {
      # initialise a new to-do list for next iteration
      new.to.do = c()
      # loop through each node in current to-do list
      for(i in to.do) {
        this.outgoing.edges = which(my.tree$edge[,1] == i)
        # loop over this node's children
        for(j in this.outgoing.edges) {
          this.child = my.tree$edge[j,2]
          this.branch.length = my.tree$edge.length[j]
          # construct state for this child based on its parent
          x[[this.child]] = x[[i]]
          # find leftmost zero in current state, and change with some probability
          # recall dynamics here are 00000 -> 10000 -> 11000 -> 11100 -> 11110 -> 11111
          ## (see first paragraph of section "Synthetic case studies")
          ref = which(x[[this.child]] == 0)[1]
          if(runif(1) < accumulation.rate*this.branch.length) { x[[this.child]][ref] = 1 } 
          # in the reversible case, allow the leftmost feature ("first feature" in the ms.:
          # second paragraph of "Synthetic case studies") to revert with some probability
          if(expt == "single.rev") {
            if(runif(1) < loss.rate*this.branch.length) { x[[this.child]][1] = 0 }
          }
          # add this child to to state list, and to next iteration's to-do
          new.to.do = c(new.to.do, this.child)
        }
      }
      # update to-do list
      to.do = new.to.do
    }
    
    my.tree$tip.label = x[1:length(my.tree$tip.label)]
    
    # if we have precise observations, construct set of tip states
    if(expt == "single" | expt == "single.rev") { # fig-1.png 1, Figure 3 of current ms.
      uncertain = FALSE
      # assign feature barcodes to tree
      # convert binary tip labels into 1-indexed decimal state refs
      tip.states = unlist(lapply(my.tree$tip.label,BinToDec))+1
      
      mk.tree = my.tree
      mk.tree$tip.label = tip.states
      
      # just for plotting
      # retrieve barcodes from (converted) 0-indexed decimal state refs
      barcodes = unlist(lapply(tip.states-1, DecToBin, L))
      my.tree2 = my.tree  
      my.tree2$tip.label = barcodes
      
      data.plot = ggtree(my.tree2, layout="circular") + geom_tiplab2(size=3)
      data.plot.nb = ggtree(my.tree2, layout="circular", branch.length="none") + geom_tiplab2(size=2)
    }
    
    # modelling uncertain observations, construct set of tip priors
    if(expt == "single.uncertain") { # fig-2.png; Figure 4 of current ms.
      uncertain = TRUE
      fraction_unknown = 0.5  # proportion of observations to make uncertain: 0.5 in the paper
      # initialise with zero probability
      tip.priors = matrix(0, nrow=length(my.tree$tip.label), ncol=2**L)
      my.tree2 = my.tree
      # loop through observations
      for(i in 1:length(my.tree$tip.label)) {
        # 0-indexed decimal state refs
        this.ref = BinToDec(my.tree$tip.label[[i]])
        # convert into 1-indexed decimal state refs for priors
        tip.priors[i,this.ref+1] = 1
        if(runif(1) < fraction_unknown) {
          # otherwise, allow another random state to be compatible with this observation
          # 0-indexed decimal state refs
          other.ref = sample(0:(2**L-1), size = 1)
          # convert into 1-indexed decimal state refs for priors, and split prior probability between original and new randomly selected state
          tip.priors[i,this.ref+1] = 0.5
          tip.priors[i,other.ref+1] = 0.5
          my.tree2$tip.label[i] = paste0(c(my.tree$tip.label[[i]], "/\n", DecToBin(other.ref, L), "?"), collapse="")
        } else {
          my.tree2$tip.label[i] = paste0(c(my.tree$tip.label[[i]]), collapse="")
        }
      }
      mk.tree = my.tree
      data.plot = ggtree(my.tree2, layout="circular") + geom_tiplab2(size=2, lineheight=0.7)
      data.plot.nb = ggtree(my.tree2, layout="circular", branch.length="none") + geom_tiplab2(size=2, lineheight=0.7)
    }
  }
  
  #### this reads tree and barcode data for the tuberculosis set
  if(expt == "TB") {
    L = 5
    cross.sectional = FALSE
    uncertain = FALSE
    my.tree = read.tree("Data/ng.2878-S2.txt")
    my.data = read.csv("Data/tuberculosis-v5-header-19-29.csv")
    
    # prune tips that don't have corresponding observations
    tips.togo = c()
    for(i in 1:length(my.tree$tip.label)) {
      ref = which(my.data$Isolate == my.tree$tip.label[i])
      if(length(ref) != 1) { tips.togo = c(tips.togo, i) }
    }
    mk.tree = ape::drop.tip(my.tree, tips.togo)
    
    # assign tip states based on data set
    tip.states = c()
    for(i in 1:length(mk.tree$tip.label)) {
      ref = which(my.data$Isolate == mk.tree$tip.label[i])
      barcode = paste(as.vector(my.data
                                [ref,2:(L+1)]), collapse="")
      # 1-indexed decimal state refs
      tip.states = c(tip.states, BinToDec(barcode)+1)
    }
    
    # duplicate tree for plotting
    my.tree2 = mk.tree
    # retrieve barcodes from (converted) 0-indexed decimal state refs
    barcodes = unlist(lapply(tip.states-1, DecToBin, L))
    my.tree2$tip.label = barcodes
    
    data.plot = ggtree(my.tree2, layout="circular") + geom_tiplab2(size=3)
    data.plot.nb = ggtree(my.tree2, layout="circular", branch.length="none") + geom_tiplab2(size=1)
  }
  
  # compile information into a return structure
  # for cross-sectional or uncertain data, we specify priors on tips
  # otherwise, we specify precise states
  if(cross.sectional == TRUE | uncertain == TRUE) {
    l.return = list(tree = mk.tree, L = L, 
                    cross.sectional = cross.sectional, uncertain = uncertain,
                    tips = tip.priors,
                    data.plot = data.plot, data.plot.nb = data.plot.nb)
  } else {
    l.return = list(tree = mk.tree, L = L, 
                    cross.sectional = cross.sectional, uncertain = uncertain,
                    tips = tip.states,
                    data.plot = data.plot, data.plot.nb = data.plot.nb)
  }
  
  return(l.return)
}

# this wrapper function for parallelisation takes an argument determining which experiment to run
# and returns a combined object containing the source data, the result of a irreversible fit, and the result of a reversible fit
parallel.fn = function(fork) {
  
  # the set of cases we will consider
  expt.set = c("cross.sectional.single", "cross.sectional.many", # these simple test cases aren't included in the manuscript but can be used for debugging -- uncomment for this
               "single", "single.rev", # fig-1.png; Figure 3 of current ms.
               "single.uncertain", # fig-2.png; Figure 4 of current ms.
               "cross.sectional.cross", # fig-2.png; Figure 4 of current ms.
               "TB", "ovarian") # fig-3.png
  
  expt = expt.set[fork]
  
  # get the data structure for this case
  dset = setup.data(expt)
  
  # for uncertain data, we specify priors on tip states rather than precise states
  # we also do this for cross-sectional data because of the way we coerce this structure into a set of binary trees
  if(dset$cross.sectional == TRUE | dset$uncertain == TRUE) { 
    use.priors = TRUE 
  } else {
    use.priors = FALSE
  }
  
  # irreversible model fit
  print("irreversible")
  mk.out.irrev = mk.inference(dset$tree, dset$L, 
                              use.priors, dset$tips, 
                              reversible = FALSE)
  # reversible model fit
  print("reversible")
  mk.out.rev = mk.inference(dset$tree, dset$L, 
                            use.priors, dset$tips, 
                            reversible = TRUE,
                            optim_max_iterations = 2000)
  
  l.return = list(dset=dset, mk.out.irrev=mk.out.irrev, mk.out.rev=mk.out.rev)
  return(l.return)
}

# prepare three-panel results figure: data, irreversible fit, reversible fit
results.fig = function(combined.obj, label="", flux.threshold.pmax = 0.01, omit.branch.lengths = FALSE) {
  
  if(combined.obj$mk.out.irrev$fitted_mk$converged != TRUE) {
    message("WARNING: irreversible fit didn't converge!")
  }
  if(combined.obj$mk.out.rev$fitted_mk$converged != TRUE) {
    message("WARNING: reversible fit didn't converge!")
  }
  
  graph.df.rev = combined.obj$mk.out.rev$mk_fluxes
  L = combined.obj$dset$L
  
  AIC.rev = combined.obj$mk.out.rev$fitted_mk$AIC
  AIC.rev.reduced = AIC.rev - 2*length(which(combined.obj$mk.out.rev$mk_fluxes$Flux==0))
  title.rev = paste0("reversible fit, simplified AIC ~ ", round(AIC.rev.reduced, digits=2), 
                     " (full ", round(AIC.rev, digits=2), ")", collapse = "")
  flux.threshold.rev = flux.threshold.pmax*max(graph.df.rev$Flux)
  g.rev = plot.hypercube2(graph.df.rev[graph.df.rev$Flux > flux.threshold.rev,], L) +
    ggtitle(title.rev)
  
  AIC.irrev = combined.obj$mk.out.irrev$fitted_mk$AIC
  AIC.irrev.reduced = AIC.irrev - 2*length(which(combined.obj$mk.out.irrev$mk_fluxes$Flux==0))
  title.irrev = paste0("irreversible fit, simplified AIC ~ ", round(AIC.irrev.reduced, digits=2), 
                       " (full ", round(AIC.irrev, digits=2), ")", collapse = "")
  
  graph.df.irrev = combined.obj$mk.out.irrev$mk_fluxes
  flux.threshold.irrev = flux.threshold.pmax*max(graph.df.irrev$Flux)
  g.irrev = plot.hypercube2(graph.df.irrev[graph.df.irrev$Flux > flux.threshold.irrev,], L) +
    ggtitle(title.irrev)
  
  if(omit.branch.lengths == FALSE) {
    g.data = combined.obj$dset$data.plot #+ ggtitle(label)
  } else {
    g.data = combined.obj$dset$data.plot.nb #+ ggtitle(label)
  }
  
  return( ggarrange(g.data, g.irrev, g.rev, nrow = 1, 
                    labels = c("A", "B", "C"), 
                    widths=c(0.8,1,1),
                    label.y=c(1, 0.1, 0.1)) )
}

##### run the set of experiments
nexpts = 8
parallelised.runs <- mcmapply(parallel.fn,
                              fork = 1:nexpts,
                              SIMPLIFY = FALSE,
                              mc.cores = min(detectCores(), nexpts))

# produce and style output plots
pmaxs = c(0.03, 0.03, 0.05, 0.02, 0.02, 0.05, 0.05, 0.05)
obls = rep(FALSE, 8); obls[7] = TRUE
fig.list = list()
sf = 2
for(i in 1:nexpts) {
  fig.list[[i]] = results.fig(parallelised.runs[[i]], omit.branch.lengths = obls[i], flux.threshold.pmax = pmaxs[i])
  png(paste0("expt-", i, ".png"), width=1000*sf, height=350*sf, res=72*sf)
  print(fig.list[[i]])
  dev.off()
}

