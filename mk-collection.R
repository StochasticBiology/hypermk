# collection of HyperMk experiments

source("mk-shared.R")

graph.df = res.df = b.df = i.df = data.frame()
data.plot = list()
for(expt in c( "single", "single.rev", "single.uncertain", "cross.sectional.single", "cross.sectional.many", "cross.sectional.cross", "TB")) {
  print(expt)
  set.seed(1)
  
  # cross-sectional data, a single sample
  if(expt == "cross.sectional.single") {
    L = 3
    
    # to deal with cross-sectional data we could in principle construct a tree with just a root and a single tip
    # but this seems to challenge the numerics of the fitting code
    # instead we construct a tree with a root and two tips, one of which is specified by our observation
    # and the other is completely ambiguous (uniform prior over all states)
    my.tree = stree(2)
    my.pruned = my.tree
    
    # the tip prior matrix
    tip.priors = matrix(0, nrow=2, ncol=2**L)
    # uniform prior for tip 2
    tip.priors[2,] = 1/(2**L)
    # deterministic prior for tip 1
    tip.priors[1,6] = 1
  }
  
  # cross-sectional data, several samples (set up for L=3)
  if(expt == "cross.sectional.many") {
    L = 3
    
    # see comment above. now we construct a list of 2-tip trees, one for each observation
    my.tree = list(stree(2), stree(2), stree(2))
    my.pruned = my.tree
    # example set of tip states
    tip.states = c(1, 3, 7)+1
    
    # initialise prior matrix for each tree with uniform prior over second tips
    zero.mat = matrix(0, nrow=2, ncol=2**L)
    zero.mat[2,] = 1/(2**L)
    tip.priors = list(zero.mat, zero.mat, zero.mat)
    # enforce deterministic prior for each cross-sectional observation
    for(i in 1:length(tip.states)) {
      tip.priors[[i]][1,tip.states[i]] = 1
    }
  }
  
  # cross-sectional data, supporting two competing pathways (set up for L=3)
  if(expt == "cross.sectional.cross") {
    L = 4
    
    # see comment above. now we construct a list of 2-tip trees, one for each observation
    my.tree = list(stree(2), stree(2), stree(2), stree(2), stree(2), stree(2))
    my.pruned = my.tree
    # example set of tip states
    tip.states = c(1, 3, 7, 8, 12, 14)+1
    
    # initialise prior matrix for each tree with uniform prior over second tips
    zero.mat = matrix(0, nrow=2, ncol=2**L)
    zero.mat[2,] = 1/(2**L)
    tip.priors = list(zero.mat, zero.mat, zero.mat, zero.mat, zero.mat, zero.mat)
    # enforce deterministic prior for each cross-sectional observation
    for(i in 1:length(tip.states)) {
      tip.priors[[i]][1,tip.states[i]] = 1
    }
  }
  
  ####### random tree with random, single-pathway dynamics
  if(expt == "single" || expt == "single.rev" || expt == "single.uncertain") {
    L = 5
    
    # parameterisation for tree construction
    tree.size = 64
    birth.rate = 1
    death.rate = 0.1
    # accumulation rate for features (and loss rate, for reversible setup)
    accumulation.rate = 2
    loss.rate = 1
    
    # create random phylogeny with 2^n nodes from birth-death process parameterised as above
    my.tree = rphylo(tree.size, birth=birth.rate, death=death.rate)
    my.tree$node.label = as.character(1:my.tree$Nnode)
    tree.labels = c(my.tree$tip.label, my.tree$node.label)
    
    my.root = getRoot(my.tree)
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
          ref = which(x[[this.child]] == 0)[1]
          if(runif(1) < accumulation.rate*this.branch.length) { x[[this.child]][ref] = 1 } 
          # in the reversible case, allow the leftmost feature to revert with some probability
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
    
    # if we have precise observations, construct set of tip states
    if(expt == "single" | expt == "single.rev") {
      # assign feature barcodes to tree   
      my.tree$tip.label = x[1:length(my.tree$tip.label)]
      tip.states = unlist(lapply(my.tree$tip.label,BinToDec))+1
      my.tree2 = my.tree
      my.tree2$tip.label = tip.states
      
      my.pruned = my.tree2
      
      data.plot[[length(data.plot)+1]] = ggtree(my.tree2, layout="circular") + geom_tiplab2(size=1)
      data.plot[[length(data.plot)+1]] = ggtree(my.tree2, layout="circular", branch.length="none") + geom_tiplab2(size=2)
    }
    
    # modelling uncertain observations, construct set of tip priors
    if(expt == "single.uncertain") {
      # initialise with zero probability
      my.tree$tip.label = x[1:length(my.tree$tip.label)]
      tip.priors = matrix(0, nrow=length(my.tree$tip.label), ncol=2**L)
      # loop through observations
      for(i in 1:length(my.tree$tip.label)) {
        this.ref = BinToDec(my.tree$tip.label[[i]])
        tip.priors[i,this.ref+1] = 1
        if(runif(1) < 0.5) {
          # otherwise, allow another random state to be compatible with this observation
          other.ref = round(runif(1, min=0, max=2**L-1))
          tip.priors[i,other.ref+1] = 1
        }
      }
      my.pruned = my.tree
    }
  }
  
  #### this reads tree and barcode data for the tuberculosis set
  if(expt == "TB") {
    L = 4
    
    my.tree = read.tree("Data/ng.2878-S2.txt")
    my.data = read.csv("Data/tuberculosis-v5-header-19-29.csv")
    
    # prune tips that don't have corresponding observations
    tips.togo = c()
    for(i in 1:length(my.tree$tip.label)) {
      ref = which(my.data$Isolate == my.tree$tip.label[i])
      if(length(ref) != 1) { tips.togo = c(tips.togo, i) }
    }
    my.pruned = drop.tip(my.tree, tips.togo)
    
    # assign tip states based on data set
    tip.states = c()
    for(i in 1:length(my.pruned$tip.label)) {
      ref = which(my.data$Isolate == my.pruned$tip.label[i])
      barcode = paste(as.vector(my.data
                                [ref,2:(L+1)]), collapse="")
      tip.states = c(tip.states, BinToDec(barcode)+1)
    }
    
    data.plot[[length(data.plot)+1]] = ggtree(my.tree2, layout="circular") + geom_tiplab2(size=1)
    data.plot[[length(data.plot)+1]] = ggtree(my.tree2, layout="circular", branch.length="none") + geom_tiplab2(size=2)
    #    ggtree(my.tree2, layout="circular") + geom_tiplab2(size=1)
    #    ggtree(my.tree2, layout="circular", branch.length="none") + geom_tiplab2(size=2)
  }
  
  ############# inference section
  # now we have tip.states (feature sets) and my.pruned (tree)
  
  # record feature sets
  barcodes = unlist(lapply(tip.states-1, DecToBin, L))
  barcodes.numeric = matrix(unlist(lapply(tip.states-1, DecToBinV, L)), ncol=L)
  
  #### irreversible transitions
  # construct matrix describing possible transitions
  index_matrix = mk_index_matrix(L, reversible=FALSE)
  
  ###### irreversible transitions
  # do the Mk model fitting
  # remember the (deterministic) prior on the root state! this is important
  
  # for cross-sectional data we use root_priors because we have the uniform prior on a dummy tip as above
  if(expt == "cross.sectional.single" | 
     expt == "cross.sectional.many" | 
     expt == "cross.sectional.cross" |
     expt == "single.uncertain"){
    print("doing irreversible model fit")
    fitted_mk.irrev = fit_mk(my.pruned, 2**L, 
                             tip_priors=tip.priors, 
                             rate_model=index_matrix, 
                             root_prior=c(1,rep(0, 2**L-1)))
    
  } else {
    print("doing irreversible model fit")
    # otherwise we have precisely specified tip states
    fitted_mk.irrev = fit_mk(my.pruned, 2**L, 
                             tip_states=tip.states, 
                             rate_model=index_matrix, 
                             root_prior=c(1,rep(0, 2**L-1)))
  }
  
  # convert inferred rate matrix into transition set
  mk_df.irrev = mk_pull_transitions(fitted_mk.irrev, reversible = FALSE)
  
  #### reversible transitions
  index_matrix = mk_index_matrix(L, reversible = TRUE)
  
  ###### reversible transitions
  # do the Mk model fitting
  # remember the (deterministic) prior on the root state! this is important
  
  # for cross-sectional data we use root_priors because we have the uniform prior on a dummy tip as above
  if(expt == "cross.sectional.single" | 
     expt == "cross.sectional.many" | 
     expt == "cross.sectional.cross" |
     expt == "single.uncertain"){
    print("doing reversible model fit")
    fitted_mk.rev = fit_mk(my.pruned, 2**L, 
                           tip_priors=tip.priors, 
                           rate_model=index_matrix, 
                           root_prior=c(1,rep(0, 2**L-1)))
    
  } else {
    print("doing reversible model fit")
    # otherwise we have precisely specified tip states
    fitted_mk.rev = fit_mk(my.pruned, 2**L, 
                           tip_states=tip.states, 
                           rate_model=index_matrix, 
                           root_prior=c(1,rep(0, 2**L-1)))
  }
  
  # convert inferred rate matrix to transition set
  mk_df.rev = mk_pull_transitions(fitted_mk.rev, reversible = TRUE)
  
  # set up data frame containing transitions and fluxes
  mk.rev.df = mk_simulate_fluxes_reversible(fitted_mk.rev)
  mk.irrev.df = mk_simulate_fluxes_irreversible(fitted_mk.irrev)
  
  mk.rev.df$Experiment = expt
  mk.rev.df$Fit = "reversible"
  graph.df = rbind(graph.df, mk.rev.df)
  mk.irrev.df$Experiment = expt
  mk.irrev.df$Fit = "irreversible"
  graph.df = rbind(graph.df, mk.rev.df)
  graph.df = rbind(graph.df, mk.irrev.df)
  
  res.df = rbind(res.df, data.frame(Experiment = expt, 
                                    Fit = "reversible", 
                                    AIC = fitted_mk.rev$AIC, 
                                    AIC.reduced = fitted_mk.rev$AIC-2*length(which(mk.rev.df$Flux==0))))
  
  res.df = rbind(res.df, data.frame(Experiment = expt, 
                                    Fit = "irreversible", 
                                    AIC = fitted_mk.irrev$AIC, 
                                    AIC.reduced = fitted_mk.irrev$AIC-2*length(which(mk.irrev.df$Flux==0))))
  
  # construct tables of observed barcodes and their decimals in the dataset  
  b.stats = as.data.frame(table(barcodes))
  i.stats = as.data.frame(table(tip.states-1))
  b.stats$Experiment = expt
  i.stats$Experiment = expt
  b.df = rbind(b.df, b.stats)
  i.df = rbind(i.df, i.stats)
}

expt.set = c( "single", "single.rev", 
              "single.uncertain", "cross.sectional.single", 
              "cross.sectional.many", "cross.sectional.cross", 
              "TB")
Lset = c(5, 5, 5, 3, 3, 4, 4)

for(i in 1:length(expt.set)) {
  expt = expt.set[i]
  L = Lset[i]
  
  # plot graph, without thresholding by flux
  mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="reversible",]
  g.rev = plot.hypercube2(graph.df[graph.df$Experiment==expt & graph.df$Fit=="reversible",], L, probability=TRUE) +
    ggtitle(paste(c(expt, ", rev fit, AIC ", round(mk.stats$AIC, digits=2), 
                    " or simplified  ", round(mk.stats$AIC.reduced, digits=2)), 
                  collapse=""))
  
  # plot graph, with thresholding by flux
  g.rev.flux = plot.hypercube2(graph.df[graph.df$Experiment==expt & graph.df$Fit=="reversible" & graph.df$Flux > 0,], L)
  
  g.rev.flux2 = g.rev.flux + 
    ggtitle(paste(c(expt, ", rev fit, AIC ", round(mk.stats$AIC, digits=2), 
                    " or simplified ", round(mk.stats$AIC.reduced, digits=2)), 
                  collapse="")) 
  
  # plot graph without pruning by flux
  mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="irreversible",]
  g.irrev = plot.hypercube2(graph.df[graph.df$Experiment==expt & graph.df$Fit=="irreversible",], L, probability=TRUE) +
    ggtitle(paste(c(expt, ", irrev fit, AIC ", round(mk.stats$AIC, digits=2), 
                    " or simplified ", round(mk.stats$AIC.reduced, digits=2)), 
                  collapse=""))
  
  # plot graph with pruning by flux 
  g.irrev.flux = plot.hypercube2(graph.df[graph.df$Experiment==expt & graph.df$Fit=="irreversible" & graph.df$Flux != 0,], L)
  
  g.irrev.flux2 = g.irrev.flux +  
    ggtitle(paste(c(expt, ", irrev fit, AIC ", round(mk.stats$AIC, digits=2), 
                    " or simplified ", round(mk.stats$AIC.reduced, digits=2)), 
                  collapse=""))
  
  # plot these as simple columns
  g.i = ggplot(i.df[i.df$Experiment==expt,], aes(x=Var1,y=Freq)) + geom_col() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  g.b = ggplot(b.df[b.df$Experiment==expt,], aes(x=barcodes,y=Freq)) + geom_col() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # output to file
  fname = paste0("mk-", expt, "-", L, ".png")
  sf = 2
  png(fname, width=800*sf, height=400*sf, res=72*sf)
  print(ggarrange(g.i, g.rev, g.irrev, g.b, g.rev.flux, g.irrev.flux, nrow=2, ncol=3))
  dev.off()
  
  # output to file
  fname = paste0("mk-graphs-", expt, "-", L, ".png")
  sf = 2
  png(fname, width=600*sf, height=400*sf, res=72*sf)
  print(ggarrange(g.rev, g.irrev, g.rev.flux, g.irrev.flux, nrow=2, ncol=2))
  dev.off()
  
  # output to file
  fname = paste0("mk-fluxes-", expt, "-", L, ".png")
  sf = 2
  png(fname, width=600*sf, height=300*sf, res=72*sf)
  print(ggarrange(g.rev.flux2, g.irrev.flux2))
  dev.off()
}

for(i in 1:length(data.plot)) {
  fname = paste0("mk-data-", i, ".png")
  png(fname, width=800*sf, height=400*sf, res=72*sf)
  print(data.plot[[i]])
  dev.off()
}

#p.tree = ggtree(my.pruned, layout="circular", branch.length="none") 
#gheatmap(p.tree, barcodes.numeric, low="white", high="black", colnames_angle=90) + theme(legend.position = "none")

