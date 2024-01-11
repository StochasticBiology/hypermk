# collection of HyperMk experiments for time benchmarking

source("mk-shared.R")

time.res.df = data.frame()
for(L in c(2, 3, 4, 5)) {
  for(tree.size in c(4, 8, 16, 32, 64)) {
    print(paste0(L, "-", tree.size))
    set.seed(1)
    
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
          if(runif(1) < loss.rate*this.branch.length) { x[[this.child]][1] = 0 }
          # add this child to to state list, and to next iteration's to-do
          new.to.do = c(new.to.do, this.child)
        }
      }
      # update to-do list
      to.do = new.to.do
    }
    
    # assign feature barcodes to tree   
    my.tree$tip.label = x[1:length(my.tree$tip.label)]
    tip.states = unlist(lapply(my.tree$tip.label,BinToDec))+1
    my.tree2 = my.tree
    my.tree2$tip.label = tip.states
    
    my.pruned = my.tree2
    
    data.plot[[length(data.plot)+1]] = ggtree(my.tree2, layout="circular") + geom_tiplab2(size=1)
    data.plot[[length(data.plot)+1]] = ggtree(my.tree2, layout="circular", branch.length="none") + geom_tiplab2(size=2)
    
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
    
    
    print("doing irreversible model fit")
    # otherwise we have precisely specified tip states
    irrev.time = system.time({
      fitted_mk.irrev = fit_mk(my.pruned, 2**L, tip_states=tip.states, rate_model=index_matrix, root_prior=c(1,rep(0, 2**L-1)))
    })
    
    #### reversible transitions
    index_matrix = mk_index_matrix(L, reversible = TRUE)
    
    ###### reversible transitions
    # do the Mk model fitting
    # remember the (deterministic) prior on the root state! this is important
    
    print("doing reversible model fit")
    # otherwise we have precisely specified tip states
    rev.time = system.time({ fitted_mk.rev = fit_mk(my.pruned, 2**L, 
                                                    tip_states=tip.states, 
                                                    rate_model=index_matrix, 
                                                    root_prior=c(1,rep(0, 2**L-1)))
    })
    
    time.res.df = rbind(time.res.df, data.frame(L=L, tree.size=tree.size, fit="reversible", time=rev.time[3]))
    time.res.df = rbind(time.res.df, data.frame(L=L, tree.size=tree.size, fit="irreversible", time=irrev.time[3]))
  }
}

# cross-sectional examples
for(L in c(3, 4, 5)) {
  for(n.states in c(4, 8, 16, 32, 64)) {
    print(paste0(L, "-", n.states))
    states = round(runif(n.states, min=0, max=2**L-1))
    mk.data = mk_cross_sectional(states, L)
    index_matrix = mk_index_matrix(L, reversible=FALSE)
    irrev.time = system.time({ fitted_mk.irrev = fit_mk(mk.data$tree, 2**L, 
                                                        tip_priors=mk.data$tips, 
                                                        rate_model=index_matrix, 
                                                        root_prior=c(1,rep(0, 2**L-1)))
    })
    index_matrix = mk_index_matrix(L, reversible=TRUE)
    rev.time = system.time({ fitted_mk.rev = fit_mk(mk.data$tree, 2**L, 
                                                    tip_priors=mk.data$tips, 
                                                    rate_model=index_matrix, 
                                                    root_prior=c(1,rep(0, 2**L-1)))
    })
    time.res.df = rbind(time.res.df, data.frame(L=L, tree.size=n.states, fit="cs.reversible", time=rev.time[3]))
    time.res.df = rbind(time.res.df, data.frame(L=L, tree.size=n.states, fit="cs.irreversible", time=irrev.time[3]))
  }
}

sf = 2
png("mk-timing.png", width=600*sf, height=300*sf, res=72*sf)
ggplot(time.res.df, aes(x=L, y=time, color=factor(tree.size))) + 
  geom_line(size=3) + facet_wrap(~fit) +
  scale_y_log10() + theme_light() + ylab("Time / s") + labs(color="Observations")
dev.off()

