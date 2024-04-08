# collection of HyperMk experiments

source("mk-shared.R")
# FIXME: "cross.sectional.single", and "cross.sectional.many" are
#        not used anywhere?
graph.df = res.df = b.df = i.df = data.frame()
data.plot = data.plot.nb = list()
for(expt in c( "single", "single.rev", # fig-1.png; Figure 3 of current ms.
              "single.uncertain", # fig-2.png; Figure 4 of current ms.
              "cross.sectional.single", # Not shown in figures
              "cross.sectional.many", # Not shown in figures
              "cross.sectional.cross", # fig-2.png; Figure 4 of current ms.
              "TB")) {
    print(expt)
    set.seed(1)
    
    # cross-sectional data, a single sample
    if(expt == "cross.sectional.single") {
      # A single cross-sectional observation with binary genotype 101
      # (0-indexed decimal: 5).
      # Not shown in figures in the paper.
      L = 3
      
      # to deal with cross-sectional data we could in principle construct a tree with just a root and a single tip
      # but this seems to challenge the numerics of the fitting code
      # instead we construct a tree with a root and two tips, one of which is specified by our observation
      # and the other is completely ambiguous (uniform prior over all states)
      my.tree = ape::stree(2)
      my.pruned = my.tree
      
      # the tip prior matrix
      tip.priors = matrix(0, nrow=2, ncol=2**L)
      # uniform prior for tip 2
      tip.priors[2,] = 1/(2**L)
      # deterministic prior for tip 1
      # 0-indexed decimal state is 5: only 6th column is a 1.
      tip.priors[1,6] = 1 
    }
    
    # cross-sectional data, several samples (set up for L=3)
    if(expt == "cross.sectional.many") {
      # Three cross-sectional observations, with binary genotypes
      # 001, 011, 111 (0-indexed decimal: 1, 3, 7, respectively)
      # Not shown in figures in the paper.
      L = 3
      
      # see comment above. now we construct a list of 2-tip trees, one for each observation
      my.tree = list(ape::stree(2), ape::stree(2), ape::stree(2))
      my.pruned = my.tree
      # example set of tip states
      # 1-indexed decimal states
      tip.states = c(1, 3, 7)+1
      
      # initialise prior matrix for each tree with uniform prior over second tips
      zero.mat = matrix(0, nrow=2, ncol=2**L)
      zero.mat[2,] = 1/(2**L)
      tip.priors = list(zero.mat, zero.mat, zero.mat)
      # enforce deterministic prior for each cross-sectional observation
      for(i in 1:length(tip.states)) {
        # 0-indexed decimal states of observations are 1, 3, 7:
        # columns 2, 4, 8, of tip.priors 1, 2, 3, are set to 1; the rest are left at 0.
        tip.priors[[i]][1,tip.states[i]] = 1
      }
      
      # record feature sets
      # get barcodes from (converted) 0-indexed decimal states
      # (barcodes are not used for analysis per se, but for plots)
      barcodes = unlist(lapply(tip.states-1, DecToBin, L))
      barcodes.numeric = matrix(unlist(lapply(tip.states-1, DecToBinV, L)), ncol=L)
      
      # construct tables of observed barcodes and their decimals in the dataset  
      b.stats = as.data.frame(table(barcodes))
      data.plot[[expt]] = ggplot(b.stats, aes(x=barcodes, y=Freq)) + geom_col() +
        theme_light() + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        xlab("Observations") + ylab("Count") +
        scale_y_continuous(breaks = seq(0, max(b.stats$Freq), by = 1)) 
      
    }
    
    # cross-sectional data, supporting two competing pathways (set up for L=3)
    if(expt == "cross.sectional.cross") {
      # fig-2.png; Figure 4 of current ms
      L = 4
      
      # see comment above. now we construct a list of 2-tip trees, one for each observation
      my.tree = list(ape::stree(2), ape::stree(2), ape::stree(2),
                     ape::stree(2), ape::stree(2), ape::stree(2))
      my.pruned = my.tree
      # example set of tip states
      # 1-indexed decimal states
      # binary genotypes are
      # 0001, 0011, 0111,  (0-indexed decimal: 1, 3, 7)
      # 1000, 1100, 1110   (0-indexed decimal: 8, 12, 14)
      tip.states = c(1, 3, 7, 8, 12, 14)+1
      
      # initialise prior matrix for each tree with uniform prior over second tips
      zero.mat = matrix(0, nrow=2, ncol=2**L)
      zero.mat[2,] = 1/(2**L)
      tip.priors = list(zero.mat, zero.mat, zero.mat, zero.mat, zero.mat, zero.mat)
      # enforce deterministic prior for each cross-sectional observation
      # (deterministic prior is placed on the first tip ---uniform prior
      # was placed above on the second tip)
      for(i in 1:length(tip.states)) {
        tip.priors[[i]][1,tip.states[i]] = 1
      }
      
      # record feature sets
      # get barcodes from (converted) 0-indexed decimal states
      # (barcodes are not used for analysis per se, but for plots)
      barcodes = unlist(lapply(tip.states-1, DecToBin, L))
      barcodes.numeric = matrix(unlist(lapply(tip.states-1, DecToBinV, L)), ncol=L)
      
      # construct tables of observed barcodes and their decimals in the dataset  
      b.stats = as.data.frame(table(barcodes))
      data.plot[[expt]] = ggplot(b.stats, aes(x=barcodes, y=Freq)) + geom_col() +
        theme_light() + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        xlab("Observations") + ylab("Count") +
        scale_y_continuous(breaks = seq(0, max(b.stats$Freq), by = 1)) 
    }
    
####### random tree with random, single-pathway dynamics
    if(expt == "single" || expt == "single.rev" || expt == "single.uncertain") {
      # single and single.rev in fig-1.png 1, Figure 3 of current ms.
      # single.uncertain in fig-2.png, Figure 4 of current ms.
      L = 5
      
      # parameterisation for tree construction
      tree.size = 64
      birth.rate = 1
      death.rate = 0.1
      # accumulation rate for features (and loss rate, for reversible setup)
      # FIXME: these rates are the same. What is the purpose of the if?
      if(expt == "single.rev") {
        accumulation.rate = 1.2
      } else {
        accumulation.rate = 1.2
      }
      loss.rate = 1
      
      # create random phylogeny with 2^n nodes from birth-death process parameterised as above
      my.tree = ape::rphylo(tree.size, birth=birth.rate, death=death.rate)
      my.tree$node.label = as.character(1:my.tree$Nnode)
      tree.labels = c(my.tree$tip.label, my.tree$node.label)
      
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
      if(expt == "single" | expt == "single.rev") { # fig-1.png 1, Figure 3 of current ms.
        # assign feature barcodes to tree   
        my.tree$tip.label = x[1:length(my.tree$tip.label)]
        # convert binary tip labels into 1-indexed decimal state refs
        tip.states = unlist(lapply(my.tree$tip.label,BinToDec))+1
        my.tree2 = my.tree
        my.pruned = my.tree
        # retrieve barcodes from (converted) 0-indexed decimal state refs
        barcodes = unlist(lapply(tip.states-1, DecToBin, L))
        my.pruned$tip.label = tip.states
        my.tree2$tip.label = barcodes
        
        data.plot[[expt]] = ggtree(my.tree2, layout="circular") + geom_tiplab2(size=3)
        data.plot.nb[[expt]] = ggtree(my.tree2, layout="circular", branch.length="none") + geom_tiplab2(size=2)
      }
      
      # modelling uncertain observations, construct set of tip priors
      if(expt == "single.uncertain") { # fig-2.png; Figure 4 of current ms.
        # initialise with zero probability
        my.tree$tip.label = x[1:length(my.tree$tip.label)]
        tip.priors = matrix(0, nrow=length(my.tree$tip.label), ncol=2**L)
        my.tree2 = my.tree
        # loop through observations
        for(i in 1:length(my.tree$tip.label)) {
          # 0-indexed decimal state refs
          this.ref = BinToDec(my.tree$tip.label[[i]])
          # convert into 1-indexed decimal state refs for priors
          tip.priors[i,this.ref+1] = 1
          if(runif(1) < 0.5) {
            # otherwise, allow another random state to be compatible with this observation
            # 0-indexed decimal state refs
            other.ref = round(runif(1, min=0, max=2**L-1))
            # convert into 1-indexed decimal state refs for priors
            tip.priors[i,other.ref+1] = 1
            my.tree2$tip.label[i] = paste0(c(my.tree$tip.label[[i]], "/\n", DecToBin(other.ref, L), "?"), collapse="")
          } else {
            my.tree2$tip.label[i] = paste0(c(my.tree$tip.label[[i]]), collapse="")
          }
        }
        my.pruned = my.tree
        data.plot[[expt]] = ggtree(my.tree2, layout="circular") + geom_tiplab2(size=2, lineheight=0.7)
        data.plot.nb[[expt]] = ggtree(my.tree2, layout="circular", branch.length="none") + geom_tiplab2(size=2, lineheight=0.7)
      }
    }
    
#### this reads tree and barcode data for the tuberculosis set
    if(expt == "TB") {
        L = 5
        
        my.tree = read.tree("Data/ng.2878-S2.txt")
        my.data = read.csv("Data/tuberculosis-v5-header-19-29.csv")
        
        # prune tips that don't have corresponding observations
        tips.togo = c()
        for(i in 1:length(my.tree$tip.label)) {
            ref = which(my.data$Isolate == my.tree$tip.label[i])
            if(length(ref) != 1) { tips.togo = c(tips.togo, i) }
        }
        my.pruned = ape::drop.tip(my.tree, tips.togo)
        
        # assign tip states based on data set
        tip.states = c()
        for(i in 1:length(my.pruned$tip.label)) {
            ref = which(my.data$Isolate == my.pruned$tip.label[i])
            barcode = paste(as.vector(my.data
                                      [ref,2:(L+1)]), collapse="")
                                        # 1-indexed decimal state refs
            tip.states = c(tip.states, BinToDec(barcode)+1)
        }
        
        my.tree2 = my.pruned
        # retrieve barcodes from (converted) 0-indexed decimal state refs
        barcodes = unlist(lapply(tip.states-1, DecToBin, L))
        my.tree2$tip.label = barcodes
        
        data.plot[[expt]] = ggtree(my.tree2, layout="circular") + geom_tiplab2(size=3)
        data.plot.nb[[expt]] = ggtree(my.tree2, layout="circular", branch.length="none") + geom_tiplab2(size=1)
    }
    
############# inference section
    # now we have tip.states (feature sets) and my.pruned (tree)
    
    #### irreversible transitions
    # construct matrix describing possible transitions
    index_matrix_irrev = mk_index_matrix(L, reversible=FALSE)
    
    ###### irreversible transitions
    # do the Mk model fitting
    # remember the (deterministic) prior on the root state! this is important

    # FIXME: I do not understand this. Isn't this the exact same thing
    # that is done in mk-cgh.R? The "standard protocol" for cross-sectional data?
    # for cross-sectional data we use root_priors because we have the uniform prior on a dummy tip as above
    if(expt == "cross.sectional.single" | 
       expt == "cross.sectional.many" | 
       expt == "cross.sectional.cross" |
       expt == "single.uncertain"){
        print("doing irreversible model fit")
        fitted_mk.irrev = castor::fit_mk(my.pruned, 2**L, 
                                         tip_priors=tip.priors, 
                                         rate_model=index_matrix_irrev, 
                                         root_prior=c(1,rep(0, 2**L-1)))
        
    } else { # single, single.rev, TB
      print("doing irreversible model fit")
      # otherwise we have precisely specified tip states
      fitted_mk.irrev = castor::fit_mk(my.pruned, 2**L, 
                                       tip_states=tip.states, 
                                       rate_model=index_matrix_irrev, 
                                       root_prior=c(1,rep(0, 2**L-1)))
    }
    
    # convert inferred rate matrix into transition set
    mk_df.irrev = mk_pull_transitions(fitted_mk.irrev, reversible = FALSE)
    
#### reversible transitions
    index_matrix_rev = mk_index_matrix(L, reversible = TRUE)
    
    ###### reversible transitions
    # do the Mk model fitting
    # remember the (deterministic) prior on the root state! this is important
    
    # for cross-sectional data we use root_priors because we have the uniform prior on a dummy tip as above
    if(expt == "cross.sectional.single" | 
         expt == "cross.sectional.many" | 
         expt == "cross.sectional.cross" |
         expt == "single.uncertain"){
      print("doing reversible model fit")
      fitted_mk.rev = castor::fit_mk(my.pruned, 2**L, 
                                     tip_priors=tip.priors, 
                                     rate_model=index_matrix_rev, 
                                     root_prior=c(1,rep(0, 2**L-1)))
      
    } else { # single, single.rev, TB
      print("doing reversible model fit")
      # otherwise we have precisely specified tip states
      fitted_mk.rev = castor::fit_mk(my.pruned, 2**L, 
                                     tip_states=tip.states, 
                                     rate_model=index_matrix_rev, 
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
    
}

######## fig-1.png; Figure 3 of current ms.
L = 5
flux.threshold.pmax = 0.01
expt = "single"
this.g.df = graph.df[graph.df$Experiment==expt & graph.df$Fit=="irreversible",]
mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="irreversible",]
flux.threshold = flux.threshold.pmax*max(this.g.df$Flux)
t.str = titlestr(expt, "irrev fit", mk.stats)
g.1 = plot.hypercube2(this.g.df[this.g.df$Flux > flux.threshold,], L) + 
  ggtitle(t.str) 

this.g.df = graph.df[graph.df$Experiment==expt & graph.df$Fit=="reversible",]
mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="reversible",]
flux.threshold = flux.threshold.pmax*max(this.g.df$Flux)
t.str = titlestr(expt, "rev fit", mk.stats)
g.2 = plot.hypercube2(this.g.df[this.g.df$Flux > flux.threshold,], L) + 
  ggtitle(t.str) 

expt = "single.rev"
this.g.df = graph.df[graph.df$Experiment==expt & graph.df$Fit=="irreversible",]
mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="irreversible",]
flux.threshold = flux.threshold.pmax*max(this.g.df$Flux)
t.str = titlestr(expt, "irrev fit", mk.stats)
g.3 = plot.hypercube2(this.g.df[this.g.df$Flux > flux.threshold,], L) + 
  ggtitle(t.str) 

this.g.df = graph.df[graph.df$Experiment==expt & graph.df$Fit=="reversible",]
mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="reversible",]
flux.threshold = flux.threshold.pmax*max(this.g.df$Flux)
t.str = titlestr(expt, "rev fit", mk.stats)
g.4 = plot.hypercube2(this.g.df[this.g.df$Flux > flux.threshold,], L) + 
  ggtitle(t.str) 

g.fig.1 = ggarrange(data.plot[["single"]], g.1, g.2,
                    data.plot[["single.rev"]], g.3, g.4, 
                    nrow=2, ncol=3, labels=c("A", "B", "C", "D", "E", "F"),
                    label.y=c(1,0.9,0.9, 1,0.9,0.9))

sf = 2
png("fig-1.png", width=800*sf, height=600*sf, res=72*sf)
print(g.fig.1)
dev.off()

######## fig-2.png; Figure 4 of current ms.
L = 5
flux.threshold.pmax = 0.01
expt = "single.uncertain"
this.g.df = graph.df[graph.df$Experiment==expt & graph.df$Fit=="irreversible",]
mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="irreversible",]
flux.threshold = flux.threshold.pmax*max(this.g.df$Flux)
t.str = titlestr(expt, "irrev fit", mk.stats)
g.1 = plot.hypercube2(this.g.df[this.g.df$Flux > flux.threshold,], L) + 
  ggtitle(t.str) 

this.g.df = graph.df[graph.df$Experiment==expt & graph.df$Fit=="reversible",]
mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="reversible",]
flux.threshold = flux.threshold.pmax*max(this.g.df$Flux)
t.str = titlestr(expt, "rev fit", mk.stats)
g.2 = plot.hypercube2(this.g.df[this.g.df$Flux > flux.threshold,], L) + 
  ggtitle(t.str) 

L = 4
expt = "cross.sectional.cross"
this.g.df = graph.df[graph.df$Experiment==expt & graph.df$Fit=="irreversible",]
mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="irreversible",]
flux.threshold = flux.threshold.pmax*max(this.g.df$Flux)
t.str = titlestr(expt, "irrev fit", mk.stats)
g.3 = plot.hypercube2(this.g.df[this.g.df$Flux > flux.threshold,], L) + 
  ggtitle(t.str) 

this.g.df = graph.df[graph.df$Experiment==expt & graph.df$Fit=="reversible",]
mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="reversible",]
flux.threshold = flux.threshold.pmax*max(this.g.df$Flux)
t.str = titlestr(expt, "rev fit", mk.stats)
g.4 = plot.hypercube2(this.g.df[this.g.df$Flux > flux.threshold,], L) + 
  ggtitle(t.str) 

g.fig.2 = ggarrange(data.plot[["single.uncertain"]], g.1, g.2,
                    data.plot[["cross.sectional.cross"]], g.3, g.4, 
                    nrow=2, ncol=3, labels=c("A", "B", "C", "D", "E", "F"),
                    label.y=c(1,0.9,0.9, 1,0.9,0.9))

sf = 2
png("fig-2.png", width=800*sf, height=600*sf, res=72*sf)
print(g.fig.2)
dev.off()

######## fig-3.png; Figure 5 of current ms.
L = 5
flux.threshold.pmax = 0.05
expt = "TB"
this.g.df = graph.df[graph.df$Experiment==expt & graph.df$Fit=="irreversible",]
mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="irreversible",]
flux.threshold = flux.threshold.pmax*max(this.g.df$Flux)
t.str = titlestr(expt, "irrev fit", mk.stats)
g.1 = plot.hypercube2(this.g.df[this.g.df$Flux > flux.threshold,], L) + 
  ggtitle(t.str) 

this.g.df = graph.df[graph.df$Experiment==expt & graph.df$Fit=="reversible",]
mk.stats = res.df[res.df$Experiment==expt & res.df$Fit=="reversible",]
flux.threshold = flux.threshold.pmax*max(this.g.df$Flux)
t.str = titlestr(expt, "rev fit", mk.stats)
g.2 = plot.hypercube2(this.g.df[this.g.df$Flux > flux.threshold,], L) + 
  ggtitle(t.str) 

g.fig.3 = ggarrange(data.plot.nb[["TB"]], g.1, g.2,
                    nrow=1, ncol=3, labels=c("A", "B", "C"),
                    label.y=c(1,0.1,0.1))

sf = 2
png("fig-3.png", width=800*sf, height=300*sf, res=72*sf)
print(g.fig.3)
dev.off()

for(i in 1:length(data.plot)) {
  fname = paste0("mk-data-", i, ".png")
  png(fname, width=800*sf, height=400*sf, res=72*sf)
  print(data.plot[[i]])
  dev.off()
}

