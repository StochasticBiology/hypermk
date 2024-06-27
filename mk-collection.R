# collection of HyperMk experiments

source("mk-specifics.R")

# the longest case study here (TB) takes several hours on a modern machine to fit one instance of the reversible model
# (that's after subsetting the full dataset)
# using more trials increases the chance we find the global optimum, but will multiply this runtime UNLESS we use more threads (below)
# minimum value is 1 -- use this if computer time is limiting
Ntrials = 1

# For a single run:
# If the computer has multiple cores, using Nthreads will decrease running time
# when Ntrials > 1 (if Nthreads = Ntrials, running time will be about the same as
## for Ntrials = 1). By default, set Nthreads to Ntrials, unless there are fewer
# cores
#Nthreads = min(Ntrials, parallel::detectCores() - 1)
# For the parallelised case studies below, we will assign each case study to
# a different core (if available), so will avoid forking threads within these
# parallel processes
Nthreads = 1

### simple demo of pruning and refitting
# simple synthetic cross-sectional dataset
dset = setup.data("cross.sectional.single")

# unpruned inference; all edges are parameters
mk.out.irrev = mk.inference(dset$tree, dset$L, 
                            use.priors=TRUE, dset$tips, 
                            reversible = FALSE,
                            Ntrials = 1)

# prune a specific edge (here, corresponding to one alternative pathway)
blanks = matrix(c(1,2), nrow=1, ncol=2, byrow = TRUE)
mk.out.irrev2 = mk.inference(dset$tree, dset$L, 
                             use.priors=TRUE, dset$tips, 
                             reversible = FALSE,
                             Ntrials = 1,
                             to.nullify=blanks)

# get a set of edges to prune, corresponding to zero-flux edges in the previous example
to.nullify = mk.out.irrev$mk_fluxes[which(mk.out.irrev$mk_fluxes$Flux==0),1:2]+1
mk.out.irrev3 = mk.inference(dset$tree, dset$L, 
                             use.priors=TRUE, dset$tips, 
                             reversible = FALSE,
                             Ntrials = 1,
                             to.nullify=to.nullify)

# graphically compare these (with AICs)
sf = 2
png("fig-pruning.png", width=600*sf, height=400*sf, res=72*sf)
ggarrange(mk.inference.plot(mk.out.irrev),
          mk.inference.plot(mk.out.irrev2),
          mk.inference.plot(mk.out.irrev3), 
          nrow=1, labels=c("A", "B", "C"), label.y=0.1)
dev.off()

##### run the set of experiments
# produce and style output plots
# 1 "cross.sectional.single", 2 "cross.sectional.many", 3 "single", 4 "single.rev", 
# 5 "single.uncertain", 6 "cross.sectional.cross",
# 7 "ovarian", 8 "TB"

# by default, we'll only run the quickest experiments (1-6); increase to 8 for the long haul
nexpts = 8
Ntrials = 1
parallelised.runs <- mcmapply(parallel.fn,
                              fork = 1:nexpts,
                              SIMPLIFY = FALSE,
                              mc.cores = min(detectCores(), nexpts))

pmaxs = c(0.03, 0.03, 0.05, 0.0, 0.02, 0.05, 0.05, 0.05)
pmaxs = rep(0,8)
obls = rep(FALSE, 8); obls[8] = TRUE
fig.list = list()
sf = 2
for(i in 1:nexpts) {    
  message(paste0("Experiment ", i, collapse=""))
  fig.list[[i]] = results.fig(parallelised.runs[[i]], 
                              omit.branch.lengths = obls[i], 
                              flux.threshold.pmax = pmaxs[i])
  png(paste0("expt-pruned-", i, ".png"), width=1000*sf, height=350*sf, res=72*sf)
  print(fig.list[[i]])
  dev.off()
}

expt = 6
results.fig(parallelised.runs[[expt]], omit.branch.lengths = obls[expt], flux.threshold.pmax = pmaxs[expt])

Ntrials = 3
tmp = parallel.fn(2)
results.fig(tmp, flux.threshold.pmax = 0)

##### If you only want to run a subset of the experiments
##   Change the indices passed to argument fork in the call to
##   parallelised runs.
##   The map between experiments and indices can be seen in expt.set in
##   parallel.fn
##   1 -> cross.sectional.single
##   2 -> cross.sectional.many
##   3 -> single
##   4 -> single.rev
##   5 -> single.uncertain
##   6 -> cross.sectional.cross
##   7 -> ovarian
##   8 -> TB

## See the mapping as:
## bpfn <- body(parallel.fn)
## expt.set.values <- eval(bpfn[[grep("expt.set =", bpfn, fixed = TRUE)]][[3]])


##  For example, to run cross.sectional.many and cross.sectional.cross do
##  (code wrapped in FALSE to prevent automatic execution)

if (FALSE) {
  ## Convenience, to pass the experiment by name, and use it in
  ## the figure
  ## Get the value of expt.set in function parallel.fn
  bpfn <- body(parallel.fn)
  expt.set.values <- eval(bpfn[[grep("expt.set =", bpfn, fixed = TRUE)]][[3]])
  the_expts <- setNames(seq_along(expt.set.values), expt.set.values)

  ## Now, specify the experiments you want by name
  these_expts <- c("cross.sectional.many", "cross.sectional.cross")


  nexpts <- length(these_expts)
  parallelised.runs <- mcmapply(parallel.fn,
                                fork = the_expts[these_expts],
                                SIMPLIFY = FALSE,
                                mc.cores = min(detectCores(), nexpts))

  pmaxs = c(0.03, 0.03, 0.05, 0.02, 0.02, 0.05, 0.05, 0.05)
  obls = rep(FALSE, 8); obls[7] = TRUE
  fig.list = list()
  sf = 2

  for(i in 1:nexpts) {
    fig.list[[i]] = results.fig(parallelised.runs[[i]], omit.branch.lengths = obls[i], flux.threshold.pmax = pmaxs[i])
    # png(paste0("expt-", i, ".png"), width=1000*sf, height=350*sf, res=72*sf)
    png(paste0("expt-", these_expts[i], ".png"), width=1000*sf, height=350*sf, res=72*sf)
    print(fig.list[[i]])
    dev.off()
  }

}
