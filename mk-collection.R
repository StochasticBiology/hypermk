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
n.sample = length(dset$x.decimal)

# unpruned inference; all edges are parameters
mk.out.irrev = mk.inference(dset$tree, dset$L, 
                            use.priors=TRUE, dset$tips, 
                            reversible = FALSE,
                            Ntrials = 1)
mk.out.irrev$myAIC = mk.out.irrev$fitted_mk$AIC
mk.out.irrev$n = n.sample
mk.out.irrev$logL = mk.out.irrev$fitted_mk$loglikelihood
mk.out.irrev$k = (mk.out.irrev$myAIC + 2*mk.out.irrev$logL)/2
mk.out.irrev$myAICc = mk.out.irrev$myAIC - (2*mk.out.irrev$k**2 + 2*mk.out.irrev$k)/(mk.out.irrev$n - mk.out.irrev$k - 1)

# prune a specific edge (here, corresponding to one alternative pathway)
blanks = matrix(c(1,2), nrow=1, ncol=2, byrow = TRUE)
mk.out.irrev2 = mk.inference(dset$tree, dset$L, 
                             use.priors=TRUE, dset$tips, 
                             reversible = FALSE,
                             Ntrials = 1,
                             to.nullify=blanks)
mk.out.irrev2$myAIC = mk.out.irrev2$fitted_mk$AIC
mk.out.irrev2$n = n.sample
mk.out.irrev2$logL = mk.out.irrev2$fitted_mk$loglikelihood
mk.out.irrev2$k = (mk.out.irrev2$myAIC + 2*mk.out.irrev2$logL)/2
mk.out.irrev2$myAICc = mk.out.irrev2$myAIC - (2*mk.out.irrev2$k**2 + 2*mk.out.irrev2$k)/(mk.out.irrev2$n - mk.out.irrev2$k - 1)

# get a set of edges to prune, corresponding to zero-flux edges in the previous example
to.nullify = mk.out.irrev$mk_fluxes[which(mk.out.irrev$mk_fluxes$Flux==0),1:2]+1
mk.out.irrev3 = mk.inference(dset$tree, dset$L, 
                             use.priors=TRUE, dset$tips, 
                             reversible = FALSE,
                             Ntrials = 1,
                             to.nullify=to.nullify)
nonzeroes = length(which(mk.out.irrev3$fitted_mk$transition_matrix!=0))
nonzeroes.diag = length(which(diag(mk.out.irrev3$fitted_mk$transition_matrix)!=0))
mk.out.irrev3$logL = mk.out.irrev3$fitted_mk$loglikelihood
mk.out.irrev3$myAIC = 2*(nonzeroes-nonzeroes.diag)-2*mk.out.irrev3$logL
mk.out.irrev3$n = n.sample
mk.out.irrev3$k = (mk.out.irrev3$myAIC + 2*mk.out.irrev3$logL)/2
mk.out.irrev3$myAICc = mk.out.irrev3$myAIC - 
  (2*mk.out.irrev3$k**2 + 2*mk.out.irrev3$k)/(mk.out.irrev3$n - mk.out.irrev3$k - 1)

# graphically compare these (with AICs)
sf = 2
png("fig-pruning.png", width=600*sf, height=400*sf, res=72*sf)
ggarrange(mk.inference.plot(mk.out.irrev, full.title=TRUE) + theme(plot.title = element_text(hjust = 0.5)),
          mk.inference.plot(mk.out.irrev2, full.title=TRUE) + theme(plot.title = element_text(hjust = 0.5)),
          mk.inference.plot(mk.out.irrev3, full.title=TRUE) + theme(plot.title = element_text(hjust = 0.5)), 
          nrow=1, labels=c("A", "B", "C"))
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

png("fig-1.png", width=1000*sf, height=700*sf, res=72*sf)
print(ggarrange(results.fig(parallelised.runs[[3]], 
                            omit.branch.lengths = obls[3], 
                            flux.threshold.pmax = pmaxs[3],
                            stats = "AICc"),
                results.fig(parallelised.runs[[4]], 
                            omit.branch.lengths = obls[4], 
                            flux.threshold.pmax = pmaxs[4], 
                            stats = "AICc",
                            labels=c("D", "E", "F")),
                nrow=2))
dev.off()

png("fig-2.png", width=1000*sf, height=700*sf, res=72*sf)
print(ggarrange(results.fig(parallelised.runs[[6]], 
                            omit.branch.lengths = obls[6], 
                            flux.threshold.pmax = pmaxs[6],
                            stats = "AICc"),
                results.fig(parallelised.runs[[5]], 
                            omit.branch.lengths = obls[5], 
                            flux.threshold.pmax = pmaxs[5], 
                            stats = "AICc",
                            labels=c("D", "E", "F")),
                nrow=2))
dev.off()

png("fig-3.png", width=1000*sf, height=350*sf, res=72*sf)
print(ggarrange(results.fig(parallelised.runs[[7]], 
                            omit.branch.lengths = obls[7], 
                            flux.threshold.pmax = pmaxs[7],
                            stats = "AICc")))
dev.off()

png("fig-4.png", width=1000*sf, height=350*sf, res=72*sf)
print(ggarrange(results.fig(parallelised.runs[[8]], 
                            omit.branch.lengths = obls[8], 
                            flux.threshold.pmax = pmaxs[8],
                            stats = "AICc")))
dev.off()

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
