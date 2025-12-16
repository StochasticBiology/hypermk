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

##### run the set of experiments
# produce and style output plots
# 1 "cross.sectional.single", 2 "cross.sectional.many", 3 "single", 4 "single.rev",
# 5 "single.uncertain", 6 "cross.sectional.cross",
# 7 "ovarian", 8 "TB"

# by default, we'll only run the quickest experiments (1-6); increase to 8 for the long haul
nexpts = 6
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

sf = 5
png("fig-1-tmp.png", width=1000*sf, height=700*sf, res=110*sf)
print(ggarrange(results.fig(parallelised.runs[[3]],
                            omit.branch.lengths = obls[3],
                            flux.threshold.pmax = pmaxs[3],
                            stats = "AIC"),
                results.fig(parallelised.runs[[4]],
                            omit.branch.lengths = obls[4],
                            flux.threshold.pmax = pmaxs[4],
                            stats = "AIC",
                            labels=c("D", "E", "F")),
                nrow=2))
dev.off()

png("fig-1.png", width=1000*sf, height=700*sf, res=72*sf)
print(ggarrange(results.fig(parallelised.runs[[3]],
                            omit.branch.lengths = obls[3],
                            flux.threshold.pmax = pmaxs[3],
                            stats = "AIC"),
                results.fig(parallelised.runs[[4]],
                            omit.branch.lengths = obls[4],
                            flux.threshold.pmax = pmaxs[4],
                            stats = "AIC",
                            labels=c("D", "E", "F")),
                nrow=2))
dev.off()

png("fig-2.png", width=1000*sf, height=700*sf, res=72*sf)
print(ggarrange(results.fig(parallelised.runs[[6]],
                            omit.branch.lengths = obls[6],
                            flux.threshold.pmax = pmaxs[6],
                            stats = "AIC"),
                results.fig(parallelised.runs[[5]],
                            omit.branch.lengths = obls[5],
                            flux.threshold.pmax = pmaxs[5],
                            stats = "AIC",
                            labels=c("D", "E", "F")),
                nrow=2))
dev.off()

png("fig-3.png", width=1000*sf, height=350*sf, res=72*sf)
print(ggarrange(results.fig(parallelised.runs[[7]],
                            omit.branch.lengths = obls[7],
                            flux.threshold.pmax = pmaxs[7],
                            stats = "AIC")))
dev.off()

png("fig-4.png", width=1000*sf, height=350*sf, res=72*sf)
print(ggarrange(results.fig(parallelised.runs[[8]],
                            omit.branch.lengths = obls[8],
                            flux.threshold.pmax = pmaxs[8],
                            stats = "AIC")))
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
