source("mk-shared.R")

L = 5
ov.table = read.csv("Data/ovarian.csv", header=FALSE)
ov.reduced = ov.table[,1:L]
for(i in 1:L) {
  val = 2**(L-i)
  ov.reduced[,i] = ov.reduced[,i]*val
}
ov.states = rowSums(ov.reduced)

ov.data = mk_cross_sectional(ov.states, L)
index_matrix = mk_index_matrix(L, reversible=FALSE)
fitted_mk.irrev = fit_mk(ov.data$tree, 2**L, 
                         tip_priors=ov.data$tips, 
                         rate_model=index_matrix, 
                         root_prior=c(1,rep(0, 2**L-1)))

index_matrix = mk_index_matrix(L, reversible=TRUE)
fitted_mk.rev = fit_mk(ov.data$tree, 2**L, 
                       tip_priors=ov.data$tips, 
                       rate_model=index_matrix, 
                       root_prior=c(1,rep(0, 2**L-1)))

ov.irrev = mk_pull_transitions(fitted_mk.irrev, reversible = FALSE)
ov.rev = mk_pull_transitions(fitted_mk.rev, reversible = TRUE)

# set up data frame containing transitions and fluxes
ov.irrev.df = mk_simulate_fluxes_irreversible(fitted_mk.irrev)
ov.rev.df = mk_simulate_fluxes_reversible(fitted_mk.rev)

ov.res.df = data.frame()
ov.res.df = rbind(ov.res.df, data.frame(Experiment = "ovarian", 
                                        Fit = "irreversible", 
                                        AIC = fitted_mk.irrev$AIC, 
                                        AIC.reduced = fitted_mk.irrev$AIC-2*length(which(ov.irrev.df$Flux==0))))
ov.res.df = rbind(ov.res.df, data.frame(Experiment = "ovarian", 
                                        Fit = "reversible", 
                                        AIC = fitted_mk.rev$AIC, 
                                        AIC.reduced = fitted_mk.rev$AIC-2*length(which(ov.rev.df$Flux==0))))


expt = "ovarian"

# plot graph, without thresholding by flux
mk.stats = ov.res.df[ov.res.df$Experiment==expt & ov.res.df$Fit=="reversible",]
g.rev = plot.hypercube2(ov.rev.df, L, rates=TRUE) +
  ggtitle(paste(c(expt, ", rev fit, AIC ", round(mk.stats$AIC, digits=2), 
                  " or simplified ", round(mk.stats$AIC.reduced, digits=2)), 
                collapse="")) 

flux.threshold.pmax = 0.1
flux.threshold = flux.threshold.pmax*max(ov.rev.df$Flux)

# plot graph, with thresholding by flux
g.rev.flux = plot.hypercube2(ov.rev.df[ov.rev.df$Flux > flux.threshold,], L) +
  ggtitle(paste(c(expt, ", rev fit, AIC ", round(mk.stats$AIC, digits=2), 
                  " or simplified ", round(mk.stats$AIC.reduced, digits=2)), 
                collapse="")) 

flux.threshold.pmax = 0.1
flux.threshold = flux.threshold.pmax*max(ov.irrev.df$Flux)

# plot graph without pruning by flux
mk.stats = ov.res.df[ov.res.df$Experiment==expt & ov.res.df$Fit=="irreversible",]
g.irrev = plot.hypercube2(ov.irrev.df, L, rates=TRUE) +
  ggtitle(paste(c(expt, ", irrev fit, AIC ", round(mk.stats$AIC, digits=2), 
                  " or simplified ", round(mk.stats$AIC.reduced, digits=2)), 
                collapse=""))

# plot graph with pruning by flux 
g.irrev.flux = plot.hypercube2(ov.irrev.df[ov.irrev.df$Flux > flux.threshold,], L) +
  ggtitle(paste(c(expt, ", irrev fit, AIC ", round(mk.stats$AIC, digits=2), 
                  " or simplified ", round(mk.stats$AIC.reduced, digits=2)), 
                collapse="")) 

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
print(ggarrange(g.rev.flux, g.irrev.flux))
dev.off()


