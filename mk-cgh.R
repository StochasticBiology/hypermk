source("mk-shared.R")

data.plot = list()

expt = "ovarian"

L = 5
ov.table = read.csv("Data/ovarian.csv", header=FALSE)
# 0-indexed decimal state representation
ov.states = apply(ov.table[, 1:L], 1, BinToDec)

# reconstruct barcodes from 0-indexed decimals
barcodes = unlist(lapply(ov.states, DecToBin, L))

# construct tables of observed barcodes and their decimals in the dataset  
b.stats = as.data.frame(table(barcodes))
data.plot[[expt]] = ggplot(b.stats, aes(x=barcodes, y=Freq)) + geom_col() +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("Observations") + ylab("Count") +
  scale_y_continuous(breaks = seq(0, max(b.stats$Freq), by = 1)) 

# pass 0-indexed state refs
ov.data = mk_cross_sectional(ov.states, L)
index_matrix_irrev = mk_index_matrix(L, reversible=FALSE)
fitted_mk.irrev = castor::fit_mk(ov.data$tree, 2**L, 
                                 tip_priors=ov.data$tips, 
                                 rate_model=index_matrix_irrev, 
                                 root_prior=c(1,rep(0, 2**L-1)))

index_matrix_rev = mk_index_matrix(L, reversible=TRUE)
fitted_mk.rev = castor::fit_mk(ov.data$tree, 2**L, 
                               tip_priors=ov.data$tips, 
                               rate_model=index_matrix_rev, 
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


# plot graph, without thresholding by flux
mk.stats = ov.res.df[ov.res.df$Experiment==expt & ov.res.df$Fit=="reversible",]
t.str = titlestr(expt, "rev fit", mk.stats)
g.rev = plot.hypercube2(ov.rev.df, L, rates=TRUE) +
  ggtitle(t.str) 

flux.threshold.pmax = 0.1
flux.threshold = flux.threshold.pmax*max(ov.rev.df$Flux)

# plot graph, with thresholding by flux
g.rev.flux = plot.hypercube2(ov.rev.df[ov.rev.df$Flux > flux.threshold,], L) +
  ggtitle(t.str) 

flux.threshold.pmax = 0.1
flux.threshold = flux.threshold.pmax*max(ov.irrev.df$Flux)

# plot graph without thresholding by flux
mk.stats = ov.res.df[ov.res.df$Experiment==expt & ov.res.df$Fit=="irreversible",]
t.str = titlestr(expt, "irrev fit", mk.stats)
g.irrev = plot.hypercube2(ov.irrev.df, L, rates=TRUE) +
  ggtitle(t.str)

# plot graph with thresholding by flux 
g.irrev.flux = plot.hypercube2(ov.irrev.df[ov.irrev.df$Flux > flux.threshold,], L) +
  ggtitle(t.str) 

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

# output to file. Figure 6 of current ms.
fname = "fig-4x.png"
sf = 2
png(fname, width=800*sf, height=600*sf, res=72*sf)
print(ggarrange(ggarrange(data.plot[[expt]], labels=c("A")),
                ggarrange(g.irrev.flux, g.rev.flux, ncol=2, 
                          labels=c("B","C"), label.y=c(0.9,0.9)), 
                nrow=2))
dev.off()

