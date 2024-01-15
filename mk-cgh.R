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
mk.rev.g = graph_from_data_frame(ov.rev.df)
mk.rev.g$barcodes = unlist(lapply(as.integer(V(mk.rev.g)$name), DecToBin, L))
mk.stats = ov.res.df[ov.res.df$Experiment==expt & ov.res.df$Fit=="reversible",]
g.rev = ggraph(mk.rev.g) + 
  geom_edge_arc(strength=0.1,aes(alpha=sqrt(Rate)), 
                arrow=arrow(length=unit(0.2, "inches"), type="closed")) + 
  geom_node_label(aes(label=name), size=2, alpha=0.4) + 
  ggtitle(paste(c(expt, ", rev fit, AIC ", round(mk.stats$AIC, digits=2), 
                  " reducable to ", round(mk.stats$AIC.reduced, digits=2)), 
                collapse=""))+ theme_void() +
  theme(legend.position="none", plot.title = element_text(size = 10)) 

# plot graph, with thresholding by flux
mk.rev.g = graph_from_data_frame(ov.rev.df[ov.rev.df$Flux > 0,])
g.rev.flux = ggraph(mk.rev.g) + 
  geom_edge_arc(strength=0.1,aes(alpha=sqrt(Flux)), 
                arrow=arrow(length=unit(0.2, "inches"), type="closed")) + 
  geom_node_label(aes(label=name), size=2, alpha=0.4) + 
  scale_edge_width(limits=c(0,NA))+ 
  ggtitle(paste(c(expt, ", rev fit, AIC ", round(mk.stats$AIC, digits=2), 
                  " reducable to ", round(mk.stats$AIC.reduced, digits=2)), 
                collapse="")) + theme_void() +
  theme(legend.position="none", plot.title = element_text(size = 10))


# plot graph without pruning by flux
mk.irrev.g = graph_from_data_frame(ov.irrev.df)
mk.irrev.g$barcodes = unlist(lapply(as.integer(V(mk.irrev.g)$name), DecToBin, L))
mk.stats = ov.res.df[ov.res.df$Experiment==expt & ov.res.df$Fit=="irreversible",]
g.irrev = ggraph(mk.irrev.g) + 
  geom_edge_arc(strength=0.1,aes(alpha=sqrt(Rate)), 
                arrow=arrow(length=unit(0.2, "inches"), type="closed")) + 
  geom_node_label(aes(label=name), size=2, alpha=0.4) +
  ggtitle(paste(c(expt, ", irrev fit, AIC ", round(mk.stats$AIC, digits=2), 
                  " or simplified ", round(mk.stats$AIC.reduced, digits=2)), 
                collapse=""))+ theme_void() +
  theme(legend.position="none", plot.title = element_text(size = 10))

# plot graph with pruning by flux 
mk.irrev.g = graph_from_data_frame(ov.irrev.df[ov.irrev.df$Flux > 0,])
g.irrev.flux = ggraph(mk.irrev.g) + 
  geom_edge_arc(strength=0.1,aes(alpha=sqrt(Flux)), 
                arrow=arrow(length=unit(0.2, "inches"), type="closed")) + 
  geom_node_label(aes(label=name), alpha=0.4, size=2) + 
  scale_edge_width(limits=c(0,NA)) + theme_void() +
  ggtitle(paste(c(expt, ", irrev fit, AIC ", round(mk.stats$AIC, digits=2), 
                  " or simplified ", round(mk.stats$AIC.reduced, digits=2)), 
                collapse="")) +
  theme(legend.position="none", plot.title = element_text(size = 10)) 

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

