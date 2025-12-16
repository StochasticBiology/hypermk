require(hypermk)
require(ape)
require(ggpubr)

# set seed for reproducibility
set.seed(1)

# create simple data matrix
m = matrix(c(0,0,1,
             0,1,1,
             1,1,1), ncol=3, nrow=3)
# create random phylogeny
tree = ape::rphylo(3, 1, 1)

# fit irreversible model to cross-sectional observations
fit.cs.irrev = mk_infer_cross_sectional(m, reversible=FALSE)
# fit reversible model to cross-sectional observations
fit.cs.rev = mk_infer_cross_sectional(m, reversible=TRUE)
# fit irreversible model to phylogenetically-linked observations
fit.phy.irrev = mk_infer_phylogenetic(m, tree, reversible=FALSE)
# fit reversible model to phylogenetically-linked observations
fit.phy.rev = mk_infer_phylogenetic(m, tree, reversible=TRUE)

# prune extraneous parameters from these fitted models
fit.cs.irrev.prune = mk_prune_model(fit.cs.irrev)
fit.cs.rev.prune = mk_prune_model(fit.cs.rev)
fit.phy.irrev.prune = mk_prune_model(fit.phy.irrev)
fit.phy.rev.prune = mk_prune_model(fit.phy.rev)

# plot comparisons and statistics
ggpubr::ggarrange(
  mk.inference.plot(fit.cs.irrev, full.title = TRUE),
  mk.inference.plot(fit.cs.rev, full.title = TRUE),
  mk.inference.plot(fit.phy.irrev, full.title = TRUE),
  mk.inference.plot(fit.phy.rev, full.title = TRUE),
  mk.inference.plot(fit.cs.irrev.prune, full.title = TRUE),
  mk.inference.plot(fit.cs.rev.prune, full.title = TRUE),
  mk.inference.plot(fit.phy.irrev.prune, full.title = TRUE),
  mk.inference.plot(fit.phy.rev.prune, full.title = TRUE),
  nrow = 2, ncol = 4
)


