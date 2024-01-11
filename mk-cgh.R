source("mk-shared.R")

L = 4
ov.table = read.table("Data/ovarian.csv", header=FALSE)
ov.reduced = ov.df[,1:L]
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
