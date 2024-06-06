## Source code and redefine mk.inference to show debug output from castor:
## set verbose and diagnostics = TRUE in call to
## castor::fit_mk


source("../mk-shared.R")
require(parallel)



## redefine, to show debug output
mk.inference = function(mk.tree, L, use.priors, tips, reversible,
                        optim_max_iterations = 200, Ntrials = 1,
                        optim_algorithm = c("optim", "nlminb"),
                        Nthreads = 1,
                        to.nullify = matrix(nrow=0, ncol=2)) {
  optim_algorithm = match.arg(optim_algorithm)
  # for cross-sectional data and uncertain data, tips = tip.priors, use.priors = TRUE
  # otherwise, tips = tip.states, use.priors = FALSE

  # construct matrix describing possible transitions
  index_matrix = mk_index_matrix(L, reversible=reversible, to.nullify=to.nullify)

  # do the Mk model fitting
  # remember the (deterministic) prior on the root state! this is important
  ## Uncomment to generate data for reproducible crashes
  ## browser()
  if(use.priors == TRUE) {
    # specify priors, rather than precise states, on the tips of the tree
    fitted_mk = castor::fit_mk(mk.tree, 2**L,
                               tip_priors=tips,
                               optim_algorithm = optim_algorithm,
                               rate_model=index_matrix,
                               root_prior=c(1,rep(0, 2**L-1)),
                               optim_max_iterations = optim_max_iterations,
                               Ntrials = Ntrials,
                               Nthreads = Nthreads,
                               verbose = TRUE, diagnostics = TRUE)
  } else {
    # specify precise states
    fitted_mk = castor::fit_mk(mk.tree, 2**L,
                               tip_states=tips,
                               optim_algorithm = optim_algorithm,
                               rate_model=index_matrix,
                               root_prior=c(1,rep(0, 2**L-1)),
                               optim_max_iterations = optim_max_iterations,
                               Ntrials = Ntrials,
                               Nthreads = Nthreads,
                               verbose = TRUE, diagnostics = TRUE)
  }

  if(fitted_mk$converged != TRUE) {
    message("WARNING: Mk model fit didn't converge!")
  }

  # convert inferred rate matrix into transition set
  mk_df = mk_pull_transitions(fitted_mk, reversible = reversible)
  # and simulate fluxes through this transition set
  mk_fluxes = mk_simulate_fluxes(fitted_mk, L, reversible = reversible)

  # return a list of useful info
  l.return = list(fitted_mk = fitted_mk,
                  reversible = reversible,
                  mk_df = mk_df,
                  mk_fluxes = mk_fluxes)

  return(l.return)
}
