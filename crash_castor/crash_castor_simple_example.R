## This is one minimal reproducible example of castor's
## crash.
## Based on crash_single.uncertain.example_2

## How to create the data and observe the crash
## 1. In crash_common.R uncomment the browser() call
## 2. When it stops save the value of arguments
## arg_trees <<- mk.tree
## arg_Nstates <<- 2**L
## arg_rate_model <<- index_matrix
## arg_root_prior <<- c(1,rep(0, 2**L-1))
## arg_optim_max_iterations <<- optim_max_iterations
## arg_Ntrials <<- Ntrials
## arg_Nthreads <<- Nthreads
## 3. Q (from the browser) and save the data:
##    save(file = "args_castor_crash.RData", list = ls(pattern = glob2rx("arg_*")))
## 4. Start a new R and run the code below
##

library(castor)
load("args_castor_crash.RData")

set.seed(50001)
fit_mk(
  trees = arg_trees,
  Nstates = 32, ## arg_Nstates,
  tip_priors = arg_tip_priors,
  optim_algorithm = "optim", ## arg_optim_algorithm,
  rate_model = arg_rate_model,
  root_prior = c(1, rep(0, 31)), ## arg_root_prior,
  optim_max_iterations = 200, ## arg_optim_max_iterations,
  Ntrials = 5, ## arg_Ntrials,
  Nthreads = 1, ## arg_Nthreads,
  verbose = TRUE, diagnostics = TRUE
)

## Will eventually fail with
##
## Trial 2, tree 1 (64 tips): Model evaluation failed:
## Error in stats::optim(start_dense_rates/rate_scale, function(x) objective_function(x *  :
## L-BFGS-B needs finite values of 'fn'
