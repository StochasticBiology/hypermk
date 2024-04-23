### What is this

## Attempts to analyze the cross-sectional data using phytools.

## geiger is not included, since it cannot be used for the phylogenetic
## examples, and its limitations would appear here too: it seems to
## only be able to model transitions involving the states actually
## observed at the tips. See file phylog-examples.R



### Limit cores/CPUs/OMP/...

## phytools, in a machine with more than 20 cores, will readily try to use
## them, but I don't see how to limit them from phytools itself.  Try
## everything I can think of (one or more things are probably unnecessary,
## but this works)

## Run this BEFORE loading phytools. Set max_threads_cores to whatever is sensible
max_threads_cores <- 30
Sys.getenv(c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS"))
Sys.setenv(OMP_NUM_THREADS = max_threads_cores)
Sys.setenv(OMP_THREAD_LIMIT = max_threads_cores)
Sys.setenv(OPENBLAS_NUM_THREADS = max_threads_cores)
library(RhpcBLASctl)
RhpcBLASctl::omp_set_num_threads(max_threads_cores)
RhpcBLASctl::blas_set_num_threads(max_threads_cores)


### Packages
library(phytools)
source("../mk-shared.R")

#### utilities for phytools

## Return the decimal data as a matrix for phytools
## o.w., phytools can't fit out models
## (see examples of need in datasets-phylog-cross-sect.R)
## Like one-hot-encoding for the decimal state.
dec_vector_to_matrix <- function(x, max_state) {
  if (max(x) > max_state) stop("largest x > max_state")
  m <- matrix(0, nrow = length(x), ncol = max_state)
  m[cbind(seq_along(x), x)] <- 1
  colnames(m) <- as.character(seq_len(max_state))
  rownames(m) <- as.character(seq_along(x))
  return(m)
}


## Given an object from phytools' fitMk, return
## the non-zero entries of the transition matrix, ordered and numbered
## as in mk_df's pull-transitions
ph_pull_transitions <- function(x) {
  m_index <- which(x$index.matrix != 0, arr.ind = TRUE)
  rate_index <- x$index.matrix[m_index]
  m_index_rate <- cbind(From = m_index[, "row"],
                        To =   m_index[, "col"],
                        Rate_index = rate_index)

  m_index_rate <- m_index_rate[order(rate_index), ]
  m_index_rate <- cbind(m_index_rate, Rate = x$rates)
  m_index_rate <- m_index_rate[m_index_rate[, "Rate"] != 0, ]
  ## Order as mk_df
  m_index_rate <- m_index_rate[order(m_index_rate[, "From"],
                                     m_index_rate[, "To"]), ]
  ## Renumber as mk_df
  m_index_rate[, "From"] <- m_index_rate[, "From"] - 1
  m_index_rate[, "To"] <- m_index_rate[, "To"] - 1
  m_index_rate <- data.frame(m_index_rate)
  row.names(m_index_rate) <- seq_len(nrow(m_index_rate))
  return(m_index_rate)
}

## Transition matrix from fitMk. Copied from phytools:::printMk
ph_matrix <- function(x) {
  Q <- matrix(NA, length(x$states), length(x$states))
  Q[] <- c(0, x$rates)[x$index.matrix + 1]
  diag(Q) <- 0
  diag(Q) <- -rowSums(Q)
  colnames(Q) <- rownames(Q) <- x$states
  return(Q)
}

## Simulate fluxes. Copied from mk_simulate_fluxes
## except for the pull_transitions code and transition matrix.
## I add remove_zero for easier comparisons
ph_simulate_fluxes <- function(x, L, reversible = TRUE, nwalker = 10000,
                               remove_zero = TRUE) {
  # to get flux matrix we'll simulate random walkers on the transition matrix
  if(reversible == TRUE) {
    Q <- ph_matrix(x)
    # set up data frame containing transitions and fluxes
    mk.rev.df = ph_pull_transitions(x)
    mk.rev.df = mk.rev.df[mk.rev.df$From != mk.rev.df$To,]
    mk.rev.df$Rate[mk.rev.df$Rate == Inf] = 10*max(mk.rev.df$Rate[mk.rev.df$Rate != Inf])
    mk.rev.df$Flux = 0
    # simulate walkers starting from 0^L
    for(walk in 1:nwalker) {
      state = 0
      for(t in 1:(2*L)) {
        # pull row of transition matrix
        trans = Q[state+1,]
        trans[state+1] = 0
        if(sum(trans)==0) {break}
        # sample the next transition
        this.out = sample(0:(2**L-1), size=1, prob=trans)

        ref = which(mk.rev.df$From == state & mk.rev.df$To == this.out)
        mk.rev.df$Flux[ref] = mk.rev.df$Flux[ref]+1
        state = this.out
      }
    }
    if (remove_zero) mk.rev.df <- mk.rev.df[mk.rev.df$Flux > 0,  ]
    return(mk.rev.df)
  }
  else {

    mk.irrev.df = ph_pull_transitions(x)
    mk.irrev.df = mk.irrev.df[mk.irrev.df$From != mk.irrev.df$To,]
    mk.irrev.df$Rate[mk.irrev.df$Rate == Inf] = 10*max(mk.irrev.df$Rate[mk.irrev.df$Rate != Inf])
    mk.irrev.df$Flux = 0
    # simulate walkers starting from 0^L
    for(walk in 1:nwalker) {
      state = 0
      for(t in 1:L) {
        # get possible out transitions
        outs = which(mk.irrev.df$From == state)
        if(length(outs) == 0) {
          break
        } else if(length(outs) == 1) {
          this.out = outs[1]
        } else {
          # sample a possible out transition
          this.out = sample(outs, size=1, prob=mk.irrev.df$Rate[outs])
        }
        mk.irrev.df$Flux[this.out] = mk.irrev.df$Flux[this.out]+1
        state = mk.irrev.df$To[this.out]
      }
    }
    if (remove_zero) mk.irrev.df <- mk.irrev.df[mk.irrev.df$Flux > 0,  ]
    return(mk.irrev.df)
  }
}



### Data sets: original + one with more replicates

if (FALSE) {
  ## The next pulls the data from the saved data. Not used
  ## now, for clarity, but this would allow users to run the other examples
  ## just loading the RData.
  ## The RData is created as per datasets-phylog-cross-sect.R
  load("setup.data_output.RData")
  ## csc are just the cross.sectional.cross data with 6 rows and two paths
  data_decimal <- d_cross.sectional.cross$x.decimal
}

## Original data. From mk-collection.R
L <- 4
m = matrix(c(0,0,0,1,
             0,0,1,1,
             0,1,1,1,
             1,0,0,0,
             1,1,0,0,
             1,1,1,0), ncol=L, byrow=TRUE)

data_decimal <- apply(m, 1, BinToDec) + 1

#### Data format for phytools
## Matrix form (one-hot encoding) for phytools
data_matrix <- dec_vector_to_matrix(data_decimal, max_state = 2^L)


## Increased sample size. For example, 14.
set.seed(1)
n_species <- 14
indices <- sample(seq_along(data_decimal), n_species, replace = TRUE)
## Enlarged data set: increase the number of rows to n_species
data_decimal_2 <- data_decimal[indices]
data_matrix_2 <- dec_vector_to_matrix(data_decimal_2, max_state = 2^L)


#### Same data, for castor-based analysis

data_castor <- mk_cross_sectional(m, L = L)

## Enlarged data set: increase the number of rows to n_species
data_castor_2 <- mk_cross_sectional(m[indices, ], L = L)


### Index matrices for transitions rates

index_matrix_irrev <- mk_index_matrix(L = L, reversible = FALSE)
index_matrix_rev <- mk_index_matrix(L = L, reversible = TRUE)

rownames(index_matrix_irrev) <- colnames(index_matrix_irrev) <-
  seq_len(nrow(index_matrix_irrev))

rownames(index_matrix_rev) <- colnames(index_matrix_rev) <-
  seq_len(nrow(index_matrix_rev))


### Phylogenies: for phytools, we need branch lenghts

## Using, in the code below, a phylogeny created as
## star <- ape::stree(6, "start", 1:6)
## or
## star <- phytools::starTree(seq_len(nrow(csc_mat)))
## breaks fitMk. Branch lengths seem to be needed.

## FIXME: branch lengths would allow us to add time of sample
## if we had it. Explore this eventually.

## FIXME: interesting that castor does not need branch lengths. What is
## it assuming or using instead?

star <- phytools::starTree(seq_len(nrow(data_matrix)),
                           rep(1, nrow(data_matrix)))

star2 <- phytools::starTree(seq_len(nrow(data_matrix_2)),
                            rep(1, nrow(data_matrix_2)))

### Fits, phytools

## FIXME: for real, might want to use pruning=TRUE (faster, as per help).

p_irr <- phytools::fitMk(tree = star,
                         x = data_matrix,
                         model = index_matrix_irrev,
                         pi = c(1, rep(0, 15)),
                         opt.method = "optim")

## To check some simple rate assumptions
p_irr_twice <-  phytools::fitMk(tree = phytools::starTree(1:6, rep(2, 6)),
                                x = data_matrix,
                                model = index_matrix_irrev,
                                pi = c(1, rep(0, 15)),
                                opt.method = "optim")

p_rev <- phytools::fitMk(tree = star,
                         x = data_matrix,
                         model = index_matrix_rev,
                         pi = c(1, rep(0, 15)),
                         opt.method = "optim")

p_2_irr <- phytools::fitMk(tree = star2,
                           x = data_matrix_2,
                           model = index_matrix_irrev,
                           pi = c(1, rep(0, 15)),
                           opt.method = "optim")

p_2_rev <- phytools::fitMk(tree = star2,
                           x = data_matrix_2,
                           model = index_matrix_rev,
                           pi = c(1, rep(0, 15)),
                           opt.method = "optim")

## Rerun all above with the (default) nlminb
p_irr_nlminb <- phytools::fitMk(tree = star,
                                x = data_matrix,
                                model = index_matrix_irrev,
                                pi = c(1, rep(0, 15)),
                                opt.method = "nlminb")

p_irr_twice_nlminb <-  phytools::fitMk(tree = phytools::starTree(1:6, rep(2, 6)),
                                       x = data_matrix,
                                       model = index_matrix_irrev,
                                       pi = c(1, rep(0, 15)),
                                       opt.method = "nlminb")


p_rev_nlminb <- phytools::fitMk(tree = star,
                                x = data_matrix,
                                model = index_matrix_rev,
                                pi = c(1, rep(0, 15)),
                                opt.method = "nlminb")

p_2_irr_nlminb <- phytools::fitMk(tree = star2,
                                  x = data_matrix_2,
                                  model = index_matrix_irrev,
                                  pi = c(1, rep(0, 15)),
                                  opt.method = "nlminb")

p_2_rev_nlminb <- phytools::fitMk(tree = star2,
                                  x = data_matrix_2,
                                  model = index_matrix_rev,
                                  pi = c(1, rep(0, 15)),
                                  opt.method = "nlminb")

save(file = "cross-sectional.phytools.RData",
     p_irr, p_irr_twice,
     p_rev, p_2_irr, p_2_rev,
     p_irr_nlminb, p_irr_twice_nlminb,
     p_rev_nlminb, p_2_irr_nlminb, p_2_rev_nlminb,
     indices)

## We can also fit the model with MCMC
## BUT: there is a bug in the mcmcMk code, and if
## geiger is available this fails. To run this, move geiger's
## library directory out of the way.
## I have reported the bug: https://github.com/liamrevell/phytools/issues/154
p_irr_mcmc <- phytools::mcmcMk(tree = star,
                               x = data_matrix,
                               model = index_matrix_irrev,
                               pi = c(1, rep(0, 15)),
                               ## Might not be large enough?
                               ngen = 50000,
                               ## plotting while running slows it; see help too
                               plot = FALSE)

### Fits, castor

castor_irr <- mk_infer_cross_sectional(m = m, reversible = FALSE)

castor_rev <- mk_infer_cross_sectional(m = m, reversible = TRUE)

castor_2_irr <- mk_infer_cross_sectional(m = m[indices, ], reversible = FALSE)

castor_2_rev <- mk_infer_cross_sectional(m = m[indices, ], reversible = TRUE)

castor_irr_nlminb <- mk_infer_cross_sectional(m = m, reversible = FALSE,
                                              optim_algorithm = "nlminb")

castor_rev_nlminb <- mk_infer_cross_sectional(m = m, reversible = TRUE,
                                              optim_algorithm = "nlminb")

castor_2_irr_nlminb <- mk_infer_cross_sectional(m = m[indices, ],
                                                reversible = FALSE,
                                                optim_algorithm = "nlminb")

castor_2_rev_nlminb <- mk_infer_cross_sectional(m = m[indices, ],
                                                reversible = TRUE,
                                                optim_algorithm = "nlminb")

save(file = "cross-sectional-castor.RData",
     castor_irr, castor_rev,
     castor_2_irr, castor_2_rev,
     castor_irr_nlminb, castor_rev_nlminb,
     castor_2_irr_nlminb, castor_2_rev_nlminb)

## Timings: castor seems faster: phytools used 30 cores here
## (multiple cores decrease phytools time, but not linearly, though)
## and times (system.time, elapsed seconds) were for phytools and castor
## 6, 14, 9, 28  [phytools]
## 4, 4,  8, 21  [castor]
## FIXME: if this is relevant, time seriously with different number
## of cores.

#### Compare castor and phytools


##### Consistency between fits on phylogenies with length 1 or 2
## Hummm... a few I do not understand
ph_pull_transitions(p_irr)
ph_pull_transitions(p_irr_twice)

ph_pull_transitions(p_irr_nlminb)
ph_pull_transitions(p_irr_twice_nlminb)


##### Transitions and fluxes for phytools and castor

## There is a scale issue. On top of that, is phytools sometimes doing
## more aggressive pruning and driving more rates to 0 with nlminb?

ph_pull_transitions(p_irr_nlminb)
castor_irr_nlminb$mk_df
ph_pull_transitions(p_irr)
castor_irr$mk_df

ph_pull_transitions(p_rev)
castor_rev$mk_df

ph_pull_transitions(p_2_irr)
castor_2_irr$mk_df

ph_pull_transitions(p_2_rev)
castor_2_rev$mk_df


##
ph_simulate_fluxes(p_irr, L = L)
castor_irr$mk_fluxes[castor_irr$mk_fluxes$Flux > 0, ]


##
ph_simulate_fluxes(p_irr_nlminb, L = L)
castor_irr_nlminb$mk_fluxes[castor_irr_nlminb$mk_fluxes$Flux > 0, ]


## Eh??!!
ph_simulate_fluxes(p_rev, L = L)
castor_rev$mk_fluxes[castor_rev$mk_fluxes$Flux > 0, ]


ph_simulate_fluxes(p_rev_nlminb, L = L)
castor_rev_nlminb$mk_fluxes[castor_rev_nlminb$mk_fluxes$Flux > 0, ]


## Very similar
ph_simulate_fluxes(p_2_irr, L = L)
castor_2_irr$mk_fluxes[castor_2_irr$mk_fluxes$Flux > 0, ]

ph_simulate_fluxes(p_2_irr_nlminb, L = L)
castor_2_irr_nlminb$mk_fluxes[castor_2_irr_nlminb$mk_fluxes$Flux > 0, ]


## Hummm...
ph_simulate_fluxes(p_2_rev, L = L)
castor_2_rev$mk_fluxes[castor_2_rev$mk_fluxes$Flux > 0, ]

ph_simulate_fluxes(p_2_rev_nlminb, L = L)
castor_2_rev_nlminb$mk_fluxes[castor_2_rev_nlminb$mk_fluxes$Flux > 0, ]




#### Can we run castor on a star phylogeny for cross-sectional data?

## This fails: Error in if (first_guess_rate == 0)
## It does not depend on the optimizer (same error with optim and nlminb)
castor_star_irr <- castor::fit_mk(trees = ape::stree(6, "star", 1:6),
                                  Nstates = 16,
                                  tip_states = data_decimal,
                                  rate_model = index_matrix_irrev,
                                  root_prior = c(1, rep(0, 15)))

## And this fails too, same error
castor_star_irr <- castor::fit_mk(trees = star,
                                  Nstates = 16,
                                  tip_states = data_decimal,
                                  rate_model = index_matrix_irrev,
                                  root_prior = c(1, rep(0, 15)))
