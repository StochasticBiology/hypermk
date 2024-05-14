### What is this

## Attempts to analyze the phylogenetic data using phytools and geiger.

#### Summary

##  - geiger cannot be used for our models as it seems it only models
##    transitions involving the states actually observed at the tips,
##    excluding other states

##  - phytools could be used at least for phylogenetic data without
##    uncertainty but: a) runs with the single data set showed lack of
##    convergence and phytools fitMk does not give the option of increasing
##    iterations; b) it does not seem possible to start the optimizer from
##    multiple places (in contrast to castor); c) I see no way to deal with
##    uncertain observations.



### Load the data
## The RData is created as per datasets-phylog-cross-sect.R
load("setup.data_output.RData")

### Limit cores/CPUs/OMP/...

## phytools, in a machine with more than 20 cores, will readily try to use
## them, but I don't see how to limit them from phytools itself.  Try
## everything I can think of (one or more things are probably unnecessary,
## but this works)

## Run this BEFORE loading phytools. Set max_threads_cores to whatever is sensible
max_threads_cores <- 4
Sys.getenv(c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS"))
Sys.setenv(OMP_NUM_THREADS = max_threads_cores)
Sys.setenv(OMP_THREAD_LIMIT = max_threads_cores)
Sys.setenv(OPENBLAS_NUM_THREADS = max_threads_cores)
library(RhpcBLASctl)
RhpcBLASctl::omp_set_num_threads(max_threads_cores)
RhpcBLASctl::blas_set_num_threads(max_threads_cores)



### Packages
library(phytools)
library(geiger)
source("../mk-shared.R")

### Phylogenetic data: single

#### Common transformations

single_index_matrix_irrev <- mk_index_matrix(L = d_single$L, reversible = FALSE)

single_index_matrix_rev <- mk_index_matrix(L = d_single$L, reversible = TRUE)

## Must add unique names to data and tree, at least for geiger
## Check first
stopifnot(identical(d_single$tree$tip.label, d_single$x.decimal))

g_s_tree <- d_single$tree
## Setting tip.label with as.character(seq_along) seems not necessary
## or even a bad idea?
g_s_tree$tip.label <- seq_along(d_single$x.decimal)
g_s_dat <- d_single$x.decimal
names(g_s_dat) <- as.character(seq_along(d_single$x.decimal))




#### phytools


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


##### Single, irreversible
## Unclear from the help if matrix(i,j) is i->j or j->i
## It seems i -> j, like for castor

## Do not set diagonal to NA, or it fails
## as there is a k <- max(rate)
mmp <- single_index_matrix_irrev
colnames(mmp) <- rownames(mmp) <-  seq_len(nrow(mmp))


## If x is a vector, this cannot work for reasons similar
## to geiger: In fitMk there is this check
## if (ncol(model) != ncol(x))
##    stop("model does not have the right number of columns")
## where ncol(x) is obtained from
## x <- to.matrix(x, sort(unique(x)))
## So use matrix shape

## p_single_irrev <- phytools::fitMk(tree = g_s_tree,
##                                   x = g_s_dat,
##                                   model = mmp)


x_mat <- dec_vector_to_matrix(g_s_dat, 32)

p_single_irrev <- phytools::fitMk(tree = g_s_tree,
                                  x = x_mat,
                                  model = mmp,
                                  pi = c(1, rep(0, 31)))

p_single_irrev

# IGJ addition -- coerce into form that can be analysed and plotted with existing functions
p_single_irrev$transition_matrix = p_single_irrev$index.matrix
for(i in 1:length(p_single_irrev$rates)) {
  p_single_irrev$transition_matrix[p_single_irrev$transition_matrix==i] = p_single_irrev$rates[i]  
}
p_single_irrev.trans = mk_pull_transitions(p_single_irrev, reversible = FALSE)
p_single_irrev.fluxes = mk_simulate_fluxes(p_single_irrev, 5, reversible=FALSE) 
plot.hypercube2(p_single_irrev.fluxes, 5)

save(file = "phytools_single_irrev.RData", p_single_irrev)


## The plot is too busy and takes long to load
## plot(p_single_irrev)

## Can also run mcmcMk
## BUT: there is a bug in the mcmcMk code, and if
## geiger is available this fails.
## I have reported the bug: https://github.com/liamrevell/phytools/issues/154
phytools::mcmcMk(tree = g_s_tree,
                 x = x_mat,
                 model = mmp,
                 pi = c(1, rep(0, 31)),
                 ngen = 100, ## yeah, ridiculous; for speed
                 plot = FALSE)

p_single_irrev_prune <- phytools::fitMk(tree = g_s_tree,
                                        x = x_mat,
                                        model = mmp,
                                        pi = c(1, rep(0, 31)),
                                        pruning = TRUE)
save(file = "phytools_single_irrev_prune.RData", p_single_irrev_prune)

p_single_irrev_prune

p_single_irrev_prune_optim <- phytools::fitMk(tree = g_s_tree,
                                              x = x_mat,
                                              model = mmp,
                                              pi = c(1, rep(0, 31)),
                                              pruning = TRUE,
                                              opt.method = "optim")
save(file = "phytools_single_irrev_prune_optim.RData", p_single_irrev_prune_optim)

p_single_irrev_prune_optim


##### Single, reversible

mmpr <- single_index_matrix_rev
colnames(mmpr) <- rownames(mmpr) <-  seq_len(nrow(mmpr))

p_single_rev <- phytools::fitMk(tree = g_s_tree,
                                x = x_mat,
                                model = mmpr,
                                pi = c(1, rep(0, 31)))

save(file = "phytools_single_rev.RData", p_single_rev)

p_single_rev

p_single_rev_prune <- phytools::fitMk(tree = g_s_tree,
                                      x = x_mat,
                                     model = mmpr,
                                     pi = c(1, rep(0, 31)),
                                     pruning = TRUE)

save(file = "phytools_single_rev_prune.RData", p_single_rev_prune)

p_single_rev_prune


p_single_rev_prune_optim <- phytools::fitMk(tree = g_s_tree,
                                            x = x_mat,
                                      model = mmpr,
                                      pi = c(1, rep(0, 31)),
                                      pruning = TRUE,
                                      opt.method = "optim")

save(file = "phytools_single_rev_prune_optim.RData", p_single_rev_prune_optim)


p_single_rev_prune_optim


##### phytools: most don't converge and there is no way to increase iterations

## The help only talks about increasing iterations for function fitHRM.
## The code of fitMk seems to have no way of passing niter/iter, nor ways
## of changing optimizers' settings.

p_single_irrev
p_single_irrev_prune
p_single_irrev_prune_optim  ## this seems to have converged

p_single_rev
p_single_rev_prune
p_single_rev_prune_optim

## We might want to compare with the output from castor, but since
## most runs do not seem to converge, this seems a moot point.




#### geiger

##### Single, irreversible

mm <- t(single_index_matrix_irrev)
colnames(mm) <- rownames(mm) <-  seq_len(nrow(mm))

## NA in diagonal? that is what is shown in the help
## but not in pp. 157 and ff. of Revell & Harmon. Shouldn't hurt.
## But it seems to hurt with fitMk
diag(mm) <- NA

## Geiger's fitDiscrete does not seem capable of directly fitting
## models where to set of all states is larger that
## the number of observed states

## See this failure
g_single_irrev <- geiger::fitDiscrete(phy = g_s_tree,
                                      dat = g_s_dat,
                                      model = mm)
## the error
## Error in .check.states(tree, states, strict.vals = 1:k)
## is due to k being 4, the number of observed states
## k <- nlevels(as.factor(dat))  in geiger:::mkn.lik

## Even this fails
mm2 <- mm
mm2 <- mm2[c(25, 29, 31, 32), c(25, 29, 31, 32)]

geiger::fitDiscrete(phy = g_s_tree,
                    dat = g_s_dat,
                     model = mm2)
## trying to pass a strict.vals = c("25", "29", "31", "32"))
## is hopeless, because the call to .check.states is
##  .check.states(tree, states, strict.vals = 1:k)

## But we can hack it this way
g2d <- g_s_dat
g2d[g2d == 25] <- 1
g2d[g2d == 29] <- 2
g2d[g2d == 31] <- 3
g2d[g2d == 32] <- 4

tmp_1 <- geiger::fitDiscrete(phy = g_s_tree,
                             dat = g2d,
                             model = mm2)

## So exchange ids of states in the big data
g3d <- g_s_dat

## exchange values
g3d[g3d == 25] <- 1001
g3d[g3d == 29] <- 1002
g3d[g3d == 31] <- 1003
g3d[g3d == 32] <- 1004

g3d[g3d == 1] <- 25
g3d[g3d == 2] <- 29
g3d[g3d == 3] <- 31
g3d[g3d == 4] <- 32

g3d[g3d == 1001] <- 1
g3d[g3d == 1002] <- 2
g3d[g3d == 1003] <- 3
g3d[g3d == 1004] <- 4

table(g3d, g_s_dat)

## exchange matrix
new_o <- c(25, 29, 31, 32, 5:24, 1, 26:28, 2, 30, 3, 4)
mm_tmp <- mm
## exchange rows and then columns
mm_tmp <- mm_tmp[new_o, ]
mm_tmp <- mm_tmp[, new_o]

## rename rows and columns
colnames(mm_tmp) <- rownames(mm_tmp) <-  seq_len(nrow(mm_tmp))

## Some checks
all.equal(mm_tmp[c(5:24, 26:28, 30), c(5:24, 26:28, 30)],
          mm[c(5:24, 26:28, 30), c(5:24, 26:28, 30)])

original <- cbind(which(mm > 0, arr.ind = TRUE),
                  mm[which(mm > 0, arr.ind = TRUE)])

transformed <- cbind(which(mm_tmp > 0, arr.ind = TRUE),
                     mm_tmp[which(mm_tmp > 0, arr.ind = TRUE)])

## I'd need to double check the indices of the parameters are the same.
## But first, verify we can even run it:
## Yes, it runs, but it ignores every other state. So this seems hopeless.
g_single_irrev <- geiger::fitDiscrete(phy = g_s_tree,
                                      dat = g3d,
                                      model = mm_tmp)
