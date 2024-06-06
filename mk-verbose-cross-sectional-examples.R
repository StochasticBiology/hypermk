## FIXME: once changes in this branch are considered definitive,
##       use the refit-prune strategy in the examples below.

## Verbose example for cross-sectional data modified from
## mk-timing.R, mk-cgh.R, mk-collection.R

## If you have data you want to analyse, go directly to the "Analyse"
## sections.

source("mk-shared.R")


### Simulate data

L <- 3 ## Number of features
sample.size <- 16 ## Number of cases or sample size

#### Simulate as matrix of subjects by features

## Data given as subject by feature, with 1 in altered feature. E.g.:
##
## |           | Gene_1 | Gene_2 | Gene_3 |
## |-----------+--------+--------+--------|
## | Subject_1 |      0 |      0 |      1 |
## | Subject_2 |      1 |      0 |      1 |
## | Subject_3 |      1 |      0 |      0 |
## | Subject_4 |      1 |      1 |      0 |
## | ...       |        |        |        |
##

## For reproducibility
set.seed(1)

subject_by_feature <- matrix(sample(c(0, 1), L * sample.size, replace = TRUE),
                             nrow = sample.size)

#### Simulate specific patterns

## Two-competing pathways: see "cross.sectional.cross" in mk-collection.R

### Analyse

##  Analysis involves these steps
##  1. Putting observed data in mk-suitable format
##  2. Creating the matrix of possible transitions for the rate_model
##     argument in castor's fit_mk.
##  3. Running castor's fit_mk
##  4. Simulating fluxes from the transition matrix
##  5. Optionally, refiting the model after pruning fluxes = 0

##  One can carry out steps 1 to 4 individually, but
##  function mk.inference provides a general wrapper. For cross-sectional
##  observations, function mk_infer_cross_sectional calls mk.inference
##  after automatically preparing the data to use mk.inference.
##  We focus on this route.  At the end, we provide code for the
##  step-by-step approach


#### Analyse: unpruned
rev_example <- mk_infer_cross_sectional(subject_by_feature,
                                        reversible = TRUE)

irrev_example <- mk_infer_cross_sectional(subject_by_feature,
                                          reversible = FALSE)

#### Analyse: prune

## Refit the model setting to 0 rates with flux of 0. In many of these small
## examples, pruning could make no difference since often no fluxes are
## 0.  With set.seed(1) some of the fluxes are 0.

to.nullify.rev <- rev_example$mk_fluxes[which(rev_example$mk_fluxes$Flux==0),1:2]+1
to.nullify.irrev <- irrev_example$mk_fluxes[which(irrev_example$mk_fluxes$Flux==0),1:2]+1

rev_example_pruned <- mk_infer_cross_sectional(subject_by_feature,
                                               reversible = TRUE,
                                               to.nullify = to.nullify.rev)

irrev_example_pruned <- mk_infer_cross_sectional(subject_by_feature,
                                                 reversible = FALSE,
                                                 to.nullify = to.nullify.irrev)


### Plot

## Obtain AICs and create titles
AIC.rev = rev_example$fitted_mk$AIC
AIC.rev.reduced = rev_example_pruned$fitted_mk$AIC
title.rev = paste0("reversible fit, simplified AIC ~ ", round(AIC.rev.reduced, digits=2),
                   " (full ", round(AIC.rev, digits=2), ")", collapse = "")

AIC.irrev = irrev_example$fitted_mk$AIC
AIC.irrev.reduced = irrev_example_pruned$fitted_mk$AIC
title.irrev = paste0("irreversible fit, simplified AIC ~ ", round(AIC.irrev.reduced, digits=2),
                     " (full ", round(AIC.irrev, digits=2), ")", collapse = "")


## Plot
plot.hypercube2(trans.f = rev_example$mk_fluxes, bigL = L, rates = FALSE) +
  ggtitle(title.rev)

plot.hypercube2(trans.f = irrev_example$mk_fluxes, bigL = L, rates = FALSE) +
  ggtitle(title.irrev)

## Plot thresholding fluxes, for example only those > 25%

plot.hypercube2(trans.f = rev_example$mk_fluxes[rev_example$mk_fluxes$Flux > 0.25 * max(rev_example$mk_fluxes$Flux), ],
                bigL = L, rates = FALSE) +
  ggtitle(title.rev)

## With the simulated data with set.seed(1) this is identical to the unthresholded figure
plot.hypercube2(trans.f = irrev_example$mk_fluxes[irrev_example$mk_fluxes$Flux > 0.25 * max(irrev_example$mk_fluxes$Flux), ],
                bigL = L, rates = FALSE) +
  ggtitle(title.irrev)



#### Appendix: step-by-step
#####  1. Put data in mk-suitable format
## Analysis requires states as 0-indexed decimal format.

## If your data are in subject_by_feature matrix format turn to 0-indexed
## as:
## states <- apply(subject_by_feature, 1, BinToDec)

## Put 0-indexed data in a forma suitable for mk-based analysis
mk.data <- mk_cross_sectional(subject_by_feature, L)

## How the input data look like:
## two lists, tree and tips, each of length sample.size
##   - tree: 2-tip tree
##   - tips: 2xL array;
##           - 1st row: 0 in all positions, except 1 for present feature;
##                      this is the deterministic prior for the observed
##                      state of the feature;
##           - 2nd row: 1/(2**L); the uniform prior over 2**L states;
(mk.data)

##### 2. Create matrix of possible transitions

## Transitions are possible only between states that differ in one feature
## (under irreversibility, only gains are possible).
## Create (2**L) by (2**L) matrix where a non-0 entry uniquely labels each
## possible transition.

###### Reversible case
index_matrix_rev <- mk_index_matrix(L, reversible = TRUE)

###### Irreversible case
index_matrix_irrev <- mk_index_matrix(L, reversible = FALSE)


##### 3. Run castor's fit_mk

fitted_mk_rev <- castor::fit_mk(trees = mk.data$tree,
                                Nstates = 2**L,
                                tip_priors = mk.data$tips,
                                rate_model = index_matrix_rev,
                                ## root prior puts all prob. in the 0 state
                                root_prior = c(1, rep(0, (2**L) - 1)))

fitted_mk_irrev <- castor::fit_mk(trees = mk.data$tree,
                                  Nstates = 2**L,
                                  tip_priors = mk.data$tips,
                                  rate_model = index_matrix_irrev,
                                  ## root prior puts all prob. in the 0 state
                                  root_prior = c(1, rep(0, (2**L) - 1)))

##### 4. Extract transitions

## From the fitted mk model, extract the estimated transitions
trans_rev <- mk_pull_transitions(fitted_mk_rev, reversible = TRUE)
trans_irrev <- mk_pull_transitions(fitted_mk_irrev, reversible = FALSE)


##### 5. Simulate fluxes by using random walkers on the transition matrix

fluxes_rev <- mk_simulate_fluxes(fitted_mk_rev, L = L, reversible = TRUE)
fluxes_irrev <- mk_simulate_fluxes(fitted_mk_irrev, L = L, reversible = FALSE)


##### 6. Plot

## reversible
plot.hypercube2(trans.f = fluxes_rev, bigL = L, rates = FALSE)
## You can threshold the fluxes, for example only those > 25%
plot.hypercube2(trans.f = fluxes_rev[fluxes_rev$Flux > 0.25 * max(fluxes_rev$Flux), ],
                bigL = L, rates = FALSE)

## irreversible
plot.hypercube2(trans.f = fluxes_irrev, bigL = L, rates = FALSE)

## You can threshold the fluxes, for example only those > 25%
plot.hypercube2(trans.f = fluxes_irrev[fluxes_irrev$Flux > 0.25 * max(fluxes_irrev$Flux), ],
                bigL = L, rates = FALSE)
