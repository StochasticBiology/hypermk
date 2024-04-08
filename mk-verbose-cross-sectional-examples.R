## Verbose example for cross-sectional data modified from
## mk-timing.R, mk-cgh.R, mk-collection.R

## If you have data you want to analyse, go directly to the "Analyse"
## sections.

source("mk-shared.R")


### Simulate data

L <- 3 ## Number of features
n.states <- 16 ## Number of cases or sample size

#### Simulate as decimal states

## States of the n.states samples, given as decimal, 0-indexed
## FIXME: clearer if we sample with replacement from the possible states?
## states <- round(runif(n.states, min=0, max=2**L-1))
possible_states <- seq(from = 0, to = (2**L) - 1)
states <- sample(possible_states, size = n.states, replace = TRUE)


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

subject_by_feature <- matrix(sample(c(0, 1), L * n.states, replace = TRUE),
                             nrow = n.states)

#### Simulate specific patterns

## Two-competing pathways: see "cross.sectional.cross" in mk-collection.R

### Analyse

##  Analysis involves three steps
##  1. Putting observed data in mk-suitable format
##  2. Creating the matrix of possible transitions for the rate_model
##     argument in castor's fit_mk.
##  3. Running castor's fit_mk

##  We illustrate the reversible and irreversible cases below. The first
##  step is common to both


#### 1. Put data in mk-suitable format
## Analysis requires states as 0-indexed decimal format.

## If your data are in subject_by_feature matrix format turn to 0-indexed
## as:
## states <- apply(subject_by_feature, 1, BinToDec)

## Put 0-indexed data in a forma suitable for mk-based analysis
mk.data <- mk_cross_sectional(states, L)

## How the input data look like:
## two lists, tree and tips, each of length n.states
##   - tree: 2-tip tree
##   - tips: 2xL array;
##           - 1st row: 0 in all positions, except 1 for present feature;
##                      this is the deterministic prior for the observed
##                      state of the feature;
##           - 2nd row: 1/(2**L); the uniform prior over 2**L states;
(mk.data)

#### 2. Create matrix of possible transitions

## Transitions are possible only between states that differ in one feature
## (under irreversibility, only gains are possible).
## Create (2**L) by (2**L) matrix where a non-0 entry uniquely labels each
## possible transition.

##### Reversible case
index_matrix_rev <- mk_index_matrix(L, reversible = TRUE)

##### Irreversible case
index_matrix_irrev <- mk_index_matrix(L, reversible = FALSE)


#### 3. Run castor's fit_mk

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

#### 4. Extract transitions

## From the fitted mk model, extract the estimated transitions
trans_rev <- mk_pull_transitions(fitted_mk_rev, reversible = TRUE)
trans_irrev <- mk_pull_transitions(fitted_mk_irrev, reversible = FALSE)

## FIXME: is there a -1 in the rate of self-transitions in reversible
## model and not irreversible? I understand why they are -1
## but why do something different for reversible and irreversible?
## See note in mk-shared, function mk_pull_transitions.

#### 5. Simulate fluxes by using random walkers on the transition matrix

fluxes_rev <- mk_simulate_fluxes_reversible(fitted_mk_rev)
fluxes_irrev <- mk_simulate_fluxes_irreversible(fitted_mk_irrev)

#### 6. Plot

## FIXME: I do not understand the logic of "reduced" or something previous:
##        fluxes that are 0 are not added to penalty, but the rate was
##        estimated. However, these fluxes == 0 seem to correspond
##        to transitions from states that are not reachable from 0
## Questions:
##  - why are there estimates of rates for impossible transitions?
##  - would we want to remove those rates from the mk_pull_transitions
##    output? (i.e., remove anything involving a non-reachable state
##    --- I think this is easy in irreversible case, 
##    and doable in the reversible---); though maybe looking at fluxes
##    is good enough? Would removing them make the simulate_fluxes code
##    faster? Would this speed increase even matter?


stats_rev <- data.frame(AIC = fitted_mk_rev$AIC,
                        AIC.reduced = fitted_mk_rev$AIC-2*length(which(fluxes_rev$Flux==0)))

title_rev <- titlestr(expt = "tiny_example", fit = "rev fit",
                      stats.df = stats_rev)

plot.hypercube2(trans.f = fluxes_rev, bigL = L, rates = FALSE) +
  ggtitle(title_rev)

## You can threshold the fluxes, for example only those > 25%
plot.hypercube2(trans.f = fluxes_rev[fluxes_rev$Flux > 0.25 * max(fluxes_rev$Flux), ],
                bigL = L, rates = FALSE) +  ggtitle(title_rev)


## Same for irreversible

stats_irrev <- data.frame(AIC = fitted_mk_irrev$AIC,
                          AIC.reduced = fitted_mk_irrev$AIC-2*length(which(fluxes_irrev$Flux==0)))

title_irrev <- titlestr(expt = "tiny_example", fit = "irrev fit",
                        stats.df = stats_irrev)

plot.hypercube2(trans.f = fluxes_irrev, bigL = L, rates = FALSE) +
  ggtitle(title_irrev)

## You can threshold the fluxes, for example only those > 25%
plot.hypercube2(trans.f = fluxes_irrev[fluxes_irrev$Flux > 0.25 * max(fluxes_irrev$Flux), ],
                bigL = L, rates = FALSE) +  ggtitle(title_irrev)
