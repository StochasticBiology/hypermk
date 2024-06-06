### Packages
require(phytools)
require(geiger)
require(corHMM)
setwd("..")
source("mk-specifics.R")
setwd("other_methods")

###############
####### Phylogenetic data: single reversible case study

data.set = setup.data("single.rev")
L = 5

# simple constructed case study
set.seed(1)
L = 3
n.tip = 16
my.tree = rtree(n.tip)

single.path = TRUE
if(single.path == TRUE) {
  obs = matrix(c(0,0,0,
                 0,0,1,
                 0,1,1,
                 1,1,1,
                 1,1,0), byrow = TRUE, ncol=3)
  obs = rbind(obs, obs, obs, 0)
  tip.df = data.frame(species=my.tree$tip.label, 
                      obs)
} else {
  tip.df = data.frame(species=my.tree$tip.label, 
                      c1=rbinom(n.tip,1,0.5), 
                      c2=rbinom(n.tip,1,0.5),
                      c3=rbinom(n.tip,1,0.5))
}
tips.dec = tip.df[,2]*4 + tip.df[,3]*2 + tip.df[,4] + 1
data.set = list()
data.set$obs = tip.df
data.set$x.decimal = data.set$tips = tips.dec
data.set$x.binary = as.matrix(tip.df[,2:(L+1)])
data.set$tree = my.tree
data.set$tree$tip.label = tips.dec
data.set$tree$node.labels
data.set$L = 3

###############
#### default (castor) core
mk.out.rev = mk.inference(data.set$tree, data.set$L, 
                          use.priors = FALSE, data.set$tips, 
                          reversible = TRUE,
                          optim_max_iterations = 2000,
                          Ntrials = 1,
                          Nthreads = 1)

g.castor = plot.hypercube2(mk.out.rev$mk_fluxes, L)

###############
#### for phytools core

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

g_s_dat = data.set$x.decimal
g_s_tree = data.set$tree
x_mat <- dec_vector_to_matrix(g_s_dat, 2**L)

mmpr <- mk_index_matrix(L, reversible=TRUE)
colnames(mmpr) <- rownames(mmpr) <-  seq_len(nrow(mmpr))

p_single_rev <- phytools::fitMk(tree = g_s_tree,
                                x = x_mat,
                                model = mmpr,
                                pi = c(1, rep(0, (2**L)-1)))

# IGJ addition -- coerce into form that can be analysed and plotted with existing functions
p_single_rev$transition_matrix = p_single_rev$index.matrix
for(i in 1:length(p_single_rev$rates)) {
  p_single_rev$transition_matrix[p_single_rev$transition_matrix==i] = p_single_rev$rates[i]  
}
p_single_rev.trans = mk_pull_transitions(p_single_rev, reversible = TRUE)
p_single_rev.fluxes = mk_simulate_fluxes(p_single_rev, L, reversible=TRUE) 

g.phytools = plot.hypercube2(p_single_rev.fluxes, L)

###############
#### for corHMM core
data.set$tree$tip.label = 1:length(data.set$tree$tip.label)
data.set$obs = data.frame(species=data.set$tree$tip.label, data.set$x.binary)

# do the fitting
fitted.mk.c = corHMM(data.set$tree, data.set$obs, rate.cat=1)

# awkwardly pull the states corresponding to each entry in our fitted transition matrix
fitted.mk.c$states = rep(0, max(fitted.mk.c$data.legend$d))
fitted.mk.c$dec.states = rep(0, max(fitted.mk.c$data.legend$d))
for(i in 1:max(fitted.mk.c$data.legend$d)) {
  this.ref = which(fitted.mk.c$data.legend$d == i)[1]
  fitted.mk.c$states[i] = paste0(fitted.mk.c$data[this.ref,2:(L+1)], collapse="")
  fitted.mk.c$dec.states[i] = BinToDec(fitted.mk.c$states[i])+1
}

# construct the full 2**L x 2**L transition matrix using the nonzero elements of the fitted solution
fitted.mk.c$transition_matrix = matrix(0, nrow=2**L, ncol=2**L)
nonzeroes = which(!is.na(fitted.mk.c$solution), arr.ind = TRUE)
for(i in 1:nrow(nonzeroes)) {
  ref.a = fitted.mk.c$dec.states[nonzeroes[i,1]]
  ref.d = fitted.mk.c$dec.states[nonzeroes[i,2]]
  val = fitted.mk.c$solution[nonzeroes[i,1], nonzeroes[i,2]]
  fitted.mk.c$transition_matrix[ref.a,ref.d] = val
}

fitted.mk.c.trans = mk_pull_transitions(fitted.mk.c, reversible = TRUE)
fitted.mk.c.fluxes = mk_simulate_fluxes(fitted.mk.c, L, reversible=TRUE) 

g.corHMM = plot.hypercube2(fitted.mk.c.fluxes, L)

fitted.mk.c$loglik
p_single_rev$logLik
mk.out.rev$fitted_mk$loglikelihood

ggarrange(g.castor, g.phytools, g.corHMM)
