### Packages
require(phytools)
require(geiger)
require(corHMM)
setwd("..")
source("mk-specifics.R")
setwd("other_methods")

# issues so far: phytools doesn't seem to capture one of the observed states (110)
# how to force corHMM to be irreversible?

###############
####### Phylogenetic data: single reversible case study

# just call this function to create an object with the standard format for a dataset
# the actual data will be overwritten by the specific case below
data.set = setup.data("single.rev")

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

mk.out.irrev = mk.inference(data.set$tree, data.set$L, 
                          use.priors = FALSE, data.set$tips, 
                          reversible = FALSE,
                          optim_max_iterations = 2000,
                          Ntrials = 1,
                          Nthreads = 1)

g.castor.rev = plot.hypercube2(mk.out.rev$mk_fluxes, L)
g.castor.irrev = plot.hypercube2(mk.out.irrev$mk_fluxes, L)

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

# irreversible case
mmp <- mk_index_matrix(L, reversible=FALSE)
colnames(mmp) <- rownames(mmp) <-  seq_len(nrow(mmp))

p_single_irrev <- phytools::fitMk(tree = g_s_tree,
                                  x = x_mat,
                                  model = mmp,
                                  pi = c(1, rep(0, 31)))

# IGJ addition -- coerce into form that can be analysed and plotted with existing functions
p_single_irrev$transition_matrix = p_single_irrev$index.matrix
for(i in 1:length(p_single_irrev$rates)) {
  p_single_irrev$transition_matrix[p_single_irrev$transition_matrix==i] = p_single_irrev$rates[i]  
}
p_single_irrev.trans = mk_pull_transitions(p_single_irrev, reversible = FALSE)
p_single_irrev.fluxes = mk_simulate_fluxes(p_single_irrev, L, reversible=FALSE) 

# reversible case
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

# plots
g.phytools.irrev = plot.hypercube2(p_single_irrev.fluxes, L)
g.phytools.rev = plot.hypercube2(p_single_rev.fluxes, L)

###############
#### for corHMM core

data.set$tree$tip.label = 1:length(data.set$tree$tip.label)
data.set$obs = data.frame(species=data.set$tree$tip.label, data.set$x.binary)

# fit without explicitly considering all binary states
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

# irreversible case
# do the fitting
fitted.mk.c.rev = corHMM(data.set$tree, data.set$obs, rate.cat=1, 
                         rate.mat = mmpr, collapse = FALSE)

# pull the transition matrix
tmp = fitted.mk.c.rev$solution
tmp[is.na(tmp)] = 0
fitted.mk.c.rev$transition_matrix = tmp

fitted.mk.c.rev.trans = mk_pull_transitions(fitted.mk.c.rev, reversible = TRUE)
fitted.mk.c.rev.fluxes = mk_simulate_fluxes(fitted.mk.c.rev, L, reversible=TRUE) 

g.corHMM.rev = plot.hypercube2(fitted.mk.c.rev.fluxes, L)

# try forcing irreversibility by changing the allowed rate matrix

fitted.mk.c.irrev = corHMM(data.set$tree, data.set$obs, rate.cat=1, 
                           rate.mat = mmp, collapse = FALSE)

# pull the transition matrix
tmp = fitted.mk.c.irrev$solution
tmp[is.na(tmp)] = 0
fitted.mk.c.irrev$transition_matrix = tmp

fitted.mk.c.irrev.trans = mk_pull_transitions(fitted.mk.c.irrev, reversible = TRUE)
fitted.mk.c.irrev.fluxes = mk_simulate_fluxes(fitted.mk.c.irrev, L, reversible=TRUE) 

g.corHMM.irrev = plot.hypercube2(fitted.mk.c.irrev.fluxes, L)

fitted.mk.c.rev$loglik
p_single_rev$logLik
mk.out.rev$fitted_mk$loglikelihood

fitted.mk.c.irrev$loglik
p_single_irrev$logLik
mk.out.irrev$fitted_mk$loglikelihood

sf = 2
png("compare-set.png", width=800*sf, height=800*sf, res=72*sf)
ggarrange(g.castor.irrev, g.corHMM.irrev, g.phytools.irrev, 
          g.castor.rev, g.corHMM, g.phytools.rev, 
          g.corHMM.rev,
          labels = c("A. castor irrev", "B. corHMM irrev", "C. phytools irrev", 
                     "D. castor rev", "E. corHMM default", "F. phytools rev", 
                     "G. corHMM rev"),
          ncol=3,nrow=3)
dev.off()
