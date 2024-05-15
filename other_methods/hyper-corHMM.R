# example using corHMM to fit the hypercubic Mk model

library(corHMM)
setwd("..")
source("mk-collection.R")
setwd("other_methods")

# very basic example to demo the syntax
n.tip = 16
my.tree = rtree(n.tip)
tip.df = data.frame(species=my.tree$tip.label, 
                    c1=rbinom(n.tip,1,0.5), 
                    c2=rbinom(n.tip,1,0.5),
                    c3=rbinom(n.tip,1,0.5))
corHMM(my.tree, tip.df, rate.cat=1)

# example using one of our synthetic cases
data.set = setup.data("single.rev")
data.set$tree$tip.label = 1:length(data.set$tree$tip.label)
data.set$obs = data.frame(species=data.set$tree$tip.label, data.set$x.binary)
fitted.mk.c = corHMM(data.set$tree, data.set$obs, rate.cat=1)

# awkwardly pull the states corresponding to each entry in our fitted transition matrix
L = 5
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
fitted.mk.c.fluxes = mk_simulate_fluxes(fitted.mk.c, 5, reversible=TRUE) 

plot.hypercube2(fitted.mk.c.fluxes, 5)
  
# what happens if our observations span more than one state?
n.tip = 16
my.tree = rtree(n.tip)
obs = matrix(0, ncol=5, nrow=n.tip)
obs[2,] = c(1, 1, 1, 1, 1)
for(i in 3:n.tip) {
  r = runif(1)
  if(r < 0.5) {
    obs[i,] = c(0, 0, 1, 1, 1)
  } else {
    obs[i,] = c(1, 1, 1, 0, 0)
  } 
} 
tip.df = data.frame(species=my.tree$tip.label, 
                    obs)
corHMM(my.tree, tip.df, rate.cat=1)
# ^ this throws an error

# larger, random state sets
# set.seed(1) generates a dataset that throws an error
# set.seed(3) generates a dataset that runs OK
set.seed(3)
tip.df = data.frame(species=my.tree$tip.label, 
                    c1=rbinom(n.tip,1,0.5), 
                    c2=rbinom(n.tip,1,0.5),
                    c3=rbinom(n.tip,1,0.5),
                    c4=rbinom(n.tip,1,0.5),
                    c5=rbinom(n.tip,1,0.5))

corHMM(my.tree, tip.df, rate.cat=1)
# ^ this throws an error for some structures

#####

# can't immediately handle a star (non-binary) phylogeny...
n.tip = 16
my.tree = stree(n.tip)
tip.df = data.frame(species=my.tree$tip.label, 
                    c1=rbinom(n.tip,1,0.5), 
                    c2=rbinom(n.tip,1,0.5),
                    c3=rbinom(n.tip,1,0.5))
corHMM(my.tree, tip.df, rate.cat=1)
# ... even after dichotomising? perhaps rootedness is an issue
my.tree = multi2di(my.tree)
corHMM(my.tree, tip.df, rate.cat=1)

# example from the vignettes
data(primates)
phy <- primates[[1]]
phy <- multi2di(phy)
data <- primates[[2]]
corHMM(phy = phy, data = data, rate.cat = 1)
