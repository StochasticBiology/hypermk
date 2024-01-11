# shared library loads and helper function for Mk model project

library(ape)
library(phangorn)
library(castor)
library(igraph)
library(ggplot2)
library(ggraph)
library(stringr)
library(ggpubr)
library(ggtree)

# binary to decimal function
BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

# decimal to binary function
DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

# decimal to binary function, returning a numerical vector
DecToBinV <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(s)
}

# this function converts a species name string from the Newick format which Common Taxonomy Tree gives us into a simpler lower-case, no quotes version comparable to Kostas' dataset
convname = function(str) {
  return(tolower(gsub("\'", "", gsub("_", " ", str))))
}

# construct matrix assigning indices to elements describing possible transitions
mk_index_matrix = function(L, reversible=TRUE) {
  index_matrix = matrix(rep(0, (2**L)**2), nrow=2**L)
  index = 1
  # go through binary strings and put indices on transitions to Hamming +1 neighbours
  for(i in 0:(2**L-1)) {
    b = as.numeric(strsplit(DecToBin(i, L),"")[[1]])
    for(j in 1:L) {
      if(b[j] == 0) {
        c = b
        c[j] = 1
        partner = BinToDec(c)
        index_matrix[i+1,partner+1] = index
        index = index+1
        if(reversible == TRUE) {
          # add reverse transition if required
          index_matrix[partner+1,i+1] = index
          index = index+1
        }
      }
    }
  }
  return(index_matrix)
}

# extract transitions as data frame from fitted Mk model
mk_pull_transitions = function(fit.mk, reversible=TRUE) {
  m = fit.mk$transition_matrix
  mk_df = data.frame()
  if(reversible == TRUE) {
    for(i in 1:(nrow(m))) {
      for(j in 1:(ncol(m))) {
        if(m[i,j] != 0) {
          mk_df = rbind(mk_df, data.frame(From=i-1, To=j-1, Rate=m[i,j]/sum(m[i,seq(1,ncol(m))[seq(1,ncol(m))!=i]])))
        }
      }
    }
  } else {
    # loop over transitions
    for(i in 1:(nrow(m)-1)) {
      for(j in (i+1):ncol(m)) {
        if(m[i,j] != 0) {
          mk_df = rbind(mk_df, data.frame(From=i-1, To=j-1, Rate=m[i,j]/sum(m[i,(i+1):ncol(m)])))
        }
      }
    }
  }
  return(mk_df)
}

mk_simulate_fluxes_reversible = function(fitted_mk.rev) {
# to get flux matrix we'll simulate random walkers on the transition matrix
nwalker = 10000
threshold = nwalker*(2*L)/10000

# set up data frame containing transitions and fluxes
mk.rev.df = mk_pull_transitions(fitted_mk.rev, reversible=TRUE)
mk.rev.df = mk.rev.df[mk.rev.df$From != mk.rev.df$To,]
mk.rev.df$Rate[mk.rev.df$Rate == Inf] = 10*max(mk.rev.df$Rate[mk.rev.df$Rate != Inf])
mk.rev.df$Flux = 0
# simulate walkers starting from 0^L
for(walk in 1:nwalker) {
  state = 0
  for(t in 1:(2*L)) {
    # pull row of transition matrix
    trans = fitted_mk.rev$transition_matrix[state+1,]
    trans[state+1] = 0
    if(sum(trans)==0) {break}
    # sample the next transition
    this.out = sample(0:(2**L-1), size=1, prob=trans)
    
    ref = which(mk.rev.df$From == state & mk.rev.df$To == this.out)
    mk.rev.df$Flux[ref] = mk.rev.df$Flux[ref]+1
    state = this.out     
  }
}
return(mk.rev.df)
}

mk_simulate_fluxes_irreversible = function(fitted_mk.irrev) {
  nwalker = 10000
  threshold = nwalker*(2*L)/10000
  
  mk.irrev.df = mk_pull_transitions(fitted_mk.irrev, reversible=FALSE)
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
  return(mk.irrev.df)
}

# adapted HyperTraPS plot for HyperHMM outputs
plot.hypercube = function(trans.p, bigL, node.labels = TRUE, use.probability = FALSE) {
  ### produce hypercube subgraph
  trans.g = graph_from_data_frame(trans.p[!is.infinite(trans.p$Probability),])
  bs = unlist(lapply(as.numeric(V(trans.g)$name), DecToBin, len=bigL))
  V(trans.g)$binname = V(trans.g)$name
  layers = str_count(bs, "1")
  this.plot =  ggraph(trans.g, layout="sugiyama", layers=layers) + 
    geom_edge_arc(aes(edge_width=Probability, edge_alpha=Probability), strength=0.1) +
    scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,0.5)) 
  if(node.labels == TRUE) {
    this.plot = this.plot + geom_node_label(aes(label=binname),angle=45,size=2) 
  }
  return(this.plot)
}
