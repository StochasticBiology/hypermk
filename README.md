# hypermk

HyperMk for reversible accumulation modelling
===

HyperMk in an R library form. Install with

`remotes::install_github("StochasticBiology/hypermk", ref = "r-library")`

You can find the demo files with

`system.file(".", package = "hypermk")`

Now with associated article in Bioinformatics: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btae737/7922554 . The original codebase was not in an R package format; you can find it here https://github.com/StochasticBiology/hypermk/tree/publication-code .

HyperMk uses a hypercubic transition matrix and the Mk (Markov k-state) model [1,2] from phylogenetics to model accumulation processes, including reversibility. The code here uses the `castor` library in R [3] to provide flexible specification and fitting of the Mk model. 

![image](https://github.com/StochasticBiology/hypermk/assets/50171196/6d03abd9-015b-4629-8a48-287a616cf65f)


Contents
---
General, in `R/`:

* `mk-shared.R` contains a set of common functions for labelling and plotting, as well as for creating the matrices for specifying parameter forms for HyperMk

Specific to research article, in `inst/`, producing the figures for the above article:

* `mk-specifics.R` contains a set of helper functions for creating and analysing a set of synthetic and real-world examples, including analysis of a cross-sectional ovarian cancer dataset and a phylogenetic TB dataset.
* `mk-collection.R` runs inference on these synthetic and real-world examples
* `mk-timing.R` scans through properties of synthetic datasets to demonstrate the time scaling of this implementation.

`inst/extdata` contains data on drug resistance evolution in tuberculosis from [4], and from the ovarian cancer CGH study [5], which are used as case studies.

Specifically:

`mk-shared.R`
----

The first few functions are higher-level and should be sufficient for application to a given dataset (see `mk-demo.R`). The others are lower-level and mainly used internally.

| Function | Description |
|----------|------------ |
| `mk_infer_cross_sectional` | simple wrapper of `mk_cross_sectional` + `mk.inference`, taking a matrix and outputting the fitted model |
| `mk_infer_phylogenetic` | simple wrapper of `mk_phylogeny_precise` + `mk.inference`, taking a matrix and tree and outputting the fitted model |
| `mk_prune_model` | prunes a fitted model by removing low-flux edges and refitting. Takes a fitted model and a threshold flux
| `mk.inference.plot` | plots a fitted model. takes a fitted model, an optional threshold flux, and an optional Boolean describing whether to use a verbose title
|----------|------------ |
| `mk.inference` | fits the hypercubic Mk model of choice to a prepared dataset, simulates fluxes on the transition graph, and returns details of the fit including the parameterised hypercube |
| `mk_index_matrix` | produces the index matrix for a given experimental design, assigning integer labels to permitted transitions between states (reversible/irreversible, specified edges removed) |
| `mk_pull_transitions` | get the transitions from a fitted Mk model|
| `mk_simulate_fluxes` | simulate walkers through this transition set to get fluxes on each edge|
| `mk_phylogeny_precise` | given a set of states and a tree linking them, prepare the data for use in `mk.inference` |
| `mk_cross_sectional` | given a set of cross-sectional states, prepare the data for use in `mk.inference` |
| `mk_infer_cross_sectional` | simple wrapper of `mk_cross_sectional` + `mk.inference`, taking a matrix and outputting the fitted model |

`mk-specifics.R`
----

| Function | Description |
|----------|------------ |
| `setup.data` | prepares data for a collection of synthetic and experimental case studies |
| `parallel.fn` | wraps an experiment for a given dataset, calling `mk.inference` for reversible and irreversible, pruned and raw fits, and producing a combined object for plotting |
| `results.fig` | plots these combined objects |


References
---

Chapters 7 and 8 of https://bio.libretexts.org/Bookshelves/Evolutionary_Developmental_Biology/Phylogenetic_Comparative_Methods_(Harmon) are excellent reading on the Mk model.

[1] Pagel, M., 1994. Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters. Proceedings of the Royal Society of London. Series B: Biological Sciences, 255(1342), pp.37-45.

[2] Lewis, P.O., 2001. A likelihood approach to estimating phylogeny from discrete morphological character data. Systematic biology, 50(6), pp.913-925.

[3] Louca, S. and Doebeli, M., 2018. Efficient comparative phylogenetics on large trees. Bioinformatics, 34(6), pp.1053-1055.

[4] Casali, N., Nikolayevskyy, V., Balabanova, Y., Harris, S.R., Ignatyeva, O., Kontsevaya, I., Corander, J., Bryant, J., Parkhill, J., Nejentsev, S. and Horstmann, R.D., 2014. Evolution and transmission of drug-resistant tuberculosis in a Russian population. Nature genetics, 46(3), pp.279-286.

[5] Knutsen, T., Gobu, V., Knaus, R., Padilla‐Nash, H., Augustus, M., Strausberg, R.L., Kirsch, I.R., Sirotkin, K. and Ried, T., 2005. The interactive online SKY/M‐FISH & CGH database and the Entrez cancer chromosomes search database: linkage of chromosomal aberrations with the genome sequence. Genes, Chromosomes and Cancer, 44(1), pp.52-64.
