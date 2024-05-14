# hypermk

HyperMk for reversible accumulation modelling
===

HyperMk uses a hypercubic transition matrix and the Mk (Markov k-state) model [1,2] from phylogenetics to model accumulation processes, including reversibility. The code here uses the `castor` library in R [3] to provide flexible specification and fitting of the Mk model. 

![image](https://github.com/StochasticBiology/hypermk/assets/50171196/2e190633-07ad-4550-9660-117d7b44b714)

Contents
---

* `mk-shared.R` contains a set of common functions for labelling and plotting, as well as for creating the matrices for specifying parameter forms for HyperMk
* `mk-collection.R` creates a set of synthetic examples and demonstrates inference using HyperMk. Includes analysis of a cross-sectional ovarian cancer dataset and a phylogenetic TB dataset.
* `mk-timing.R` scans through properties of synthetic datasets to demonstrate the time scaling of this implementation.
* `mk-verbose-cross-sectional-examples.R` contains small, heavily commented
  cross-sectional examples.

`Data` contains data on drug resistance evolution in tuberculosis from [4], and from the ovarian cancer CGH study [5], which are used as case studies.

Specifically:

`mk-shared.R`
----

| Function | Description |
|----------|------------ |
| `mk.inference` | fits the hypercubic Mk model of choice to a prepared dataset, simulates fluxes on the transition graph, and returns details of the fit including the parameterised hypercube |
| `mk_index_matrix` | produces the index matrix for a given experimental design, assigning integer labels to permitted transitions between states (reversible/irreversible, specified edges removed) |
| `mk_pull_transitions` | get the transitions from a fitted Mk model|
| `mk_simulate_fluxes` | simulate walkers through this transition set to get fluxes on each edge|
| `mk_phylogeny_precise` | given a set of states and a tree linking them, prepare the data for use in `mk.inference` |
| `mk_cross_sectional` | given a set of cross-sectional states, prepare the data for use in `mk.inference` |
| `mk_infer_cross_sectional` | simple wrapper of `mk_cross_sectional` + `mk.inference`, taking a matrix and outputting the fitted model |
| `mk.inference.plot` | plots the output of `mk.inference` |

`mk-collection.R`
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
