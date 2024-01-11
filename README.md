# hypermk
HyperMk for reversible accumulation modelling

HyperMk uses a hypercubic transition matrix and the Mk (Markov k-state) model [1,2] from phylogenetics to model accumulation processes, including reversibility. The code here uses the `castor` library in R [3] to provide flexible specification and fitting of the Mk model. `mk-collection.R` creates a set of synthetic examples and demonstrates inference using HyperMk. `mk-timing.R` scans through properties of synthetic datasets to demonstrate the time scaling of this implementation. `Data` contains data on drug resistance evolution in tuberculosis from [4], which is used as a case study.

TO DO: add details of cancer case study

References
---

Chapters 7 and 8 of https://bio.libretexts.org/Bookshelves/Evolutionary_Developmental_Biology/Phylogenetic_Comparative_Methods_(Harmon) are excellent reading on the Mk model.

[1] Pagel, M., 1994. Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters. Proceedings of the Royal Society of London. Series B: Biological Sciences, 255(1342), pp.37-45.

[2] Lewis, P.O., 2001. A likelihood approach to estimating phylogeny from discrete morphological character data. Systematic biology, 50(6), pp.913-925.

[3] Louca, S. and Doebeli, M., 2018. Efficient comparative phylogenetics on large trees. Bioinformatics, 34(6), pp.1053-1055.

[4] Casali, N., Nikolayevskyy, V., Balabanova, Y., Harris, S.R., Ignatyeva, O., Kontsevaya, I., Corander, J., Bryant, J., Parkhill, J., Nejentsev, S. and Horstmann, R.D., 2014. Evolution and transmission of drug-resistant tuberculosis in a Russian population. Nature genetics, 46(3), pp.279-286.
