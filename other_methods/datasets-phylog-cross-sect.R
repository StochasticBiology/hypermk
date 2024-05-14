### Generating/reading data sets as in the main repo
##  This includes the three simulated cross sectional,
##  ovarian, and the three simulated phylogenetic

#### How to generate

## 1. Source up to line 269 of mk-collection.R (i.e.,
##    from line 1 to the end of definition of function setup.data).
##    Line numbers: using mk-collection.R, branch igj-wip-rdu,
##    commit a873202
## 2. Then run the following code

d_cross.sectional.single <- setup.data("cross.sectional.single")
d_cross.sectional.many <- setup.data("cross.sectional.many")
d_cross.sectional.cross <- setup.data("cross.sectional.cross")
d_ovarian <- setup.data("ovarian")
d_single <- setup.data("single")
d_single.rev <- setup.data("single.rev")
d_single.uncertain <- setup.data("single.uncertain")
d_TB <-  setup.data("TB")

save(file = "setup.data_output.RData", list = ls(pattern = glob2rx("^d_*")))


## For the sake of it, verify one assertion
identical(d_single$x.decimal, d_single$tips)
identical(d_single.rev$x.decimal, d_single.rev$tips)
identical(d_TB$x.decimal, d_TB$tips)
