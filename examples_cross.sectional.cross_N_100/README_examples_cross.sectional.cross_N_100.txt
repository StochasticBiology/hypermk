What
====

R script and png files from running cross.sectional.cross example but with larger sample size: N = 100, obtained by adding rows to matrix "m" from the original "m" (by resampling with replacement).  This is done repeatedly, so there are a total of 18 pngs.



Input files and output
======================

Two different sets of runs:

-  Using optim (not nlminb) for speed, and Ntrials = 1 to avoid crashes.

   - R file: mk-collection-cross.sectional.cross-N-100-optim-one-trial.R

    - pngs: expt-cross.sectional.cross-N-100_nrun_[1-10]optim_single_trial.png

- Using nlminb, and Ntrials = 20

   - R file: mk-collection-cross.sectional.cross-N-100.R

   - pngs:
     expt-cross.sectional.cross-N-100_nrun_[1-8].png
     (why not 9 and 10? because I accidentally killed to running R process)


- According to castor, all of the runs converged.

- The datasets used are always different (as we resample randomly).



Why?
====

- To compare nlminb with optim, and see the stability of recovering the single true paths, when using an N that could be representative of real data sets.
