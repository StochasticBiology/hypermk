What
====

R scripts and png files from running cross.sectional.cross example but with larger sample size: N = 50, obtained by sampling with replacement from the original 6 rows of matrix m. This is done 10 times (10 different random seeds), and the data are analysed with "optim" and "nlminb" optimizers.  Thus, there are a total of 20 pngs.



Input files and output
======================

Two different sets of runs:

-  Using optim and Ntrials = 1 to avoid crashes.

   - R file: mk-collection-cross.sectional.cross-N-50-optim-one-trial.R

    - pngs: expt-cross.sectional.cross-N-50_nrun_[1-10]optim_single_trial.png

- Using nlminb, and Ntrials = 20

   - R file: mk-collection-cross.sectional.cross-N-50.R

   - pngs:
     expt-cross.sectional.cross-N-50_nrun_[1-10].png

- According to castor, all of the runs converged.


Why?
====

- To compare nlminb with optim, and see the stability of recovering the single true paths, when using an N that could be representative of real data sets.
