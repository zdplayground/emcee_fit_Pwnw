# emcee_fit_Pwnw
We apply emcee package (http://dfm.io/emcee/current/) to fit matter (galaxy) power spectrum with the Baryon Acoustic Oscillation (BAO) feature. This work is referenced largely from Florian Beutler's MCMC code. 

We can test the code by 
```
mpirun -n 4 python pre_mcmc.py input_params1.yaml
```
where we run code parallelly with mpirun, and the input_params1.yaml contains the input parameters which can be modified.
 
We use Python3 version (better in anaconda environment). mpi4py and emcee packages are needed in the python directory. We calcuate the value of fitting model and logarithm of prior of Bayesian statistics in Fortran code. Here we set uniform prior distribution. f2py is needed to compile the Fortran routines which are called in python code.  

Here we include two types of fitting models. One is called SBRS model, based on the paper (arXiv:1511.00663) Seo et al. (2016). The other is called EFT model referenced from (arXiv:1708.01297) Ding et al. (2018). Please see more details about the models in these papers.

For each type of models, we include cases fitting both pre- and post-reconstruction power spectrum. pre_mcmc.py fits pre-reconstruction, calling module from pre_log_prob_rsd_EFT.f95. If you need to use SBRS model, just modify function lnlike in pre_mcmc.py. Similarly, post_mcmc_EFT(SBRS).py fits post-reconstruction power spectrum with EFT(SBRS) model, and it uses model from post_log_prob_rsd_EFT(SBRS).f95. Input parameter files input_params1(2).yaml are for mcmc_funs.py, and input_params_post_EFT(SBRS).yaml are for post_mcmc_EFT(SBRS).py. 

In the input folder, there are linear power spectra with/without BAO signal, calculated from CAMB. kaves.W(N)10_sub_1.0000ga.dat is the simulated (observed) matter (subsampled) power spectrum with/without BAO signal at redshift z=0. In the subfolder z1mcut34, we include simulated galaxy power spectra at z=1.0 for both pre- and post-reconstruction.
Folder output_files contain output files from these code.
