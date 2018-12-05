# emcee_fit_Pwnw
We apply emcee package (http://dfm.io/emcee/current/) to fit matter (galaxy) power spectrum with the Baryon Acoustic Oscillation (BAO) feature. This work is referenced largely from Florian Beutler' MCMC code. 

We can test the code by 
```
mpirun -n 4 python pre_mcmc.py input_params1.yaml
```
where we run code parallelly with mpirun, and the input_params1.yaml contains the input parameters which you can modify.
 
We use Python3 version (better in anaconda environment). mpi4py and emcee are required in the python directory. We calcuate the value of fitting model and logarithm of prior of Bayesian statistics in Fortran code. Here we set uniform prior distribution. f2py is needed to compile the Fortran routines which are called by python code.  

Here we include two types of fitting models. One is called SBRS model, based on the paper (arXiv:1511.00663) Seo et al. (2016). The other is called EFT model referenced from (arXiv:1708.01297) Ding et al. (2018). Please more details about the models in these papers.

For each type of models, we show cases for both pre- and post-reconstruction power spectrum, which are run by pre_mcmc.py and post_mcmc.py respectively. 
