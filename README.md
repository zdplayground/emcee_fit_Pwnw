# emcee_fit_Pwnw
We apply emcee package (http://dfm.io/emcee/current/) to fit matter (galaxy) power spectrum with the Baryon Acoustic Oscillation (BAO) feature. This work is referenced largely from Florian Beutler' MCMC code. 

We can test the code by 
```
mpirun -n 4 python test_mcmc.py input_params.yaml
```
where we run code parallelly with mpirun, and the input_params.yaml contains the input parameters which you can modify them correspondingly.
 
Python code is python3 version (better in anaconda environment). We need to have mpi4py and emcee installed in the python directory. We calcuate the value of fitting model and logarithm of prior of Bayesian statistics in Fortran code. Here we set uniform prior distribution. f2py is needed to compile the Fortran routines which are called by python code.  
