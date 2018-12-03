#!/Users/ding/anaconda3/bin/python
# Try to compact mcmc fitting routine. --11/30/2018
#
import os, sys
import time
import numpy as np
import scipy.optimize as op
from scipy import interpolate
import emcee
from emcee.utils import MPIPool
#import corner
#import matplotlib.pyplot as plt
#from matplotlib.ticker import MaxNLocator
import yaml
from mcmc_funs import growth_factor, gelman_rubin_convergence, set_params, lnlike, lnprob, write_params

def chi2(*args):
    return -2 * lnlike(*args)

# MCMC routine
def mcmc_routine(ndim, N_walkers, theta, params_T, params_indices, fix_params, k_range, mu_range, Pk_wnow_obs, ivar_Pk_wnow, tck_Pk_linw, tck_Pk_sm, norm_gf, params_name, pool):
    ti = time.time()

    Nchains = 4
    minlength = 800
    epsilon = 0.01
    ichaincheck = 50
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    #result = op.fmin_powell(chi2, theta, args=(params_indices, fix_params, k_range, mu_range, Pk_wnow_obs, ivar_Pk_wnow, tck_Pk_linw, tck_Pk_sm, norm_gf))
    # result = op.minimize(chi2, theta, args=(params_indices, fix_params, k_range, mu_range, Pk_wnow_obs, ivar_Pk_wnow, tck_Pk_linw, tck_Pk_sm, norm_gf), method='Powell')
    # theta_optimize = result["x"]
    # print("Parameters from Powell optimization: ", theta_optimize) # only output parameters which are free to change

    theta_optimize = theta
    num_alpha = params_indices[0] + params_indices[1]
    num_Sigma = sum(params_indices[2: 2+1])     # we have 1 parameter for \Sigma.
    print("# of Sigma params: ", num_Sigma)
    # in case the initial sigma for MCMC is negative, fix it to be positive.
    if params_indices[2] == 1: #
        for i in range(num_alpha, num_alpha+num_Sigma):
            theta_optimize[i] = abs(theta_optimize[i])

    print("Initial parameters for MCMC: ", theta_optimize)

    pos = []
    sampler = []
    rstate = np.random.get_state()
    # Set up the sampler.
    for jj in range(Nchains):
        pos.append([theta_optimize + params_T*np.random.uniform(-1.0,1.0, ndim) for i in range(N_walkers)])

        sampler.append(emcee.EnsembleSampler(N_walkers, ndim, lnprob, a=2.0, args=(params_indices, fix_params, k_range, mu_range, Pk_wnow_obs, ivar_Pk_wnow, tck_Pk_linw, tck_Pk_sm, norm_gf), pool=pool))
    print(type(sampler))

    # Clear and run the production chain.
    print("Running MCMC...")

    withinchainvar = np.zeros((Nchains,ndim))
    meanchain = np.zeros((Nchains,ndim))
    scalereduction = np.arange(ndim,dtype=np.float)
    for jj in range(0, ndim):
        scalereduction[jj] = 2.

    itercounter = 0
    chainstep = minlength
    loopcriteria = 1
    num_iteration = 1
    while loopcriteria and num_iteration < 50:  # If iteration is larger than 50, the fitting couldn't coverge exactly to reach elsilon<0.01.
        itercounter = itercounter + chainstep
        print("chain length =",itercounter," minlength =",minlength)

        for jj in range(Nchains):
            # Since we write the chain to a file we could put storechain=False, but in that case
            # the function sampler.get_autocorr_time() below will give an error
            for result in sampler[jj].sample(pos[jj], iterations=chainstep, rstate0=np.random.get_state(), storechain=True, thin=1):
                pos[jj] = result[0]
                #print(pos)
                chainchi2 = -2.*result[1]
                rstate = result[2]

            # we do the convergence test on the second half of the current chain (itercounter/2)
            chainsamples = sampler[jj].chain[:, itercounter//2:, :].reshape((-1, ndim))
            #print("len chain = ", chainsamples.shape)
            withinchainvar[jj] = np.var(chainsamples, axis=0)
            meanchain[jj] = np.mean(chainsamples, axis=0)

        scalereduction = gelman_rubin_convergence(withinchainvar, meanchain, itercounter//2, Nchains, ndim)
        print("scalereduction = ", scalereduction)

        loopcriteria = 0
        for jj in range(0, ndim):
            if np.absolute(1.0-scalereduction[jj]) > epsilon:
                loopcriteria = 1

        chainstep = ichaincheck
        num_iteration = num_iteration + 1
    print("Done.")

    # Print out the mean acceptance fraction. In general, acceptance_fraction
    # has an entry for each walker so, in this case, it is a 250-dimensional vector.
    for jj in range(0, Nchains):
        print("Mean acceptance fraction for chain ", jj,": ", np.mean(sampler[jj].acceptance_fraction))
        # Estimate the integrated autocorrelation time for the time series in each parameter.
        #print("Autocorrelation time for chain ", jj,": ", sampler[jj].get_autocorr_time())

    ###################################
    ## Compute the quantiles ##########
    ###################################

    mergedsamples=[]
    for jj in range(0, Nchains):
        mergedsamples.extend(sampler[jj].chain[:, itercounter//2:, :].reshape((-1, ndim)))
    print("length of merged chain = ", sum(map(len,mergedsamples))//ndim)

    theta_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(mergedsamples, [15.86555, 50, 84.13445], axis=0)))
    theta_mcmc = list(theta_mcmc)

    print("MCMC result: ")
    for i in range(len(theta)):
        print("{0}={1[0]}+{1[1]}-{1[2]}".format(params_name[i], theta_mcmc[i]))

    del sampler
    tf = time.time()
    print("One mcmc running set time: ", tf-ti)
    return np.array(theta_mcmc)


def main():
    # firstly, read one file and get k bins we want for the fitting range; It's dark matter power spectrum in redshift space
    parameter_file = sys.argv[1]
    with open(parameter_file, 'r') as fr:
        input_params = yaml.load(fr)

    data_m = np.genfromtxt(input_params['Pwig_ifile'], dtype='f8', comments='#', delimiter='', skip_header = 11) # skip the first data row, the first 10 rows are comments.
    #print(data_m)
    num_kbin = np.size(data_m, axis=0)
    kk = data_m[:, 0]
    # get indices based on the ascending k
    indices_sort = [i[0] for i in sorted(enumerate(kk), key=lambda x: x[1])]
    # sort out the indices whose k<=0.3 h/Mpc
    for i in range(num_kbin):
        if kk[indices_sort[i]] > 0.3:
            break
    print(indices_sort[i-1])
    indices_p = indices_sort[0: i]
    k_p = kk[indices_p]
    N_fitbin = len(k_p)

    mu_p, Pwig = data_m[indices_p, 1], data_m[indices_p, 2]
    print(k_p, N_fitbin)

    # input Pnow, note the (k, mu) indices have the same order as those of Pwig data file
    data_m = np.genfromtxt(input_params['Pnow_ifile'], dtype='f8', comments='#', delimiter='', skip_header = 11)
    Pnow = data_m[indices_p, 2]
    Pwnw_diff_obs = Pwig - Pnow

    # input diagonal terms of covariance matrix of (Pwig-Pnow)
    diag_Cov_Pwnw = np.loadtxt(input_params['diag_Cov_Pwnw_ifile'], dtype='f8', comments='#', usecols=(2,))
    ivar_Pk_wnow = 1.0/diag_Cov_Pwnw
    #print(ivar_Pk_wnow)

    # input (theoretical) linear power spectrum
    k_wiggle, Pk_wiggle = np.loadtxt(input_params['Pwig_linear'], dtype='f8', comments='#', unpack=True)
    tck_Pk_linw = interpolate.splrep(k_wiggle, Pk_wiggle)

    k_smooth, Pk_smooth = np.loadtxt(input_params['Pnow_linear'], dtype='f8', comments='#', unpack=True)
    tck_Pk_sm = interpolate.splrep(k_smooth, Pk_smooth)

    all_params = list(input_params['init_params'].values())
    params_indices = input_params['params_indices']
    params_name = list(input_params['init_params'].keys())
    all_temperature = input_params['all_temperature']

    N_params, theta, fix_params, params_T, params_name = set_params(all_params, params_indices, params_name, all_temperature)
    #print(N_params, theta, fix_params, params_T, params_name)
    sim_z = input_params['sim_z']    # redshift of the simulated power spectrum
    N_walkers = input_params['N_walkers']
    Omega_m = input_params['Omega_m']

    G_0 = growth_factor(0.0, Omega_m) # G_0 at z=0, normalization factor
    norm_gf = growth_factor(sim_z, Omega_m)/G_0

    pool = MPIPool(loadbalance=True)
    np.random.seed(1)    # set random seed for random number generator
    params_mcmc = mcmc_routine(N_params, N_walkers, theta, params_T, params_indices, fix_params, k_p, mu_p, Pwnw_diff_obs, ivar_Pk_wnow, tck_Pk_linw, tck_Pk_sm, norm_gf, params_name, pool)
    print(params_mcmc)
    chi_square = chi2(params_mcmc[:, 0], params_indices, fix_params, k_p, mu_p, Pwnw_diff_obs, ivar_Pk_wnow, tck_Pk_linw, tck_Pk_sm, norm_gf)
    dof = N_fitbin-N_params
    reduced_chi2 = chi_square/dof
    odir = './output_files/'
    if not os.path.exists(odir):
        os.makedirs(odir)
    ofile_params = odir + input_params['ofile_name'].format(sim_z, ''.join(map(str, params_indices)))
    write_params(ofile_params, params_mcmc, params_name, reduced_chi2, fix_params, dof)
    pool.close()


if __name__ == '__main__':
    main()
