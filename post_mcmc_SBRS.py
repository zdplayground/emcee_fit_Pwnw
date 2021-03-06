#!/Users/ding/anaconda3/bin/python
# We can set parameters alpha_1(alpha_perp), alpha_2(alpha_para), Sigma_xy, Sigma_z, Sigma_fog, f, b and b_scale as free or fixed by setting params_indices.
# The fitting model is
# (P_wig - P_now)_obs = (1+b_scale *k^2)G^2 *(b_0+ f* mu^2*(1-S(k)))^2 /(1+k^2 *mu^2*sigma_fog^2)^2*(Plin - Psm) C_G^2, where both b_0, b_scale, f, and sigma_fog and Smoothing factor are the fitting parameters.
# Compared with the real space, there are three more parameters f, Sigma_sm counting Kaiser effect and sigma_fog counting FoG effect, respectively.
import os, sys
import time
import numpy as np
import scipy.optimize as op
from scipy import interpolate, integrate
import emcee
from emcee.utils import MPIPool
#import corner
#import matplotlib.pyplot as plt
#from matplotlib.ticker import MaxNLocator
import yaml
from mcmc_funs import growth_factor, gelman_rubin_convergence, set_params, write_params
from post_lnprob_module_rsd_SBRS import match_params, cal_pk_model, lnprior

def Sigma2_dd_integrand(k, tck_Pk_linw, R_bao, Sigma2_sm):
    Pk_lin_0 = interpolate.splev(k, tck_Pk_linw, der=0)
    #print(k, Sigma2_sm)
    Sm_kernel = np.exp(-0.25*k*k * Sigma2_sm)
    return Pk_lin_0*(1.0- np.sin(k*R_bao)/(k*R_bao))*(1.0-Sm_kernel)**2.0

# input Pk_obs = \hat{P}_wig - \hat{P}_now; sigma_fog represents \Sigma_s in the damping term of finger-of-god
def lnlike(theta, params_indices, fix_params, k_p, mu_p, Pk_obs, ivar, tck_Pk_linw, tck_Pk_sm, norm_gf):
    alpha_1, alpha_2, sigma_xy, sigma_z, sigma_sm, sigma_fog, f, b_0, b_scale = match_params(theta, params_indices, fix_params)
    # set alpha_1=alpha_2 and sigma_xy = sigma_z
    #alpha_1 = alpha_2  # be careful to comment it if both alpha_1 and alpha_2 are free parameters.
    ##sigma_xy = sigma_z
    coeff = 1.0/alpha_1*(1.0+mu_p**2.0*(pow(alpha_1/alpha_2, 2.0)-1.0))**0.5
    k_t = k_p*coeff
    mu_t = mu_p/(alpha_2*coeff)
    Pk_linw = interpolate.splev(k_t, tck_Pk_linw, der=0)
    Pk_sm = interpolate.splev(k_t, tck_Pk_sm, der=0)

    Pk_model = cal_pk_model(Pk_linw, Pk_sm, k_t, mu_t, sigma_xy, sigma_z, sigma_sm, sigma_fog, f, b_0, b_scale, norm_gf)
    diff = Pk_model - Pk_obs
    return -0.5* np.sum(diff**2.0 *ivar)

def lnprob(theta, params_indices, fix_params, k_p, mu_p, Pk_obs, ivar, tck_Pk_linw, tck_Pk_sm, norm_gf):
    lp = lnprior(theta, params_indices, fix_params)
    if (lp < -1.e20):
        return -np.inf
    return lp + lnlike(theta, params_indices, fix_params, k_p, mu_p, Pk_obs, ivar, tck_Pk_linw, tck_Pk_sm, norm_gf)

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

    #all_params = list(input_params['init_params'].values())
    params_indices = input_params['params_indices']
    params_name = list(input_params['init_params'].keys())
    all_temperature = input_params['all_temperature']

    #print(N_params, theta, fix_params, params_T, params_name)
    sim_z = input_params['sim_z']    # redshift of the simulated power spectrum
    N_walkers = input_params['N_walkers']
    Omega_m = input_params['Omega_m']

    G_0 = growth_factor(0.0, Omega_m) # G_0 at z=0, normalization factor
    norm_gf = growth_factor(sim_z, Omega_m)/G_0
    # if Sigma_sm is fixed by theoretical value, we calculate other Sigma parameters theoretically in case they are fixed as well
    q_max = 110.0  # Mpc/h, BAO radius
    const = 1.0/(6.0*np.pi**2.0) * norm_gf**2.0
    if params_indices[4] == 0:
        Sigma2_sm = input_params['init_params']['Sigma_sm']**2.0
        Sigma2_dd = const * integrate.quad(Sigma2_dd_integrand, k_wiggle[0], 100.0, args=(tck_Pk_linw, q_max, Sigma2_sm), epsabs=1.e-4, epsrel=1.e-4)[0]
        sigma_xy = (2.0*Sigma2_dd)**0.5         # There is factor 2 due to different expression bewteen two models.
        sigma_z = (1.0+input_params['f_theory'])*sigma_xy
        input_params['init_params']['Sigma_xy'] = sigma_xy
        input_params['init_params']['Sigma_z'] = sigma_z

    all_params = list(input_params['init_params'].values())
    N_params, theta, fix_params, params_T, params_name = set_params(all_params, params_indices, params_name, all_temperature)

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
