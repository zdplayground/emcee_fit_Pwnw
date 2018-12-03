#!/Users/ding/anaconda3/bin/python
import numpy as np
from scipy import special, interpolate
from pre_lnprob_module_rsd_ZV_lagrange import match_params, cal_pk_model, lnprior

# define growth factor G(z)
def growth_factor(z, Omega_m):
    a = 1.0/(1.0+z)
    v = (1.0+z)*(Omega_m/(1.0-Omega_m))**(1.0/3.0)
    phi = np.arccos((v+1.0-3.0**0.5)/(v+1.0+3.0**0.5))
    m = (np.sin(75.0/180.0* np.pi))**2.0
    part1c = 3.0**0.25 * (1.0+ v**3.0)**0.5
# first elliptic integral
    F_elliptic = special.ellipkinc(phi, m)
# second elliptic integral
    Se_elliptic = special.ellipeinc(phi, m)
    part1 = part1c * ( Se_elliptic - 1.0/(3.0+3.0**0.5)*F_elliptic)
    part2 = (1.0 - (3.0**0.5 + 1.0)*v*v)/(v+1.0+3.0**0.5)
    d_1 = 5.0/3.0*v*(part1 + part2)
# if a goes to 0, use d_11, when z=1100, d_1 is close to d_11
#    d_11 = 1.0 - 2.0/11.0/v**3.0 + 16.0/187.0/v**6.0
    return a*d_1


# Gelman&Rubin convergence criterion, referenced from Florian's code RSDfit_challenge_hex_steps_fc_hoppper.py
def gelman_rubin_convergence(withinchainvar, meanchain, n, Nchains, ndim):

    # Calculate Gelman & Rubin diagnostic
    # 1. Remove the first half of the current chains
    # 2. Calculate the within chain and between chain variances
    # 3. estimate your variance from the within chain and between chain variance
    # 4. Calculate the potential scale reduction parameter

    meanall = np.mean(meanchain, axis=0)
    W = np.mean(withinchainvar, axis=0)
    B = np.arange(ndim,dtype=np.float)
    for jj in range(0, ndim):
        B[jj] = 0.
    for jj in range(0, Nchains):
        B = B + n*(meanall - meanchain[jj])**2/(Nchains-1.)
    estvar = (1. - 1./n)*W + B/n
    scalereduction = np.sqrt(estvar/W)

    return scalereduction

# write parameters fitted in files
def write_params(filename, params_mcmc, params_name, reduced_chi2, fix_params, dof):
    header_line = '# The fitted parameters {} (by row) with their upward and downward one sigma error, and the reduced \chi^2.\n'.format(params_name[:])
    header_line += '# And the other parameters are fixed (0. may mean unfixed) with value: {}. \n'.format(fix_params[:])
    with open(filename, 'w') as fwriter:
        fwriter.write(header_line)
        for i in range(len(params_name)):
            fwriter.write("#{0}: {1:.7f} {2:.7f} {3:.7f}\n".format(params_name[i], params_mcmc[i][0], params_mcmc[i][1], params_mcmc[i][2]))
        fwriter.write("#reduced_chi2: {0:.7f}\n".format(reduced_chi2))
        fwriter.write("#dof: {0}".format(dof))

# Define a function to set parameters which are free and which are fixed.
def set_params(all_params, params_indices, all_names, all_temperature):
    fix_params = np.array([], dtype=np.float)
    theta = np.array([], dtype=np.float)
    params_T = np.array([], dtype=np.float)
    params_name = []
    N_params = 0
    count = 0
    for i in params_indices:
        if i == 1:
            fix_params = np.append(fix_params, 0.)
            theta = np.append(theta, all_params[count])
            params_T = np.append(params_T, all_temperature[count])
            params_name.append(all_names[count])
            N_params += 1
        else:
            fix_params = np.append(fix_params, all_params[count])
        count += 1
    print(theta, params_name, N_params)
    print("fixed params: ", fix_params)
    return N_params, theta, fix_params, params_T, params_name


def lnlike(theta, params_indices, fix_params, k_p, mu_p, Pk_obs, ivar, tck_Pk_linw, tck_Pk_sm, norm_gf):
    alpha_1, alpha_2, sigma, f, b_0, b_scale = match_params(theta, params_indices, fix_params)
    # set alpha_1=alpha_2 and sigma_xy = sigma_z
    ##alpha_1 = alpha_2  # be careful to comment it if both alpha_1 and alpha_2 are free parameters.
    ##sigma_xy = sigma_z
    coeff = 1.0/alpha_1*(1.0+mu_p**2.0*(pow(alpha_1/alpha_2, 2.0)-1.0))**0.5
    k_t = k_p*coeff
    mu_t = mu_p/(alpha_2*coeff)
    Pk_linw = interpolate.splev(k_t, tck_Pk_linw, der=0)
    Pk_sm = interpolate.splev(k_t, tck_Pk_sm, der=0)

    Pk_model = cal_pk_model(Pk_linw, Pk_sm, k_t, mu_t, sigma, f, b_0, b_scale, norm_gf)
    diff = Pk_model - Pk_obs
    return -0.5* np.sum(diff**2.0 *ivar)

def lnprob(theta, params_indices, fix_params, k_p, mu_p, Pk_obs, ivar, tck_Pk_linw, tck_Pk_sm, norm_gf):
    lp = lnprior(theta, params_indices, fix_params)
    if (lp < -1.e20):
        return -np.inf
    return lp + lnlike(theta, params_indices, fix_params, k_p, mu_p, Pk_obs, ivar, tck_Pk_linw, tck_Pk_sm, norm_gf)
