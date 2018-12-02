! Fit pre-reconstruction power spectrum in redshift space.
! The likelihood function is referenced from Zvonimir Vlah. -- 11/10/2016.
! To match the exponential term in Seo et al. (2015), we change \Sigma to \Sigma/sqrt(2) in Zvonimir model, i.e.,
! (Pwig-Pnow)_obs = (b1^2+2*f*b1*nu^2+f^2*nu^4+(b1*bp+f*bp*nu^2)*(k/k_L)^2)*(Pwig-Pnow)_linear*exp(-k^2*(1+f*(2+f)*nu^2)*Sigma^2/2). --11/21/2016.
! b1: galaxy bias (scale independent)
! bp: higher order galaxy bias
! f: structure growth factor
! nu: cosine of the angle between k vector and the line-of-sight
! Sigma: BAO wiggle damping factor
!
!----------------- Need to recompile the code after change is made in the code. ----------------------------
! Use f2py -m pre_lnprob_module_rsd_ZV_lagrange -h pre_lnprob_module_rsd_ZV_lagrange.pyf pre_log_prob_rsd_ZV_lagrange.f95 --overwrite-signature to generate
! signature file.
! Use f2py -c --fcompiler=gnu95 pre_lnprob_module_rsd_ZV_lagrange.pyf pre_log_prob_rsd_ZV_lagrange.f95 to generate the module
!
subroutine match_params(theta, params_indices, fix_params, params_array, dim_theta, dim_params)
    implicit none
    integer:: dim_theta, dim_params, count, i
    double precision:: theta(dim_theta), params_indices(dim_params), fix_params(dim_params)
!f2py intent(in):: theta, params_indices, fix_params
    double precision:: params_array(dim_params)
!f2py intent(out):: params_array

    count = 1  ! be very careful that the starting value is different from Python's. It's 1 in fortran!
    do i=1, dim_params
        if (params_indices(i) == 1) then
            params_array(i) = theta(count)
            count = count + 1
        else
            params_array(i) = fix_params(i)
        endif
    end do
    return

end subroutine

subroutine cal_Pk_model(Pk_linw, Pk_sm, k_t, mu_t, sigma, f, b_0, b_scale, norm_gf, Pk_model, dim_kt)
    implicit none
    integer:: dim_kt, i
    double precision, dimension(dim_kt):: Pk_linw, Pk_sm, k_t, mu_t
    double precision:: sigma, f, b_0, b_scale, norm_gf
    double precision:: mu, k, t1, t2, Pk_model(dim_kt)
!f2py intent(in):: Pk_linw, Pk_sm, k_t, mu_t, sigma, f, b_0, b_scale, norm_gf
!f2py intent(out):: Pk_model
    do i=1, dim_kt
        mu = mu_t(i)
        k = k_t(i)
        t1 = b_0*b_0 + 2.d0*f*b_0*mu*mu + f*f*mu**4.d0 + &
             (b_0*b_scale + f*b_scale*mu*mu)*k*k
        t2 = exp(-(k*sigma)**2.d0/2.d0 * (1.d0+f*(2.d0 + f)*mu*mu))
        Pk_model(i) = (Pk_linw(i)-Pk_sm(i))*t1*t2* norm_gf**2.d0
    enddo
    return
end subroutine

subroutine lnprior(theta, params_indices, fix_params, lp, dim_theta, dim_params)
    implicit none
    integer:: dim_theta, dim_params
    double precision:: theta(dim_theta), params_indices(dim_params), fix_params(dim_params)
    double precision:: lp
    double precision:: params_array(dim_params), alpha_1, alpha_2, sigma, f, b_0, b_scale

!f2py intent(in):: theta, params_indices, fix_params
!f2py intent(out):: lp
    call match_params(theta, params_indices, fix_params, params_array, dim_theta, dim_params)
    alpha_1 = params_array(1)
    alpha_2 = params_array(2)
    sigma = params_array(3)
    f = params_array(4)
    b_0 = params_array(5)
    b_scale = params_array(6)

    if (alpha_1>-1.d-7 .and. alpha_1<1.1 .and. alpha_2>-1.d-7 .and. alpha_2<1.1 .and. sigma>-1.d-7 .and. sigma<100.0 &
        .and. f > -1.d-7 .and. f < 4.0 .and. b_0>0. .and. b_0<6.0 .and. b_scale>-1000.0 .and. b_scale<1000.0) then
        lp = 0.d0
    else
        lp = -1.d30  ! return a negative infinitely large number
    endif
    !print*, lp
    return
end subroutine
