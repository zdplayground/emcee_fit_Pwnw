! Copy the code from log_prob.f90, modify it to get subroutine of lnlike function in redshift space.
! The likelihood function is referenced from Seo et al. 2015. -- 10/06/2016.
! 1. Modify it to include one more parameter sigma_fog, --10/07/2016.
! 2. Realize that for pre-reconstruction, sigma_sm=\inf, which leads S(k)=0. Here we don't include S(k) term. -- 10/13/2016
! 3. Rename the file to pre-log_prob_rsd_HS.f90
! Use f2py -m pre_lnprob_module_rsd_SBRS -h pre_lnprob_module_rsd_SBRS.pyf pre_log_prob_rsd_SBRS.f95 --overwrite-signature to generate
! signature file.
! Use f2py -c --fcompiler=gnu95 pre_lnprob_module_rsd_SBRS.pyf pre_log_prob_rsd_SBRS.f95 to generate the module
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

end subroutine match_params

subroutine cal_Pk_model(Pk_linw, Pk_sm, k_t, mu_t, sigma_xy, sigma_z, sigma_fog, f, b_0, &
                       b_scale, norm_gf, Pk_model, dim_kt)
    implicit none
    integer:: dim_kt, i
    double precision, dimension(dim_kt):: Pk_linw, Pk_sm, k_t, mu_t, Pk_model
    double precision:: sigma_xy, sigma_z, sigma_fog, f, b_0, b_scale, norm_gf
    double precision:: rsd_fog_term, exp_term
!f2py intent(in):: Pk_linw, Pk_sm, k_t, mu_t, sigma_xy, sigma_z, sigma_fog, f, b_0, b_scale, norm_gf
!f2py intent(out):: Pk_model
    do i=1, dim_kt
!       smooth_term = exp(-(sigma_sm*k_t(i))**2.d0/4.d0)  ! reference from Eq.(5) in Seo et al. 2015, for post-reconstruction
!       for pre-construction, smooth_term = 0.d0
!        rsd_fog_term = (b_0+ f *mu_t(i)**2.d0*(1.d0-smooth_term)) /(1.d0+(k_t(i)*mu_t(i)* sigma_fog)**2.d0 /2.d0)
        rsd_fog_term = (b_0+ f *mu_t(i)**2.d0) /(1.d0+(k_t(i)*mu_t(i)* sigma_fog)**2.d0 /2.d0)
        exp_term = exp(k_t(i)*k_t(i) *(mu_t(i)*mu_t(i)*(sigma_xy+sigma_z)*(sigma_xy-sigma_z)-sigma_xy**2.d0)/2.d0)
        Pk_model(i) = (Pk_linw(i)-Pk_sm(i))* exp_term * (1.d0+b_scale * k_t(i)**2.d0) * (rsd_fog_term * norm_gf)**2.d0
    enddo
    return
end subroutine cal_Pk_model

subroutine lnprior(theta, params_indices, fix_params, lp, dim_theta, dim_params)
    implicit none
    integer:: dim_theta, dim_params
    double precision, dimension(dim_theta):: theta
    double precision, dimension(dim_params):: params_indices, fix_params, params_array
    double precision:: lp
    double precision:: alpha_1, alpha_2, sigma_xy, sigma_z, sigma_fog, f, b_0, b_scale
!f2py intent(in):: theta, params_indices, fix_params
!f2py intent(out):: lp
    call match_params(theta, params_indices, fix_params, params_array, dim_theta, dim_params)
    alpha_1 = params_array(1)
    alpha_2 = params_array(2)
    sigma_xy = params_array(3)
    sigma_z = params_array(4)  ! remove abs() and set sigma_xy and sigma_z positive in params_array before passed by.--10/06/2016
    sigma_fog = params_array(5)
    f = params_array(6)        ! introduce two more parameters f and sigma_fog in redshift space.
    b_0 = params_array(7)
    b_scale = params_array(8)

    if (alpha_1>-1.d-7 .and. alpha_1<1.1 .and. alpha_2>-1.d-7 .and. alpha_2<1.1 .and. sigma_xy>-1.d-7 .and. sigma_xy<15.0 &
        .and. sigma_z>-1.d-7 .and. sigma_z<15.0 .and. b_0>0. .and. b_0<6.0 .and. b_scale>-100.0 .and. b_scale<100.0 &
        .and. sigma_fog > -1.d-7 .and. sigma_fog < 20.0 .and. f > -1.d-7 .and. f < 2.0) then
        lp = 0.d0
    else
        lp = -1.d30  ! return a negative infinitely large number
    endif
    !print*, lp
    return
end subroutine lnprior
