!
! Use f2py -m post_lnprob_module_rsd_EFT -h post_lnprob_module_rsd_EFT.pyf post_log_prob_rsd_EFT.f95 --overwrite-signature
! to generate signature file.
! Use f2py -c --fcompiler=gnu95 post_lnprob_module_rsd_EFT.pyf post_log_prob_rsd_EFT.f95 to generate the module
!
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

subroutine cal_Pk_model(Pk_linw, Pk_sm, k_t, mu_t, sigma2_sm, sigma2_dd, sigma2_sd, sigma2_ss, f, b_0, &
                        b_scale, norm_gf, Pk_model, dim_kt)
    implicit none
    integer:: dim_kt, i
    double precision, dimension(dim_kt):: Pk_linw, Pk_sm, k_t, mu_t
    double precision:: sigma2_sm, sigma2_dd, sigma2_sd, sigma2_ss, f, b_0, b_scale, norm_gf
    double precision:: mu, k, k_2, sm_k, t1, t2, x1, delta_P_dd, delta_P_sd, delta_P_ss, Pk_model(dim_kt)
!f2py intent(in):: Pk_linw, Pk_sm, k_t, mu_t, sigma2_sm, sigma2_dd, sigma2_sd, sigma2_ss, f, b_0, b_scale, norm_gf
!f2py intent(out):: Pk_model
    do i=1, dim_kt
        mu = mu_t(i)
        k = k_t(i)
        k_2 = k*k
        sm_k = exp(-0.25d0 * k_2 * sigma2_sm)
        !!sm_k = 0.d0 ! this is only for testing pre-reconstruction case.
        t1 = f* mu*mu
        t2 = 1.d0 - sm_k
        x1 = b_0 - sm_k + t1*t2
        !call Sigma2_machine(sigma2_dd, sigma2_sd, sigma2_ss)

        delta_P_dd = exp(-k_2*(1.d0+ t1*(2.d0+f))*sigma2_dd)* x1*(x1+b_scale*k_2)
        delta_P_sd = -exp(-k_2*(1.d0+ t1)*sigma2_sd) * (x1 + 0.5d0*b_scale*k_2)*sm_k
        delta_P_ss = exp(-k_2*sigma2_ss) * sm_k**2.d0

        Pk_model(i) = (delta_P_dd - 2.d0*delta_P_sd + delta_P_ss)*(Pk_linw(i)-Pk_sm(i))* norm_gf**2.d0
    enddo
    return
end subroutine

subroutine lnprior(theta, params_indices, fix_params, lp, dim_theta, dim_params)
    implicit none
    integer:: dim_theta, dim_params
    double precision:: theta(dim_theta), params_indices(dim_params), fix_params(dim_params)
    double precision:: lp, params_array(dim_params)
    double precision:: alpha_1, alpha_2, Sigma2_sm, Sigma2_dd, Sigma2_sd, Sigma2_ss, f, b_0, b_scale

!f2py intent(in):: theta, params_indices, fix_params
!f2py intent(out):: lp
    call match_params(theta, params_indices, fix_params, params_array, dim_theta, dim_params)
    alpha_1 = params_array(1)
    alpha_2 = params_array(2)
    Sigma2_sm = params_array(3)
    Sigma2_dd = params_array(4)
    Sigma2_sd = params_array(5)
    Sigma2_ss = params_array(6)
    f = params_array(7)
    b_0 = params_array(8)
    b_scale = params_array(9)

    if (alpha_1>-1.d-7 .and. alpha_1<1.1 .and. alpha_2>-1.d-7 .and. alpha_2<1.1 .and. Sigma2_sm>-1.d-7 .and. Sigma2_sm<600.0 &
       .and. f > -1.d-7 .and. f < 2.0 .and. b_0>0.5 .and. b_0<6.0 .and. b_scale>-1000.0 .and. b_scale<1000.0) then
        lp = 0.d0
    else
        lp = -1.d30  ! return a negative infinitely large number
    endif
    !print*, lp
    return
end subroutine
