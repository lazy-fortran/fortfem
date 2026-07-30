module fortfem_spherical_maxwell_dtn
    !! Exact outgoing Maxwell capacity operator on a spherical boundary.
    !!
    !! The map is T(E_t) = n x curl(E). In the tangential vector spherical
    !! harmonic basis, TE and TM modes of degree l >= 1 are diagonal. Writing
    !! d_l = k h_l'(kR)/h_l(kR), their eigenvalues are
    !!
    !!   T_TE = -(d_l + 1/R),
    !!   T_TM = k^2/(d_l + 1/R).
    !!
    !! These follow from the M/N vector spherical waves and the
    !! Riccati-Hankel derivative (x h_l)'/h_l. The scalar logarithmic
    !! derivative is shared with the spherical Helmholtz DtN implementation.
    use fortfem_kinds, only: dp
    use fortfem_spherical_helmholtz_dtn, only: &
        spherical_helmholtz_dtn_eigenvalue, &
        spherical_helmholtz_dtn_eigenvalue_jvp
    implicit none

    private

    public :: apply_spherical_maxwell_dtn
    public :: apply_spherical_maxwell_dtn_jvp
    public :: apply_spherical_maxwell_dtn_vjp
    public :: spherical_maxwell_dtn_eigenvalues
    public :: spherical_maxwell_dtn_eigenvalues_jvp
    public :: spherical_maxwell_dtn_eigenvalues_vjp

contains

    subroutine apply_spherical_maxwell_dtn( &
            degrees, te_trace, tm_trace, wavenumber, radius, te_result, &
            tm_result, status)
        integer, intent(in) :: degrees(:)
        complex(dp), intent(in) :: te_trace(:), tm_trace(:)
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(out) :: te_result(:), tm_result(:)
        integer, intent(out) :: status

        complex(dp) :: te_eigenvalue, tm_eigenvalue
        integer :: coefficient

        te_result = cmplx(0.0_dp, 0.0_dp, dp)
        tm_result = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (size(degrees) /= size(te_trace)) return
        if (size(tm_trace) /= size(te_trace)) return
        if (size(te_result) /= size(te_trace)) return
        if (size(tm_result) /= size(te_trace)) return
        do coefficient = 1, size(degrees)
            call spherical_maxwell_dtn_eigenvalues( &
                degrees(coefficient), wavenumber, radius, te_eigenvalue, &
                tm_eigenvalue, status)
            if (status /= 0) return
            te_result(coefficient) = te_eigenvalue*te_trace(coefficient)
            tm_result(coefficient) = tm_eigenvalue*tm_trace(coefficient)
        end do
        status = 0
    end subroutine apply_spherical_maxwell_dtn

    subroutine apply_spherical_maxwell_dtn_jvp( &
            degrees, te_trace, tm_trace, wavenumber, radius, te_trace_dot, &
            tm_trace_dot, wavenumber_dot, radius_dot, te_result_dot, &
            tm_result_dot, status)
        integer, intent(in) :: degrees(:)
        complex(dp), intent(in) :: te_trace(:), tm_trace(:)
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(in) :: te_trace_dot(:), tm_trace_dot(:)
        real(dp), intent(in) :: wavenumber_dot, radius_dot
        complex(dp), intent(out) :: te_result_dot(:), tm_result_dot(:)
        integer, intent(out) :: status

        complex(dp) :: te_eigenvalue, te_eigenvalue_dot
        complex(dp) :: tm_eigenvalue, tm_eigenvalue_dot
        integer :: coefficient

        te_result_dot = cmplx(0.0_dp, 0.0_dp, dp)
        tm_result_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call validate_product_inputs( &
            degrees, te_trace, tm_trace, te_trace_dot, tm_trace_dot, &
            te_result_dot, tm_result_dot, status)
        if (status /= 0) return
        do coefficient = 1, size(degrees)
            call spherical_maxwell_dtn_eigenvalues_jvp( &
                degrees(coefficient), wavenumber, radius, wavenumber_dot, &
                radius_dot, te_eigenvalue, tm_eigenvalue, te_eigenvalue_dot, &
                tm_eigenvalue_dot, status)
            if (status /= 0) return
            te_result_dot(coefficient) = &
                te_eigenvalue*te_trace_dot(coefficient) + &
                te_eigenvalue_dot*te_trace(coefficient)
            tm_result_dot(coefficient) = &
                tm_eigenvalue*tm_trace_dot(coefficient) + &
                tm_eigenvalue_dot*tm_trace(coefficient)
        end do
        status = 0
    end subroutine apply_spherical_maxwell_dtn_jvp

    subroutine apply_spherical_maxwell_dtn_vjp( &
            degrees, te_trace, tm_trace, wavenumber, radius, te_result_bar, &
            tm_result_bar, te_trace_bar, tm_trace_bar, wavenumber_bar, &
            radius_bar, status)
        integer, intent(in) :: degrees(:)
        complex(dp), intent(in) :: te_trace(:), tm_trace(:)
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(in) :: te_result_bar(:), tm_result_bar(:)
        complex(dp), intent(out) :: te_trace_bar(:), tm_trace_bar(:)
        real(dp), intent(out) :: wavenumber_bar, radius_bar
        integer, intent(out) :: status

        complex(dp) :: te_eigenvalue, te_eigenvalue_bar
        complex(dp) :: tm_eigenvalue, tm_eigenvalue_bar
        real(dp) :: local_radius_bar, local_wavenumber_bar
        integer :: coefficient

        te_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        tm_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wavenumber_bar = 0.0_dp
        radius_bar = 0.0_dp
        call validate_product_inputs( &
            degrees, te_trace, tm_trace, te_result_bar, tm_result_bar, &
            te_trace_bar, tm_trace_bar, status)
        if (status /= 0) return
        do coefficient = 1, size(degrees)
            call spherical_maxwell_dtn_eigenvalues( &
                degrees(coefficient), wavenumber, radius, te_eigenvalue, &
                tm_eigenvalue, status)
            if (status /= 0) return
            te_trace_bar(coefficient) = &
                conjg(te_eigenvalue)*te_result_bar(coefficient)
            tm_trace_bar(coefficient) = &
                conjg(tm_eigenvalue)*tm_result_bar(coefficient)
            te_eigenvalue_bar = &
                conjg(te_trace(coefficient))*te_result_bar(coefficient)
            tm_eigenvalue_bar = &
                conjg(tm_trace(coefficient))*tm_result_bar(coefficient)
            call spherical_maxwell_dtn_eigenvalues_vjp( &
                degrees(coefficient), wavenumber, radius, te_eigenvalue_bar, &
                tm_eigenvalue_bar, local_wavenumber_bar, local_radius_bar, &
                status)
            if (status /= 0) return
            wavenumber_bar = wavenumber_bar + local_wavenumber_bar
            radius_bar = radius_bar + local_radius_bar
        end do
        status = 0
    end subroutine apply_spherical_maxwell_dtn_vjp

    subroutine spherical_maxwell_dtn_eigenvalues( &
            degree, wavenumber, radius, te_eigenvalue, tm_eigenvalue, status)
        integer, intent(in) :: degree
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(out) :: te_eigenvalue, tm_eigenvalue
        integer, intent(out) :: status

        complex(dp) :: logarithmic_derivative, riccati_derivative

        te_eigenvalue = cmplx(0.0_dp, 0.0_dp, dp)
        tm_eigenvalue = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (degree < 1 .or. wavenumber <= 0.0_dp .or. radius <= 0.0_dp) return
        call spherical_helmholtz_dtn_eigenvalue( &
            degree, wavenumber, radius, logarithmic_derivative, status)
        if (status /= 0) return
        riccati_derivative = logarithmic_derivative + 1.0_dp/radius
        if (abs(riccati_derivative) <= tiny(1.0_dp)) return
        te_eigenvalue = -riccati_derivative
        tm_eigenvalue = wavenumber**2/riccati_derivative
        status = 0
    end subroutine spherical_maxwell_dtn_eigenvalues

    subroutine spherical_maxwell_dtn_eigenvalues_jvp( &
            degree, wavenumber, radius, wavenumber_dot, radius_dot, &
            te_eigenvalue, tm_eigenvalue, te_eigenvalue_dot, &
            tm_eigenvalue_dot, status)
        integer, intent(in) :: degree
        real(dp), intent(in) :: wavenumber, radius
        real(dp), intent(in) :: wavenumber_dot, radius_dot
        complex(dp), intent(out) :: te_eigenvalue, tm_eigenvalue
        complex(dp), intent(out) :: te_eigenvalue_dot, tm_eigenvalue_dot
        integer, intent(out) :: status

        complex(dp) :: logarithmic_derivative, logarithmic_derivative_dot
        complex(dp) :: riccati_derivative, riccati_derivative_dot

        te_eigenvalue_dot = cmplx(0.0_dp, 0.0_dp, dp)
        tm_eigenvalue_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call spherical_maxwell_dtn_eigenvalues( &
            degree, wavenumber, radius, te_eigenvalue, tm_eigenvalue, status)
        if (status /= 0) return
        call spherical_helmholtz_dtn_eigenvalue_jvp( &
            degree, wavenumber, radius, wavenumber_dot, radius_dot, &
            logarithmic_derivative, logarithmic_derivative_dot, status)
        if (status /= 0) return
        riccati_derivative = logarithmic_derivative + 1.0_dp/radius
        riccati_derivative_dot = &
            logarithmic_derivative_dot - radius_dot/radius**2
        te_eigenvalue_dot = -riccati_derivative_dot
        tm_eigenvalue_dot = &
            2.0_dp*wavenumber*wavenumber_dot/riccati_derivative - &
            wavenumber**2*riccati_derivative_dot/riccati_derivative**2
        status = 0
    end subroutine spherical_maxwell_dtn_eigenvalues_jvp

    subroutine spherical_maxwell_dtn_eigenvalues_vjp( &
            degree, wavenumber, radius, te_eigenvalue_bar, tm_eigenvalue_bar, &
            wavenumber_bar, radius_bar, status)
        integer, intent(in) :: degree
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(in) :: te_eigenvalue_bar, tm_eigenvalue_bar
        real(dp), intent(out) :: wavenumber_bar, radius_bar
        integer, intent(out) :: status

        complex(dp) :: te_eigenvalue, tm_eigenvalue
        complex(dp) :: te_eigenvalue_dot, tm_eigenvalue_dot

        call spherical_maxwell_dtn_eigenvalues_jvp( &
            degree, wavenumber, radius, 1.0_dp, 0.0_dp, te_eigenvalue, &
            tm_eigenvalue, te_eigenvalue_dot, tm_eigenvalue_dot, status)
        if (status /= 0) then
            wavenumber_bar = 0.0_dp
            radius_bar = 0.0_dp
            return
        end if
        wavenumber_bar = real( &
            conjg(te_eigenvalue_bar)*te_eigenvalue_dot + &
            conjg(tm_eigenvalue_bar)*tm_eigenvalue_dot, dp)
        call spherical_maxwell_dtn_eigenvalues_jvp( &
            degree, wavenumber, radius, 0.0_dp, 1.0_dp, te_eigenvalue, &
            tm_eigenvalue, te_eigenvalue_dot, tm_eigenvalue_dot, status)
        if (status /= 0) then
            radius_bar = 0.0_dp
            return
        end if
        radius_bar = real( &
            conjg(te_eigenvalue_bar)*te_eigenvalue_dot + &
            conjg(tm_eigenvalue_bar)*tm_eigenvalue_dot, dp)
    end subroutine spherical_maxwell_dtn_eigenvalues_vjp

    subroutine validate_product_inputs( &
            degrees, te_trace, tm_trace, te_seed, tm_seed, te_result, &
            tm_result, status)
        integer, intent(in) :: degrees(:)
        complex(dp), intent(in) :: te_trace(:), tm_trace(:)
        complex(dp), intent(in) :: te_seed(:), tm_seed(:)
        complex(dp), intent(in) :: te_result(:), tm_result(:)
        integer, intent(out) :: status

        status = 1
        if (size(degrees) /= size(te_trace)) return
        if (size(tm_trace) /= size(te_trace)) return
        if (size(te_seed) /= size(te_trace)) return
        if (size(tm_seed) /= size(te_trace)) return
        if (size(te_result) /= size(te_trace)) return
        if (size(tm_result) /= size(te_trace)) return
        status = 0
    end subroutine validate_product_inputs

end module fortfem_spherical_maxwell_dtn
