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
        spherical_helmholtz_dtn_eigenvalue
    implicit none

    private

    public :: apply_spherical_maxwell_dtn
    public :: spherical_maxwell_dtn_eigenvalues

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

end module fortfem_spherical_maxwell_dtn
