module fortfem_spherical_helmholtz_dtn
    !! Exact outgoing scalar Helmholtz DtN map on a spherical boundary.
    !!
    !! For each spherical-harmonic degree l, the exterior radial factor is
    !! h_l^(1)(k r), hence T_l = k h_l^(1)'(k R)/h_l^(1)(k R).
    !! Recurrences follow NIST DLMF 10.51.1--10.51.2:
    !! https://dlmf.nist.gov/10.51
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: apply_spherical_helmholtz_dtn
    public :: spherical_helmholtz_dtn_eigenvalue

contains

    subroutine apply_spherical_helmholtz_dtn( &
            degrees, trace, wavenumber, radius, normal_derivative, status)
        integer, intent(in) :: degrees(:)
        complex(dp), intent(in) :: trace(:)
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(out) :: normal_derivative(:)
        integer, intent(out) :: status

        complex(dp) :: eigenvalue
        integer :: coefficient

        status = 0
        normal_derivative = cmplx(0.0_dp, 0.0_dp, dp)
        if (size(degrees) /= size(trace) .or. &
            size(normal_derivative) /= size(trace)) then
            status = 1
            return
        end if

        do coefficient = 1, size(trace)
            call spherical_helmholtz_dtn_eigenvalue( &
                degrees(coefficient), wavenumber, radius, eigenvalue, status)
            if (status /= 0) return
            normal_derivative(coefficient) = eigenvalue*trace(coefficient)
        end do
    end subroutine apply_spherical_helmholtz_dtn

    subroutine spherical_helmholtz_dtn_eigenvalue( &
            degree, wavenumber, radius, eigenvalue, status)
        integer, intent(in) :: degree
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(out) :: eigenvalue
        integer, intent(out) :: status

        complex(dp) :: derivative, hankel, hankel_minus
        real(dp) :: argument

        eigenvalue = cmplx(0.0_dp, 0.0_dp, dp)
        status = 0
        if (degree < 0 .or. wavenumber <= 0.0_dp .or. radius <= 0.0_dp) then
            status = 1
            return
        end if

        argument = wavenumber*radius
        call spherical_hankel_h1(degree, argument, hankel)
        if (degree == 0) then
            call spherical_hankel_h1(1, argument, derivative)
            derivative = -derivative
        else
            call spherical_hankel_h1(degree - 1, argument, hankel_minus)
            derivative = hankel_minus - real(degree + 1, dp)*hankel/argument
        end if
        eigenvalue = wavenumber*derivative/hankel
    end subroutine spherical_helmholtz_dtn_eigenvalue

    pure subroutine spherical_hankel_h1(degree, argument, value)
        integer, intent(in) :: degree
        real(dp), intent(in) :: argument
        complex(dp), intent(out) :: value

        complex(dp) :: previous, current, following
        integer :: order

        previous = cmplx(sin(argument), -cos(argument), dp)/argument
        if (degree == 0) then
            value = previous
            return
        end if

        current = cmplx( &
            sin(argument)/argument**2 - cos(argument)/argument, &
            -cos(argument)/argument**2 - sin(argument)/argument, dp)
        do order = 1, degree - 1
            following = real(2*order + 1, dp)*current/argument - previous
            previous = current
            current = following
        end do
        value = current
    end subroutine spherical_hankel_h1

end module fortfem_spherical_helmholtz_dtn
