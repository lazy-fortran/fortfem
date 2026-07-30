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
    public :: apply_spherical_helmholtz_dtn_jvp
    public :: apply_spherical_helmholtz_dtn_vjp
    public :: spherical_helmholtz_dtn_eigenvalue
    public :: spherical_helmholtz_dtn_eigenvalue_jvp
    public :: spherical_helmholtz_dtn_eigenvalue_vjp

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

    subroutine apply_spherical_helmholtz_dtn_jvp( &
            degrees, trace, wavenumber, radius, trace_dot, wavenumber_dot, &
            radius_dot, normal_derivative_dot, status)
        integer, intent(in) :: degrees(:)
        complex(dp), intent(in) :: trace(:), trace_dot(:)
        real(dp), intent(in) :: wavenumber, radius
        real(dp), intent(in) :: wavenumber_dot, radius_dot
        complex(dp), intent(out) :: normal_derivative_dot(:)
        integer, intent(out) :: status

        complex(dp) :: eigenvalue, eigenvalue_dot
        integer :: coefficient

        normal_derivative_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call validate_product_inputs( &
            degrees, trace, trace_dot, normal_derivative_dot, status)
        if (status /= 0) return
        do coefficient = 1, size(trace)
            call spherical_helmholtz_dtn_eigenvalue_jvp( &
                degrees(coefficient), wavenumber, radius, wavenumber_dot, &
                radius_dot, eigenvalue, eigenvalue_dot, status)
            if (status /= 0) return
            normal_derivative_dot(coefficient) = &
                eigenvalue*trace_dot(coefficient) + &
                eigenvalue_dot*trace(coefficient)
        end do
        status = 0
    end subroutine apply_spherical_helmholtz_dtn_jvp

    subroutine apply_spherical_helmholtz_dtn_vjp( &
            degrees, trace, wavenumber, radius, normal_derivative_bar, &
            trace_bar, wavenumber_bar, radius_bar, status)
        integer, intent(in) :: degrees(:)
        complex(dp), intent(in) :: trace(:), normal_derivative_bar(:)
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(out) :: trace_bar(:)
        real(dp), intent(out) :: wavenumber_bar, radius_bar
        integer, intent(out) :: status

        complex(dp) :: eigenvalue, eigenvalue_bar
        real(dp) :: local_radius_bar, local_wavenumber_bar
        integer :: coefficient

        trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wavenumber_bar = 0.0_dp
        radius_bar = 0.0_dp
        call validate_product_inputs( &
            degrees, trace, normal_derivative_bar, trace_bar, status)
        if (status /= 0) return
        do coefficient = 1, size(trace)
            call spherical_helmholtz_dtn_eigenvalue( &
                degrees(coefficient), wavenumber, radius, eigenvalue, status)
            if (status /= 0) return
            trace_bar(coefficient) = &
                conjg(eigenvalue)*normal_derivative_bar(coefficient)
            eigenvalue_bar = &
                conjg(trace(coefficient))*normal_derivative_bar(coefficient)
            call spherical_helmholtz_dtn_eigenvalue_vjp( &
                degrees(coefficient), wavenumber, radius, eigenvalue_bar, &
                local_wavenumber_bar, local_radius_bar, status)
            if (status /= 0) return
            wavenumber_bar = wavenumber_bar + local_wavenumber_bar
            radius_bar = radius_bar + local_radius_bar
        end do
        status = 0
    end subroutine apply_spherical_helmholtz_dtn_vjp

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

    subroutine spherical_helmholtz_dtn_eigenvalue_jvp( &
            degree, wavenumber, radius, wavenumber_dot, radius_dot, &
            eigenvalue, eigenvalue_dot, status)
        integer, intent(in) :: degree
        real(dp), intent(in) :: wavenumber, radius
        real(dp), intent(in) :: wavenumber_dot, radius_dot
        complex(dp), intent(out) :: eigenvalue, eigenvalue_dot
        integer, intent(out) :: status

        complex(dp) :: ratio, ratio_argument_derivative
        real(dp) :: argument, argument_dot

        eigenvalue_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call spherical_helmholtz_dtn_eigenvalue( &
            degree, wavenumber, radius, eigenvalue, status)
        if (status /= 0) return
        argument = wavenumber*radius
        argument_dot = wavenumber_dot*radius + wavenumber*radius_dot
        ratio = eigenvalue/wavenumber
        ratio_argument_derivative = &
            -2.0_dp*ratio/argument - 1.0_dp + &
            real(degree*(degree + 1), dp)/argument**2 - ratio**2
        eigenvalue_dot = &
            wavenumber_dot*ratio + &
            wavenumber*ratio_argument_derivative*argument_dot
        status = 0
    end subroutine spherical_helmholtz_dtn_eigenvalue_jvp

    subroutine spherical_helmholtz_dtn_eigenvalue_vjp( &
            degree, wavenumber, radius, eigenvalue_bar, wavenumber_bar, &
            radius_bar, status)
        integer, intent(in) :: degree
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(in) :: eigenvalue_bar
        real(dp), intent(out) :: wavenumber_bar, radius_bar
        integer, intent(out) :: status

        complex(dp) :: eigenvalue, eigenvalue_dot

        call spherical_helmholtz_dtn_eigenvalue_jvp( &
            degree, wavenumber, radius, 1.0_dp, 0.0_dp, eigenvalue, &
            eigenvalue_dot, status)
        if (status /= 0) then
            wavenumber_bar = 0.0_dp
            radius_bar = 0.0_dp
            return
        end if
        wavenumber_bar = real(conjg(eigenvalue_bar)*eigenvalue_dot, dp)
        call spherical_helmholtz_dtn_eigenvalue_jvp( &
            degree, wavenumber, radius, 0.0_dp, 1.0_dp, eigenvalue, &
            eigenvalue_dot, status)
        if (status /= 0) then
            radius_bar = 0.0_dp
            return
        end if
        radius_bar = real(conjg(eigenvalue_bar)*eigenvalue_dot, dp)
    end subroutine spherical_helmholtz_dtn_eigenvalue_vjp

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

    subroutine validate_product_inputs( &
            degrees, trace, seed, result, status)
        integer, intent(in) :: degrees(:)
        complex(dp), intent(in) :: trace(:), seed(:), result(:)
        integer, intent(out) :: status

        status = 1
        if (size(degrees) /= size(trace)) return
        if (size(seed) /= size(trace)) return
        if (size(result) /= size(trace)) return
        status = 0
    end subroutine validate_product_inputs

end module fortfem_spherical_helmholtz_dtn
