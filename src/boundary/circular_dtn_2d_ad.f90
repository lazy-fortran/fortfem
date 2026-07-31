module fortfem_circular_dtn_2d_ad
    !! Analytical products for the circular exterior Helmholtz DtN map.
    !!
    !! The modal eigenvalue is
    !!
    !!   lambda_n = k H_n^(1)'(k R) / H_n^(1)(k R).
    !!
    !! Its parameter products use the cylindrical Bessel equation (DLMF
    !! 10.2.1), avoiding finite differences of the special-function kernel.
    use fortfem_circular_dtn_2d, only: circular_helmholtz_dtn_eigenvalue
    use fortfem_kinds, only: dp
    use fortnum_fft, only: fft_c2c, fft_c2c_vjp
    use fortnum_special_complex_bessel, only: hankel_h1_real
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    public :: apply_circular_helmholtz_dtn_jvp
    public :: apply_circular_helmholtz_dtn_vjp
    public :: circular_helmholtz_dtn_eigenvalue_jvp
    public :: circular_helmholtz_dtn_eigenvalue_vjp

contains

    subroutine circular_helmholtz_dtn_eigenvalue_jvp( &
            mode, wavenumber, radius, wavenumber_dot, radius_dot, &
            eigenvalue_dot, status)
        integer, intent(in) :: mode
        real(dp), intent(in) :: wavenumber, radius
        real(dp), intent(in) :: wavenumber_dot, radius_dot
        complex(dp), intent(out) :: eigenvalue_dot
        integer, intent(out) :: status

        complex(dp) :: ratio, ratio_dot_argument
        real(dp) :: argument, argument_dot
        integer :: order

        eigenvalue_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call circular_helmholtz_dtn_eigenvalue( &
            mode, wavenumber, radius, ratio, status)
        if (status /= 0) return
        order = abs(mode)
        argument = wavenumber*radius
        call hankel_log_derivative(order, argument, ratio_dot_argument, status)
        if (status /= 0) return
        argument_dot = wavenumber_dot*radius + wavenumber*radius_dot
        eigenvalue_dot = wavenumber_dot*(ratio/wavenumber) + &
            wavenumber*ratio_dot_argument*argument_dot
    end subroutine circular_helmholtz_dtn_eigenvalue_jvp

    subroutine circular_helmholtz_dtn_eigenvalue_vjp( &
            mode, wavenumber, radius, eigenvalue_bar, wavenumber_bar, &
            radius_bar, status)
        integer, intent(in) :: mode
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(in) :: eigenvalue_bar
        real(dp), intent(out) :: wavenumber_bar, radius_bar
        integer, intent(out) :: status

        complex(dp) :: eigenvalue_dot

        wavenumber_bar = 0.0_dp
        radius_bar = 0.0_dp
        call circular_helmholtz_dtn_eigenvalue_jvp( &
            mode, wavenumber, radius, 1.0_dp, 0.0_dp, eigenvalue_dot, status)
        if (status /= 0) return
        wavenumber_bar = real(conjg(eigenvalue_bar)*eigenvalue_dot, dp)
        call circular_helmholtz_dtn_eigenvalue_jvp( &
            mode, wavenumber, radius, 0.0_dp, 1.0_dp, eigenvalue_dot, status)
        if (status /= 0) then
            wavenumber_bar = 0.0_dp
            return
        end if
        radius_bar = real(conjg(eigenvalue_bar)*eigenvalue_dot, dp)
    end subroutine circular_helmholtz_dtn_eigenvalue_vjp

    subroutine apply_circular_helmholtz_dtn_jvp( &
            trace, wavenumber, radius, maximum_mode, trace_dot, &
            wavenumber_dot, radius_dot, normal_derivative_dot, &
            discarded_relative_norm_dot, status)
        complex(dp), intent(in) :: trace(:), trace_dot(:)
        real(dp), intent(in) :: wavenumber, radius
        integer, intent(in) :: maximum_mode
        real(dp), intent(in) :: wavenumber_dot, radius_dot
        complex(dp), intent(out) :: normal_derivative_dot(:)
        real(dp), intent(out) :: discarded_relative_norm_dot
        integer, intent(out) :: status

        complex(dp), allocatable :: spectrum(:), spectrum_dot(:)
        complex(dp) :: eigenvalue, eigenvalue_dot
        real(dp) :: discarded_energy, discarded_energy_dot
        real(dp) :: total_energy, total_energy_dot
        real(dp) :: discarded_relative_norm
        integer :: bin, mode, mode_status, point_count

        normal_derivative_dot = cmplx(0.0_dp, 0.0_dp, dp)
        discarded_relative_norm_dot = 0.0_dp
        status = -1
        point_count = size(trace)
        if (point_count < 1 .or. size(trace_dot) /= point_count .or. &
            size(normal_derivative_dot) /= point_count) return
        if (wavenumber <= 0.0_dp .or. radius <= 0.0_dp) return
        if (maximum_mode < 0) return

        allocate(spectrum(point_count), spectrum_dot(point_count))
        spectrum = trace
        spectrum_dot = trace_dot
        call fft_c2c(spectrum, -1)
        call fft_c2c(spectrum_dot, -1)
        total_energy = sum(abs(spectrum)**2)
        total_energy_dot = 2.0_dp*real( &
            sum(conjg(spectrum)*spectrum_dot), dp)
        discarded_energy = 0.0_dp
        discarded_energy_dot = 0.0_dp

        do bin = 0, point_count - 1
            mode = signed_mode(bin, point_count)
            if (abs(mode) > maximum_mode) then
                discarded_energy = discarded_energy + abs(spectrum(bin + 1))**2
                discarded_energy_dot = discarded_energy_dot + 2.0_dp*real( &
                    conjg(spectrum(bin + 1))*spectrum_dot(bin + 1), dp)
                spectrum_dot(bin + 1) = cmplx(0.0_dp, 0.0_dp, dp)
                cycle
            end if
            call circular_helmholtz_dtn_eigenvalue( &
                mode, wavenumber, radius, eigenvalue, mode_status)
            if (mode_status /= 0) then
                status = -2
                return
            end if
            call circular_helmholtz_dtn_eigenvalue_jvp( &
                mode, wavenumber, radius, wavenumber_dot, radius_dot, &
                eigenvalue_dot, mode_status)
            if (mode_status /= 0) then
                status = -2
                return
            end if
            spectrum_dot(bin + 1) = eigenvalue*spectrum_dot(bin + 1) + &
                eigenvalue_dot*spectrum(bin + 1)
            spectrum(bin + 1) = eigenvalue*spectrum(bin + 1)
        end do

        call fft_c2c(spectrum_dot, 1)
        normal_derivative_dot = &
            spectrum_dot/real(point_count, dp)
        if (total_energy > tiny(1.0_dp)) then
            discarded_relative_norm = sqrt(discarded_energy/total_energy)
            if (discarded_relative_norm > tiny(1.0_dp)) then
                discarded_relative_norm_dot = 0.5_dp*discarded_relative_norm*( &
                    discarded_energy_dot/discarded_energy - &
                    total_energy_dot/total_energy)
            end if
        end if
        status = 0
    end subroutine apply_circular_helmholtz_dtn_jvp

    subroutine apply_circular_helmholtz_dtn_vjp( &
            trace, wavenumber, radius, maximum_mode, normal_derivative, &
            normal_derivative_bar, discarded_relative_norm, discarded_bar, &
            trace_bar, wavenumber_bar, radius_bar, status)
        complex(dp), intent(in) :: trace(:), normal_derivative(:)
        complex(dp), intent(in) :: normal_derivative_bar(:)
        real(dp), intent(in) :: wavenumber, radius
        integer, intent(in) :: maximum_mode
        real(dp), intent(in) :: discarded_relative_norm, discarded_bar
        complex(dp), intent(out) :: trace_bar(:)
        real(dp), intent(out) :: wavenumber_bar, radius_bar
        integer, intent(out) :: status

        complex(dp), allocatable :: spectrum(:), spectrum_bar(:)
        complex(dp), allocatable :: normal_spectrum_bar(:)
        complex(dp) :: eigenvalue, eigenvalue_bar
        real(dp) :: discarded_energy, total_energy
        real(dp) :: alpha_discarded, alpha_total, discarded_relative_norm_check
        real(dp) :: local_wavenumber_bar, local_radius_bar
        integer :: bin, mode, mode_status, point_count

        trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wavenumber_bar = 0.0_dp
        radius_bar = 0.0_dp
        status = -1
        point_count = size(trace)
        if (point_count < 1 .or. size(normal_derivative) /= point_count .or. &
            size(normal_derivative_bar) /= point_count .or. &
            size(trace_bar) /= point_count) return
        if (wavenumber <= 0.0_dp .or. radius <= 0.0_dp) return
        if (maximum_mode < 0 .or. discarded_relative_norm < 0.0_dp) return

        allocate( &
            spectrum(point_count), spectrum_bar(point_count), &
            normal_spectrum_bar(point_count))
        spectrum = trace
        call fft_c2c(spectrum, -1)
        total_energy = sum(abs(spectrum)**2)
        discarded_energy = 0.0_dp
        do bin = 0, point_count - 1
            mode = signed_mode(bin, point_count)
            if (abs(mode) > maximum_mode) then
                discarded_energy = discarded_energy + abs(spectrum(bin + 1))**2
            end if
        end do

        spectrum_bar = normal_derivative_bar/real(point_count, dp)
        call fft_c2c_vjp(spectrum_bar, 1)
        normal_spectrum_bar = spectrum_bar
        spectrum_bar = cmplx(0.0_dp, 0.0_dp, dp)
        do bin = 0, point_count - 1
            mode = signed_mode(bin, point_count)
            if (abs(mode) > maximum_mode) cycle
            call circular_helmholtz_dtn_eigenvalue( &
                mode, wavenumber, radius, eigenvalue, mode_status)
            if (mode_status /= 0) then
                status = -2
                return
            end if
            eigenvalue_bar = &
                conjg(spectrum(bin + 1))*normal_spectrum_bar(bin + 1)
            spectrum_bar(bin + 1) = &
                conjg(eigenvalue)*normal_spectrum_bar(bin + 1)
            call circular_helmholtz_dtn_eigenvalue_vjp( &
                mode, wavenumber, radius, eigenvalue_bar, &
                local_wavenumber_bar, local_radius_bar, mode_status)
            if (mode_status /= 0) then
                status = -2
                return
            end if
            wavenumber_bar = wavenumber_bar + local_wavenumber_bar
            radius_bar = radius_bar + local_radius_bar
        end do
        if (total_energy > tiny(1.0_dp) .and. discarded_bar /= 0.0_dp) then
            discarded_relative_norm_check = sqrt(discarded_energy/total_energy)
            if (discarded_relative_norm_check > tiny(1.0_dp)) then
                alpha_discarded = discarded_bar/( &
                    2.0_dp*discarded_relative_norm_check*total_energy)
                alpha_total = -discarded_bar*discarded_energy/( &
                    2.0_dp*discarded_relative_norm_check*total_energy**2)
                do bin = 0, point_count - 1
                    mode = signed_mode(bin, point_count)
                    spectrum_bar(bin + 1) = spectrum_bar(bin + 1) + &
                        2.0_dp*alpha_total*spectrum(bin + 1)
                    if (abs(mode) > maximum_mode) then
                        spectrum_bar(bin + 1) = spectrum_bar(bin + 1) + &
                            2.0_dp*alpha_discarded*spectrum(bin + 1)
                    end if
                end do
            end if
        end if

        call fft_c2c_vjp(spectrum_bar, -1)
        trace_bar = spectrum_bar
        status = 0
    end subroutine apply_circular_helmholtz_dtn_vjp

    subroutine hankel_log_derivative( &
            order, argument, ratio_derivative, status)
        integer, intent(in) :: order
        real(dp), intent(in) :: argument
        complex(dp), intent(out) :: ratio_derivative
        integer, intent(out) :: status

        complex(dp) :: derivative, hankel, hankel_minus, hankel_plus
        type(fortnum_status_t) :: special_status

        ratio_derivative = cmplx(0.0_dp, 0.0_dp, dp)
        status = -1
        if (order < 0 .or. argument <= 0.0_dp) return
        call hankel_h1_real(order, argument, hankel, special_status)
        if (.not. status_ok(special_status)) return
        if (order == 0) then
            call hankel_h1_real(1, argument, hankel_plus, special_status)
            derivative = -hankel_plus
        else
            call hankel_h1_real(order - 1, argument, hankel_minus, special_status)
            call hankel_h1_real(order + 1, argument, hankel_plus, special_status)
            derivative = 0.5_dp*(hankel_minus - hankel_plus)
        end if
        if (.not. status_ok(special_status)) return
        ratio_derivative = -derivative/hankel/argument - 1.0_dp + &
            real(order*order, dp)/argument**2 - (derivative/hankel)**2
        status = 0
    end subroutine hankel_log_derivative

    pure integer function signed_mode(bin, point_count) result(mode)
        integer, intent(in) :: bin, point_count

        if (bin <= point_count/2) then
            mode = bin
        else
            mode = bin - point_count
        end if
    end function signed_mode

end module fortfem_circular_dtn_2d_ad
