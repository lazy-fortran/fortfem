module fortfem_circular_dtn_2d
    use fortfem_kinds, only: dp
    use fortnum_fft, only: fft_c2c
    use fortnum_special_complex_bessel, only: hankel_h1_real
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    public :: apply_circular_helmholtz_dtn
    public :: circular_helmholtz_dtn_eigenvalue

contains

    subroutine apply_circular_helmholtz_dtn( &
            trace, wavenumber, radius, maximum_mode, normal_derivative, &
            discarded_relative_norm, status)
        complex(dp), intent(in) :: trace(:)
        real(dp), intent(in) :: wavenumber, radius
        integer, intent(in) :: maximum_mode
        complex(dp), intent(out) :: normal_derivative(:)
        real(dp), intent(out) :: discarded_relative_norm
        integer, intent(out) :: status

        complex(dp), allocatable :: spectrum(:)
        complex(dp) :: eigenvalue
        real(dp) :: discarded_energy, total_energy
        integer :: bin, mode, mode_status, point_count

        normal_derivative = (0.0_dp, 0.0_dp)
        discarded_relative_norm = 0.0_dp
        status = -1
        point_count = size(trace)
        if (point_count < 1) return
        if (size(normal_derivative) /= point_count) return
        if (wavenumber <= 0.0_dp .or. radius <= 0.0_dp) return
        if (maximum_mode < 0) return

        allocate(spectrum(point_count))
        spectrum = trace
        call fft_c2c(spectrum, -1)
        total_energy = sum(abs(spectrum)**2)
        discarded_energy = 0.0_dp

        do bin = 0, point_count - 1
            if (bin <= point_count / 2) then
                mode = bin
            else
                mode = bin - point_count
            end if
            if (abs(mode) > maximum_mode) then
                discarded_energy = &
                    discarded_energy + abs(spectrum(bin + 1))**2
                spectrum(bin + 1) = (0.0_dp, 0.0_dp)
                cycle
            end if
            call circular_helmholtz_dtn_eigenvalue( &
                mode, wavenumber, radius, eigenvalue, mode_status)
            if (mode_status /= 0) then
                status = -2
                return
            end if
            spectrum(bin + 1) = eigenvalue * spectrum(bin + 1)
        end do

        if (total_energy > 0.0_dp) then
            discarded_relative_norm = sqrt(discarded_energy / total_energy)
        end if
        call fft_c2c(spectrum, 1)
        normal_derivative = spectrum / real(point_count, dp)
        status = 0
    end subroutine apply_circular_helmholtz_dtn

    subroutine circular_helmholtz_dtn_eigenvalue( &
            mode, wavenumber, radius, eigenvalue, status)
        integer, intent(in) :: mode
        real(dp), intent(in) :: wavenumber, radius
        complex(dp), intent(out) :: eigenvalue
        integer, intent(out) :: status

        complex(dp) :: derivative, hankel, hankel_minus, hankel_plus
        type(fortnum_status_t) :: special_status
        real(dp) :: argument
        integer :: order

        eigenvalue = (0.0_dp, 0.0_dp)
        status = -1
        if (wavenumber <= 0.0_dp .or. radius <= 0.0_dp) return

        order = abs(mode)
        argument = wavenumber * radius
        call hankel_h1_real(order, argument, hankel, special_status)
        if (.not. status_ok(special_status)) then
            status = -2
            return
        end if
        if (order == 0) then
            call hankel_h1_real(1, argument, hankel_plus, special_status)
            derivative = -hankel_plus
        else
            call hankel_h1_real( &
                order - 1, argument, hankel_minus, special_status)
            call hankel_h1_real( &
                order + 1, argument, hankel_plus, special_status)
            derivative = 0.5_dp * (hankel_minus - hankel_plus)
        end if
        if (.not. status_ok(special_status)) then
            status = -2
            return
        end if

        eigenvalue = wavenumber * derivative / hankel
        status = 0
    end subroutine circular_helmholtz_dtn_eigenvalue

end module fortfem_circular_dtn_2d
