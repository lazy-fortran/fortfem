module fortfem_planar_acoustic_displacement_dtn
    use fortfem_kinds, only: dp
    use fortnum_fft, only: fft_c2c
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: apply_planar_acoustic_displacement_dtn
    public :: assemble_planar_acoustic_displacement_dtn_form

contains

    subroutine assemble_planar_acoustic_displacement_dtn_form( &
            sample_count, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, form, status)
        integer, intent(in) :: sample_count, maximum_mode
        real(dp), intent(in) :: angular_frequency, sound_speed
        real(dp), intent(in) :: fluid_density, period
        complex(dp), intent(out) :: form(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: displacement(:), pressure(:)
        complex(dp), allocatable :: strong_operator(:, :)
        real(dp) :: spacing
        integer :: column, next, operator_status, previous, row

        form = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (sample_count < 3) return
        if (size(form, 1) /= sample_count .or. &
            size(form, 2) /= sample_count) return
        allocate(displacement(sample_count), pressure(sample_count))
        allocate(strong_operator(sample_count, sample_count))
        do column = 1, sample_count
            displacement = cmplx(0.0_dp, 0.0_dp, dp)
            displacement(column) = cmplx(1.0_dp, 0.0_dp, dp)
            call apply_planar_acoustic_displacement_dtn( &
                displacement, angular_frequency, sound_speed, fluid_density, &
                period, maximum_mode, pressure, operator_status)
            if (operator_status /= 0) return
            strong_operator(:, column) = pressure
        end do
        spacing = period/real(sample_count, dp)
        do row = 1, sample_count
            previous = modulo(row - 2, sample_count) + 1
            next = modulo(row, sample_count) + 1
            form(row, :) = spacing/6.0_dp*( &
                strong_operator(previous, :) + &
                4.0_dp*strong_operator(row, :) + &
                strong_operator(next, :))
        end do
        status = 0
    end subroutine assemble_planar_acoustic_displacement_dtn_form

    subroutine apply_planar_acoustic_displacement_dtn( &
            displacement, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, pressure, status)
        complex(dp), intent(in) :: displacement(:)
        real(dp), intent(in) :: angular_frequency, sound_speed
        real(dp), intent(in) :: fluid_density, period
        integer, intent(in) :: maximum_mode
        complex(dp), intent(out) :: pressure(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: spectrum(:)
        complex(dp) :: normal_wavenumber
        real(dp) :: acoustic_wavenumber, radicand, tangential_wavenumber
        integer :: bin, mode, sample_count

        pressure = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        sample_count = size(displacement)
        if (sample_count < 1 .or. size(pressure) /= sample_count) return
        if (angular_frequency <= 0.0_dp .or. sound_speed <= 0.0_dp .or. &
            fluid_density <= 0.0_dp .or. period <= 0.0_dp) return
        if (maximum_mode < 0 .or. maximum_mode > sample_count/2) return

        allocate(spectrum(sample_count))
        spectrum = displacement
        call fft_c2c(spectrum, -1)
        acoustic_wavenumber = angular_frequency/sound_speed
        do bin = 0, sample_count - 1
            if (bin <= sample_count/2) then
                mode = bin
            else
                mode = bin - sample_count
            end if
            if (abs(mode) > maximum_mode) then
                spectrum(bin + 1) = cmplx(0.0_dp, 0.0_dp, dp)
                cycle
            end if
            tangential_wavenumber = 2.0_dp*pi*real(mode, dp)/period
            radicand = acoustic_wavenumber**2 - tangential_wavenumber**2
            if (abs(radicand) <= &
                32.0_dp*epsilon(1.0_dp)*acoustic_wavenumber**2) then
                if (abs(spectrum(bin + 1)) > &
                    64.0_dp*epsilon(1.0_dp)*max( &
                    1.0_dp, maxval(abs(spectrum)))) return
                spectrum(bin + 1) = cmplx(0.0_dp, 0.0_dp, dp)
                cycle
            end if
            if (radicand > 0.0_dp) then
                normal_wavenumber = cmplx(sqrt(radicand), 0.0_dp, dp)
            else
                normal_wavenumber = cmplx(0.0_dp, sqrt(-radicand), dp)
            end if
            spectrum(bin + 1) = &
                -cmplx(0.0_dp, 1.0_dp, dp)*fluid_density* &
                angular_frequency**2/normal_wavenumber*spectrum(bin + 1)
        end do
        call fft_c2c(spectrum, 1)
        pressure = spectrum/real(sample_count, dp)
        status = 0
    end subroutine apply_planar_acoustic_displacement_dtn
end module fortfem_planar_acoustic_displacement_dtn
