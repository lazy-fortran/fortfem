module fortfem_planar_acoustic_displacement_dtn
    !! Planar acoustic displacement-to-pressure map: a Fourier multiplier on
    !! equispaced periodic samples (generated mode multiplier), its circulant
    !! Galerkin form with the periodic P1 mass, and the form's JVP and VJP in the
    !! angular frequency, the complex sound speed and the fluid density.
    use fortfem_generated_periodic_p1_mass, only: generated_periodic_p1_mass
    use fortfem_generated_planar_acoustic_dtn_multiplier, only: &
        generated_planar_acoustic_dtn_multiplier
    use fortfem_generated_planar_acoustic_dtn_multiplier_jvp, only: &
        generated_planar_acoustic_dtn_multiplier_jvp
    use fortfem_generated_planar_acoustic_dtn_multiplier_vjp, only: &
        generated_planar_acoustic_dtn_multiplier_vjp
    use fortfem_generated_planar_acoustic_dtn_radicand, only: &
        generated_planar_acoustic_dtn_radicand
    use fortfem_kinds, only: dp
    use fortnum_fft, only: fft_c2c
    implicit none
    private

    public :: apply_planar_acoustic_displacement_dtn
    public :: assemble_planar_acoustic_displacement_dtn_form
    public :: assemble_planar_acoustic_displacement_dtn_form_jvp
    public :: assemble_planar_acoustic_displacement_dtn_form_vjp

    interface apply_planar_acoustic_displacement_dtn
        module procedure apply_planar_acoustic_displacement_dtn_real
        module procedure apply_planar_acoustic_displacement_dtn_complex
    end interface

    interface assemble_planar_acoustic_displacement_dtn_form
        module procedure assemble_planar_acoustic_dtn_form_real
        module procedure assemble_planar_acoustic_dtn_form_complex
    end interface

contains

    subroutine assemble_planar_acoustic_dtn_form_real( &
            sample_count, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, form, status)
        integer, intent(in) :: sample_count, maximum_mode
        real(dp), intent(in) :: angular_frequency, sound_speed
        real(dp), intent(in) :: fluid_density, period
        complex(dp), intent(out) :: form(:, :)
        integer, intent(out) :: status

        call assemble_planar_acoustic_dtn_form_complex( &
            sample_count, angular_frequency, &
            cmplx(sound_speed, 0.0_dp, dp), fluid_density, period, &
            maximum_mode, form, status)
    end subroutine assemble_planar_acoustic_dtn_form_real

    subroutine assemble_planar_acoustic_dtn_form_complex( &
            sample_count, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, form, status)
        integer, intent(in) :: sample_count, maximum_mode
        real(dp), intent(in) :: angular_frequency
        complex(dp), intent(in) :: sound_speed
        real(dp), intent(in) :: fluid_density, period
        complex(dp), intent(out) :: form(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: displacement(:), pressure(:)
        integer :: operator_status

        form = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (sample_count < 3) return
        if (size(form, 1) /= sample_count .or. &
            size(form, 2) /= sample_count) return
        allocate(displacement(sample_count), pressure(sample_count))
        ! The map is a Fourier multiplier and the periodic P1 mass a circulant
        ! stencil, so the form is circulant: one column of the strong operator
        ! and one mass application determine it.
        displacement = cmplx(0.0_dp, 0.0_dp, dp)
        displacement(1) = cmplx(1.0_dp, 0.0_dp, dp)
        call apply_planar_acoustic_displacement_dtn( &
            displacement, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, pressure, operator_status)
        if (operator_status /= 0) return
        call circulant_form(pressure, period, form)
        status = 0
    end subroutine assemble_planar_acoustic_dtn_form_complex

    subroutine assemble_planar_acoustic_displacement_dtn_form_jvp( &
            sample_count, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, angular_frequency_dot, sound_speed_dot, &
            fluid_density_dot, form_dot, status)
        !! Directional derivative of the complex form; status as for the form.
        integer, intent(in) :: sample_count, maximum_mode
        real(dp), intent(in) :: angular_frequency
        complex(dp), intent(in) :: sound_speed
        real(dp), intent(in) :: fluid_density, period
        real(dp), intent(in) :: angular_frequency_dot
        complex(dp), intent(in) :: sound_speed_dot
        real(dp), intent(in) :: fluid_density_dot
        complex(dp), intent(out) :: form_dot(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: spectrum(:)
        integer :: bin
        real(dp) :: mode

        form_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_form_inputs(sample_count, angular_frequency, &
            sound_speed, fluid_density, period, maximum_mode, form_dot)) return
        ! The impulse e_1 has a flat spectrum, so the column's spectrum is the
        ! multiplier itself, and its derivative is the multiplier's.
        allocate(spectrum(sample_count))
        spectrum = cmplx(0.0_dp, 0.0_dp, dp)
        do bin = 0, sample_count - 1
            if (.not. retained_mode(bin, sample_count, maximum_mode, mode)) cycle
            if (at_cutoff(angular_frequency, sound_speed, period, mode)) return
            call generated_planar_acoustic_dtn_multiplier_jvp(angular_frequency, &
                sound_speed, fluid_density, period, mode, angular_frequency_dot, &
                sound_speed_dot, fluid_density_dot, 0.0_dp, spectrum(bin + 1))
        end do
        call fft_c2c(spectrum, 1)
        call circulant_form(spectrum/real(sample_count, dp), period, form_dot)
        status = 0
    end subroutine assemble_planar_acoustic_displacement_dtn_form_jvp

    subroutine assemble_planar_acoustic_displacement_dtn_form_vjp( &
            sample_count, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, form_bar, angular_frequency_bar, &
            sound_speed_bar, fluid_density_bar, status)
        !! Cotangents for form_bar such that Re<form_bar, form_dot> equals
        !! omega_bar omega_dot + Re(conj(c_bar) c_dot) + rho_bar rho_dot.
        integer, intent(in) :: sample_count, maximum_mode
        real(dp), intent(in) :: angular_frequency
        complex(dp), intent(in) :: sound_speed
        real(dp), intent(in) :: fluid_density, period
        complex(dp), intent(in) :: form_bar(:, :)
        real(dp), intent(out) :: angular_frequency_bar
        complex(dp), intent(out) :: sound_speed_bar
        real(dp), intent(out) :: fluid_density_bar
        integer, intent(out) :: status

        complex(dp), allocatable :: column_bar(:), spectrum_bar(:)
        complex(dp) :: omega_bar, c_bar, rho_bar, period_bar
        integer :: bin, offset, row, previous, next
        real(dp) :: mode

        angular_frequency_bar = 0.0_dp
        sound_speed_bar = cmplx(0.0_dp, 0.0_dp, dp)
        fluid_density_bar = 0.0_dp
        status = 1
        if (.not. valid_form_inputs(sample_count, angular_frequency, &
            sound_speed, fluid_density, period, maximum_mode, form_bar)) return
        ! Adjoint of the circulant layout form(a, o + 1) = column(a - o).
        allocate(column_bar(sample_count), spectrum_bar(sample_count))
        column_bar = cmplx(0.0_dp, 0.0_dp, dp)
        do offset = 0, sample_count - 1
            column_bar = column_bar + cshift(form_bar(:, offset + 1), offset)
        end do
        ! The periodic P1 mass is a real symmetric stencil, so it is self-adjoint.
        do row = 1, sample_count
            previous = modulo(row - 2, sample_count) + 1
            next = modulo(row, sample_count) + 1
            call generated_periodic_p1_mass(period, real(sample_count, dp), &
                column_bar(previous), column_bar(row), column_bar(next), &
                spectrum_bar(row))
        end do
        ! The adjoint of the unnormalized inverse DFT is the forward DFT.
        call fft_c2c(spectrum_bar, -1)
        spectrum_bar = spectrum_bar/real(sample_count, dp)
        do bin = 0, sample_count - 1
            if (.not. retained_mode(bin, sample_count, maximum_mode, mode)) cycle
            if (at_cutoff(angular_frequency, sound_speed, period, mode)) return
            ! The generated VJP is the plain transpose: feed conj(bar) and
            ! conjugate the result back (Re pairing for the real inputs).
            call generated_planar_acoustic_dtn_multiplier_vjp(angular_frequency, &
                sound_speed, fluid_density, period, mode, &
                conjg(spectrum_bar(bin + 1)), omega_bar, c_bar, rho_bar, &
                period_bar)
            angular_frequency_bar = angular_frequency_bar + real(omega_bar, dp)
            sound_speed_bar = sound_speed_bar + conjg(c_bar)
            fluid_density_bar = fluid_density_bar + real(rho_bar, dp)
        end do
        status = 0
    end subroutine assemble_planar_acoustic_displacement_dtn_form_vjp

    subroutine circulant_form(pressure, period, form)
        !! Circulant form whose first column is the periodic P1 mass of pressure.
        complex(dp), intent(in) :: pressure(:)
        real(dp), intent(in) :: period
        complex(dp), intent(out) :: form(:, :)

        complex(dp), allocatable :: column(:)
        integer :: next, offset, previous, row, sample_count

        sample_count = size(pressure)
        allocate(column(sample_count))
        do row = 1, sample_count
            previous = modulo(row - 2, sample_count) + 1
            next = modulo(row, sample_count) + 1
            call generated_periodic_p1_mass(period, real(sample_count, dp), &
                pressure(previous), pressure(row), pressure(next), column(row))
        end do
        do offset = 0, sample_count - 1
            form(:, offset + 1) = cshift(column, -offset)
        end do
    end subroutine circulant_form

    logical function valid_form_inputs(sample_count, angular_frequency, &
            sound_speed, fluid_density, period, maximum_mode, form) result(valid)
        integer, intent(in) :: sample_count, maximum_mode
        real(dp), intent(in) :: angular_frequency
        complex(dp), intent(in) :: sound_speed
        real(dp), intent(in) :: fluid_density, period
        complex(dp), intent(in) :: form(:, :)

        valid = .false.
        if (sample_count < 3) return
        if (size(form, 1) /= sample_count .or. &
            size(form, 2) /= sample_count) return
        if (angular_frequency <= 0.0_dp .or. abs(sound_speed) <= 0.0_dp .or. &
            fluid_density <= 0.0_dp .or. period <= 0.0_dp) return
        if (maximum_mode < 0 .or. maximum_mode > sample_count/2) return
        valid = .true.
    end function valid_form_inputs

    logical function retained_mode(bin, sample_count, maximum_mode, mode)
        !! Signed mode of a DFT bin, and whether it is within maximum_mode.
        integer, intent(in) :: bin, sample_count, maximum_mode
        real(dp), intent(out) :: mode

        if (bin <= sample_count/2) then
            mode = real(bin, dp)
        else
            mode = real(bin - sample_count, dp)
        end if
        retained_mode = abs(mode) <= real(maximum_mode, dp)
    end function retained_mode

    logical function at_cutoff(angular_frequency, sound_speed, period, mode)
        !! The normal wavenumber vanishes to roundoff: the multiplier is singular.
        real(dp), intent(in) :: angular_frequency, period, mode
        complex(dp), intent(in) :: sound_speed

        complex(dp) :: radicand

        call generated_planar_acoustic_dtn_radicand(angular_frequency, &
            sound_speed, period, mode, radicand)
        at_cutoff = abs(radicand) <= &
            32.0_dp*epsilon(1.0_dp)*abs(angular_frequency/sound_speed)**2
    end function at_cutoff

    subroutine apply_planar_acoustic_displacement_dtn_real( &
            displacement, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, pressure, status)
        complex(dp), intent(in) :: displacement(:)
        real(dp), intent(in) :: angular_frequency, sound_speed
        real(dp), intent(in) :: fluid_density, period
        integer, intent(in) :: maximum_mode
        complex(dp), intent(out) :: pressure(:)
        integer, intent(out) :: status

        call apply_planar_acoustic_displacement_dtn_complex( &
            displacement, angular_frequency, &
            cmplx(sound_speed, 0.0_dp, dp), fluid_density, period, &
            maximum_mode, pressure, status)
    end subroutine apply_planar_acoustic_displacement_dtn_real

    subroutine apply_planar_acoustic_displacement_dtn_complex( &
            displacement, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, pressure, status)
        complex(dp), intent(in) :: displacement(:)
        real(dp), intent(in) :: angular_frequency
        complex(dp), intent(in) :: sound_speed
        real(dp), intent(in) :: fluid_density, period
        integer, intent(in) :: maximum_mode
        complex(dp), intent(out) :: pressure(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: spectrum(:)
        complex(dp) :: multiplier
        real(dp) :: mode
        integer :: bin, sample_count

        pressure = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        sample_count = size(displacement)
        if (sample_count < 1 .or. size(pressure) /= sample_count) return
        if (angular_frequency <= 0.0_dp .or. abs(sound_speed) <= 0.0_dp .or. &
            fluid_density <= 0.0_dp .or. period <= 0.0_dp) return
        if (maximum_mode < 0 .or. maximum_mode > sample_count/2) return

        allocate(spectrum(sample_count))
        spectrum = displacement
        call fft_c2c(spectrum, -1)
        do bin = 0, sample_count - 1
            if (.not. retained_mode(bin, sample_count, maximum_mode, mode)) then
                spectrum(bin + 1) = cmplx(0.0_dp, 0.0_dp, dp)
                cycle
            end if
            if (at_cutoff(angular_frequency, sound_speed, period, mode)) then
                if (abs(spectrum(bin + 1)) > &
                    64.0_dp*epsilon(1.0_dp)*max( &
                    1.0_dp, maxval(abs(spectrum)))) return
                spectrum(bin + 1) = cmplx(0.0_dp, 0.0_dp, dp)
                cycle
            end if
            call generated_planar_acoustic_dtn_multiplier(angular_frequency, &
                sound_speed, fluid_density, period, mode, multiplier)
            spectrum(bin + 1) = multiplier*spectrum(bin + 1)
        end do
        call fft_c2c(spectrum, 1)
        pressure = spectrum/real(sample_count, dp)
        status = 0
    end subroutine apply_planar_acoustic_displacement_dtn_complex
end module fortfem_planar_acoustic_displacement_dtn
