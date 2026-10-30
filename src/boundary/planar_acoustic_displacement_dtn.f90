module fortfem_planar_acoustic_displacement_dtn
    !! Planar acoustic displacement-to-pressure map: a Fourier multiplier on
    !! equispaced periodic samples (generated mode multiplier), its circulant
    !! Galerkin form with the periodic P1 mass, and the form's JVP and VJP in the
    !! angular frequency, the complex sound speed and the fluid density.
    !!
    !! With period_sample_count = M > sample_count = m, the period holds M
    !! samples of which the first m are the interface and the other M - m lie on
    !! a hard baffle (zero normal displacement): the displacement is
    !! zero-extended over the period, and the form is the leading m x m block
    !! of the M x M circulant form, a Toeplitz matrix. This is the classic
    !! baffled DtN; it couples the interface ends only through the fluid.
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
            period, maximum_mode, form, status, period_sample_count)
        integer, intent(in) :: sample_count, maximum_mode
        real(dp), intent(in) :: angular_frequency, sound_speed
        real(dp), intent(in) :: fluid_density, period
        complex(dp), intent(out) :: form(:, :)
        integer, intent(out) :: status
        integer, intent(in), optional :: period_sample_count

        call assemble_planar_acoustic_dtn_form_complex( &
            sample_count, angular_frequency, &
            cmplx(sound_speed, 0.0_dp, dp), fluid_density, period, &
            maximum_mode, form, status, period_sample_count)
    end subroutine assemble_planar_acoustic_dtn_form_real

    subroutine assemble_planar_acoustic_dtn_form_complex( &
            sample_count, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, form, status, period_sample_count)
        !! period and maximum_mode refer to the full period of
        !! period_sample_count samples (default sample_count).
        integer, intent(in) :: sample_count, maximum_mode
        real(dp), intent(in) :: angular_frequency
        complex(dp), intent(in) :: sound_speed
        real(dp), intent(in) :: fluid_density, period
        complex(dp), intent(out) :: form(:, :)
        integer, intent(out) :: status
        integer, intent(in), optional :: period_sample_count

        complex(dp), allocatable :: displacement(:), pressure(:)
        integer :: operator_status, total

        form = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        total = total_samples(sample_count, period_sample_count)
        if (sample_count < 3 .or. total < sample_count) return
        if (size(form, 1) /= sample_count .or. &
            size(form, 2) /= sample_count) return
        allocate(displacement(total), pressure(total))
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
            fluid_density_dot, form_dot, status, period_sample_count)
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
        integer, intent(in), optional :: period_sample_count

        complex(dp), allocatable :: spectrum(:)
        integer :: bin, total
        real(dp) :: mode

        form_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        total = total_samples(sample_count, period_sample_count)
        if (.not. valid_form_inputs(sample_count, total, angular_frequency, &
            sound_speed, fluid_density, period, maximum_mode, form_dot)) return
        ! The impulse e_1 has a flat spectrum, so the column's spectrum is the
        ! multiplier itself, and its derivative is the multiplier's.
        allocate(spectrum(total))
        spectrum = cmplx(0.0_dp, 0.0_dp, dp)
        do bin = 0, total - 1
            if (.not. retained_mode(bin, total, maximum_mode, mode)) cycle
            if (at_cutoff(angular_frequency, sound_speed, period, mode)) return
            call generated_planar_acoustic_dtn_multiplier_jvp(angular_frequency, &
                sound_speed, fluid_density, period, mode, angular_frequency_dot, &
                sound_speed_dot, fluid_density_dot, 0.0_dp, spectrum(bin + 1))
        end do
        call fft_c2c(spectrum, 1)
        call circulant_form(spectrum/real(total, dp), period, form_dot)
        status = 0
    end subroutine assemble_planar_acoustic_displacement_dtn_form_jvp

    subroutine assemble_planar_acoustic_displacement_dtn_form_vjp( &
            sample_count, angular_frequency, sound_speed, fluid_density, &
            period, maximum_mode, form_bar, angular_frequency_bar, &
            sound_speed_bar, fluid_density_bar, status, period_sample_count)
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
        integer, intent(in), optional :: period_sample_count

        complex(dp), allocatable :: column_bar(:), spectrum_bar(:)
        complex(dp) :: omega_bar, c_bar, rho_bar, period_bar
        integer :: a, b, bin, row, previous, next, total
        real(dp) :: mode

        angular_frequency_bar = 0.0_dp
        sound_speed_bar = cmplx(0.0_dp, 0.0_dp, dp)
        fluid_density_bar = 0.0_dp
        status = 1
        total = total_samples(sample_count, period_sample_count)
        if (.not. valid_form_inputs(sample_count, total, angular_frequency, &
            sound_speed, fluid_density, period, maximum_mode, form_bar)) return
        ! Adjoint of the layout form(a, b) = column(a - b mod total).
        allocate(column_bar(total), spectrum_bar(total))
        column_bar = cmplx(0.0_dp, 0.0_dp, dp)
        do b = 1, sample_count
            do a = 1, sample_count
                row = modulo(a - b, total) + 1
                column_bar(row) = column_bar(row) + form_bar(a, b)
            end do
        end do
        ! The periodic P1 mass is a real symmetric stencil, so it is self-adjoint.
        do row = 1, total
            previous = modulo(row - 2, total) + 1
            next = modulo(row, total) + 1
            call generated_periodic_p1_mass(period, real(total, dp), &
                column_bar(previous), column_bar(row), column_bar(next), &
                spectrum_bar(row))
        end do
        ! The adjoint of the unnormalized inverse DFT is the forward DFT.
        call fft_c2c(spectrum_bar, -1)
        spectrum_bar = spectrum_bar/real(total, dp)
        do bin = 0, total - 1
            if (.not. retained_mode(bin, total, maximum_mode, mode)) cycle
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
        !! Leading block of the circulant form whose first column is the
        !! periodic P1 mass of pressure: form(a, b) = column(a - b mod total).
        complex(dp), intent(in) :: pressure(:)
        real(dp), intent(in) :: period
        complex(dp), intent(out) :: form(:, :)

        complex(dp), allocatable :: column(:)
        integer :: a, b, next, previous, row, total

        total = size(pressure)
        allocate(column(total))
        do row = 1, total
            previous = modulo(row - 2, total) + 1
            next = modulo(row, total) + 1
            call generated_periodic_p1_mass(period, real(total, dp), &
                pressure(previous), pressure(row), pressure(next), column(row))
        end do
        do b = 1, size(form, 2)
            do a = 1, size(form, 1)
                form(a, b) = column(modulo(a - b, total) + 1)
            end do
        end do
    end subroutine circulant_form

    pure integer function total_samples(sample_count, period_sample_count)
        integer, intent(in) :: sample_count
        integer, intent(in), optional :: period_sample_count

        total_samples = sample_count
        if (present(period_sample_count)) total_samples = period_sample_count
    end function total_samples

    logical function valid_form_inputs(sample_count, total, angular_frequency, &
            sound_speed, fluid_density, period, maximum_mode, form) result(valid)
        integer, intent(in) :: sample_count, total, maximum_mode
        real(dp), intent(in) :: angular_frequency
        complex(dp), intent(in) :: sound_speed
        real(dp), intent(in) :: fluid_density, period
        complex(dp), intent(in) :: form(:, :)

        valid = .false.
        if (sample_count < 3 .or. total < sample_count) return
        if (size(form, 1) /= sample_count .or. &
            size(form, 2) /= sample_count) return
        if (angular_frequency <= 0.0_dp .or. abs(sound_speed) <= 0.0_dp .or. &
            fluid_density <= 0.0_dp .or. period <= 0.0_dp) return
        if (maximum_mode < 0 .or. maximum_mode > total/2) return
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
