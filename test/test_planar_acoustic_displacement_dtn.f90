program test_planar_acoustic_displacement_dtn
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: apply_planar_acoustic_displacement_dtn, &
        assemble_planar_acoustic_displacement_dtn_form
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    logical :: all_passed

    all_passed = .true.
    call check_mode(0, cmplx(0.0_dp, -6.0_dp, dp), &
        "Uniform displacement produces outgoing plane-wave pressure")
    call check_mode(2, cmplx(0.0_dp, -18.0_dp/sqrt(5.0_dp), dp), &
        "Propagating oblique displacement uses the outgoing branch")
    call check_mode(4, cmplx(-18.0_dp/sqrt(7.0_dp), 0.0_dp, dp), &
        "Evanescent displacement produces a real decaying pressure")
    call check_truncation()
    call check_grazing_rejection()
    call check_complex_sound_speed()
    call check_weak_form()

    call check_summary("Planar acoustic displacement DtN")
    if (.not. all_passed) error stop 1

contains

    subroutine check_mode(mode, multiplier, description)
        integer, intent(in) :: mode
        complex(dp), intent(in) :: multiplier
        character(len=*), intent(in) :: description

        integer, parameter :: sample_count = 11
        complex(dp) :: displacement(sample_count), pressure(sample_count)
        real(dp) :: coordinate
        integer :: point, status

        do point = 1, sample_count
            coordinate = 2.0_dp*pi*real(point - 1, dp)/ &
                real(sample_count, dp)
            displacement(point) = exp(cmplx( &
                0.0_dp, real(mode, dp)*coordinate, dp))
        end do
        call apply_planar_acoustic_displacement_dtn( &
            displacement, 3.0_dp, 1.0_dp, 2.0_dp, 2.0_dp*pi, 5, &
            pressure, status)
        call record_condition(status == 0 .and. &
            maxval(abs(pressure - multiplier*displacement)) < 2.0e-12_dp, &
            description)
    end subroutine check_mode

    subroutine check_truncation()
        integer, parameter :: sample_count = 11
        complex(dp) :: displacement(sample_count), pressure(sample_count)
        real(dp) :: coordinate
        integer :: point, status

        do point = 1, sample_count
            coordinate = 2.0_dp*pi*real(point - 1, dp)/ &
                real(sample_count, dp)
            displacement(point) = 1.0_dp + exp(cmplx(0.0_dp, &
                4.0_dp*coordinate, dp))
        end do
        call apply_planar_acoustic_displacement_dtn( &
            displacement, 3.0_dp, 1.0_dp, 2.0_dp, 2.0_dp*pi, 2, &
            pressure, status)
        call record_condition(status == 0 .and. &
            maxval(abs(pressure - cmplx(0.0_dp, -6.0_dp, dp))) < &
            2.0e-12_dp, "Acoustic DtN discards modes above its cutoff")
    end subroutine check_truncation

    subroutine check_grazing_rejection()
        complex(dp) :: displacement(7), pressure(7)
        real(dp) :: coordinate
        integer :: point, status

        do point = 1, 7
            coordinate = 2.0_dp*pi*real(point - 1, dp)/7.0_dp
            displacement(point) = exp(cmplx(0.0_dp, &
                3.0_dp*coordinate, dp))
        end do
        call apply_planar_acoustic_displacement_dtn( &
            displacement, 3.0_dp, 1.0_dp, 2.0_dp, 2.0_dp*pi, 3, &
            pressure, status)
        call record_condition(status /= 0, &
            "Acoustic DtN reports a singular grazing mode")
    end subroutine check_grazing_rejection

    subroutine check_complex_sound_speed()
        complex(dp) :: displacement(7), pressure(7)
        integer :: status

        displacement = cmplx(1.0_dp, 0.0_dp, dp)
        call apply_planar_acoustic_displacement_dtn( &
            displacement, 3.0_dp, cmplx(1.0_dp, 0.1_dp, dp), 2.0_dp, &
            2.0_dp*pi, 3, pressure, status)
        call record_condition(status == 0 .and. maxval(abs(pressure - &
            cmplx(0.6_dp, -6.0_dp, dp))) < 2.0e-12_dp, &
            "Acoustic DtN supports a complex lossy sound speed")
    end subroutine check_complex_sound_speed

    subroutine check_weak_form()
        integer, parameter :: sample_count = 11
        complex(dp) :: form(sample_count, sample_count)
        complex(dp) :: mode_values(sample_count), result(sample_count)
        complex(dp) :: pressure_eigenvalue
        real(dp) :: coordinate, mass_eigenvalue
        integer :: mode, point, status

        do mode = 0, 4, 2
            do point = 1, sample_count
                coordinate = 2.0_dp*pi*real(point - 1, dp)/ &
                    real(sample_count, dp)
                mode_values(point) = exp(cmplx( &
                    0.0_dp, real(mode, dp)*coordinate, dp))
            end do
            call assemble_planar_acoustic_displacement_dtn_form( &
                sample_count, 3.5_dp, 1.0_dp, 2.0_dp, 2.0_dp*pi, 5, &
                form, status)
            if (mode == 0) then
                pressure_eigenvalue = cmplx(0.0_dp, -7.0_dp, dp)
            else if (mode == 2) then
                pressure_eigenvalue = &
                    cmplx(0.0_dp, -24.5_dp/sqrt(8.25_dp), dp)
            else
                pressure_eigenvalue = &
                    cmplx(-24.5_dp/sqrt(3.75_dp), 0.0_dp, dp)
            end if
            mass_eigenvalue = 2.0_dp*pi/real(sample_count, dp)/3.0_dp*( &
                2.0_dp + cos(2.0_dp*pi*real(mode, dp)/ &
                real(sample_count, dp)))
            result = matmul(form, mode_values)
            call record_condition(status == 0 .and. maxval(abs(result - &
                mass_eigenvalue*pressure_eigenvalue*mode_values)) < &
                2.0e-12_dp, &
                "Acoustic displacement DtN form preserves Fourier spectrum")
        end do
    end subroutine check_weak_form

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition
end program test_planar_acoustic_displacement_dtn
