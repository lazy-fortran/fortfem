program test_curved_acoustic_displacement_ntd_2d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_curved_acoustic_displacement_ntd_2d, &
        circular_helmholtz_dtn_eigenvalue
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: mesh_sizes(3) = [24, 48, 96]
    real(dp) :: errors(size(mesh_sizes))
    integer :: mesh
    logical :: all_passed

    all_passed = .true.
    do mesh = 1, size(mesh_sizes)
        call circle_mode_error(mesh_sizes(mesh), errors(mesh))
    end do
    write (*, '(a,3(es12.4,1x))') "curved acoustic NtD errors: ", errors

    call record_condition( &
        errors(2) < 0.35_dp*errors(1) .and. &
        errors(3) < 0.35_dp*errors(2), &
        "Curved acoustic NtD circle mode converges at second order")
    call record_condition(errors(3) < 2.0e-3_dp, &
        "Curved acoustic NtD matches the outgoing Hankel mode")
    call check_summary("Curved acoustic displacement NtD")
    if (.not. all_passed) error stop 1

contains

    subroutine circle_mode_error(panel_count, relative_error)
        integer, intent(in) :: panel_count
        real(dp), intent(out) :: relative_error

        real(dp), parameter :: density = 1.7_dp
        real(dp), parameter :: angular_frequency = 2.3_dp
        real(dp), parameter :: radius = 0.8_dp
        real(dp), parameter :: wavenumber = 1.4_dp
        complex(dp), allocatable :: displacement(:), exact(:), pressure(:)
        integer, allocatable :: panel_nodes(:, :)
        real(dp), allocatable :: panel_end(:, :), panel_start(:, :)
        complex(dp) :: dtn_eigenvalue, pressure_amplitude
        real(dp) :: angle, pi
        integer :: panel, status

        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        allocate(panel_nodes(2, panel_count))
        allocate(displacement(panel_count), exact(panel_count))
        allocate(pressure(panel_count))
        pi = acos(-1.0_dp)
        do panel = 1, panel_count
            angle = 2.0_dp*pi*real(panel - 1, dp)/real(panel_count, dp)
            panel_start(:, panel) = radius*[cos(angle), sin(angle)]
            angle = 2.0_dp*pi*real(panel, dp)/real(panel_count, dp)
            panel_end(:, panel) = radius*[cos(angle), sin(angle)]
            panel_nodes(:, panel) = [panel, modulo(panel, panel_count) + 1]
            angle = 2.0_dp*pi*(real(panel, dp) - 0.5_dp)/ &
                real(panel_count, dp)
            displacement(panel) = exp(cmplx(0.0_dp, angle, dp))
        end do

        call circular_helmholtz_dtn_eigenvalue( &
            1, wavenumber, radius, dtn_eigenvalue, status)
        if (status /= 0) error stop "analytical circular DtN failed"
        pressure_amplitude = density*angular_frequency**2/dtn_eigenvalue
        do panel = 1, panel_count
            angle = 2.0_dp*pi*real(panel - 1, dp)/real(panel_count, dp)
            exact(panel) = pressure_amplitude* &
                exp(cmplx(0.0_dp, angle, dp))
        end do

        call apply_curved_acoustic_displacement_ntd_2d( &
            panel_start, panel_end, panel_nodes, wavenumber, density, &
            angular_frequency, 16, displacement, pressure, status)
        call record_condition(status == 0, &
            "Curved acoustic displacement NtD solve succeeds")
        relative_error = sqrt(sum(abs(pressure - exact)**2)/sum(abs(exact)**2))
    end subroutine circle_mode_error

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_curved_acoustic_displacement_ntd_2d
