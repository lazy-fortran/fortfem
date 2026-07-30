program test_curved_acoustic_displacement_form_2d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_curved_acoustic_displacement_ntd_form_2d, &
        circular_helmholtz_dtn_eigenvalue
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: mesh_sizes(3) = [16, 32, 64]
    real(dp) :: errors(size(mesh_sizes))
    integer :: mesh
    logical :: all_passed

    all_passed = .true.
    do mesh = 1, size(mesh_sizes)
        call circle_form_error(mesh_sizes(mesh), errors(mesh))
    end do
    write (*, '(a,3(es12.4,1x))') "curved acoustic form errors: ", errors

    call record_condition( &
        errors(2) < 0.4_dp*errors(1) .and. &
        errors(3) < 0.4_dp*errors(2), &
        "Curved acoustic weak form converges at second order")
    call record_condition(errors(3) < 1.0e-3_dp, &
        "Curved acoustic weak form matches the outgoing Hankel traction")
    call check_summary("Curved acoustic displacement weak form")
    if (.not. all_passed) error stop 1

contains

    subroutine circle_form_error(panel_count, relative_error)
        integer, intent(in) :: panel_count
        real(dp), intent(out) :: relative_error

        real(dp), parameter :: angular_frequency = 2.3_dp
        real(dp), parameter :: density = 1.7_dp
        real(dp), parameter :: radius = 0.8_dp
        real(dp), parameter :: wavenumber = 1.4_dp
        complex(dp), allocatable :: exact_load(:), form(:, :), pressure(:)
        complex(dp), allocatable :: vector_displacement(:)
        integer, allocatable :: panel_nodes(:, :)
        real(dp), allocatable :: panel_end(:, :), panel_start(:, :)
        complex(dp) :: dtn_eigenvalue, pressure_amplitude
        real(dp) :: angle, half_step, length, normal(2), pi
        integer :: component, first, panel, second, status

        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        allocate(panel_nodes(2, panel_count))
        allocate(form(2*panel_count, 2*panel_count))
        allocate(vector_displacement(2*panel_count))
        allocate(exact_load(2*panel_count), pressure(panel_count))
        pi = acos(-1.0_dp)
        half_step = pi/real(panel_count, dp)
        do panel = 1, panel_count
            angle = 2.0_dp*pi*real(panel - 1, dp)/real(panel_count, dp)
            panel_start(1, panel) = radius*cos(angle)
            panel_start(2, panel) = radius*sin(angle)
            panel_nodes(1, panel) = panel
            vector_displacement(2*panel - 1) = cos(angle)* &
                exp(cmplx(0.0_dp, angle, dp))
            vector_displacement(2*panel) = sin(angle)* &
                exp(cmplx(0.0_dp, angle, dp))
            angle = 2.0_dp*pi*real(panel, dp)/real(panel_count, dp)
            panel_end(1, panel) = radius*cos(angle)
            panel_end(2, panel) = radius*sin(angle)
            panel_nodes(2, panel) = modulo(panel, panel_count) + 1
        end do

        call circular_helmholtz_dtn_eigenvalue( &
            1, wavenumber, radius, dtn_eigenvalue, status)
        if (status /= 0) error stop "analytical circular DtN failed"
        pressure_amplitude = density*angular_frequency**2* &
            cos(half_step)**2/dtn_eigenvalue
        do panel = 1, panel_count
            angle = 2.0_dp*pi*real(panel - 1, dp)/real(panel_count, dp)
            pressure(panel) = pressure_amplitude* &
                exp(cmplx(0.0_dp, angle, dp))
        end do

        exact_load = (0.0_dp, 0.0_dp)
        do panel = 1, panel_count
            length = norm2(panel_end(:, panel) - panel_start(:, panel))
            normal = [ &
                panel_end(2, panel) - panel_start(2, panel), &
                panel_start(1, panel) - panel_end(1, panel)]/length
            do first = 1, 2
                do second = 1, 2
                    do component = 1, 2
                        exact_load( &
                            2*panel_nodes(first, panel) - 2 + component) = &
                            exact_load( &
                            2*panel_nodes(first, panel) - 2 + component) + &
                            normal(component)*length* &
                            merge(2.0_dp, 1.0_dp, first == second)/6.0_dp* &
                            pressure(panel_nodes(second, panel))
                    end do
                end do
            end do
        end do

        call assemble_curved_acoustic_displacement_ntd_form_2d( &
            panel_start, panel_end, panel_nodes, wavenumber, density, &
            angular_frequency, 16, form, status)
        call record_condition(status == 0, &
            "Curved acoustic displacement weak form assembles")
        relative_error = sqrt(sum(abs( &
            matmul(form, vector_displacement) - exact_load)**2)/ &
            sum(abs(exact_load)**2))
    end subroutine circle_form_error

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_curved_acoustic_displacement_form_2d
