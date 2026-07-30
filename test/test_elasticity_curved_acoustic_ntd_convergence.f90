program test_elasticity_curved_acoustic_ntd_convergence
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        circular_helmholtz_dtn_eigenvalue, &
        solve_elasticity_curved_acoustic_ntd_p1
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: mesh_sizes(3) = [12, 24, 48]
    real(dp) :: errors(size(mesh_sizes))
    integer :: mesh
    logical :: all_passed

    all_passed = .true.
    do mesh = 1, size(mesh_sizes)
        call rigid_translation_error(mesh_sizes(mesh), errors(mesh))
    end do
    write (*, '(a,3(es12.4,1x))') "curved elastic-acoustic errors: ", errors
    call record_condition( &
        errors(2) < 0.45_dp*errors(1) .and. &
        errors(3) < 0.45_dp*errors(2), &
        "Curved elastic-acoustic solve converges at second order")
    call record_condition(errors(3) < 2.0e-3_dp, &
        "Curved elastic-acoustic solve reproduces rigid translation")
    call check_summary("Curved elastic-acoustic coupling")
    if (.not. all_passed) error stop 1

contains

    subroutine rigid_translation_error(panel_count, relative_error)
        integer, intent(in) :: panel_count
        real(dp), intent(out) :: relative_error

        real(dp), parameter :: angular_frequency = 1.1_dp
        real(dp), parameter :: fluid_density = 1.7_dp
        real(dp), parameter :: poisson_ratio = 0.25_dp
        real(dp), parameter :: radius = 0.8_dp
        real(dp), parameter :: solid_density = 2.0_dp
        real(dp), parameter :: sound_speed = angular_frequency/1.4_dp
        complex(dp), parameter :: young_modulus = (10.0_dp, 0.2_dp)
        complex(dp), allocatable :: dirichlet_values(:), incident_pressure(:)
        complex(dp), allocatable :: load(:), solution(:)
        integer, allocatable :: dirichlet_dofs(:), interface_nodes(:)
        integer, allocatable :: panel_nodes(:, :), triangles(:, :)
        real(dp), allocatable :: vertices(:, :)
        complex(dp) :: dtn_eigenvalue, pressure_amplitude
        real(dp) :: angle, area, determinant, length, normal(2), pi
        integer :: component, first, node, panel, second, status

        allocate(vertices(2, panel_count + 1), triangles(3, panel_count))
        allocate(interface_nodes(panel_count), panel_nodes(2, panel_count))
        vertices(:, 1) = 0.0_dp
        pi = acos(-1.0_dp)
        do panel = 1, panel_count
            angle = 2.0_dp*pi*real(panel - 1, dp)/real(panel_count, dp)
            vertices(1, panel + 1) = radius*cos(angle)
            vertices(2, panel + 1) = radius*sin(angle)
            interface_nodes(panel) = panel + 1
            panel_nodes(1, panel) = panel
            panel_nodes(2, panel) = modulo(panel, panel_count) + 1
            triangles(:, panel) = [ &
                1, panel + 1, modulo(panel, panel_count) + 2]
        end do

        allocate(load(2*(panel_count + 1)), solution(2*(panel_count + 1)))
        load = (0.0_dp, 0.0_dp)
        do panel = 1, panel_count
            determinant = (vertices(1, triangles(2, panel)) - &
                vertices(1, 1))*(vertices(2, triangles(3, panel)) - &
                vertices(2, 1)) - &
                (vertices(1, triangles(3, panel)) - vertices(1, 1))* &
                (vertices(2, triangles(2, panel)) - vertices(2, 1))
            area = 0.5_dp*abs(determinant)
            do node = 1, 3
                load(2*triangles(node, panel) - 1) = &
                    load(2*triangles(node, panel) - 1) - &
                    solid_density*angular_frequency**2*area/3.0_dp
            end do
        end do

        call circular_helmholtz_dtn_eigenvalue( &
            1, angular_frequency/sound_speed, radius, dtn_eigenvalue, status)
        if (status /= 0) error stop "analytical circular DtN failed"
        pressure_amplitude = fluid_density*angular_frequency**2/dtn_eigenvalue
        do panel = 1, panel_count
            length = norm2( &
                vertices(:, interface_nodes(panel_nodes(2, panel))) - &
                vertices(:, interface_nodes(panel_nodes(1, panel))))
            normal = [ &
                vertices(2, interface_nodes(panel_nodes(2, panel))) - &
                vertices(2, interface_nodes(panel_nodes(1, panel))), &
                vertices(1, interface_nodes(panel_nodes(1, panel))) - &
                vertices(1, interface_nodes(panel_nodes(2, panel)))]/length
            do first = 1, 2
                do second = 1, 2
                    angle = 2.0_dp*pi* &
                        real(panel_nodes(second, panel) - 1, dp)/ &
                        real(panel_count, dp)
                    do component = 1, 2
                        node = 2*interface_nodes( &
                            panel_nodes(first, panel)) - 2 + component
                        load(node) = load(node) + normal(component)*length* &
                            merge(2.0_dp, 1.0_dp, first == second)/6.0_dp* &
                            pressure_amplitude*cos(angle)
                    end do
                end do
            end do
        end do

        allocate(incident_pressure(panel_count))
        incident_pressure = (0.0_dp, 0.0_dp)
        allocate(dirichlet_dofs(0), dirichlet_values(0))
        call solve_elasticity_curved_acoustic_ntd_p1( &
            vertices, triangles, interface_nodes, panel_nodes, &
            angular_frequency, sound_speed, fluid_density, young_modulus, &
            poisson_ratio, solid_density, load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, solution, status)
        call record_condition(status == 0, &
            "Curved elastic-acoustic solve succeeds")
        relative_error = 0.0_dp
        do node = 1, panel_count + 1
            relative_error = relative_error + &
                abs(solution(2*node - 1) - 1.0_dp)**2 + &
                abs(solution(2*node))**2
        end do
        relative_error = sqrt(relative_error/real(panel_count + 1, dp))
    end subroutine rigid_translation_error

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_elasticity_curved_acoustic_ntd_convergence
