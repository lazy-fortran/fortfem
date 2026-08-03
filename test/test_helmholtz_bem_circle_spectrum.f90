program test_helmholtz_bem_circle_spectrum
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_helmholtz_double_layer_constant, &
        assemble_helmholtz_hypersingular_linear, &
        assemble_helmholtz_single_layer_constant
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: mesh_sizes(3) = [12, 24, 48]
    complex(dp), parameter :: exact_single_layer(2) = [ &
        (0.4497818998782933_dp, 0.4280549908588178_dp), &
        (0.3249907729053004_dp, 0.052619748735733766_dp)]
    complex(dp), parameter :: exact_double_layer(2) = [ &
        (-0.25522568497710474_dp, 0.2329503859714022_dp), &
        (0.05502603904479714_dp, 0.08986510741594829_dp)]
    complex(dp), parameter :: exact_hypersingular(2) = [ &
        (0.41099886362254484_dp, -0.1267731564473764_dp), &
        (0.7599358370060741_dp, -0.15347350994467426_dp)]
    real(dp) :: errors(3)
    integer :: mesh_id
    logical :: all_passed

    all_passed = .true.
    do mesh_id = 1, size(mesh_sizes)
        call circle_spectral_error(mesh_sizes(mesh_id), errors(mesh_id))
    end do

    call record_condition(errors(1) > 3.5_dp * errors(2) .and. &
        errors(2) > 3.5_dp * errors(3), &
        "Outgoing Helmholtz boundary spectra converge at second order")
    call record_condition(errors(3) < 2.0e-2_dp, &
        "Helmholtz Calderon operators resolve the first two circle modes")

    call check_summary("Helmholtz BEM circle spectrum")
    if (.not. all_passed) error stop 1

contains

    subroutine circle_spectral_error(panel_count, maximum_error)
        integer, intent(in) :: panel_count
        real(dp), intent(out) :: maximum_error

        complex(dp), allocatable :: double_layer(:, :), hypersingular(:, :)
        complex(dp), allocatable :: single_layer(:, :)
        integer, allocatable :: panel_nodes(:, :)
        real(dp), allocatable :: nodal_mode(:), panel_end(:, :)
        real(dp), allocatable :: panel_mode(:), panel_start(:, :)
        complex(dp) :: numerical
        real(dp) :: angle, length, mass_norm, pi, relative_error
        integer :: mode, panel, status

        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        allocate(panel_nodes(2, panel_count))
        allocate(single_layer(panel_count, panel_count))
        allocate(double_layer(panel_count, panel_count))
        allocate(hypersingular(panel_count, panel_count))
        allocate(panel_mode(panel_count), nodal_mode(panel_count))

        pi = acos(-1.0_dp)
        do panel = 1, panel_count
            angle = 2.0_dp * pi * real(panel - 1, dp) / real(panel_count, dp)
            panel_start(:, panel) = [cos(angle), sin(angle)]
            angle = 2.0_dp * pi * real(panel, dp) / real(panel_count, dp)
            panel_end(:, panel) = [cos(angle), sin(angle)]
            panel_nodes(1, panel) = panel
            panel_nodes(2, panel) = modulo(panel, panel_count) + 1
        end do
        length = norm2(panel_end(:, 1) - panel_start(:, 1))

        call assemble_helmholtz_single_layer_constant( &
            panel_start, panel_end, 1.3_dp, 16, single_layer, status)
        call record_condition(status == 0, &
            "Circle Helmholtz single-layer assembly succeeds")
        call assemble_helmholtz_double_layer_constant( &
            panel_start, panel_end, 1.3_dp, 16, double_layer, status)
        call record_condition(status == 0, &
            "Circle Helmholtz double-layer assembly succeeds")
        call assemble_helmholtz_hypersingular_linear( &
            panel_start, panel_end, panel_nodes, panel_count, 1.3_dp, 16, &
            hypersingular, status)
        call record_condition(status == 0, &
            "Circle Helmholtz hypersingular assembly succeeds")

        maximum_error = 0.0_dp
        do mode = 1, 2
            do panel = 1, panel_count
                angle = 2.0_dp * pi * &
                    (real(panel, dp) - 0.5_dp) / real(panel_count, dp)
                panel_mode(panel) = cos(real(mode, dp) * angle)
                angle = 2.0_dp * pi * real(panel - 1, dp) / &
                    real(panel_count, dp)
                nodal_mode(panel) = cos(real(mode, dp) * angle)
            end do

            mass_norm = length * dot_product(panel_mode, panel_mode)
            numerical = dot_product( &
                panel_mode, matmul(single_layer, panel_mode)) / mass_norm
            relative_error = abs(numerical - exact_single_layer(mode)) / &
                abs(exact_single_layer(mode))
            maximum_error = max(maximum_error, relative_error)

            numerical = dot_product( &
                panel_mode, matmul(double_layer, panel_mode)) / mass_norm
            relative_error = abs(numerical - exact_double_layer(mode)) / &
                abs(exact_double_layer(mode))
            maximum_error = max(maximum_error, relative_error)

            mass_norm = linear_mass_norm(nodal_mode, length)
            numerical = dot_product( &
                nodal_mode, matmul(hypersingular, nodal_mode)) / mass_norm
            relative_error = abs(numerical - exact_hypersingular(mode)) / &
                abs(exact_hypersingular(mode))
            maximum_error = max(maximum_error, relative_error)
        end do
    end subroutine circle_spectral_error

    pure function linear_mass_norm(values, panel_length) result(norm)
        real(dp), intent(in) :: values(:), panel_length
        real(dp) :: norm

        integer :: panel, successor

        norm = 0.0_dp
        do panel = 1, size(values)
            successor = modulo(panel, size(values)) + 1
            norm = norm + panel_length / 3.0_dp * ( &
                values(panel)**2 + values(panel) * values(successor) + &
                values(successor)**2)
        end do
    end function linear_mass_norm

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_helmholtz_bem_circle_spectrum
