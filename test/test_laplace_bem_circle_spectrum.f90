program test_laplace_bem_circle_spectrum
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_laplace_hypersingular_linear, &
        assemble_laplace_single_layer_constant
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: mesh_sizes(3) = [12, 24, 48]
    real(dp) :: errors(3)
    integer :: mesh_id
    logical :: all_passed

    all_passed = .true.
    do mesh_id = 1, size(mesh_sizes)
        call circle_spectral_error(mesh_sizes(mesh_id), errors(mesh_id))
    end do

    call record_condition(errors(2) < errors(1) .and. errors(3) < errors(2), &
        "Laplace boundary spectra converge under circle refinement")
    call record_condition(errors(3) < 1.0e-2_dp, &
        "Laplace boundary spectra resolve the first two circle modes")

    call check_summary("Laplace BEM circle spectrum")
    if (.not. all_passed) error stop 1

contains

    subroutine circle_spectral_error(panel_count, maximum_error)
        integer, intent(in) :: panel_count
        real(dp), intent(out) :: maximum_error

        integer, allocatable :: panel_nodes(:, :)
        real(dp), allocatable :: hypersingular(:, :), nodal_mode(:)
        real(dp), allocatable :: panel_end(:, :), panel_mode(:)
        real(dp), allocatable :: panel_start(:, :), single_layer(:, :)
        real(dp) :: angle, length, mass_norm, numerical, pi
        real(dp) :: relative_error
        integer :: mode, panel, status

        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        allocate(panel_nodes(2, panel_count))
        allocate(single_layer(panel_count, panel_count))
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

        call assemble_laplace_single_layer_constant( &
            panel_start, panel_end, 20, single_layer, status)
        call record_condition(status == 0, &
            "Circle single-layer assembly succeeds")
        call assemble_laplace_hypersingular_linear( &
            panel_start, panel_end, panel_nodes, panel_count, 20, &
            hypersingular, status)
        call record_condition(status == 0, &
            "Circle hypersingular assembly succeeds")

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
            relative_error = abs(numerical - 0.5_dp / real(mode, dp)) / &
                (0.5_dp / real(mode, dp))
            maximum_error = max(maximum_error, relative_error)

            mass_norm = linear_mass_norm(nodal_mode, length)
            numerical = dot_product( &
                nodal_mode, matmul(hypersingular, nodal_mode)) / mass_norm
            relative_error = abs(numerical - 0.5_dp * real(mode, dp)) / &
                (0.5_dp * real(mode, dp))
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

end program test_laplace_bem_circle_spectrum
