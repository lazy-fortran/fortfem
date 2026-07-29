program test_mixed_poisson_rt0_convergence
    use check, only: check_condition, check_summary
    use fortfem_api, only: mesh_t, rectangle_mesh, solve_mixed_poisson_rt0
    use fortfem_gauss_quadrature_2d, only: &
        gauss_quadrature_triangle_t, get_gauss_quadrature_triangle
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_rt_field_2d, only: evaluate_rt_field_2d
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp) :: coarse_errors(2), fine_errors(2), first_rates(2)
    real(dp) :: medium_errors(2), second_rates(2)
    logical :: all_passed

    all_passed = .true.
    call solve_manufactured_problem(9, coarse_errors)
    call solve_manufactured_problem(17, medium_errors)
    call solve_manufactured_problem(33, fine_errors)
    first_rates = log(coarse_errors / medium_errors) / log(2.0_dp)
    second_rates = log(medium_errors / fine_errors) / log(2.0_dp)

    call record_condition( &
        minval(first_rates) > 0.8_dp .and. minval(second_rates) > 0.8_dp, &
        "RT0-DG0 flux and pressure attain first-order convergence")
    call record_condition(fine_errors(1) < 1.2e-1_dp, &
        "Mixed Poisson flux reaches the analytical gradient")
    call record_condition(fine_errors(2) < 2.5e-2_dp, &
        "Mixed Poisson pressure reaches the analytical solution")

    call check_summary("RT0-DG0 mixed Poisson convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine solve_manufactured_problem(vertex_count, errors)
        integer, intent(in) :: vertex_count
        real(dp), intent(out) :: errors(2)

        type(fortsparse_status_t) :: status
        type(mesh_t) :: mesh
        real(dp), allocatable :: cell_pressure(:), cell_source(:), flux_dofs(:)

        mesh = rectangle_mesh( &
            vertex_count, vertex_count, &
            [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
        call project_cell_source(mesh%data, cell_source)
        call solve_mixed_poisson_rt0( &
            mesh%data, cell_source, flux_dofs, cell_pressure, status)
        if (status%code /= 0) error stop "Manufactured mixed solve failed"
        call compute_errors( &
            mesh%data, flux_dofs, cell_pressure, errors)
    end subroutine solve_manufactured_problem

    subroutine project_cell_source(mesh, cell_source)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), allocatable, intent(out) :: cell_source(:)

        type(gauss_quadrature_triangle_t) :: quadrature
        real(dp) :: point(2), vertex_a(2), vertex_b(2), vertex_c(2)
        integer :: point_index, triangle

        quadrature = get_gauss_quadrature_triangle(7)
        allocate(cell_source(mesh%n_triangles))
        do triangle = 1, mesh%n_triangles
            vertex_a = mesh%vertices(:, mesh%triangles(1, triangle))
            vertex_b = mesh%vertices(:, mesh%triangles(2, triangle))
            vertex_c = mesh%vertices(:, mesh%triangles(3, triangle))
            cell_source(triangle) = 0.0_dp
            do point_index = 1, quadrature%n_points
                point = vertex_a + &
                    quadrature%xi(point_index) * (vertex_b - vertex_a) + &
                    quadrature%eta(point_index) * (vertex_c - vertex_a)
                cell_source(triangle) = cell_source(triangle) + &
                    2.0_dp * quadrature%weights(point_index) * &
                    exact_source(point)
            end do
        end do
        call quadrature%destroy()
    end subroutine project_cell_source

    subroutine compute_errors(mesh, flux_dofs, cell_pressure, errors)
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(in) :: flux_dofs(:), cell_pressure(:)
        real(dp), intent(out) :: errors(2)

        type(gauss_quadrature_triangle_t) :: quadrature
        complex(dp), allocatable :: complex_flux(:)
        complex(dp) :: numerical_flux(2)
        real(dp) :: determinant, exact_flux(2), point(2), weight
        real(dp) :: vertex_a(2), vertex_b(2), vertex_c(2)
        integer :: point_index, triangle

        quadrature = get_gauss_quadrature_triangle(7)
        allocate(complex_flux(size(flux_dofs)))
        complex_flux = cmplx(flux_dofs, 0.0_dp, dp)
        errors = 0.0_dp
        do triangle = 1, mesh%n_triangles
            vertex_a = mesh%vertices(:, mesh%triangles(1, triangle))
            vertex_b = mesh%vertices(:, mesh%triangles(2, triangle))
            vertex_c = mesh%vertices(:, mesh%triangles(3, triangle))
            determinant = &
                (vertex_b(1) - vertex_a(1)) * &
                (vertex_c(2) - vertex_a(2)) - &
                (vertex_b(2) - vertex_a(2)) * &
                (vertex_c(1) - vertex_a(1))
            do point_index = 1, quadrature%n_points
                point = vertex_a + &
                    quadrature%xi(point_index) * (vertex_b - vertex_a) + &
                    quadrature%eta(point_index) * (vertex_c - vertex_a)
                weight = determinant * quadrature%weights(point_index)
                call evaluate_rt_field_2d( &
                    mesh, triangle, point(1), point(2), complex_flux, &
                    numerical_flux)
                exact_flux = exact_darcy_flux(point)
                errors(1) = errors(1) + weight * &
                    sum((real(numerical_flux, dp) - exact_flux)**2)
                errors(2) = errors(2) + weight * &
                    (cell_pressure(triangle) - exact_pressure(point))**2
            end do
        end do
        call quadrature%destroy()
        errors = sqrt(errors)
    end subroutine compute_errors

    pure real(dp) function exact_pressure(point) result(value)
        real(dp), intent(in) :: point(2)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = sin(pi * point(1)) * sin(pi * point(2))
    end function exact_pressure

    pure function exact_darcy_flux(point) result(value)
        real(dp), intent(in) :: point(2)
        real(dp) :: value(2)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value(1) = -pi * cos(pi * point(1)) * sin(pi * point(2))
        value(2) = -pi * sin(pi * point(1)) * cos(pi * point(2))
    end function exact_darcy_flux

    pure real(dp) function exact_source(point) result(value)
        real(dp), intent(in) :: point(2)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = 2.0_dp * pi**2 * exact_pressure(point)
    end function exact_source

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_mixed_poisson_rt0_convergence
