program test_mixed_poisson_rt1_convergence
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_triangle_discontinuous_dof_map, &
        build_triangle_trimmed_dof_map, evaluate_triangle_lagrange_basis, &
        evaluate_triangle_rt_interpolant, initialize_triangle_lagrange_basis, &
        initialize_triangle_raviart_thomas, mesh_t, rectangle_mesh, &
        solve_mixed_poisson_rt, triangle_duffy_quadrature, &
        triangle_lagrange_basis_t, triangle_rt_basis_t
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp) :: coarse_errors(2), fine_errors(2), first_rates(2)
    real(dp) :: medium_errors(2), second_rates(2)
    logical :: all_passed

    all_passed = .true.
    call solve_and_measure(5, coarse_errors)
    call solve_and_measure(9, medium_errors)
    call solve_and_measure(17, fine_errors)
    first_rates = log(coarse_errors / medium_errors) / log(2.0_dp)
    second_rates = log(medium_errors / fine_errors) / log(2.0_dp)

    call record_condition( &
        minval(first_rates) > 1.7_dp .and. minval(second_rates) > 1.7_dp, &
        "RT1-DG1 flux and pressure attain second-order convergence")
    call record_condition(fine_errors(1) < 2.0e-2_dp, &
        "RT1 flux reaches the analytical Darcy field")
    call record_condition(fine_errors(2) < 3.0e-3_dp, &
        "DG1 pressure reaches the analytical solution")

    call check_summary("RT1-DG1 mixed Poisson convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine solve_and_measure(vertex_count, errors)
        integer, intent(in) :: vertex_count
        real(dp), intent(out) :: errors(2)

        type(fortsparse_status_t) :: status
        type(mesh_t) :: mesh
        real(dp), allocatable :: flux_dofs(:), pressure_dofs(:)

        mesh = rectangle_mesh( &
            vertex_count, vertex_count, &
            [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
        call solve_mixed_poisson_rt( &
            mesh%data, 1, 8, exact_source, flux_dofs, pressure_dofs, status)
        if (status%code /= 0) error stop "RT1-DG1 mixed solve failed"
        call compute_errors( &
            mesh%data, 1, flux_dofs, pressure_dofs, errors)
    end subroutine solve_and_measure

    subroutine compute_errors( &
            mesh, degree, flux_dofs, pressure_dofs, errors)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree
        real(dp), intent(in) :: flux_dofs(:), pressure_dofs(:)
        real(dp), intent(out) :: errors(2)

        type(triangle_lagrange_basis_t) :: pressure_basis
        type(triangle_rt_basis_t) :: flux_basis
        integer, allocatable :: flux_map(:, :), pressure_map(:, :)
        integer, allocatable :: transforms(:, :)
        real(dp), allocatable :: eta(:), local_flux(:), local_pressure(:)
        real(dp), allocatable :: pressure_gradients(:, :)
        real(dp), allocatable :: pressure_values(:), weights(:), xi(:)
        real(dp) :: determinant, divergence, exact_flux(2)
        real(dp) :: numerical_flux(2), numerical_pressure, point(2)
        real(dp) :: vertices(2, 3)
        integer :: flux_count, local_dof, local_status, point_index
        integer :: pressure_count, triangle

        call build_triangle_trimmed_dof_map( &
            mesh, degree + 1, flux_map, transforms, flux_count, local_status)
        if (local_status /= 0) error stop "RT1 flux map failed"
        call build_triangle_discontinuous_dof_map( &
            mesh, degree, pressure_map, pressure_count, local_status)
        if (local_status /= 0) error stop "DG1 pressure map failed"
        call initialize_triangle_raviart_thomas( &
            degree, flux_basis, local_status)
        if (local_status /= 0) error stop "RT1 basis failed"
        call initialize_triangle_lagrange_basis( &
            degree, pressure_basis, local_status)
        if (local_status /= 0) error stop "DG1 basis failed"
        call triangle_duffy_quadrature(8, xi, eta, weights, local_status)
        if (local_status /= 0) error stop "Mixed error quadrature failed"

        allocate(local_flux(size(flux_map, 1)))
        allocate(local_pressure(size(pressure_map, 1)))
        allocate(pressure_values(size(pressure_map, 1)))
        allocate(pressure_gradients(2, size(pressure_map, 1)))
        errors = 0.0_dp
        do triangle = 1, mesh%n_triangles
            do local_dof = 1, size(flux_map, 1)
                local_flux(local_dof) = &
                    real(transforms(local_dof, triangle), dp) * &
                    flux_dofs(flux_map(local_dof, triangle))
            end do
            local_pressure = pressure_dofs(pressure_map(:, triangle))
            vertices(:, 1) = &
                mesh%vertices(:, mesh%triangles(1, triangle))
            vertices(:, 2) = &
                mesh%vertices(:, mesh%triangles(2, triangle))
            vertices(:, 3) = &
                mesh%vertices(:, mesh%triangles(3, triangle))
            determinant = &
                (vertices(1, 2) - vertices(1, 1)) * &
                (vertices(2, 3) - vertices(2, 1)) - &
                (vertices(2, 2) - vertices(2, 1)) * &
                (vertices(1, 3) - vertices(1, 1))
            do point_index = 1, size(weights)
                point = vertices(:, 1) + &
                    xi(point_index) * (vertices(:, 2) - vertices(:, 1)) + &
                    eta(point_index) * (vertices(:, 3) - vertices(:, 1))
                call evaluate_triangle_rt_interpolant( &
                    vertices, flux_basis, local_flux, xi(point_index), &
                    eta(point_index), numerical_flux, divergence, local_status)
                call evaluate_triangle_lagrange_basis( &
                    pressure_basis, xi(point_index), eta(point_index), &
                    pressure_values, pressure_gradients, local_status)
                numerical_pressure = dot_product( &
                    pressure_values, local_pressure)
                exact_flux = exact_darcy_flux(point)
                errors(1) = errors(1) + determinant * weights(point_index) * &
                    sum((numerical_flux - exact_flux)**2)
                errors(2) = errors(2) + determinant * weights(point_index) * &
                    (numerical_pressure - exact_pressure(point))**2
            end do
        end do
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

    pure real(dp) function exact_source(x, y) result(value)
        real(dp), intent(in) :: x, y
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = 2.0_dp * pi**2 * sin(pi * x) * sin(pi * y)
    end function exact_source

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_mixed_poisson_rt1_convergence
