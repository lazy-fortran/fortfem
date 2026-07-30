program test_maxwell_surface_rt_basis
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_maxwell_surface_rt_dof_map, &
        evaluate_maxwell_surface_rt_basis, &
        initialize_triangle_raviart_thomas, triangle_rt_basis_t
    use fortfem_kinds, only: dp
    implicit none

    type(triangle_rt_basis_t) :: basis
    real(dp) :: divergences(3), flux_moments(3, 3), values(3, 3)
    real(dp), parameter :: edge_midpoints(2, 3) = reshape([ &
        0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp], [2, 3])
    real(dp), parameter :: outward_conormals(3, 3) = reshape([ &
        0.0_dp, -1.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, &
        -1.0_dp, 0.0_dp, 0.0_dp], [3, 3])
    integer, allocatable :: edge_vertices(:, :), global_dofs(:, :)
    integer, allocatable :: transforms(:, :)
    real(dp), allocatable :: first_divergences(:), first_values(:, :)
    real(dp), allocatable :: second_divergences(:), second_values(:, :)
    integer, parameter :: triangles(3, 2) = reshape([1, 2, 3, 1, 3, 4], [3, 2])
    real(dp), parameter :: vertices(3, 4) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], [3, 4])
    real(dp) :: first_flux, normal_jump, second_flux, trace_coordinate
    integer :: edge, moment, sample, shared_dof, status

    call initialize_triangle_raviart_thomas(0, basis, status)
    if (status /= 0) error stop "surface RT basis initialization failed"
    do edge = 1, 3
        call evaluate_maxwell_surface_rt_basis( &
            basis, edge_midpoints(1, edge), edge_midpoints(2, edge), &
            [1.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], &
            1.0_dp, values, divergences, status)
        if (status /= 0) error stop "surface RT basis evaluation failed"
        flux_moments(edge, :) = matmul( &
            transpose(values), outward_conormals(:, edge))
    end do
    write (*, "(a,9(1x,es10.2))") "surface RT0 edge moments:", flux_moments
    call check_condition( &
        maxval(abs(flux_moments - identity3())) < 2.0e-12_dp, &
        "Surface RT0 basis has Kronecker normal-edge fluxes")

    call build_maxwell_surface_rt_dof_map( &
        triangles, 0, edge_vertices, global_dofs, transforms, status)
    if (status /= 0) error stop "surface RT global map failed"
    shared_dof = global_dofs(3, 1)
    call check_condition(shared_dof == global_dofs(1, 2), &
        "Adjacent panels share the surface RT edge degree of freedom")
    call evaluate_maxwell_surface_rt_basis( &
        basis, 0.0_dp, 0.5_dp, vertices(:, 2) - vertices(:, 1), &
        vertices(:, 3) - vertices(:, 1), 1.0_dp, values, divergences, status)
    first_flux = real(transforms(3, 1), dp)*dot_product( &
        values(:, 3), [-1.0_dp, 1.0_dp, 0.0_dp]/sqrt(2.0_dp))
    call evaluate_maxwell_surface_rt_basis( &
        basis, 0.5_dp, 0.0_dp, vertices(:, 3) - vertices(:, 1), &
        vertices(:, 4) - vertices(:, 1), 1.0_dp, values, divergences, status)
    second_flux = real(transforms(1, 2), dp)*dot_product( &
        values(:, 1), [1.0_dp, -1.0_dp, 0.0_dp]/sqrt(2.0_dp))
    call check_condition(abs(first_flux + second_flux) < 2.0e-12_dp, &
        "Surface RT global orientation preserves normal continuity")

    call initialize_triangle_raviart_thomas(3, basis, status)
    call build_maxwell_surface_rt_dof_map( &
        triangles, 3, edge_vertices, global_dofs, transforms, status)
    allocate( &
        first_values(3, 24), first_divergences(24), &
        second_values(3, 24), second_divergences(24))
    normal_jump = 0.0_dp
    do sample = 1, 4
        trace_coordinate = real(sample, dp)/5.0_dp
        call evaluate_maxwell_surface_rt_basis( &
            basis, 0.0_dp, trace_coordinate, &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1), 1.0_dp, first_values, &
            first_divergences, status)
        call evaluate_maxwell_surface_rt_basis( &
            basis, trace_coordinate, 0.0_dp, &
            vertices(:, 3) - vertices(:, 1), &
            vertices(:, 4) - vertices(:, 1), 1.0_dp, second_values, &
            second_divergences, status)
        do moment = 1, 4
            first_flux = real(transforms(8 + moment, 1), dp)*dot_product( &
                first_values(:, 8 + moment), &
                [-1.0_dp, 1.0_dp, 0.0_dp]/sqrt(2.0_dp))
            second_flux = real(transforms(moment, 2), dp)*dot_product( &
                second_values(:, moment), &
                [1.0_dp, -1.0_dp, 0.0_dp]/sqrt(2.0_dp))
            normal_jump = max(normal_jump, abs(first_flux + second_flux))
        end do
    end do
    write (*, "(a,es12.4)") "degree-three surface RT normal jump: ", normal_jump
    call check_condition(normal_jump < 2.0e-11_dp, &
        "Degree-three surface RT traces are pointwise normal-continuous")
    call check_summary("Higher-order Maxwell surface RT basis")

contains

    pure function identity3() result(matrix)
        real(dp) :: matrix(3, 3)

        matrix = 0.0_dp
        matrix(1, 1) = 1.0_dp
        matrix(2, 2) = 1.0_dp
        matrix(3, 3) = 1.0_dp
    end function identity3

end program test_maxwell_surface_rt_basis
