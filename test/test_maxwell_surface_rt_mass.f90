program test_maxwell_surface_rt_mass
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_maxwell_rwg_mass_matrix, &
        assemble_maxwell_surface_rt_mass_matrix, &
        build_maxwell_rwg_surface_space, &
        build_maxwell_surface_rt_dof_map, evaluate_maxwell_rwg_basis, &
        evaluate_maxwell_surface_rt_global_basis
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: triangles(3, 4) = reshape([ &
        1, 3, 2, 1, 2, 4, 2, 3, 4, 3, 1, 4], [3, 4])
    real(dp), parameter :: vertices(3, 4) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    integer, allocatable :: edge_vertices(:, :), global_dofs(:, :)
    integer, allocatable :: rwg_edge_triangles(:, :), rwg_edges(:, :)
    integer, allocatable :: transforms(:, :)
    real(dp), allocatable :: coefficients(:), high_mass(:, :)
    real(dp), allocatable :: rt0_mass(:, :), rwg_mass(:, :)
    real(dp), allocatable :: scaled_mass(:, :), scales(:)
    real(dp) :: basis_error, centroid(3), divergence, quadratic_form
    real(dp) :: rt0_error, rt_divergence, rt_value(3), rwg_value(3)
    real(dp) :: symmetry_error
    integer :: dof, panel, status

    call assemble_maxwell_rwg_mass_matrix( &
        vertices, triangles, 8, rwg_mass, status)
    if (status /= 0) error stop "RWG mass oracle failed"
    call assemble_maxwell_surface_rt_mass_matrix( &
        vertices, triangles, 0, 8, rt0_mass, status)
    if (status /= 0) error stop "surface RT0 mass assembly failed"
    call build_maxwell_surface_rt_dof_map( &
        triangles, 0, edge_vertices, global_dofs, transforms, status)
    allocate(scales(size(edge_vertices, 2)))
    do dof = 1, size(scales)
        scales(dof) = norm2( &
            vertices(:, edge_vertices(2, dof)) - &
            vertices(:, edge_vertices(1, dof)))
    end do
    scaled_mass = spread(scales, 2, size(scales))*rt0_mass* &
        spread(scales, 1, size(scales))
    rt0_error = maxval(abs(scaled_mass - rwg_mass))
    write (*, "(a,es12.4)") "surface RT0/RWG mass error: ", rt0_error
    call check_condition(rt0_error < 3.0e-12_dp, &
        "Normalized surface RT0 mass equals edge-scaled RWG assembly")
    call build_maxwell_rwg_surface_space( &
        vertices, triangles, rwg_edges, rwg_edge_triangles, status)
    basis_error = 0.0_dp
    do dof = 1, size(rwg_edges, 2)
        do panel = 1, 2
            centroid = sum( &
                vertices(:, triangles(:, rwg_edge_triangles(panel, dof))), &
                dim=2)/3.0_dp
            call evaluate_maxwell_rwg_basis( &
                vertices, triangles, rwg_edges, rwg_edge_triangles, dof, &
                rwg_edge_triangles(panel, dof), centroid, rwg_value, &
                divergence, status)
            call evaluate_maxwell_surface_rt_global_basis( &
                vertices, triangles, 0, rwg_edge_triangles(panel, dof), dof, &
                centroid, rt_value, rt_divergence, status)
            basis_error = max( &
                basis_error, maxval(abs( &
                rt_value - rwg_value/scales(dof))))
            basis_error = max( &
                basis_error, abs(rt_divergence - divergence/scales(dof)))
        end do
    end do
    write (*, "(a,es12.4)") "surface RT0/RWG basis error: ", basis_error
    call check_condition(basis_error < 3.0e-12_dp, &
        "Global surface RT0 evaluator equals normalized RWG traces")

    call assemble_maxwell_surface_rt_mass_matrix( &
        vertices, triangles, 3, 12, high_mass, status)
    if (status /= 0) error stop "higher-order surface RT mass assembly failed"
    call check_condition(all(shape(high_mass) == [72, 72]), &
        "Degree-three surface RT mass includes edge and panel moments")
    symmetry_error = maxval(abs(high_mass - transpose(high_mass)))
    allocate(coefficients(size(high_mass, 1)))
    do dof = 1, size(coefficients)
        coefficients(dof) = sin(real(3*dof + 1, dp))
    end do
    quadratic_form = dot_product(coefficients, matmul(high_mass, coefficients))
    write (*, "(a,es12.4)") &
        "degree-three surface RT mass symmetry error: ", symmetry_error
    write (*, "(a,es12.4)") &
        "degree-three surface RT mass quadratic form: ", quadratic_form
    call check_condition( &
        symmetry_error < 3.0e-12_dp .and. quadratic_form > 1.0e-8_dp .and. &
        all([(high_mass(dof, dof) > 0.0_dp, dof=1, size(high_mass, 1))]), &
        "Degree-three surface RT mass is symmetric positive definite")
    call check_summary("Higher-order Maxwell surface RT mass")

end program test_maxwell_surface_rt_mass
