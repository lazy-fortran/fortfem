program test_maxwell_surface_rt_efie
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_efie_rwg_3d, assemble_maxwell_surface_rt_efie_3d, &
        build_maxwell_surface_rt_dof_map
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: triangles(3, 4) = reshape([ &
        1, 3, 2, 1, 2, 4, 2, 3, 4, 3, 1, 4], [3, 4])
    real(dp), parameter :: vertices(3, 4) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    complex(dp), allocatable :: rt0(:, :), rt1(:, :), rwg(:, :), scaled(:, :)
    integer, allocatable :: edges(:, :), global_dofs(:, :), transforms(:, :)
    real(dp), allocatable :: scales(:)
    real(dp) :: equivalence_error, symmetry_error
    integer :: edge, status

    call assemble_maxwell_efie_rwg_3d( &
        vertices, triangles, 0.7_dp, 1.3_dp, 4, 2.0e-5_dp, 2, rwg, status)
    if (status /= 0) error stop "RWG EFIE oracle failed"
    call assemble_maxwell_surface_rt_efie_3d( &
        vertices, triangles, 0, 0.7_dp, 1.3_dp, 4, 2.0e-5_dp, 2, rt0, status)
    if (status /= 0) error stop "surface RT0 EFIE failed"
    call build_maxwell_surface_rt_dof_map( &
        triangles, 0, edges, global_dofs, transforms, status)
    allocate(scales(size(edges, 2)))
    do edge = 1, size(edges, 2)
        scales(edge) = norm2( &
            vertices(:, edges(2, edge)) - vertices(:, edges(1, edge)))
    end do
    scaled = spread(scales, 2, size(scales))*rt0* &
        spread(scales, 1, size(scales))
    equivalence_error = maxval(abs(scaled - rwg))
    write (*, "(a,es12.4)") "surface RT0/RWG EFIE error: ", equivalence_error
    call check_condition(equivalence_error < 3.0e-4_dp, &
        "Normalized surface RT0 EFIE agrees with independent RWG assembly")

    call assemble_maxwell_surface_rt_efie_3d( &
        vertices, triangles, 1, 0.7_dp, 1.3_dp, 5, 2.0e-5_dp, 2, rt1, status)
    if (status /= 0) error stop "surface RT1 EFIE failed"
    symmetry_error = maxval(abs(rt1 - transpose(rt1)))
    write (*, "(a,es12.4)") "surface RT1 EFIE symmetry error: ", symmetry_error
    call check_condition( &
        all(shape(rt1) == [20, 20]) .and. symmetry_error < 3.0e-11_dp, &
        "Surface RT1 EFIE is a complex-symmetric edge-and-panel operator")
    call check_summary("Higher-order Maxwell surface RT EFIE")

end program test_maxwell_surface_rt_efie
