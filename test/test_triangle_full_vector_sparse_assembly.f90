program test_triangle_full_vector_sparse_assembly
    use check, only: check_condition, check_summary
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    use fortfem_feec, only: &
        assemble_triangle_bdm_div_mass_csc, &
        assemble_triangle_nedelec_second_curl_mass_csc, &
        build_triangle_full_vector_dof_map, interpolate_triangle_bdm, &
        interpolate_triangle_nedelec_second_kind
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none

    type(mesh_2d_t) :: mesh
    type(csc_t) :: bdm_matrix, nedelec_matrix
    type(fortsparse_status_t) :: sparse_status
    real(dp), allocatable :: bdm_dofs(:), local_dofs(:)
    real(dp), allocatable :: matrix_times_dofs(:), nedelec_dofs(:)
    integer, allocatable :: global_dofs(:, :), transforms(:, :)
    logical, allocatable :: assigned(:)
    real(dp) :: exact_energy, vertices(2, 3)
    integer :: active_degree, degree, dof, global_dof_count
    integer :: local_status, triangle
    logical :: all_passed

    all_passed = .true.
    mesh%n_vertices = 4
    mesh%n_triangles = 2
    mesh%has_triangles = .true.
    allocate(mesh%vertices(2, 4), mesh%triangles(3, 2))
    mesh%vertices(:, 1) = [0.0_dp, 0.0_dp]
    mesh%vertices(:, 2) = [1.0_dp, 0.0_dp]
    mesh%vertices(:, 3) = [1.0_dp, 1.0_dp]
    mesh%vertices(:, 4) = [0.0_dp, 1.0_dp]
    mesh%triangles(:, 1) = [1, 2, 3]
    mesh%triangles(:, 2) = [1, 3, 4]

    do degree = 1, 4
        active_degree = degree
        call build_triangle_full_vector_dof_map( &
            mesh, degree, global_dofs, transforms, global_dof_count, &
            local_status)
        allocate(nedelec_dofs(global_dof_count), bdm_dofs(global_dof_count))
        allocate(assigned(global_dof_count))
        nedelec_dofs = 0.0_dp
        assigned = .false.
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call interpolate_triangle_nedelec_second_kind( &
                vertices, degree, 2 * degree + 2, nedelec_field, &
                local_dofs, local_status)
            call insert_local_dofs( &
                local_dofs, triangle, global_dofs, transforms, assigned, &
                nedelec_dofs)
            deallocate(local_dofs)
        end do

        assigned = .false.
        bdm_dofs = 0.0_dp
        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call interpolate_triangle_bdm( &
                vertices, degree, 2 * degree + 2, bdm_field, local_dofs, &
                local_status)
            call insert_local_dofs( &
                local_dofs, triangle, global_dofs, transforms, assigned, &
                bdm_dofs)
            deallocate(local_dofs)
        end do

        exact_energy = real(degree * degree, dp) / &
            real(2 * degree - 1, dp) + 1.0_dp / real(2 * degree + 1, dp)
        call assemble_triangle_nedelec_second_curl_mass_csc( &
            mesh, degree, 2 * degree, nedelec_matrix, sparse_status)
        allocate(matrix_times_dofs(global_dof_count))
        matrix_times_dofs = csc_matvec(nedelec_matrix, nedelec_dofs)
        call record_condition(sparse_status%code == 0 .and. &
            nedelec_matrix%nrow == global_dof_count, &
            "Second-kind sparse operator has exact full-family dimensions")
        call record_condition(abs( &
            dot_product(nedelec_dofs, matrix_times_dofs) - exact_energy) < &
            2.0e-7_dp, &
            "Second-kind sparse operator reproduces exact curl-mass energy")
        deallocate(matrix_times_dofs)

        call assemble_triangle_bdm_div_mass_csc( &
            mesh, degree, 2 * degree, bdm_matrix, sparse_status)
        allocate(matrix_times_dofs(global_dof_count))
        matrix_times_dofs = csc_matvec(bdm_matrix, bdm_dofs)
        call record_condition(sparse_status%code == 0 .and. &
            bdm_matrix%nrow == global_dof_count, &
            "BDM sparse operator has exact full-family dimensions")
        call record_condition(abs( &
            dot_product(bdm_dofs, matrix_times_dofs) - exact_energy) < &
            2.0e-7_dp, &
            "BDM sparse operator reproduces exact div-mass energy")

        deallocate( &
            assigned, bdm_dofs, global_dofs, matrix_times_dofs, &
            nedelec_dofs, transforms)
    end do

    call assemble_triangle_bdm_div_mass_csc( &
        mesh, 0, 2, bdm_matrix, sparse_status)
    call record_condition(sparse_status%code /= 0, &
        "Full-family sparse assembly rejects degree zero")

    call check_summary("Full polynomial vector sparse assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine insert_local_dofs( &
            local_dofs, triangle, global_dofs, transforms, assigned, dofs)
        real(dp), intent(in) :: local_dofs(:)
        integer, intent(in) :: triangle
        integer, intent(in) :: global_dofs(:, :), transforms(:, :)
        logical, intent(inout) :: assigned(:)
        real(dp), intent(inout) :: dofs(:)

        integer :: local_dof

        do local_dof = 1, size(local_dofs)
            if (.not. assigned(global_dofs(local_dof, triangle))) then
                dofs(global_dofs(local_dof, triangle)) = &
                    real(transforms(local_dof, triangle), dp) * &
                    local_dofs(local_dof)
                assigned(global_dofs(local_dof, triangle)) = .true.
            end if
        end do
    end subroutine insert_local_dofs

    subroutine nedelec_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [0.0_dp, x**active_degree]
        associate (unused_y => y)
            if (kind(unused_y) /= dp) error stop
        end associate
    end subroutine nedelec_field

    subroutine bdm_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [x**active_degree, 0.0_dp]
        associate (unused_y => y)
            if (kind(unused_y) /= dp) error stop
        end associate
    end subroutine bdm_field

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_full_vector_sparse_assembly
