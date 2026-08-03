program test_curl_curl_solution_field
    use check, only: check_condition, check_summary
    use fortfem_core, only: invert_tetra_affine_map
    use fortfem_feec, only: &
        evaluate_tetra_nedelec_interpolant_at_point, &
        initialize_tetra_nedelec_first_kind, solve_tetra_nedelec_curl_mass, &
        tetra_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    use fortfem_tetra_nedelec_global_dof_map, only: &
        build_tetra_nedelec_basis_transform, build_tetra_nedelec_dof_map
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp) :: vertices(3, 8), point(3), reference(3)
    integer :: tetrahedra(4, 6)
    integer, allocatable :: edges(:, :), faces(:, :), global_dofs(:, :)
    integer, allocatable :: edge_orientations(:, :)
    integer, allocatable :: face_permutations(:, :, :)
    real(dp), allocatable :: solution(:), basis_transform(:, :), local_dofs(:)
    real(dp) :: value(3), curl(3), expected(3)
    type(tetra_nedelec_first_kind_t) :: basis
    type(fortsparse_status_t) :: sparse_status
    integer :: cell, status
    logical :: found

    call build_cube_mesh(vertices, tetrahedra)
    call solve_tetra_nedelec_curl_mass( &
        vertices, tetrahedra, 5, manufactured_source, 1.0_dp, 1.0_dp, &
        solution, sparse_status, zero_tangential_boundary=.true.)
    if (sparse_status%code /= 0) error stop "curl-curl solve oracle failed"
    call check_condition( &
        sparse_status%code == 0, "curl-curl manufactured solve succeeds")

    call initialize_tetra_nedelec_first_kind(5, basis, status)
    if (status /= 0) error stop "curl-curl basis oracle failed"
    call check_condition(status == 0, "curl-curl basis initializes")
    call build_tetra_nedelec_dof_map( &
        5, tetrahedra, edges, faces, global_dofs, edge_orientations, &
        face_permutations, status)
    if (status /= 0) error stop "curl-curl map oracle failed"
    call check_condition(status == 0, "curl-curl global map builds")
    allocate( &
        basis_transform(size(global_dofs, 1), size(global_dofs, 1)), &
        local_dofs(size(global_dofs, 1)))

    point = [0.31_dp, 0.27_dp, 0.42_dp]
    found = .false.
    do cell = 1, size(tetrahedra, 2)
        call reference_coordinates( &
            vertices(:, tetrahedra(:, cell)), point, reference, status)
        if (status /= 0) cycle
        if (minval(reference) < -1.0e-12_dp .or. &
            sum(reference) > 1.0_dp + 1.0e-12_dp) cycle
        call build_tetra_nedelec_basis_transform( &
            5, edge_orientations(:, cell), face_permutations(:, :, cell), &
            basis_transform, status)
        if (status /= 0) error stop "curl-curl orientation oracle failed"
        call check_condition( &
            status == 0, "curl-curl local orientation map builds")
        local_dofs = matmul( &
            basis_transform, solution(global_dofs(:, cell)))
        call evaluate_tetra_nedelec_interpolant_at_point( &
            vertices(:, tetrahedra(:, cell)), basis, local_dofs, point, &
            value, curl, status)
        if (status /= 0) error stop "curl-curl field oracle failed"
        call check_condition(status == 0, "curl-curl solved field evaluates")
        expected = exact_field(point)
        call check_condition( &
            maxval(abs(value - expected)) < 0.2_dp, &
            "curl-curl solved field follows manufactured solution")
        if (maxval(abs(value - expected)) >= 0.2_dp) &
            error stop "curl-curl solved field oracle mismatch"
        found = .true.
        exit
    end do
    call check_condition(found, "curl-curl oracle point lies in the mesh")
    if (.not. found) error stop "curl-curl oracle point is outside the mesh"
    call check_summary("curl-curl solved field oracle")

contains

    pure subroutine manufactured_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        associate(unused => z)
        end associate
        value = [0.0_dp, 0.0_dp, &
            (2.0_dp*pi*pi + 1.0_dp)*sin(pi*x)*sin(pi*y)]
    end subroutine manufactured_source

    pure function exact_field(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = [0.0_dp, 0.0_dp, &
            sin(pi*point(1))*sin(pi*point(2))]
    end function exact_field

    subroutine reference_coordinates(cell_vertices, physical_point, &
            local_reference, local_status)
        real(dp), intent(in) :: cell_vertices(3, 4), physical_point(3)
        real(dp), intent(out) :: local_reference(3)
        integer, intent(out) :: local_status

        call invert_tetra_affine_map( &
            cell_vertices, physical_point, local_reference, local_status)
    end subroutine reference_coordinates

    subroutine build_cube_mesh(mesh_vertices, cells)
        use fortnum_linalg, only: det3

        real(dp), intent(out) :: mesh_vertices(3, 8)
        integer, intent(out) :: cells(4, 6)
        real(dp) :: jacobian(3, 3)
        integer :: cell, temporary

        mesh_vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        mesh_vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        mesh_vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        mesh_vertices(:, 4) = [1.0_dp, 1.0_dp, 0.0_dp]
        mesh_vertices(:, 5) = [0.0_dp, 0.0_dp, 1.0_dp]
        mesh_vertices(:, 6) = [1.0_dp, 0.0_dp, 1.0_dp]
        mesh_vertices(:, 7) = [0.0_dp, 1.0_dp, 1.0_dp]
        mesh_vertices(:, 8) = [1.0_dp, 1.0_dp, 1.0_dp]
        cells(:, 1) = [1, 2, 4, 8]
        cells(:, 2) = [1, 2, 6, 8]
        cells(:, 3) = [1, 3, 4, 8]
        cells(:, 4) = [1, 3, 7, 8]
        cells(:, 5) = [1, 5, 6, 8]
        cells(:, 6) = [1, 5, 7, 8]
        do cell = 1, size(cells, 2)
            jacobian(:, 1) = &
                mesh_vertices(:, cells(2, cell)) - &
                mesh_vertices(:, cells(1, cell))
            jacobian(:, 2) = &
                mesh_vertices(:, cells(3, cell)) - &
                mesh_vertices(:, cells(1, cell))
            jacobian(:, 3) = &
                mesh_vertices(:, cells(4, cell)) - &
                mesh_vertices(:, cells(1, cell))
            if (det3(jacobian) < 0.0_dp) then
                temporary = cells(3, cell)
                cells(3, cell) = cells(4, cell)
                cells(4, cell) = temporary
            end if
        end do
    end subroutine build_cube_mesh

end program test_curl_curl_solution_field
