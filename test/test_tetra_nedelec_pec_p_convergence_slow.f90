program test_tetra_nedelec_pec_p_convergence_slow
    use check, only: check_condition, check_summary
    use fortfem_tetra_nedelec_solver_3d, only: solve_tetra_nedelec_weighted_curl_mass
    use fortfem_kinds, only: dp
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, tetra_nedelec_first_kind_t
    use fortfem_tetra_nedelec_global_dof_map, only: &
        build_tetra_nedelec_basis_transform, build_tetra_nedelec_dof_map
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortnum_linalg, only: det3
    use fortsparse, only: fortsparse_status_t
    implicit none

    type(fortsparse_status_t) :: status
    real(dp), allocatable :: solution(:)
    real(dp) :: errors(2, 5), vertices(3, 8)
    integer :: order, tetrahedra(4, 6)

    call build_cube_mesh(vertices, tetrahedra)
    do order = 1, 5
        call solve_tetra_nedelec_weighted_curl_mass( &
            vertices, tetrahedra, order, paper_reluctivity, magnetic_source, &
            1.0_dp, solution, status, &
            .true.)
        if (status%code /= 0) error stop "PEC tetra Hcurl solve failed"
        call measure_field_error( &
            vertices, tetrahedra, order, solution, errors(:, order))
        write(*, '(a,i0,a,2(es12.4,1x))') &
            "PEC Nedelec order ", order, " field/curl errors ", &
            errors(:, order)
        deallocate(solution)
    end do

    call check_condition(all(errors(:, 2:5) < errors(:, 1:4)), &
        "PEC magnetic field and curl errors decrease with polynomial order")
    call check_condition(errors(1, 5) < 1.0e-2_dp, &
        "Order-five PEC solve reaches the analytical vector field")
    call check_condition(errors(2, 5) < 1.5e-1_dp, &
        "Order-five PEC solve reaches the analytical magnetic curl")
    call check_summary("Tetrahedral PEC magnetic p-convergence")

contains

    pure subroutine paper_reluctivity(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3, 3)

        associate(unused => x + y)
        end associate
        value = 0.0_dp
        value(1, 1) = 1.0_dp
        value(2, 2) = 1.0_dp
        value(3, 3) = z
    end subroutine paper_reluctivity

    pure subroutine magnetic_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        associate(unused => z)
        end associate
        value = [0.0_dp, 0.0_dp, &
            (2.0_dp*pi*pi + 1.0_dp)*sin(pi*x)*sin(pi*y)]
    end subroutine magnetic_source

    pure function exact_field(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = [0.0_dp, 0.0_dp, &
            sin(pi*point(1))*sin(pi*point(2))]
    end function exact_field

    pure function exact_curl(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = [ &
            pi*sin(pi*point(1))*cos(pi*point(2)), &
            -pi*cos(pi*point(1))*sin(pi*point(2)), 0.0_dp]
    end function exact_curl

    subroutine build_cube_mesh(mesh_vertices, cells)
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

    subroutine measure_field_error( &
            mesh_vertices, cells, order, coefficients, error)
        real(dp), intent(in) :: mesh_vertices(:, :), coefficients(:)
        integer, intent(in) :: cells(:, :), order
        real(dp), intent(out) :: error(2)

        type(tetra_nedelec_first_kind_t) :: basis
        integer, allocatable :: edge_orientations(:, :), edges(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_transform(:, :), local_dofs(:)
        real(dp), allocatable :: physical_curls(:, :), physical_values(:, :)
        real(dp), allocatable :: reference_curls(:, :)
        real(dp), allocatable :: reference_values(:, :)
        real(dp) :: curl_value(3), jacobian(3, 3), point(3)
        real(dp) :: reference_point(3), value(3)
        integer :: dof_count, local_status, tetrahedron

        call build_tetra_nedelec_dof_map( &
            order, cells, edges, faces, global_dofs, edge_orientations, &
            face_permutations, local_status)
        if (local_status /= 0) error stop "PEC Hcurl map failed"
        call initialize_tetra_nedelec_first_kind(order, basis, local_status)
        if (local_status /= 0) error stop "PEC Hcurl basis failed"
        dof_count = size(global_dofs, 1)
        allocate(basis_transform(dof_count, dof_count))
        allocate(local_dofs(dof_count))
        allocate(reference_values(3, dof_count), reference_curls(3, dof_count))
        allocate(physical_values(3, dof_count), physical_curls(3, dof_count))
        reference_point = [0.21_dp, 0.24_dp, 0.19_dp]
        error = 0.0_dp
        do tetrahedron = 1, size(cells, 2)
            call build_tetra_nedelec_basis_transform( &
                order, edge_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) error stop "PEC Hcurl transform failed"
            local_dofs = matmul( &
                basis_transform, coefficients(global_dofs(:, tetrahedron)))
            jacobian(:, 1) = mesh_vertices(:, cells(2, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            jacobian(:, 2) = mesh_vertices(:, cells(3, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            jacobian(:, 3) = mesh_vertices(:, cells(4, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            call evaluate_tetra_nedelec_first_kind( &
                basis, reference_point, reference_values, reference_curls, &
                local_status)
            if (local_status /= 0) error stop "PEC Hcurl evaluation failed"
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, local_status)
            if (local_status /= 0) error stop "PEC Hcurl Piola map failed"
            point = mesh_vertices(:, cells(1, tetrahedron)) + &
                matmul(jacobian, reference_point)
            value = matmul(physical_values, local_dofs)
            curl_value = matmul(physical_curls, local_dofs)
            error(1) = max( &
                error(1), maxval(abs(value - exact_field(point))))
            error(2) = max( &
                error(2), maxval(abs(curl_value - exact_curl(point))))
        end do
    end subroutine measure_field_error

end program test_tetra_nedelec_pec_p_convergence_slow
