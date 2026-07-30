program curl_curl_example
    !! Verified arbitrary-order tetrahedral Nedelec curl-curl solve.
    use fortfem_api, only: solve_tetra_nedelec_curl_mass
    use fortfem_kinds, only: dp
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, tetra_nedelec_first_kind_t
    use fortfem_tetra_nedelec_global_dof_map, only: &
        build_tetra_nedelec_basis_transform, build_tetra_nedelec_dof_map
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortnum_linalg, only: det3
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: maximum_order = 5
    character(*), parameter :: output_directory = "output/example/curl_curl"
    type(fortsparse_status_t) :: status
    real(dp), allocatable :: solution(:)
    real(dp) :: errors(2, maximum_order), orders(maximum_order)
    real(dp) :: vertices(3, 8)
    integer :: order, tetrahedra(4, 6), unit

    call execute_command_line("mkdir -p "//output_directory)
    call build_cube_mesh(vertices, tetrahedra)
    do order = 1, maximum_order
        call solve_tetra_nedelec_curl_mass( &
            vertices, tetrahedra, order, manufactured_source, 1.0_dp, &
            1.0_dp, solution, status, zero_tangential_boundary=.true.)
        if (status%code /= 0) error stop "Nedelec curl-curl solve failed"
        call measure_error( &
            vertices, tetrahedra, order, solution, errors(:, order))
        orders(order) = real(order, dp)
        write (*, "(a,i0,a,2(es12.4,1x))") &
            "Nedelec order ", order, " field/curl errors ", errors(:, order)
        deallocate(solution)
    end do
    if (.not. all(errors(:, 2:) < errors(:, :maximum_order - 1))) &
        error stop "curl-curl p-convergence regression"

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot(orders, errors(1, :), label="field error", marker="o")
    call plot(orders, errors(2, :), label="curl error", marker="s")
    call xlabel("Nedelec polynomial order")
    call ylabel("maximum sampled error")
    call title("Manufactured three-dimensional curl-curl convergence")
    call legend()
    call savefig(output_directory//"/curl_curl_p_convergence.png")

    open (newunit=unit, file=output_directory//"/convergence.csv", &
        status="replace", action="write")
    write (unit, "(a)") "order,field_error,curl_error"
    do order = 1, maximum_order
        write (unit, "(*(es24.16,:,','))") orders(order), errors(:, order)
    end do
    close (unit)

contains

    pure subroutine manufactured_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        associate(unused => z)
        end associate
        value = [ &
            0.0_dp, 0.0_dp, &
            (2.0_dp*pi*pi + 1.0_dp)*sin(pi*x)*sin(pi*y)]
    end subroutine manufactured_source

    pure function exact_field(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = [ &
            0.0_dp, 0.0_dp, &
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

    subroutine measure_error( &
            mesh_vertices, cells, basis_order, coefficients, error)
        real(dp), intent(in) :: mesh_vertices(:, :), coefficients(:)
        integer, intent(in) :: cells(:, :), basis_order
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
            basis_order, cells, edges, faces, global_dofs, edge_orientations, &
            face_permutations, local_status)
        if (local_status /= 0) error stop "Nedelec DOF map failed"
        call initialize_tetra_nedelec_first_kind( &
            basis_order, basis, local_status)
        if (local_status /= 0) error stop "Nedelec basis failed"
        dof_count = size(global_dofs, 1)
        allocate( &
            basis_transform(dof_count, dof_count), local_dofs(dof_count), &
            reference_values(3, dof_count), reference_curls(3, dof_count), &
            physical_values(3, dof_count), physical_curls(3, dof_count))
        reference_point = [0.21_dp, 0.24_dp, 0.19_dp]
        error = 0.0_dp
        do tetrahedron = 1, size(cells, 2)
            call build_tetra_nedelec_basis_transform( &
                basis_order, edge_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) error stop "Nedelec transform failed"
            local_dofs = matmul( &
                basis_transform, coefficients(global_dofs(:, tetrahedron)))
            jacobian(:, 1) = &
                mesh_vertices(:, cells(2, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            jacobian(:, 2) = &
                mesh_vertices(:, cells(3, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            jacobian(:, 3) = &
                mesh_vertices(:, cells(4, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            point = mesh_vertices(:, cells(1, tetrahedron)) + &
                matmul(jacobian, reference_point)
            call evaluate_tetra_nedelec_first_kind( &
                basis, reference_point, reference_values, reference_curls, &
                local_status)
            if (local_status /= 0) error stop "Nedelec evaluation failed"
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, local_status)
            if (local_status /= 0) error stop "Nedelec Piola map failed"
            value = matmul(physical_values, local_dofs)
            curl_value = matmul(physical_curls, local_dofs)
            error(1) = max(error(1), maxval(abs(value - exact_field(point))))
            error(2) = max( &
                error(2), maxval(abs(curl_value - exact_curl(point))))
        end do
    end subroutine measure_error

end program curl_curl_example
