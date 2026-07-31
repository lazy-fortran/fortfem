program test_tetra_nedelec_sparse_assembly_slow
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_tetra_nedelec_curl_mass_csc, &
        build_tetra_edge_dof_map, build_tetra_nedelec_dof_map
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    integer, allocatable :: edges(:, :), global_dofs(:, :), orientations(:, :)
    integer, allocatable :: face_permutations(:, :, :), faces(:, :)
    integer :: edge, polynomial_order, status, tetrahedra(4, 2)
    real(dp), allocatable :: dofs(:), matrix_times_dofs(:)
    real(dp) :: constant_field(3), edge_midpoint(3), edge_vector(3), energy
    real(dp) :: curl_field(3), vertices(3, 5)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    vertices(:, 5) = [0.0_dp, 0.0_dp, -1.0_dp]
    tetrahedra(:, 1) = [1, 2, 3, 4]
    tetrahedra(:, 2) = [1, 3, 2, 5]
    call build_tetra_edge_dof_map( &
        tetrahedra, edges, global_dofs, orientations, status)
    allocate(dofs(size(edges, 2)), matrix_times_dofs(size(edges, 2)))

    constant_field = [1.0_dp, 2.0_dp, -1.0_dp]
    do edge = 1, size(edges, 2)
        edge_vector = vertices(:, edges(2, edge)) - &
            vertices(:, edges(1, edge))
        dofs(edge) = dot_product(constant_field, edge_vector)
    end do
    call assemble_tetra_nedelec_curl_mass_csc( &
        vertices, tetrahedra, matrix, sparse_status, 4.0_dp, 1.5_dp)
    matrix_times_dofs = csc_matvec(matrix, dofs)
    energy = dot_product(dofs, matrix_times_dofs)
    call record_condition(sparse_status%code == 0 .and. &
        abs(energy - 3.0_dp) < 2.0e-12_dp, &
        "Tetrahedral mass assembly has exact constant-field energy")

    curl_field = [1.0_dp, -2.0_dp, 3.0_dp]
    do edge = 1, size(edges, 2)
        edge_vector = vertices(:, edges(2, edge)) - &
            vertices(:, edges(1, edge))
        edge_midpoint = 0.5_dp * ( &
            vertices(:, edges(1, edge)) + vertices(:, edges(2, edge)))
        dofs(edge) = dot_product( &
            0.5_dp * cross_product(curl_field, edge_midpoint), edge_vector)
    end do
    call assemble_tetra_nedelec_curl_mass_csc( &
        vertices, tetrahedra, matrix, sparse_status, 2.0_dp, 0.0_dp)
    matrix_times_dofs = csc_matvec(matrix, dofs)
    energy = dot_product(dofs, matrix_times_dofs)
    call record_condition(abs(energy - 28.0_dp / 3.0_dp) < 3.0e-12_dp, &
        "Tetrahedral curl assembly has exact rotational-field energy")

    deallocate(dofs, matrix_times_dofs)
    constant_field = [1.0_dp, 2.0_dp, -1.0_dp]
    do polynomial_order = 2, 6
        call build_tetra_nedelec_dof_map( &
            polynomial_order, tetrahedra, edges, faces, global_dofs, &
            orientations, face_permutations, status)
        allocate( &
            dofs(maxval(global_dofs)), &
            matrix_times_dofs(maxval(global_dofs)))
        call set_constant_field_dofs( &
            polynomial_order, vertices, tetrahedra, edges, faces, &
            global_dofs, constant_field, dofs)
        call assemble_tetra_nedelec_curl_mass_csc( &
            vertices, tetrahedra, matrix, sparse_status, 4.0_dp, 1.5_dp, &
            order=polynomial_order)
        matrix_times_dofs = csc_matvec(matrix, dofs)
        energy = dot_product(dofs, matrix_times_dofs)
        print "(a,i0,a,es12.4)", &
            "Nedelec order ", polynomial_order, " constant energy error: ", &
            abs(energy - 3.0_dp)
        call record_condition(sparse_status%code == 0 .and. &
            abs(energy - 3.0_dp) < 1.0e-8_dp, &
            "Higher-order sparse mass assembly has exact constant-field energy")
        deallocate(dofs, matrix_times_dofs)
    end do

    call check_summary("Sparse tetrahedral Nedelec assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine set_constant_field_dofs( &
            order, mesh_vertices, cells, mesh_edges, mesh_faces, &
            cell_dofs, field, dofs)
        integer, intent(in) :: order, cells(:, :)
        real(dp), intent(in) :: mesh_vertices(:, :), field(3)
        integer, intent(in) :: mesh_edges(:, :), mesh_faces(:, :)
        integer, intent(in) :: cell_dofs(:, :)
        real(dp), intent(out) :: dofs(:)

        real(dp) :: jacobian(3, 3), reference_field(3), tangent(3)
        integer :: cell, component, degree, dof, edge_index, face
        integer :: face_dof_count, face_offset, local_dof, monomial
        integer :: x_degree, y_degree, z_degree

        dofs = 0.0_dp
        do edge_index = 1, size(mesh_edges, 2)
            tangent = mesh_vertices(:, mesh_edges(2, edge_index)) - &
                mesh_vertices(:, mesh_edges(1, edge_index))
            dofs((edge_index - 1) * order + 1) = dot_product(field, tangent)
        end do

        face_dof_count = order * (order - 1)
        face_offset = order * size(mesh_edges, 2)
        do face = 1, size(mesh_faces, 2)
            do component = 1, 2
                tangent = mesh_vertices( &
                    :, mesh_faces(component + 1, face)) - &
                    mesh_vertices(:, mesh_faces(1, face))
                monomial = 0
                do degree = 0, order - 2
                    do x_degree = 0, degree
                        y_degree = degree - x_degree
                        monomial = monomial + 1
                        dof = face_offset + (face - 1) * face_dof_count + &
                            (component - 1) * face_dof_count / 2 + monomial
                        if (order <= 5) then
                            dofs(dof) = dot_product(field, tangent)* &
                                triangle_monomial_integral( &
                                x_degree, y_degree)
                        else
                            dofs(dof) = merge( &
                                0.5_dp*dot_product(field, tangent), &
                                0.0_dp, degree == 0)
                        end if
                    end do
                end do
            end do
        end do

        local_dof = 6 * order + 4 * face_dof_count
        do cell = 1, size(cells, 2)
            jacobian(:, 1) = mesh_vertices(:, cells(2, cell)) - &
                mesh_vertices(:, cells(1, cell))
            jacobian(:, 2) = mesh_vertices(:, cells(3, cell)) - &
                mesh_vertices(:, cells(1, cell))
            jacobian(:, 3) = mesh_vertices(:, cells(4, cell)) - &
                mesh_vertices(:, cells(1, cell))
            reference_field = matmul(transpose(jacobian), field)
            do component = 1, 3
                do degree = 0, order - 3
                    do x_degree = 0, degree
                        do y_degree = 0, degree - x_degree
                            z_degree = degree - x_degree - y_degree
                            local_dof = local_dof + 1
                            if (order <= 5) then
                                dofs(cell_dofs(local_dof, cell)) = &
                                    reference_field(component)* &
                                    tetra_monomial_integral( &
                                    x_degree, y_degree, z_degree)
                            else
                                dofs(cell_dofs(local_dof, cell)) = merge( &
                                    reference_field(component)/6.0_dp, &
                                    0.0_dp, degree == 0)
                            end if
                        end do
                    end do
                end do
            end do
            local_dof = 6 * order + 4 * face_dof_count
        end do
    end subroutine set_constant_field_dofs

    pure function triangle_monomial_integral( &
            x_degree, y_degree) result(value)
        integer, intent(in) :: x_degree, y_degree
        real(dp) :: value

        value = real( &
            factorial(x_degree) * factorial(y_degree), dp) / &
            real(factorial(x_degree + y_degree + 2), dp)
    end function triangle_monomial_integral

    pure function tetra_monomial_integral( &
            x_degree, y_degree, z_degree) result(value)
        integer, intent(in) :: x_degree, y_degree, z_degree
        real(dp) :: value

        value = real( &
            factorial(x_degree) * factorial(y_degree) * &
            factorial(z_degree), dp) / &
            real(factorial(x_degree + y_degree + z_degree + 3), dp)
    end function tetra_monomial_integral

    pure function factorial(argument) result(value)
        integer, intent(in) :: argument
        integer :: value

        integer :: factor

        value = 1
        do factor = 2, argument
            value = value * factor
        end do
    end function factorial

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product(1) = first(2) * second(3) - first(3) * second(2)
        product(2) = first(3) * second(1) - first(1) * second(3)
        product(3) = first(1) * second(2) - first(2) * second(1)
    end function cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_sparse_assembly_slow
