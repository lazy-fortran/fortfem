program test_triangle_rt_arbitrary_order_sparse
    use check, only: check_condition, check_summary
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    use fortfem_api, only: &
        assemble_triangle_rt_div_mass_csc, build_triangle_trimmed_dof_map, &
        triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    type(mesh_2d_t) :: mesh
    type(csc_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    real(dp), allocatable :: dofs(:), local_dofs(:), matrix_times_dofs(:)
    integer, allocatable :: global_dofs(:, :), transforms(:, :)
    logical, allocatable :: assigned(:)
    real(dp) :: exact_energy, vertices(2, 3)
    integer :: degree, dof, global_dof_count, local_status, triangle
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

    do degree = 0, 3
        call build_triangle_trimmed_dof_map( &
            mesh, degree + 1, global_dofs, transforms, global_dof_count, &
            local_status)
        allocate(dofs(global_dof_count), assigned(global_dof_count))
        allocate(local_dofs(size(global_dofs, 1)))
        dofs = 0.0_dp
        assigned = .false.

        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call physical_polynomial_rt_dofs(vertices, degree, local_dofs)
            do dof = 1, size(local_dofs)
                if (assigned(global_dofs(dof, triangle))) then
                    call record_condition(abs( &
                        dofs(global_dofs(dof, triangle)) - &
                        real(transforms(dof, triangle), dp) * &
                        local_dofs(dof)) < 2.0e-11_dp, &
                        "Neighbouring RT elements agree on shared flux moments")
                else
                    dofs(global_dofs(dof, triangle)) = &
                        real(transforms(dof, triangle), dp) * local_dofs(dof)
                    assigned(global_dofs(dof, triangle)) = .true.
                end if
            end do
        end do

        call assemble_triangle_rt_div_mass_csc( &
            mesh, degree, 2 * degree + 2, matrix, sparse_status)
        allocate(matrix_times_dofs(global_dof_count))
        matrix_times_dofs = csc_matvec(matrix, dofs)
        exact_energy = 1.0_dp / real(2 * degree + 1, dp)
        if (degree > 0) then
            exact_energy = exact_energy + &
                real(degree * degree, dp) / real(2 * degree - 1, dp)
        end if
        call record_condition(sparse_status%code == 0 .and. &
            matrix%nrow == global_dof_count .and. &
            matrix%ncol == global_dof_count, &
            "Sparse arbitrary-degree RT operator has exact global dimensions")
        call record_condition(abs( &
            dot_product(dofs, matrix_times_dofs) - exact_energy) < &
            8.0e-8_dp, &
            "Sparse RT operator reproduces exact square polynomial energy")

        deallocate( &
            assigned, dofs, global_dofs, local_dofs, matrix_times_dofs, &
            transforms)
    end do

    call assemble_triangle_rt_div_mass_csc( &
        mesh, -1, 2, matrix, sparse_status)
    call record_condition(sparse_status%code /= 0, &
        "Sparse arbitrary-degree RT assembly rejects negative degree")

    call check_summary("Arbitrary-order triangle RT sparse assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine physical_polynomial_rt_dofs(vertices, degree, dofs)
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: degree
        real(dp), intent(out) :: dofs(:)

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: eta(:), triangle_weights(:), xi(:)
        real(dp) :: edge_point(2), jacobian(2, 2), normal(2)
        real(dp) :: physical_x, reference_field(2)
        integer :: component, edge, exponent, moment, node, status
        integer :: total_degree, x_degree, y_degree

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        allocate(edge_nodes(degree + 3), edge_weights(degree + 3))
        call gauss_legendre_ab( &
            degree + 3, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        dofs = 0.0_dp
        moment = 0
        do edge = 1, 3
            do exponent = 0, degree
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, normal)
                    call reference_polynomial_field( &
                        vertices, jacobian, edge_point(1), edge_point(2), &
                        degree, physical_x, reference_field)
                    dofs(moment) = dofs(moment) + edge_weights(node) * &
                        shifted_legendre(exponent, edge_nodes(node)) * &
                        dot_product(reference_field, normal)
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * degree + 1, xi, eta, triangle_weights, status)
        do component = 1, 2
            do total_degree = 0, degree - 1
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    moment = moment + 1
                    do node = 1, size(xi)
                        call reference_polynomial_field( &
                            vertices, jacobian, xi(node), eta(node), &
                            degree, physical_x, reference_field)
                        dofs(moment) = dofs(moment) + &
                            triangle_weights(node) * &
                            reference_field(component) * &
                            xi(node)**x_degree * eta(node)**y_degree
                    end do
                end do
            end do
        end do
    end subroutine physical_polynomial_rt_dofs

    pure subroutine reference_polynomial_field( &
            vertices, jacobian, xi, eta, degree, physical_x, reference_field)
        real(dp), intent(in) :: vertices(2, 3), jacobian(2, 2)
        real(dp), intent(in) :: xi, eta
        integer, intent(in) :: degree
        real(dp), intent(out) :: physical_x, reference_field(2)

        real(dp) :: field_value

        physical_x = vertices(1, 1) + &
            jacobian(1, 1) * xi + jacobian(1, 2) * eta
        field_value = physical_x**degree
        reference_field(1) = jacobian(2, 2) * field_value
        reference_field(2) = -jacobian(2, 1) * field_value
    end subroutine reference_polynomial_field

    pure subroutine reference_edge(edge, parameter, point, normal)
        integer, intent(in) :: edge
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: point(2), normal(2)

        select case (edge)
        case (1)
            point = [parameter, 0.0_dp]
            normal = [0.0_dp, -1.0_dp]
        case (2)
            point = [1.0_dp - parameter, parameter]
            normal = [1.0_dp, 1.0_dp]
        case (3)
            point = [0.0_dp, 1.0_dp - parameter]
            normal = [-1.0_dp, 0.0_dp]
        end select
    end subroutine reference_edge

    pure function shifted_legendre(degree, parameter) result(value)
        integer, intent(in) :: degree
        real(dp), intent(in) :: parameter
        real(dp) :: value

        real(dp) :: current, previous, coordinate
        integer :: polynomial_degree

        coordinate = 2.0_dp * parameter - 1.0_dp
        if (degree == 0) then
            value = 1.0_dp
            return
        end if
        previous = 1.0_dp
        current = coordinate
        do polynomial_degree = 1, degree - 1
            value = ( &
                real(2 * polynomial_degree + 1, dp) * coordinate * current - &
                real(polynomial_degree, dp) * previous) / &
                real(polynomial_degree + 1, dp)
            previous = current
            current = value
        end do
        value = current
    end function shifted_legendre

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_rt_arbitrary_order_sparse
