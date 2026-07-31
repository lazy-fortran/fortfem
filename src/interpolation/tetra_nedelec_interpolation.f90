module fortfem_tetra_nedelec_interpolation
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t
    use fortfem_tetra_vector_samples, only: &
        tetra_vector_sample_gradients_t, tetra_vector_samples_t, &
        zero_tetra_vector_samples_like
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: det3
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_special_jacobi, only: tetrahedron_koornwinder, &
        triangle_dubiner
    implicit none

    private

    public :: interpolate_reference_tetra_nedelec
    public :: interpolate_physical_tetra_nedelec
    public :: interpolate_sampled_physical_tetra_nedelec
    public :: interpolate_sampled_physical_tetra_nedelec_jvp
    public :: interpolate_sampled_physical_tetra_nedelec_vjp
    public :: tetra_nedelec_interpolation_points

    abstract interface
        pure subroutine physical_vector_field_3d(x, y, z, value)
            import :: dp
            real(dp), intent(in) :: x, y, z
            real(dp), intent(out) :: value(3)
        end subroutine physical_vector_field_3d

        pure subroutine reference_vector_field_3d(point, value)
            import :: dp
            real(dp), intent(in) :: point(3)
            real(dp), intent(out) :: value(3)
        end subroutine reference_vector_field_3d
    end interface

contains

    subroutine interpolate_physical_tetra_nedelec( &
            basis, vertices, field, dofs, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        procedure(physical_vector_field_3d) :: field
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        real(dp) :: determinant, jacobian(3, 3)

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        determinant = det3(jacobian)
        status = 1
        if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        call interpolate_reference_tetra_nedelec( &
            basis, pulled_back_field, dofs, status)

    contains

        pure subroutine pulled_back_field(reference_point, value)
            real(dp), intent(in) :: reference_point(3)
            real(dp), intent(out) :: value(3)

            real(dp) :: physical_point(3), physical_value(3)

            physical_point = vertices(:, 1) + &
                matmul(jacobian, reference_point)
            call field( &
                physical_point(1), physical_point(2), physical_point(3), &
                physical_value)
            value = matmul(transpose(jacobian), physical_value)
        end subroutine pulled_back_field

    end subroutine interpolate_physical_tetra_nedelec

    subroutine interpolate_reference_tetra_nedelec( &
            basis, field, dofs, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        procedure(reference_vector_field_3d) :: field
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: tetra_weights(:), triangle_weights(:)
        real(dp), allocatable :: x(:), y(:), z(:)
        real(dp) :: field_value(3), point(3), tangent(3), tangents(3, 2)
        integer :: component, edge, exponent, face, moment, node, order
        integer :: total_degree, x_degree, y_degree, z_degree

        status = 1
        order = order_from_dof_count(tetra_nedelec_dof_count(basis))
        if (order == 0) return
        allocate(dofs(tetra_nedelec_dof_count(basis)))
        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        dofs = 0.0_dp
        moment = 0
        do edge = 1, 6
            do exponent = 0, order - 1
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), point, tangent)
                    call field(point, field_value)
                    dofs(moment) = dofs(moment) + edge_weights(node) * &
                        shifted_legendre(exponent, edge_nodes(node)) * &
                        dot_product(field_value, tangent)
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * order + 4, x, y, triangle_weights, status)
        if (status /= 0) return
        do face = 1, 4
            do component = 1, 2
                do total_degree = 0, order - 2
                    do x_degree = 0, total_degree
                        y_degree = total_degree - x_degree
                        moment = moment + 1
                        do node = 1, size(x)
                            call reference_face( &
                                face, x(node), y(node), point, tangents)
                            call field(point, field_value)
                            dofs(moment) = dofs(moment) + &
                                triangle_weights(node) * &
                                face_moment_value( &
                                order, x_degree, y_degree, &
                                x(node), y(node))* &
                                dot_product( &
                                field_value, tangents(:, component))
                        end do
                    end do
                end do
            end do
        end do

        call tetra_duffy_quadrature( &
            2 * order + 4, x, y, z, tetra_weights, status)
        if (status /= 0) return
        do component = 1, 3
            do total_degree = 0, order - 3
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(x)
                            point = [x(node), y(node), z(node)]
                            call field(point, field_value)
                            dofs(moment) = dofs(moment) + &
                                tetra_weights(node)*volume_moment_value( &
                                order, x_degree, y_degree, z_degree, &
                                point)* &
                                field_value(component)
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= size(dofs)) return
        status = 0
    end subroutine interpolate_reference_tetra_nedelec

    subroutine tetra_nedelec_interpolation_points( &
            basis, vertices, points, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        type(tetra_vector_samples_t), intent(out) :: points
        integer, intent(out) :: status

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: jacobian(3, 3), reference_point(3), tangent(3)
        real(dp) :: tangents(3, 2)
        integer :: edge, face, node, order

        status = 1
        order = order_from_dof_count(tetra_nedelec_dof_count(basis))
        if (order == 0) return
        call tetrahedron_jacobian(vertices, jacobian)
        if (.not. valid_jacobian(jacobian)) return
        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        allocate(points%edge(3, size(edge_nodes), 6))
        do edge = 1, 6
            do node = 1, size(edge_nodes)
                call reference_edge( &
                    edge, edge_nodes(node), reference_point, tangent)
                points%edge(:, node, edge) = vertices(:, 1) + &
                    matmul(jacobian, reference_point)
            end do
        end do
        call triangle_duffy_quadrature( &
            2*order + 4, x, y, weights, status)
        if (status /= 0) return
        allocate(points%face(3, size(weights), 4))
        do face = 1, 4
            do node = 1, size(weights)
                call reference_face( &
                    face, x(node), y(node), reference_point, tangents)
                points%face(:, node, face) = vertices(:, 1) + &
                    matmul(jacobian, reference_point)
            end do
        end do
        call tetra_duffy_quadrature( &
            2*order + 4, x, y, z, weights, status)
        if (status /= 0) return
        allocate(points%volume(3, size(weights)))
        do node = 1, size(weights)
            reference_point = [x(node), y(node), z(node)]
            points%volume(:, node) = vertices(:, 1) + &
                matmul(jacobian, reference_point)
        end do
        status = 0
    end subroutine tetra_nedelec_interpolation_points

    subroutine interpolate_sampled_physical_tetra_nedelec( &
            basis, vertices, samples, dofs, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        type(tetra_vector_samples_t), intent(in) :: samples
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        type(tetra_vector_samples_t) :: points, pulled
        real(dp) :: jacobian(3, 3)
        integer :: edge, face, node

        call tetra_nedelec_interpolation_points( &
            basis, vertices, points, status)
        if (status /= 0) return
        if (.not. same_sample_shape(samples, points)) then
            status = 1
            return
        end if
        call tetrahedron_jacobian(vertices, jacobian)
        call zero_tetra_vector_samples_like(samples, pulled)
        do edge = 1, 6
            do node = 1, size(samples%edge, 2)
                pulled%edge(:, node, edge) = &
                    matmul(transpose(jacobian), samples%edge(:, node, edge))
            end do
        end do
        do face = 1, 4
            do node = 1, size(samples%face, 2)
                pulled%face(:, node, face) = &
                    matmul(transpose(jacobian), samples%face(:, node, face))
            end do
        end do
        do node = 1, size(samples%volume, 2)
            pulled%volume(:, node) = &
                matmul(transpose(jacobian), samples%volume(:, node))
        end do
        call interpolate_reference_tetra_nedelec_samples( &
            basis, pulled, dofs, status)
    end subroutine interpolate_sampled_physical_tetra_nedelec

    subroutine interpolate_sampled_physical_tetra_nedelec_jvp( &
            basis, vertices, samples, gradients, vertices_dot, parameter_dot, &
            dofs_dot, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4), vertices_dot(3, 4)
        type(tetra_vector_samples_t), intent(in) :: samples, parameter_dot
        type(tetra_vector_sample_gradients_t), intent(in) :: gradients
        real(dp), intent(out) :: dofs_dot(:)
        integer, intent(out) :: status

        type(tetra_vector_samples_t) :: points, pulled_dot
        real(dp) :: field_dot(3), jacobian(3, 3), jacobian_dot(3, 3)
        real(dp) :: point_dot(3), reference_point(3), tangent(3)
        real(dp) :: tangents(3, 2)
        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        integer :: edge, face, node, order

        dofs_dot = 0.0_dp
        call tetra_nedelec_interpolation_points( &
            basis, vertices, points, status)
        if (status /= 0) return
        if (.not. valid_product_shapes( &
            samples, gradients, parameter_dot, points)) then
            status = 1
            return
        end if
        if (size(dofs_dot) /= tetra_nedelec_dof_count(basis)) then
            status = 1
            return
        end if
        call tetrahedron_jacobian(vertices, jacobian)
        call tetrahedron_jacobian(vertices_dot, jacobian_dot)
        call zero_tetra_vector_samples_like(samples, pulled_dot)
        order = order_from_dof_count(tetra_nedelec_dof_count(basis))
        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        do edge = 1, 6
            do node = 1, size(edge_nodes)
                call reference_edge( &
                    edge, edge_nodes(node), reference_point, tangent)
                point_dot = vertices_dot(:, 1) + &
                    matmul(jacobian_dot, reference_point)
                field_dot = parameter_dot%edge(:, node, edge) + &
                    matmul(gradients%edge(:, :, node, edge), point_dot)
                pulled_dot%edge(:, node, edge) = &
                    matmul(transpose(jacobian_dot), &
                    samples%edge(:, node, edge)) + &
                    matmul(transpose(jacobian), field_dot)
            end do
        end do
        call triangle_duffy_quadrature( &
            2*order + 4, x, y, weights, status)
        if (status /= 0) return
        do face = 1, 4
            do node = 1, size(weights)
                call reference_face( &
                    face, x(node), y(node), reference_point, tangents)
                point_dot = vertices_dot(:, 1) + &
                    matmul(jacobian_dot, reference_point)
                field_dot = parameter_dot%face(:, node, face) + &
                    matmul(gradients%face(:, :, node, face), point_dot)
                pulled_dot%face(:, node, face) = &
                    matmul(transpose(jacobian_dot), &
                    samples%face(:, node, face)) + &
                    matmul(transpose(jacobian), field_dot)
            end do
        end do
        call tetra_duffy_quadrature( &
            2*order + 4, x, y, z, weights, status)
        if (status /= 0) return
        do node = 1, size(weights)
            reference_point = [x(node), y(node), z(node)]
            point_dot = vertices_dot(:, 1) + &
                matmul(jacobian_dot, reference_point)
            field_dot = parameter_dot%volume(:, node) + &
                matmul(gradients%volume(:, :, node), point_dot)
            pulled_dot%volume(:, node) = &
                matmul(transpose(jacobian_dot), samples%volume(:, node)) + &
                matmul(transpose(jacobian), field_dot)
        end do
        call interpolate_reference_tetra_nedelec_samples_into( &
            basis, pulled_dot, dofs_dot, status)
    end subroutine interpolate_sampled_physical_tetra_nedelec_jvp

    subroutine interpolate_sampled_physical_tetra_nedelec_vjp( &
            basis, vertices, samples, gradients, dofs_bar, vertices_bar, &
            samples_bar, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        type(tetra_vector_samples_t), intent(in) :: samples
        type(tetra_vector_sample_gradients_t), intent(in) :: gradients
        real(dp), intent(in) :: dofs_bar(:)
        real(dp), intent(out) :: vertices_bar(3, 4)
        type(tetra_vector_samples_t), intent(out) :: samples_bar
        integer, intent(out) :: status

        type(tetra_vector_samples_t) :: parameter_zero, points, pulled_bar
        real(dp) :: field_bar(3), jacobian(3, 3), jacobian_bar(3, 3)
        real(dp) :: point_bar(3), reference_point(3), tangent(3)
        real(dp) :: tangents(3, 2)
        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        integer :: edge, face, node, order

        vertices_bar = 0.0_dp
        call tetra_nedelec_interpolation_points( &
            basis, vertices, points, status)
        if (status /= 0) return
        call zero_tetra_vector_samples_like(samples, parameter_zero)
        if (.not. valid_product_shapes( &
            samples, gradients, parameter_zero, points)) then
            status = 1
            return
        end if
        if (size(dofs_bar) /= tetra_nedelec_dof_count(basis)) then
            status = 1
            return
        end if
        call interpolate_reference_tetra_nedelec_samples_vjp( &
            basis, samples, dofs_bar, pulled_bar, status)
        if (status /= 0) return
        call zero_tetra_vector_samples_like(samples, samples_bar)
        call tetrahedron_jacobian(vertices, jacobian)
        jacobian_bar = 0.0_dp
        order = order_from_dof_count(tetra_nedelec_dof_count(basis))
        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        do edge = 1, 6
            do node = 1, size(edge_nodes)
                call reference_edge( &
                    edge, edge_nodes(node), reference_point, tangent)
                call reverse_sample( &
                    jacobian, reference_point, samples%edge(:, node, edge), &
                    gradients%edge(:, :, node, edge), &
                    pulled_bar%edge(:, node, edge), field_bar, point_bar, &
                    jacobian_bar)
                samples_bar%edge(:, node, edge) = field_bar
                vertices_bar(:, 1) = vertices_bar(:, 1) + point_bar
            end do
        end do
        call triangle_duffy_quadrature( &
            2*order + 4, x, y, weights, status)
        if (status /= 0) return
        do face = 1, 4
            do node = 1, size(weights)
                call reference_face( &
                    face, x(node), y(node), reference_point, tangents)
                call reverse_sample( &
                    jacobian, reference_point, samples%face(:, node, face), &
                    gradients%face(:, :, node, face), &
                    pulled_bar%face(:, node, face), field_bar, point_bar, &
                    jacobian_bar)
                samples_bar%face(:, node, face) = field_bar
                vertices_bar(:, 1) = vertices_bar(:, 1) + point_bar
            end do
        end do
        call tetra_duffy_quadrature( &
            2*order + 4, x, y, z, weights, status)
        if (status /= 0) return
        do node = 1, size(weights)
            reference_point = [x(node), y(node), z(node)]
            call reverse_sample( &
                jacobian, reference_point, samples%volume(:, node), &
                gradients%volume(:, :, node), pulled_bar%volume(:, node), &
                field_bar, point_bar, jacobian_bar)
            samples_bar%volume(:, node) = field_bar
            vertices_bar(:, 1) = vertices_bar(:, 1) + point_bar
        end do
        call scatter_jacobian_bar(jacobian_bar, vertices_bar)
        status = 0
    end subroutine interpolate_sampled_physical_tetra_nedelec_vjp

    subroutine interpolate_reference_tetra_nedelec_samples( &
            basis, samples, dofs, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        type(tetra_vector_samples_t), intent(in) :: samples
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        allocate(dofs(tetra_nedelec_dof_count(basis)))
        call interpolate_reference_tetra_nedelec_samples_into( &
            basis, samples, dofs, status)
    end subroutine interpolate_reference_tetra_nedelec_samples

    subroutine interpolate_reference_tetra_nedelec_samples_into( &
            basis, samples, dofs, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        type(tetra_vector_samples_t), intent(in) :: samples
        real(dp), intent(out) :: dofs(:)
        integer, intent(out) :: status

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: point(3), tangent(3), tangents(3, 2)
        integer :: component, edge, exponent, face, moment, node, order
        integer :: total_degree, x_degree, y_degree, z_degree

        dofs = 0.0_dp
        status = 1
        order = order_from_dof_count(tetra_nedelec_dof_count(basis))
        if (order == 0) return
        if (size(dofs) /= tetra_nedelec_dof_count(basis)) return
        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        moment = 0
        do edge = 1, 6
            do exponent = 0, order - 1
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), point, tangent)
                    dofs(moment) = dofs(moment) + edge_weights(node)* &
                        shifted_legendre(exponent, edge_nodes(node))* &
                        dot_product(samples%edge(:, node, edge), tangent)
                end do
            end do
        end do
        call triangle_duffy_quadrature( &
            2*order + 4, x, y, weights, status)
        if (status /= 0) return
        do face = 1, 4
            do component = 1, 2
                do total_degree = 0, order - 2
                    do x_degree = 0, total_degree
                        y_degree = total_degree - x_degree
                        moment = moment + 1
                        do node = 1, size(weights)
                            call reference_face( &
                                face, x(node), y(node), point, tangents)
                            dofs(moment) = dofs(moment) + weights(node)* &
                                face_moment_value( &
                                order, x_degree, y_degree, x(node), y(node))* &
                                dot_product( &
                                samples%face(:, node, face), &
                                tangents(:, component))
                        end do
                    end do
                end do
            end do
        end do
        call tetra_duffy_quadrature( &
            2*order + 4, x, y, z, weights, status)
        if (status /= 0) return
        do component = 1, 3
            do total_degree = 0, order - 3
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(weights)
                            point = [x(node), y(node), z(node)]
                            dofs(moment) = dofs(moment) + weights(node)* &
                                volume_moment_value( &
                                order, x_degree, y_degree, z_degree, point)* &
                                samples%volume(component, node)
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= size(dofs)) return
        status = 0
    end subroutine interpolate_reference_tetra_nedelec_samples_into

    subroutine interpolate_reference_tetra_nedelec_samples_vjp( &
            basis, samples, dofs_bar, samples_bar, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        type(tetra_vector_samples_t), intent(in) :: samples
        real(dp), intent(in) :: dofs_bar(:)
        type(tetra_vector_samples_t), intent(out) :: samples_bar
        integer, intent(out) :: status

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: point(3), seed, tangent(3), tangents(3, 2)
        integer :: component, edge, exponent, face, moment, node, order
        integer :: total_degree, x_degree, y_degree, z_degree

        call zero_tetra_vector_samples_like(samples, samples_bar)
        status = 1
        order = order_from_dof_count(tetra_nedelec_dof_count(basis))
        if (order == 0) return
        if (size(dofs_bar) /= tetra_nedelec_dof_count(basis)) return
        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        moment = 0
        do edge = 1, 6
            do exponent = 0, order - 1
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), point, tangent)
                    seed = dofs_bar(moment)*edge_weights(node)* &
                        shifted_legendre(exponent, edge_nodes(node))
                    samples_bar%edge(:, node, edge) = &
                        samples_bar%edge(:, node, edge) + seed*tangent
                end do
            end do
        end do
        call triangle_duffy_quadrature( &
            2*order + 4, x, y, weights, status)
        if (status /= 0) return
        do face = 1, 4
            do component = 1, 2
                do total_degree = 0, order - 2
                    do x_degree = 0, total_degree
                        y_degree = total_degree - x_degree
                        moment = moment + 1
                        do node = 1, size(weights)
                            call reference_face( &
                                face, x(node), y(node), point, tangents)
                            seed = dofs_bar(moment)*weights(node)* &
                                face_moment_value( &
                                order, x_degree, y_degree, x(node), y(node))
                            samples_bar%face(:, node, face) = &
                                samples_bar%face(:, node, face) + &
                                seed*tangents(:, component)
                        end do
                    end do
                end do
            end do
        end do
        call tetra_duffy_quadrature( &
            2*order + 4, x, y, z, weights, status)
        if (status /= 0) return
        do component = 1, 3
            do total_degree = 0, order - 3
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(weights)
                            point = [x(node), y(node), z(node)]
                            seed = dofs_bar(moment)*weights(node)* &
                                volume_moment_value( &
                                order, x_degree, y_degree, z_degree, point)
                            samples_bar%volume(component, node) = &
                                samples_bar%volume(component, node) + seed
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= size(dofs_bar)) return
        status = 0
    end subroutine interpolate_reference_tetra_nedelec_samples_vjp

    pure subroutine reverse_sample( &
            jacobian, reference_point, sample, gradient, pulled_bar, &
            sample_bar, point_bar, jacobian_bar)
        real(dp), intent(in) :: jacobian(3, 3), reference_point(3)
        real(dp), intent(in) :: sample(3), gradient(3, 3), pulled_bar(3)
        real(dp), intent(out) :: sample_bar(3), point_bar(3)
        real(dp), intent(inout) :: jacobian_bar(3, 3)

        sample_bar = matmul(jacobian, pulled_bar)
        point_bar = matmul(transpose(gradient), sample_bar)
        jacobian_bar = jacobian_bar + &
            spread(sample, 2, 3)*spread(pulled_bar, 1, 3) + &
            spread(point_bar, 2, 3)*spread(reference_point, 1, 3)
    end subroutine reverse_sample

    pure subroutine scatter_jacobian_bar(jacobian_bar, vertices_bar)
        real(dp), intent(in) :: jacobian_bar(3, 3)
        real(dp), intent(inout) :: vertices_bar(3, 4)

        vertices_bar(:, 1) = vertices_bar(:, 1) - jacobian_bar(:, 1) - &
            jacobian_bar(:, 2) - jacobian_bar(:, 3)
        vertices_bar(:, 2) = vertices_bar(:, 2) + jacobian_bar(:, 1)
        vertices_bar(:, 3) = vertices_bar(:, 3) + jacobian_bar(:, 2)
        vertices_bar(:, 4) = vertices_bar(:, 4) + jacobian_bar(:, 3)
    end subroutine scatter_jacobian_bar

    pure subroutine tetrahedron_jacobian(vertices, jacobian)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3)

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
    end subroutine tetrahedron_jacobian

    pure logical function valid_jacobian(jacobian) result(valid)
        real(dp), intent(in) :: jacobian(3, 3)
        real(dp) :: determinant

        determinant = det3(jacobian)
        valid = determinant > 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)
    end function valid_jacobian

    pure logical function same_sample_shape(first, second) result(valid)
        type(tetra_vector_samples_t), intent(in) :: first, second

        valid = allocated(first%edge)
        if (.not. valid) return
        valid = allocated(first%face)
        if (.not. valid) return
        valid = allocated(first%volume)
        if (.not. valid) return
        valid = allocated(second%edge)
        if (.not. valid) return
        valid = allocated(second%face)
        if (.not. valid) return
        valid = allocated(second%volume)
        if (.not. valid) return
        valid = all(shape(first%edge) == shape(second%edge))
        if (.not. valid) return
        valid = all(shape(first%face) == shape(second%face))
        if (.not. valid) return
        valid = all(shape(first%volume) == shape(second%volume))
    end function same_sample_shape

    pure logical function valid_product_shapes( &
            samples, gradients, parameter_dot, points) result(valid)
        type(tetra_vector_samples_t), intent(in) :: samples, parameter_dot
        type(tetra_vector_samples_t), intent(in) :: points
        type(tetra_vector_sample_gradients_t), intent(in) :: gradients

        valid = same_sample_shape(samples, points)
        if (.not. valid) return
        valid = same_sample_shape(parameter_dot, points)
        if (.not. valid) return
        valid = allocated(gradients%edge)
        if (.not. valid) return
        valid = allocated(gradients%face)
        if (.not. valid) return
        valid = allocated(gradients%volume)
        if (.not. valid) return
        valid = all(shape(gradients%edge) == [ &
            3, 3, size(points%edge, 2), size(points%edge, 3)])
        if (.not. valid) return
        valid = all(shape(gradients%face) == [ &
            3, 3, size(points%face, 2), size(points%face, 3)])
        if (.not. valid) return
        valid = all(shape(gradients%volume) == [3, 3, size(points%volume, 2)])
    end function valid_product_shapes

    pure function order_from_dof_count(dof_count) result(order)
        integer, intent(in) :: dof_count
        integer :: order

        do order = 1, 1000
            if (dof_count == order * (order + 2) * (order + 3) / 2) return
            if (order * (order + 2) * (order + 3) / 2 > dof_count) exit
        end do
        order = 0
    end function order_from_dof_count

    pure real(dp) function face_moment_value( &
            order, first_degree, second_degree, x, y) result(value)
        integer, intent(in) :: order, first_degree, second_degree
        real(dp), intent(in) :: x, y

        if (order <= 5) then
            value = x**first_degree*y**second_degree
        else
            value = triangle_dubiner(first_degree, second_degree, x, y)
        end if
    end function face_moment_value

    pure real(dp) function volume_moment_value( &
            order, first_degree, second_degree, third_degree, point) &
            result(value)
        integer, intent(in) :: order
        integer, intent(in) :: first_degree, second_degree, third_degree
        real(dp), intent(in) :: point(3)

        if (order <= 5) then
            value = point(1)**first_degree*point(2)**second_degree* &
                point(3)**third_degree
        else
            value = tetrahedron_koornwinder( &
                first_degree, second_degree, third_degree, &
                point(1), point(2), point(3))
        end if
    end function volume_moment_value

    pure subroutine reference_edge(edge, parameter, point, tangent)
        integer, intent(in) :: edge
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: point(3), tangent(3)

        real(dp) :: vertices(3, 4)
        integer :: edge_vertices(2, 6), first, second

        call reference_topology(vertices, edge_vertices)
        first = edge_vertices(1, edge)
        second = edge_vertices(2, edge)
        point = (1.0_dp - parameter) * vertices(:, first) + &
            parameter * vertices(:, second)
        tangent = vertices(:, second) - vertices(:, first)
    end subroutine reference_edge

    pure subroutine reference_face(face, u, v, point, tangents)
        integer, intent(in) :: face
        real(dp), intent(in) :: u, v
        real(dp), intent(out) :: point(3), tangents(3, 2)

        real(dp) :: vertices(3, 4)
        integer :: edge_vertices(2, 6), face_vertices(3, 4)
        integer :: first, second, third

        call reference_topology(vertices, edge_vertices)
        face_vertices(:, 1) = [1, 2, 3]
        face_vertices(:, 2) = [1, 2, 4]
        face_vertices(:, 3) = [1, 3, 4]
        face_vertices(:, 4) = [2, 3, 4]
        first = face_vertices(1, face)
        second = face_vertices(2, face)
        third = face_vertices(3, face)
        tangents(:, 1) = vertices(:, second) - vertices(:, first)
        tangents(:, 2) = vertices(:, third) - vertices(:, first)
        point = vertices(:, first) + u * tangents(:, 1) + &
            v * tangents(:, 2)
    end subroutine reference_face

    pure subroutine reference_topology(vertices, edge_vertices)
        real(dp), intent(out) :: vertices(3, 4)
        integer, intent(out) :: edge_vertices(2, 6)

        vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
        edge_vertices(:, 1) = [1, 2]
        edge_vertices(:, 2) = [1, 3]
        edge_vertices(:, 3) = [1, 4]
        edge_vertices(:, 4) = [2, 3]
        edge_vertices(:, 5) = [2, 4]
        edge_vertices(:, 6) = [3, 4]
    end subroutine reference_topology

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

end module fortfem_tetra_nedelec_interpolation
