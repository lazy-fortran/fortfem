module fortfem_tetra_nedelec_interpolation
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: det3
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    private

    public :: interpolate_reference_tetra_nedelec
    public :: interpolate_physical_tetra_nedelec

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
                                x(node)**x_degree * y(node)**y_degree * &
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
                                tetra_weights(node) * x(node)**x_degree * &
                                y(node)**y_degree * z(node)**z_degree * &
                                field_value(component)
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= size(dofs)) return
        status = 0
    end subroutine interpolate_reference_tetra_nedelec

    pure function order_from_dof_count(dof_count) result(order)
        integer, intent(in) :: dof_count
        integer :: order

        do order = 1, 5
            if (dof_count == order * (order + 2) * (order + 3) / 2) return
        end do
        order = 0
    end function order_from_dof_count

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
