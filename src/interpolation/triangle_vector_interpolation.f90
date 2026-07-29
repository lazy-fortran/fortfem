module fortfem_triangle_vector_interpolation
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_nedelec_arbitrary_order, only: &
        evaluate_triangle_nedelec_first_kind, triangle_nedelec_first_kind_t
    use fortfem_triangle_piola_maps, only: &
        map_triangle_nedelec_covariant, map_triangle_rt_contravariant
    use fortfem_triangle_rt_arbitrary_order, only: &
        evaluate_triangle_raviart_thomas, triangle_rt_basis_t
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    private

    public :: interpolate_triangle_nedelec
    public :: evaluate_triangle_nedelec_interpolant
    public :: evaluate_triangle_rt_interpolant
    public :: interpolate_triangle_rt

    abstract interface
        subroutine vector_field_2d(x, y, value)
            import :: dp
            real(dp), intent(in) :: x, y
            real(dp), intent(out) :: value(2)
        end subroutine vector_field_2d
    end interface

contains

    subroutine evaluate_triangle_nedelec_interpolant( &
            vertices, basis, dofs, xi, eta, value, curl, status)
        real(dp), intent(in) :: vertices(2, 3)
        type(triangle_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: dofs(:), xi, eta
        real(dp), intent(out) :: value(2), curl
        integer, intent(out) :: status

        real(dp), allocatable :: physical_curls(:), physical_values(:, :)
        real(dp), allocatable :: reference_curls(:), reference_values(:, :)
        real(dp) :: jacobian(2, 2)
        integer :: dof_count

        value = 0.0_dp
        curl = 0.0_dp
        status = 1
        dof_count = size(dofs)
        if (dof_count < 1) return
        allocate(reference_values(2, dof_count), reference_curls(dof_count))
        allocate(physical_values(2, dof_count), physical_curls(dof_count))
        call evaluate_triangle_nedelec_first_kind( &
            basis, xi, eta, reference_values, reference_curls, status)
        if (status /= 0) return
        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        call map_triangle_nedelec_covariant( &
            jacobian, reference_values, reference_curls, physical_values, &
            physical_curls, status)
        if (status /= 0) return
        value = matmul(physical_values, dofs)
        curl = dot_product(physical_curls, dofs)
        status = 0
    end subroutine evaluate_triangle_nedelec_interpolant

    subroutine evaluate_triangle_rt_interpolant( &
            vertices, basis, dofs, xi, eta, value, divergence, status)
        real(dp), intent(in) :: vertices(2, 3)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: dofs(:), xi, eta
        real(dp), intent(out) :: value(2), divergence
        integer, intent(out) :: status

        real(dp), allocatable :: physical_divergences(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :)
        real(dp) :: jacobian(2, 2)
        integer :: dof_count

        value = 0.0_dp
        divergence = 0.0_dp
        status = 1
        dof_count = size(dofs)
        if (dof_count < 1) return
        allocate(reference_values(2, dof_count))
        allocate(reference_divergences(dof_count))
        allocate(physical_values(2, dof_count))
        allocate(physical_divergences(dof_count))
        call evaluate_triangle_raviart_thomas( &
            basis, xi, eta, reference_values, reference_divergences, status)
        if (status /= 0) return
        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        call map_triangle_rt_contravariant( &
            jacobian, reference_values, reference_divergences, &
            physical_values, physical_divergences, status)
        if (status /= 0) return
        value = matmul(physical_values, dofs)
        divergence = dot_product(physical_divergences, dofs)
        status = 0
    end subroutine evaluate_triangle_rt_interpolant

    subroutine interpolate_triangle_nedelec( &
            vertices, order, quadrature_degree, field, dofs, status)
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: order, quadrature_degree
        procedure(vector_field_2d) :: field
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: eta(:), triangle_weights(:), xi(:)
        real(dp) :: determinant, edge_point(2), jacobian(2, 2)
        real(dp) :: physical_field(2), physical_point(2), polynomial
        real(dp) :: reference_field(2), tangent(2)
        integer :: component, edge, exponent, moment, node, node_count
        integer :: total_degree, x_degree, y_degree

        status = 1
        if (order < 1 .or. quadrature_degree < 0) return
        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        determinant = jacobian(1, 1) * jacobian(2, 2) - &
            jacobian(1, 2) * jacobian(2, 1)
        if (determinant <= 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, maxval(abs(jacobian))**2)) return

        allocate(dofs(order * (order + 2)))
        dofs = 0.0_dp
        node_count = max(order + 1, (quadrature_degree + 2) / 2)
        allocate(edge_nodes(node_count), edge_weights(node_count))
        call gauss_legendre_ab( &
            node_count, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)

        moment = 0
        do edge = 1, 3
            do exponent = 0, order - 1
                moment = moment + 1
                do node = 1, node_count
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, tangent)
                    physical_point = vertices(:, 1) + &
                        matmul(jacobian, edge_point)
                    call field( &
                        physical_point(1), physical_point(2), physical_field)
                    reference_field = &
                        matmul(transpose(jacobian), physical_field)
                    polynomial = &
                        shifted_legendre(exponent, edge_nodes(node))
                    dofs(moment) = dofs(moment) + edge_weights(node) * &
                        polynomial * dot_product(reference_field, tangent)
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, triangle_weights, status)
        if (status /= 0) return
        do component = 1, 2
            do total_degree = 0, order - 2
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    moment = moment + 1
                    do node = 1, size(xi)
                        physical_point = vertices(:, 1) + &
                            matmul(jacobian, [xi(node), eta(node)])
                        call field( &
                            physical_point(1), physical_point(2), &
                            physical_field)
                        reference_field = &
                            matmul(transpose(jacobian), physical_field)
                        dofs(moment) = dofs(moment) + &
                            triangle_weights(node) * &
                            xi(node)**x_degree * eta(node)**y_degree * &
                            reference_field(component)
                    end do
                end do
            end do
        end do
        if (moment /= size(dofs)) then
            status = 2
            return
        end if
        status = 0
    end subroutine interpolate_triangle_nedelec

    subroutine interpolate_triangle_rt( &
            vertices, degree, quadrature_degree, field, dofs, status)
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: degree, quadrature_degree
        procedure(vector_field_2d) :: field
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: eta(:), triangle_weights(:), xi(:)
        real(dp) :: determinant, edge_point(2), jacobian(2, 2), normal(2)
        real(dp) :: physical_field(2), physical_point(2), polynomial
        real(dp) :: reference_field(2)
        integer :: component, edge, exponent, moment, node, node_count
        integer :: total_degree, x_degree, y_degree

        status = 1
        if (degree < 0 .or. quadrature_degree < 0) return
        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        determinant = jacobian(1, 1) * jacobian(2, 2) - &
            jacobian(1, 2) * jacobian(2, 1)
        if (determinant <= 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, maxval(abs(jacobian))**2)) return

        allocate(dofs((degree + 1) * (degree + 3)))
        dofs = 0.0_dp
        node_count = max(degree + 2, (quadrature_degree + 2) / 2)
        allocate(edge_nodes(node_count), edge_weights(node_count))
        call gauss_legendre_ab( &
            node_count, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)

        moment = 0
        do edge = 1, 3
            do exponent = 0, degree
                moment = moment + 1
                do node = 1, node_count
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, normal)
                    normal = [normal(2), -normal(1)]
                    physical_point = vertices(:, 1) + &
                        matmul(jacobian, edge_point)
                    call field( &
                        physical_point(1), physical_point(2), physical_field)
                    call pull_back_rt( &
                        jacobian, physical_field, reference_field)
                    polynomial = &
                        shifted_legendre(exponent, edge_nodes(node))
                    dofs(moment) = dofs(moment) + edge_weights(node) * &
                        polynomial * dot_product(reference_field, normal)
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, triangle_weights, status)
        if (status /= 0) return
        do component = 1, 2
            do total_degree = 0, degree - 1
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    moment = moment + 1
                    do node = 1, size(xi)
                        physical_point = vertices(:, 1) + &
                            matmul(jacobian, [xi(node), eta(node)])
                        call field( &
                            physical_point(1), physical_point(2), &
                            physical_field)
                        call pull_back_rt( &
                            jacobian, physical_field, reference_field)
                        dofs(moment) = dofs(moment) + &
                            triangle_weights(node) * &
                            xi(node)**x_degree * eta(node)**y_degree * &
                            reference_field(component)
                    end do
                end do
            end do
        end do
        if (moment /= size(dofs)) then
            status = 2
            return
        end if
        status = 0
    end subroutine interpolate_triangle_rt

    pure subroutine pull_back_rt(jacobian, physical_field, reference_field)
        real(dp), intent(in) :: jacobian(2, 2), physical_field(2)
        real(dp), intent(out) :: reference_field(2)

        reference_field(1) = &
            jacobian(2, 2) * physical_field(1) - &
            jacobian(1, 2) * physical_field(2)
        reference_field(2) = &
            -jacobian(2, 1) * physical_field(1) + &
            jacobian(1, 1) * physical_field(2)
    end subroutine pull_back_rt

    pure subroutine reference_edge(edge, parameter, point, tangent)
        integer, intent(in) :: edge
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: point(2), tangent(2)

        select case (edge)
        case (1)
            point = [parameter, 0.0_dp]
            tangent = [1.0_dp, 0.0_dp]
        case (2)
            point = [1.0_dp - parameter, parameter]
            tangent = [-1.0_dp, 1.0_dp]
        case (3)
            point = [0.0_dp, 1.0_dp - parameter]
            tangent = [0.0_dp, -1.0_dp]
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

end module fortfem_triangle_vector_interpolation
