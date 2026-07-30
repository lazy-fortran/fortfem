module fortfem_triangle_nedelec_second_kind
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_rt_arbitrary_order, only: &
        evaluate_triangle_raviart_thomas, initialize_triangle_raviart_thomas, &
        triangle_rt_basis_t, triangle_rt_dof_count
    use fortnum_linalg, only: dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    private

    public :: assignment(=)
    public :: evaluate_triangle_nedelec_second_kind
    public :: initialize_triangle_nedelec_second_kind
    public :: triangle_nedelec_second_kind_dof_count
    public :: triangle_nedelec_second_kind_t

    type :: triangle_nedelec_second_kind_t
        private
        integer :: degree = 0
        integer :: dof_count = 0
        integer, allocatable :: candidate_component(:)
        integer, allocatable :: powers(:, :)
        real(dp), allocatable :: coefficients(:, :)
    end type triangle_nedelec_second_kind_t

    interface assignment(=)
        module procedure assign_triangle_nedelec_second_kind
    end interface

contains

    subroutine initialize_triangle_nedelec_second_kind( &
            degree, basis, status)
        integer, intent(in) :: degree
        type(triangle_nedelec_second_kind_t), intent(out) :: basis
        integer, intent(out) :: status

        real(dp), allocatable :: inverse(:, :), moment_matrix(:, :)
        integer :: candidate, component, info, total_degree
        integer :: x_degree, y_degree

        status = 1
        if (degree < 1) return
        basis%degree = degree
        basis%dof_count = (degree + 1) * (degree + 2)
        allocate(basis%candidate_component(basis%dof_count))
        allocate(basis%powers(2, basis%dof_count))
        allocate(basis%coefficients(basis%dof_count, basis%dof_count))
        allocate(moment_matrix(basis%dof_count, basis%dof_count))
        allocate(inverse(basis%dof_count, basis%dof_count))

        candidate = 0
        do component = 1, 2
            do total_degree = 0, degree
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    candidate = candidate + 1
                    basis%candidate_component(candidate) = component
                    basis%powers(1, candidate) = x_degree
                    basis%powers(2, candidate) = y_degree
                end do
            end do
        end do
        if (candidate /= basis%dof_count) then
            status = 2
            return
        end if

        call build_second_kind_moment_matrix( &
            basis, moment_matrix, status)
        if (status /= 0) return
        basis%coefficients = 0.0_dp
        do candidate = 1, basis%dof_count
            basis%coefficients(candidate, candidate) = 1.0_dp
        end do
        call dense_solve(moment_matrix, basis%coefficients, inverse, info)
        if (info /= 0) then
            status = 2
            return
        end if
        basis%coefficients = inverse
        status = 0
    end subroutine initialize_triangle_nedelec_second_kind

    subroutine evaluate_triangle_nedelec_second_kind( &
            basis, xi, eta, values, curls, status)
        type(triangle_nedelec_second_kind_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(:, :), curls(:)
        integer, intent(out) :: status

        real(dp), allocatable :: candidate_curls(:), candidate_values(:, :)
        integer :: basis_dof, candidate

        values = 0.0_dp
        curls = 0.0_dp
        status = 1
        if (basis%degree < 1 .or. basis%dof_count < 1) return
        if (size(values, 1) /= 2) return
        if (size(values, 2) /= basis%dof_count) return
        if (size(curls) /= basis%dof_count) return
        if (xi < -64.0_dp * epsilon(1.0_dp) .or. &
            eta < -64.0_dp * epsilon(1.0_dp)) return
        if (xi + eta > 1.0_dp + 64.0_dp * epsilon(1.0_dp)) return

        allocate(candidate_values(2, basis%dof_count))
        allocate(candidate_curls(basis%dof_count))
        call evaluate_second_kind_candidates( &
            basis, xi, eta, candidate_values, candidate_curls)
        do basis_dof = 1, basis%dof_count
            do candidate = 1, basis%dof_count
                values(:, basis_dof) = values(:, basis_dof) + &
                    basis%coefficients(candidate, basis_dof) * &
                    candidate_values(:, candidate)
                curls(basis_dof) = curls(basis_dof) + &
                    basis%coefficients(candidate, basis_dof) * &
                    candidate_curls(candidate)
            end do
        end do
        status = 0
    end subroutine evaluate_triangle_nedelec_second_kind

    pure function triangle_nedelec_second_kind_dof_count( &
            basis) result(dof_count)
        type(triangle_nedelec_second_kind_t), intent(in) :: basis
        integer :: dof_count

        dof_count = basis%dof_count
    end function triangle_nedelec_second_kind_dof_count

    subroutine build_second_kind_moment_matrix(basis, matrix, status)
        type(triangle_nedelec_second_kind_t), intent(in) :: basis
        real(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        type(triangle_rt_basis_t) :: rt_basis
        real(dp), allocatable :: candidate_curls(:), candidate_values(:, :)
        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: eta(:), rt_divergences(:), rt_values(:, :)
        real(dp), allocatable :: triangle_weights(:), xi(:)
        real(dp) :: edge_point(2), polynomial, tangent(2)
        integer :: candidate, edge, exponent, moment, node
        integer :: rt_dof, rt_dof_count

        matrix = 0.0_dp
        status = 1
        allocate(candidate_values(2, basis%dof_count))
        allocate(candidate_curls(basis%dof_count))
        allocate(edge_nodes(basis%degree + 2))
        allocate(edge_weights(basis%degree + 2))
        call gauss_legendre_ab( &
            basis%degree + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)

        moment = 0
        do edge = 1, 3
            do exponent = 0, basis%degree
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, tangent)
                    call evaluate_second_kind_candidates( &
                        basis, edge_point(1), edge_point(2), candidate_values, &
                        candidate_curls)
                    polynomial = shifted_legendre(exponent, edge_nodes(node))
                    do candidate = 1, basis%dof_count
                        matrix(moment, candidate) = matrix(moment, candidate) + &
                            edge_weights(node) * polynomial * dot_product( &
                            candidate_values(:, candidate), tangent)
                    end do
                end do
            end do
        end do

        if (basis%degree >= 2) then
            call initialize_triangle_raviart_thomas( &
                basis%degree - 2, rt_basis, status)
            if (status /= 0) return
            rt_dof_count = triangle_rt_dof_count(rt_basis)
            allocate(rt_values(2, rt_dof_count))
            allocate(rt_divergences(rt_dof_count))
            call triangle_duffy_quadrature( &
                2 * basis%degree, xi, eta, triangle_weights, status)
            if (status /= 0) return
            do rt_dof = 1, rt_dof_count
                moment = moment + 1
                do node = 1, size(xi)
                    call evaluate_second_kind_candidates( &
                        basis, xi(node), eta(node), candidate_values, &
                        candidate_curls)
                    call evaluate_triangle_raviart_thomas( &
                        rt_basis, xi(node), eta(node), rt_values, &
                        rt_divergences, status)
                    if (status /= 0) return
                    do candidate = 1, basis%dof_count
                        matrix(moment, candidate) = &
                            matrix(moment, candidate) + &
                            triangle_weights(node) * dot_product( &
                            candidate_values(:, candidate), &
                            rt_values(:, rt_dof))
                    end do
                end do
            end do
        end if
        if (moment /= basis%dof_count) then
            status = 2
            return
        end if
        status = 0
    end subroutine build_second_kind_moment_matrix

    pure subroutine evaluate_second_kind_candidates( &
            basis, xi, eta, values, curls)
        type(triangle_nedelec_second_kind_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(:, :), curls(:)

        real(dp) :: monomial
        integer :: candidate, component, x_degree, y_degree

        values = 0.0_dp
        curls = 0.0_dp
        do candidate = 1, basis%dof_count
            component = basis%candidate_component(candidate)
            x_degree = basis%powers(1, candidate)
            y_degree = basis%powers(2, candidate)
            monomial = integer_power(xi, x_degree) * &
                integer_power(eta, y_degree)
            values(component, candidate) = monomial
            if (component == 1 .and. y_degree > 0) then
                curls(candidate) = -real(y_degree, dp) * &
                    integer_power(xi, x_degree) * &
                    integer_power(eta, y_degree - 1)
            else if (component == 2 .and. x_degree > 0) then
                curls(candidate) = real(x_degree, dp) * &
                    integer_power(xi, x_degree - 1) * &
                    integer_power(eta, y_degree)
            end if
        end do
    end subroutine evaluate_second_kind_candidates

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

    pure function integer_power(base, exponent) result(value)
        real(dp), intent(in) :: base
        integer, intent(in) :: exponent
        real(dp) :: value

        integer :: factor

        value = 1.0_dp
        do factor = 1, exponent
            value = value * base
        end do
    end function integer_power

    subroutine assign_triangle_nedelec_second_kind(left, right)
        type(triangle_nedelec_second_kind_t), intent(out) :: left
        type(triangle_nedelec_second_kind_t), intent(in) :: right

        left%degree = right%degree
        left%dof_count = right%dof_count
        if (allocated(right%candidate_component)) then
            allocate(left%candidate_component(size(right%candidate_component)))
            left%candidate_component = right%candidate_component
        end if
        if (allocated(right%powers)) then
            allocate(left%powers(size(right%powers, 1), size(right%powers, 2)))
            left%powers = right%powers
        end if
        if (allocated(right%coefficients)) then
            allocate(left%coefficients( &
                size(right%coefficients, 1), size(right%coefficients, 2)))
            left%coefficients = right%coefficients
        end if
    end subroutine assign_triangle_nedelec_second_kind

end module fortfem_triangle_nedelec_second_kind
