program test_triangle_nedelec_arbitrary_order
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_triangle_nedelec_first_kind, &
        initialize_triangle_nedelec_first_kind, triangle_duffy_quadrature, &
        triangle_nedelec_dof_count, triangle_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    type(triangle_nedelec_first_kind_t) :: basis, copied_basis
    real(dp), allocatable :: curls(:), dofs(:), triangle_weights(:)
    real(dp), allocatable :: values(:, :), xi(:), eta(:)
    real(dp) :: expected(2), point(2), reconstructed(2), reconstructed_curl
    integer :: basis_dof, dof_count, order, status
    logical :: all_passed

    all_passed = .true.
    point = [0.27_dp, 0.19_dp]
    do order = 1, 4
        call initialize_triangle_nedelec_first_kind(order, basis, status)
        dof_count = triangle_nedelec_dof_count(basis)
        call record_condition(status == 0 .and. &
            dof_count == order * (order + 2), &
            "Arbitrary-order first-kind Nedelec basis has the exact dimension")
        allocate(values(2, dof_count), curls(dof_count), dofs(dof_count))

        do basis_dof = 1, dof_count
            call check_basis_moments( &
                basis, order, basis_dof, dof_count)
        end do

        call polynomial_gradient_dofs(order, dofs)
        call evaluate_triangle_nedelec_first_kind( &
            basis, point(1), point(2), values, curls, status)
        reconstructed = matmul(values, dofs)
        reconstructed_curl = dot_product(curls, dofs)
        expected = [real(order, dp) * point(1)**(order - 1), 0.0_dp]
        call record_condition(status == 0 .and. &
            maxval(abs(reconstructed - expected)) < 2.0e-10_dp, &
            "Nedelec interpolation reproduces a polynomial gradient")
        call record_condition(abs(reconstructed_curl) < 3.0e-10_dp, &
            "Nedelec interpolation preserves zero curl of a gradient")

        copied_basis = basis
        call initialize_triangle_nedelec_first_kind(1, basis, status)
        call evaluate_triangle_nedelec_first_kind( &
            copied_basis, point(1), point(2), values, curls, status)
        call record_condition(status == 0, &
            "Assigned Nedelec basis retains an independent deep copy")

        deallocate(values, curls, dofs)
    end do

    call initialize_triangle_nedelec_first_kind(0, basis, status)
    call record_condition(status /= 0, &
        "First-kind Nedelec basis rejects order zero")

    call check_summary("Arbitrary-order triangle first-kind Nedelec basis")
    if (.not. all_passed) error stop 1

contains

    subroutine check_basis_moments(basis, order, basis_dof, dof_count)
        type(triangle_nedelec_first_kind_t), intent(in) :: basis
        integer, intent(in) :: order, basis_dof, dof_count

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: local_curls(:), local_values(:, :)
        real(dp) :: edge_point(2), tangent(2), moment, polynomial
        integer :: component, edge, exponent, moment_dof, node
        integer :: status, total_degree, x_degree, y_degree

        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        allocate(local_values(2, dof_count), local_curls(dof_count))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        moment_dof = 0
        do edge = 1, 3
            do exponent = 0, order - 1
                moment_dof = moment_dof + 1
                moment = 0.0_dp
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, tangent)
                    call evaluate_triangle_nedelec_first_kind( &
                        basis, edge_point(1), edge_point(2), local_values, &
                        local_curls, status)
                    polynomial = shifted_legendre(exponent, edge_nodes(node))
                    moment = moment + edge_weights(node) * polynomial * &
                        dot_product(local_values(:, basis_dof), tangent)
                end do
                call record_condition(abs(moment - kronecker_delta( &
                    moment_dof, basis_dof)) < 3.0e-10_dp, &
                    "Nedelec edge moments form a Kronecker basis")
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * order, xi, eta, triangle_weights, status)
        do component = 1, 2
            do total_degree = 0, order - 2
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    moment_dof = moment_dof + 1
                    moment = 0.0_dp
                    do node = 1, size(xi)
                        call evaluate_triangle_nedelec_first_kind( &
                            basis, xi(node), eta(node), local_values, &
                            local_curls, status)
                        moment = moment + triangle_weights(node) * &
                            local_values(component, basis_dof) * &
                            xi(node)**x_degree * eta(node)**y_degree
                    end do
                    call record_condition(abs(moment - kronecker_delta( &
                        moment_dof, basis_dof)) < 3.0e-10_dp, &
                        "Nedelec cell moments form a Kronecker basis")
                end do
            end do
        end do
        deallocate(edge_nodes, edge_weights, local_values, local_curls)
        if (allocated(xi)) deallocate(xi)
        if (allocated(eta)) deallocate(eta)
        if (allocated(triangle_weights)) deallocate(triangle_weights)
    end subroutine check_basis_moments

    subroutine polynomial_gradient_dofs(order, dofs)
        integer, intent(in) :: order
        real(dp), intent(out) :: dofs(:)

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp) :: edge_point(2), field(2), tangent(2)
        integer :: component, edge, exponent, moment_dof, node
        integer :: status, total_degree, x_degree, y_degree

        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        dofs = 0.0_dp
        moment_dof = 0
        do edge = 1, 3
            do exponent = 0, order - 1
                moment_dof = moment_dof + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, tangent)
                    field = [ &
                        real(order, dp) * edge_point(1)**(order - 1), 0.0_dp]
                    dofs(moment_dof) = dofs(moment_dof) + &
                        edge_weights(node) * &
                        shifted_legendre(exponent, edge_nodes(node)) * &
                        dot_product(field, tangent)
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * order, xi, eta, triangle_weights, status)
        do component = 1, 2
            do total_degree = 0, order - 2
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    moment_dof = moment_dof + 1
                    if (component == 1) then
                        dofs(moment_dof) = sum(triangle_weights * &
                            real(order, dp) * xi**(order - 1 + x_degree) * &
                            eta**y_degree)
                    end if
                end do
            end do
        end do
        deallocate(edge_nodes, edge_weights)
        if (allocated(xi)) deallocate(xi)
        if (allocated(eta)) deallocate(eta)
        if (allocated(triangle_weights)) deallocate(triangle_weights)
    end subroutine polynomial_gradient_dofs

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

    pure function kronecker_delta(first, second) result(value)
        integer, intent(in) :: first, second
        real(dp) :: value

        if (first == second) then
            value = 1.0_dp
        else
            value = 0.0_dp
        end if
    end function kronecker_delta

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_nedelec_arbitrary_order
