module fortfem_tetra_lagrange_arbitrary_order
    use fortfem_kinds, only: dp
    implicit none

    private

    type :: tetra_lagrange_t
        integer :: degree = -1
        integer, allocatable :: barycentric_indices(:, :)
        real(dp), allocatable :: nodes(:, :)
    end type tetra_lagrange_t

    interface assignment(=)
        module procedure assign_tetra_lagrange
    end interface

    public :: assignment(=)
    public :: evaluate_tetra_lagrange
    public :: evaluate_tetra_lagrange_jvp
    public :: evaluate_tetra_lagrange_vjp
    public :: initialize_tetra_lagrange
    public :: tetra_lagrange_dof_count
    public :: tetra_lagrange_barycentric_indices
    public :: tetra_lagrange_nodes
    public :: tetra_lagrange_t

contains

    subroutine initialize_tetra_lagrange(degree, basis, status)
        integer, intent(in) :: degree
        type(tetra_lagrange_t), intent(out) :: basis
        integer, intent(out) :: status

        integer :: dof, first, fourth, second, third

        basis%degree = -1
        if (allocated(basis%barycentric_indices)) then
            deallocate(basis%barycentric_indices)
        end if
        if (allocated(basis%nodes)) deallocate(basis%nodes)
        status = 1
        if (degree < 0) return
        allocate(basis%barycentric_indices(4, &
            (degree + 1)*(degree + 2)*(degree + 3)/6))
        allocate(basis%nodes(3, size(basis%barycentric_indices, 2)))
        if (degree == 0) then
            basis%barycentric_indices(:, 1) = 0
            basis%nodes(:, 1) = 0.25_dp
        else
            dof = 0
            do first = 0, degree
                do second = 0, degree - first
                    do third = 0, degree - first - second
                        fourth = degree - first - second - third
                        dof = dof + 1
                        basis%barycentric_indices(:, dof) = &
                            [first, second, third, fourth]
                        basis%nodes(:, dof) = real( &
                            [second, third, fourth], dp)/real(degree, dp)
                    end do
                end do
            end do
        end if
        basis%degree = degree
        status = 0
    end subroutine initialize_tetra_lagrange

    pure subroutine evaluate_tetra_lagrange( &
            basis, point, values, gradients, status)
        type(tetra_lagrange_t), intent(in) :: basis
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(:), gradients(:, :)
        integer, intent(out) :: status

        real(dp) :: barycentric(4), derivatives(4), factors(4)
        real(dp) :: second_derivatives(4)
        real(dp) :: lambda_derivatives(4)
        integer :: basis_id, component, coordinate

        values = 0.0_dp
        gradients = 0.0_dp
        status = 1
        if (basis%degree < 0) return
        if (.not. allocated(basis%barycentric_indices)) return
        if (size(values) /= size(basis%barycentric_indices, 2)) return
        if (size(gradients, 1) /= 3) return
        if (size(gradients, 2) /= size(values)) return
        barycentric = [1.0_dp - sum(point), point]
        if (any(barycentric < -64.0_dp*epsilon(1.0_dp))) return

        do basis_id = 1, size(values)
            do component = 1, 4
                call cardinal_factor( &
                    basis%barycentric_indices(component, basis_id), &
                    basis%degree, barycentric(component), factors(component), &
                    derivatives(component), second_derivatives(component))
            end do
            values(basis_id) = product(factors)
            lambda_derivatives = 0.0_dp
            do component = 1, 4
                lambda_derivatives(component) = derivatives(component)
                do coordinate = 1, 4
                    if (coordinate /= component) then
                        lambda_derivatives(component) = &
                            lambda_derivatives(component)*factors(coordinate)
                    end if
                end do
            end do
            do coordinate = 1, 3
                gradients(coordinate, basis_id) = &
                    lambda_derivatives(coordinate + 1) - &
                    lambda_derivatives(1)
            end do
        end do
        status = 0
    end subroutine evaluate_tetra_lagrange

    pure subroutine evaluate_tetra_lagrange_jvp( &
            basis, point, point_dot, values_dot, gradients_dot, status)
        type(tetra_lagrange_t), intent(in) :: basis
        real(dp), intent(in) :: point(3), point_dot(3)
        real(dp), intent(out) :: values_dot(:), gradients_dot(:, :)
        integer, intent(out) :: status

        real(dp) :: barycentric(4), barycentric_dot(4)
        real(dp) :: derivatives(4), factors(4), hessian(4, 4)
        real(dp) :: second_derivatives(4)
        integer :: basis_id, component, coordinate, factor

        values_dot = 0.0_dp
        gradients_dot = 0.0_dp
        status = 1
        if (basis%degree < 0) return
        if (.not. allocated(basis%barycentric_indices)) return
        if (size(values_dot) /= size(basis%barycentric_indices, 2)) return
        if (size(gradients_dot, 1) /= 3) return
        if (size(gradients_dot, 2) /= size(values_dot)) return
        barycentric = [1.0_dp - sum(point), point]
        barycentric_dot = [-sum(point_dot), point_dot]
        if (any(barycentric < -64.0_dp*epsilon(1.0_dp))) return

        do basis_id = 1, size(values_dot)
            do component = 1, 4
                call cardinal_factor( &
                    basis%barycentric_indices(component, basis_id), &
                    basis%degree, barycentric(component), factors(component), &
                    derivatives(component), second_derivatives(component))
            end do
            values_dot(basis_id) = 0.0_dp
            hessian = 0.0_dp
            do component = 1, 4
                hessian(component, component) = second_derivatives(component)
                do factor = 1, 4
                    if (factor /= component) then
                        hessian(component, component) = &
                            hessian(component, component)*factors(factor)
                    end if
                end do
                values_dot(basis_id) = values_dot(basis_id) + &
                    derivatives(component)*barycentric_dot(component)* &
                    product_except(factors, component)
                do coordinate = 1, 4
                    if (coordinate == component) cycle
                    hessian(component, coordinate) = &
                        derivatives(component)*derivatives(coordinate)* &
                        product_except_two(factors, component, coordinate)
                end do
            end do
            do coordinate = 1, 3
                gradients_dot(coordinate, basis_id) = dot_product( &
                    hessian(coordinate + 1, :) - hessian(1, :), &
                    barycentric_dot)
            end do
        end do
        status = 0
    end subroutine evaluate_tetra_lagrange_jvp

    pure subroutine evaluate_tetra_lagrange_vjp( &
            basis, point, values_bar, gradients_bar, point_bar, status)
        type(tetra_lagrange_t), intent(in) :: basis
        real(dp), intent(in) :: point(3), values_bar(:), gradients_bar(:, :)
        real(dp), intent(out) :: point_bar(3)
        integer, intent(out) :: status

        real(dp) :: point_dot(3)
        real(dp), allocatable :: gradients_dot(:, :), values_dot(:)
        integer :: direction

        point_bar = 0.0_dp
        status = 1
        if (size(values_bar) /= tetra_lagrange_dof_count(basis)) return
        if (size(gradients_bar, 1) /= 3 .or. &
            size(gradients_bar, 2) /= size(values_bar)) return
        allocate(values_dot(size(values_bar)))
        allocate(gradients_dot(3, size(values_bar)))
        do direction = 1, 3
            point_dot = 0.0_dp
            point_dot(direction) = 1.0_dp
            call evaluate_tetra_lagrange_jvp( &
                basis, point, point_dot, values_dot, gradients_dot, status)
            if (status /= 0) return
            point_bar(direction) = dot_product(values_bar, values_dot) + &
                sum(gradients_bar*gradients_dot)
        end do
    end subroutine evaluate_tetra_lagrange_vjp

    pure subroutine cardinal_factor( &
            index, degree, lambda, value, derivative, second_derivative)
        integer, intent(in) :: index, degree
        real(dp), intent(in) :: lambda
        real(dp), intent(out) :: value, derivative, second_derivative

        real(dp) :: term
        integer :: differentiated_factor, factor, second_factor

        value = 1.0_dp
        do factor = 0, index - 1
            value = value*(real(degree, dp)*lambda - real(factor, dp))
        end do
        value = value/factorial(index)
        derivative = 0.0_dp
        do differentiated_factor = 0, index - 1
            term = real(degree, dp)
            do factor = 0, index - 1
                if (factor /= differentiated_factor) then
                    term = term*( &
                        real(degree, dp)*lambda - real(factor, dp))
                end if
            end do
            derivative = derivative + term
        end do
        derivative = derivative/factorial(index)
        second_derivative = 0.0_dp
        do differentiated_factor = 0, index - 1
            do second_factor = 0, index - 1
                if (second_factor == differentiated_factor) cycle
                term = real(degree*degree, dp)
                do factor = 0, index - 1
                    if (factor /= differentiated_factor .and. &
                        factor /= second_factor) then
                        term = term*( &
                            real(degree, dp)*lambda - real(factor, dp))
                    end if
                end do
                second_derivative = second_derivative + term
            end do
        end do
        second_derivative = second_derivative/factorial(index)
    end subroutine cardinal_factor

    pure real(dp) function product_except(values, excluded) result(value)
        real(dp), intent(in) :: values(4)
        integer, intent(in) :: excluded
        integer :: component

        value = 1.0_dp
        do component = 1, 4
            if (component /= excluded) value = value*values(component)
        end do
    end function product_except

    pure real(dp) function product_except_two( &
            values, first_excluded, second_excluded) result(value)
        real(dp), intent(in) :: values(4)
        integer, intent(in) :: first_excluded, second_excluded
        integer :: component

        value = 1.0_dp
        do component = 1, 4
            if (component /= first_excluded .and. &
                component /= second_excluded) value = value*values(component)
        end do
    end function product_except_two

    pure function factorial(argument) result(value)
        integer, intent(in) :: argument
        real(dp) :: value

        integer :: factor

        value = 1.0_dp
        do factor = 2, argument
            value = value*real(factor, dp)
        end do
    end function factorial

    pure integer function tetra_lagrange_dof_count(basis) result(dof_count)
        type(tetra_lagrange_t), intent(in) :: basis

        dof_count = 0
        if (allocated(basis%nodes)) dof_count = size(basis%nodes, 2)
    end function tetra_lagrange_dof_count

    subroutine tetra_lagrange_barycentric_indices(basis, indices)
        type(tetra_lagrange_t), intent(in) :: basis
        integer, allocatable, intent(out) :: indices(:, :)

        allocate(indices(4, tetra_lagrange_dof_count(basis)))
        indices = basis%barycentric_indices
    end subroutine tetra_lagrange_barycentric_indices

    subroutine tetra_lagrange_nodes(basis, nodes)
        type(tetra_lagrange_t), intent(in) :: basis
        real(dp), allocatable, intent(out) :: nodes(:, :)

        allocate(nodes(3, tetra_lagrange_dof_count(basis)))
        nodes = basis%nodes
    end subroutine tetra_lagrange_nodes

    subroutine assign_tetra_lagrange(left, right)
        type(tetra_lagrange_t), intent(out) :: left
        type(tetra_lagrange_t), intent(in) :: right

        left%degree = right%degree
        if (allocated(right%barycentric_indices)) then
            allocate( &
                left%barycentric_indices, source=right%barycentric_indices)
        end if
        if (allocated(right%nodes)) then
            allocate(left%nodes, source=right%nodes)
        end if
    end subroutine assign_tetra_lagrange

end module fortfem_tetra_lagrange_arbitrary_order
