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
    public :: initialize_tetra_lagrange
    public :: tetra_lagrange_dof_count
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
        if (degree < 0 .or. degree > 4) return
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
                    derivatives(component))
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

    pure subroutine cardinal_factor( &
            index, degree, lambda, value, derivative)
        integer, intent(in) :: index, degree
        real(dp), intent(in) :: lambda
        real(dp), intent(out) :: value, derivative

        real(dp) :: term
        integer :: differentiated_factor, factor

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
    end subroutine cardinal_factor

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
