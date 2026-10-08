module fortfem_tetra_lagrange_arbitrary_order
    use fortfem_kinds, only: dp
    use fortfem_generated_barycentric_jet3, only: generated_barycentric_jet3
    use fortfem_generated_cardinal_product_jet1_order1, only: &
        generated_cardinal_product_jet1_order1
    use fortfem_generated_cardinal_product_jet1_order2, only: &
        generated_cardinal_product_jet1_order2
    use fortfem_generated_factor_product_jet3_order1, only: &
        generated_factor_product_jet3_order1
    use fortfem_generated_factor_product_jet3_order2, only: &
        generated_factor_product_jet3_order2
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

        real(dp) :: barycentric(4), barycentric_gradients(3, 4)
        real(dp) :: factors(4), derivatives(4)
        integer :: basis_id, component

        values = 0.0_dp
        gradients = 0.0_dp
        status = 1
        if (basis%degree < 0) return
        if (.not. allocated(basis%barycentric_indices)) return
        if (size(values) /= size(basis%barycentric_indices, 2)) return
        if (size(gradients, 1) /= 3) return
        if (size(gradients, 2) /= size(values)) return
        call reference_barycentric(point, barycentric, barycentric_gradients)
        if (any(barycentric < -64.0_dp*epsilon(1.0_dp))) return

        do basis_id = 1, size(values)
            do component = 1, 4
                call cardinal_factor_jet(basis%barycentric_indices(component, &
                    basis_id), basis%degree, barycentric(component), &
                    factors(component), derivatives(component))
            end do
            call generated_factor_product_jet3_order1(factors(1), derivatives(1), &
                factors(2), derivatives(2), factors(3), derivatives(3), &
                factors(4), derivatives(4), values(basis_id), &
                gradients(1, basis_id), gradients(2, basis_id), &
                gradients(3, basis_id))
        end do
        status = 0
    end subroutine evaluate_tetra_lagrange

    pure subroutine evaluate_tetra_lagrange_jvp( &
            basis, point, point_dot, values_dot, gradients_dot, status)
        type(tetra_lagrange_t), intent(in) :: basis
        real(dp), intent(in) :: point(3), point_dot(3)
        real(dp), intent(out) :: values_dot(:), gradients_dot(:, :)
        integer, intent(out) :: status

        real(dp) :: barycentric(4), barycentric_gradients(3, 4)
        real(dp) :: factors(4), derivatives(4), second_derivatives(4)
        integer :: basis_id, component

        values_dot = 0.0_dp
        gradients_dot = 0.0_dp
        status = 1
        if (basis%degree < 0) return
        if (.not. allocated(basis%barycentric_indices)) return
        if (size(values_dot) /= size(basis%barycentric_indices, 2)) return
        if (size(gradients_dot, 1) /= 3) return
        if (size(gradients_dot, 2) /= size(values_dot)) return
        call reference_barycentric(point, barycentric, barycentric_gradients)
        if (any(barycentric < -64.0_dp*epsilon(1.0_dp))) return

        do basis_id = 1, size(values_dot)
            do component = 1, 4
                call cardinal_factor_jet(basis%barycentric_indices(component, &
                    basis_id), basis%degree, barycentric(component), &
                    factors(component), derivatives(component), &
                    second_derivatives(component))
            end do
            call generated_factor_product_jet3_order2( &
                factors(1), derivatives(1), second_derivatives(1), &
                factors(2), derivatives(2), second_derivatives(2), &
                factors(3), derivatives(3), second_derivatives(3), &
                factors(4), derivatives(4), second_derivatives(4), &
                point_dot(1), point_dot(2), point_dot(3), values_dot(basis_id), &
                gradients_dot(1, basis_id), gradients_dot(2, basis_id), &
                gradients_dot(3, basis_id))
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

    pure subroutine reference_barycentric(point, barycentric, gradients)
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: barycentric(4), gradients(3, 4)
        call generated_barycentric_jet3(point(1), point(2), point(3), &
            barycentric(1), gradients(1, 1), gradients(2, 1), gradients(3, 1), &
            barycentric(2), gradients(1, 2), gradients(2, 2), gradients(3, 2), &
            barycentric(3), gradients(1, 3), gradients(2, 3), gradients(3, 3), &
            barycentric(4), gradients(1, 4), gradients(2, 4), gradients(3, 4))
    end subroutine reference_barycentric

    pure subroutine cardinal_factor_jet( &
            index, degree, lambda, value, derivative, second_derivative)
        integer, intent(in) :: index, degree
        real(dp), intent(in) :: lambda
        real(dp), intent(out) :: value, derivative
        real(dp), intent(out), optional :: second_derivative
        real(dp) :: next_value, next_derivative, next_second_derivative, normalization
        integer :: factor

        value = 1.0_dp
        derivative = 0.0_dp
        normalization = 1.0_dp
        if (present(second_derivative)) second_derivative = 0.0_dp
        do factor = 0, index - 1
            normalization = normalization*real(factor + 1, dp)
            if (present(second_derivative)) then
                call generated_cardinal_product_jet1_order2(value, derivative, &
                    second_derivative, real(degree, dp), real(factor, dp), &
                    lambda, next_value, next_derivative, next_second_derivative)
                second_derivative = next_second_derivative
            else
                call generated_cardinal_product_jet1_order1(value, derivative, &
                    real(degree, dp), real(factor, dp), lambda, &
                    next_value, next_derivative)
            end if
            value = next_value
            derivative = next_derivative
        end do
        ! The common factorial normalizes the whole generated jet.
        value = value/normalization
        derivative = derivative/normalization
        if (present(second_derivative)) then
            second_derivative = second_derivative/normalization
        end if
    end subroutine cardinal_factor_jet

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
