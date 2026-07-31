module fortfem_tetra_discontinuous_arbitrary_order
    use fortfem_kinds, only: dp
    use fortnum_special_jacobi, only: tetrahedron_koornwinder, &
        tetrahedron_koornwinder_gradient
    implicit none
    private

    type :: tetra_discontinuous_t
        integer :: degree = -1
        integer, allocatable :: exponents(:, :)
    end type tetra_discontinuous_t

    interface assignment(=)
        module procedure assign_tetra_discontinuous
    end interface

    public :: assignment(=)
    public :: evaluate_tetra_discontinuous
    public :: evaluate_tetra_discontinuous_jvp
    public :: evaluate_tetra_discontinuous_vjp
    public :: initialize_tetra_discontinuous
    public :: tetra_discontinuous_dof_count
    public :: tetra_discontinuous_t

contains

    subroutine initialize_tetra_discontinuous(degree, basis, status)
        integer, intent(in) :: degree
        type(tetra_discontinuous_t), intent(out) :: basis
        integer, intent(out) :: status

        integer :: first, second, entry, third, total

        basis%degree = -1
        if (allocated(basis%exponents)) deallocate(basis%exponents)
        status = 1
        if (degree < 0) return
        allocate(basis%exponents(3, (degree + 1)*(degree + 2)* &
            (degree + 3)/6))
        entry = 0
        do total = 0, degree
            do first = 0, total
                do second = 0, total - first
                    third = total - first - second
                    entry = entry + 1
                    basis%exponents(:, entry) = [first, second, third]
                end do
            end do
        end do
        basis%degree = degree
        status = 0
    end subroutine initialize_tetra_discontinuous

    pure subroutine evaluate_tetra_discontinuous( &
            basis, x, y, z, values, status)
        type(tetra_discontinuous_t), intent(in) :: basis
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: status

        integer :: basis_id

        values = 0.0_dp
        status = 1
        if (.not. allocated(basis%exponents)) return
        if (size(values) /= size(basis%exponents, 2)) return
        do basis_id = 1, size(values)
            if (basis%degree <= 5) then
                values(basis_id) = &
                    x**basis%exponents(1, basis_id)* &
                    y**basis%exponents(2, basis_id)* &
                    z**basis%exponents(3, basis_id)
            else
                values(basis_id) = tetrahedron_koornwinder( &
                    basis%exponents(1, basis_id), &
                    basis%exponents(2, basis_id), &
                    basis%exponents(3, basis_id), x, y, z)
            end if
        end do
        status = 0
    end subroutine evaluate_tetra_discontinuous

    pure subroutine evaluate_tetra_discontinuous_jvp( &
            basis, point, point_dot, values_dot, status)
        type(tetra_discontinuous_t), intent(in) :: basis
        real(dp), intent(in) :: point(3), point_dot(3)
        real(dp), intent(out) :: values_dot(:)
        integer, intent(out) :: status

        real(dp) :: gradient(3)
        integer :: basis_id

        values_dot = 0.0_dp
        status = 1
        if (.not. allocated(basis%exponents)) return
        if (size(values_dot) /= size(basis%exponents, 2)) return
        do basis_id = 1, size(values_dot)
            if (basis%degree <= 5) then
                call monomial_gradient( &
                    point, basis%exponents(:, basis_id), gradient)
            else
                call tetrahedron_koornwinder_gradient( &
                    basis%exponents(1, basis_id), &
                    basis%exponents(2, basis_id), &
                    basis%exponents(3, basis_id), &
                    point(1), point(2), point(3), gradient)
            end if
            values_dot(basis_id) = dot_product(gradient, point_dot)
        end do
        status = 0
    end subroutine evaluate_tetra_discontinuous_jvp

    pure subroutine evaluate_tetra_discontinuous_vjp( &
            basis, point, values_bar, point_bar, status)
        type(tetra_discontinuous_t), intent(in) :: basis
        real(dp), intent(in) :: point(3), values_bar(:)
        real(dp), intent(out) :: point_bar(3)
        integer, intent(out) :: status

        real(dp), allocatable :: values_dot(:)
        real(dp) :: point_dot(3)
        integer :: direction

        point_bar = 0.0_dp
        status = 1
        if (size(values_bar) /= tetra_discontinuous_dof_count(basis)) return
        allocate(values_dot(size(values_bar)))
        do direction = 1, 3
            point_dot = 0.0_dp
            point_dot(direction) = 1.0_dp
            call evaluate_tetra_discontinuous_jvp( &
                basis, point, point_dot, values_dot, status)
            if (status /= 0) return
            point_bar(direction) = dot_product(values_bar, values_dot)
        end do
    end subroutine evaluate_tetra_discontinuous_vjp

    pure subroutine monomial_gradient(point, exponents, gradient)
        real(dp), intent(in) :: point(3)
        integer, intent(in) :: exponents(3)
        real(dp), intent(out) :: gradient(3)

        integer :: direction, reduced(3)

        gradient = 0.0_dp
        do direction = 1, 3
            if (exponents(direction) == 0) cycle
            reduced = exponents
            reduced(direction) = reduced(direction) - 1
            gradient(direction) = real(exponents(direction), dp)* &
                point(1)**reduced(1)*point(2)**reduced(2)* &
                point(3)**reduced(3)
        end do
    end subroutine monomial_gradient

    pure integer function tetra_discontinuous_dof_count(basis) &
            result(dof_count)
        type(tetra_discontinuous_t), intent(in) :: basis

        dof_count = 0
        if (allocated(basis%exponents)) dof_count = size(basis%exponents, 2)
    end function tetra_discontinuous_dof_count

    subroutine assign_tetra_discontinuous(left, right)
        type(tetra_discontinuous_t), intent(out) :: left
        type(tetra_discontinuous_t), intent(in) :: right

        left%degree = right%degree
        if (allocated(left%exponents)) deallocate(left%exponents)
        if (allocated(right%exponents)) then
            allocate(left%exponents, source=right%exponents)
        end if
    end subroutine assign_tetra_discontinuous
end module fortfem_tetra_discontinuous_arbitrary_order
