module fortfem_tetra_discontinuous_arbitrary_order
    use fortfem_kinds, only: dp
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
        if (degree < 0 .or. degree > 4) return
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
            values(basis_id) = &
                x**basis%exponents(1, basis_id)* &
                y**basis%exponents(2, basis_id)* &
                z**basis%exponents(3, basis_id)
        end do
        status = 0
    end subroutine evaluate_tetra_discontinuous

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
