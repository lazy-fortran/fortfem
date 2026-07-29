module fortfem_triangle_bdm_arbitrary_order
    use fortfem_kinds, only: dp
    use fortfem_triangle_nedelec_second_kind, only: &
        assignment(=), evaluate_triangle_nedelec_second_kind, &
        initialize_triangle_nedelec_second_kind, &
        triangle_nedelec_second_kind_dof_count, &
        triangle_nedelec_second_kind_t
    implicit none

    private

    public :: assignment(=)
    public :: evaluate_triangle_bdm
    public :: initialize_triangle_bdm
    public :: triangle_bdm_basis_t
    public :: triangle_bdm_dof_count

    type :: triangle_bdm_basis_t
        private
        integer :: degree = 0
        integer :: dof_count = 0
        type(triangle_nedelec_second_kind_t) :: rotated_basis
    end type triangle_bdm_basis_t

    interface assignment(=)
        module procedure assign_triangle_bdm_basis
    end interface

contains

    subroutine initialize_triangle_bdm(degree, basis, status)
        integer, intent(in) :: degree
        type(triangle_bdm_basis_t), intent(out) :: basis
        integer, intent(out) :: status

        status = 1
        if (degree < 1) return
        basis%degree = degree
        call initialize_triangle_nedelec_second_kind( &
            degree, basis%rotated_basis, status)
        if (status /= 0) return
        basis%dof_count = &
            triangle_nedelec_second_kind_dof_count(basis%rotated_basis)
    end subroutine initialize_triangle_bdm

    subroutine evaluate_triangle_bdm( &
            basis, xi, eta, values, divergences, status)
        type(triangle_bdm_basis_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(:, :), divergences(:)
        integer, intent(out) :: status

        real(dp), allocatable :: curls(:), nedelec_values(:, :)
        integer :: basis_dof

        values = 0.0_dp
        divergences = 0.0_dp
        status = 1
        if (basis%degree < 1 .or. basis%dof_count < 1) return
        if (size(values, 1) /= 2) return
        if (size(values, 2) /= basis%dof_count) return
        if (size(divergences) /= basis%dof_count) return

        allocate(curls(basis%dof_count))
        allocate(nedelec_values(2, basis%dof_count))
        call evaluate_triangle_nedelec_second_kind( &
            basis%rotated_basis, xi, eta, nedelec_values, curls, status)
        if (status /= 0) return
        do basis_dof = 1, basis%dof_count
            values(1, basis_dof) = nedelec_values(2, basis_dof)
            values(2, basis_dof) = -nedelec_values(1, basis_dof)
        end do
        divergences = curls
        status = 0
    end subroutine evaluate_triangle_bdm

    pure function triangle_bdm_dof_count(basis) result(dof_count)
        type(triangle_bdm_basis_t), intent(in) :: basis
        integer :: dof_count

        dof_count = basis%dof_count
    end function triangle_bdm_dof_count

    subroutine assign_triangle_bdm_basis(left, right)
        type(triangle_bdm_basis_t), intent(out) :: left
        type(triangle_bdm_basis_t), intent(in) :: right

        left%degree = right%degree
        left%dof_count = right%dof_count
        left%rotated_basis = right%rotated_basis
    end subroutine assign_triangle_bdm_basis

end module fortfem_triangle_bdm_arbitrary_order
