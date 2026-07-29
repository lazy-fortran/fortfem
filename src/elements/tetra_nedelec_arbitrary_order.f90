module fortfem_tetra_nedelec_arbitrary_order
    use fortfem_generated_tetra_nedelec_candidates_order_1, only: &
        evaluate_candidates_order_1
    use fortfem_generated_tetra_nedelec_candidates_order_2, only: &
        evaluate_candidates_order_2
    use fortfem_generated_tetra_nedelec_candidates_order_3, only: &
        evaluate_candidates_order_3
    use fortfem_generated_tetra_nedelec_candidates_order_4, only: &
        evaluate_candidates_order_4
    use fortfem_generated_tetra_nedelec_coefficients, only: &
        load_tetra_nedelec_coefficients
    use fortfem_kinds, only: dp
    implicit none

    private

    public :: assignment(=)
    public :: evaluate_tetra_nedelec_first_kind
    public :: initialize_tetra_nedelec_first_kind
    public :: tetra_nedelec_dof_count
    public :: tetra_nedelec_first_kind_t

    type :: tetra_nedelec_first_kind_t
        private
        integer :: order = 0
        integer :: dof_count = 0
        real(dp), allocatable :: coefficients(:, :)
    end type tetra_nedelec_first_kind_t

    interface assignment(=)
        module procedure assign_tetra_nedelec_first_kind
    end interface

contains

    subroutine initialize_tetra_nedelec_first_kind(order, basis, status)
        integer, intent(in) :: order
        type(tetra_nedelec_first_kind_t), intent(out) :: basis
        integer, intent(out) :: status

        status = 1
        if (order < 1 .or. order > 4) return

        basis%order = order
        basis%dof_count = order * (order + 2) * (order + 3) / 2
        call load_tetra_nedelec_coefficients( &
            order, basis%coefficients, status)
        if (status /= 0) return
        if (size(basis%coefficients, 1) /= basis%dof_count .or. &
            size(basis%coefficients, 2) /= basis%dof_count) then
            status = 2
            return
        end if
        status = 0
    end subroutine initialize_tetra_nedelec_first_kind

    subroutine evaluate_tetra_nedelec_first_kind( &
            basis, point, values, curls, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(:, :), curls(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: candidate_curls(:, :)
        real(dp), allocatable :: candidate_values(:, :)
        real(dp) :: tolerance

        values = 0.0_dp
        curls = 0.0_dp
        status = 1
        if (basis%order < 1 .or. basis%dof_count < 1) return
        if (size(values, 1) /= 3 .or. &
            size(values, 2) /= basis%dof_count) return
        if (size(curls, 1) /= 3 .or. &
            size(curls, 2) /= basis%dof_count) return
        tolerance = 64.0_dp * epsilon(1.0_dp)
        if (any(point < -tolerance)) return
        if (sum(point) > 1.0_dp + tolerance) return

        allocate( &
            candidate_values(3, basis%dof_count), &
            candidate_curls(3, basis%dof_count))
        call evaluate_nedelec_candidates( &
            basis, point, candidate_values, candidate_curls)
        values = matmul(candidate_values, basis%coefficients)
        curls = matmul(candidate_curls, basis%coefficients)
        status = 0
    end subroutine evaluate_tetra_nedelec_first_kind

    pure function tetra_nedelec_dof_count(basis) result(dof_count)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        integer :: dof_count

        dof_count = basis%dof_count
    end function tetra_nedelec_dof_count

    pure subroutine evaluate_nedelec_candidates( &
            basis, point, values, curls)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(:, :), curls(:, :)

        select case (basis%order)
        case (1)
            call evaluate_candidates_order_1( &
                point(1), point(2), point(3), values, curls)
        case (2)
            call evaluate_candidates_order_2( &
                point(1), point(2), point(3), values, curls)
        case (3)
            call evaluate_candidates_order_3( &
                point(1), point(2), point(3), values, curls)
        case (4)
            call evaluate_candidates_order_4( &
                point(1), point(2), point(3), values, curls)
        case default
            values = 0.0_dp
            curls = 0.0_dp
        end select
    end subroutine evaluate_nedelec_candidates

    subroutine assign_tetra_nedelec_first_kind(left, right)
        type(tetra_nedelec_first_kind_t), intent(out) :: left
        type(tetra_nedelec_first_kind_t), intent(in) :: right

        left%order = right%order
        left%dof_count = right%dof_count
        if (allocated(right%coefficients)) then
            allocate(left%coefficients( &
                size(right%coefficients, 1), size(right%coefficients, 2)))
            left%coefficients = right%coefficients
        end if
    end subroutine assign_tetra_nedelec_first_kind

end module fortfem_tetra_nedelec_arbitrary_order
