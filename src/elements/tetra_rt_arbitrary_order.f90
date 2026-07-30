module fortfem_tetra_rt_arbitrary_order
    use fortfem_kinds, only: dp
    use fortfem_generated_tetra_rt_candidates_degree_0, only: &
        evaluate_rt_candidates_degree_0
    use fortfem_generated_tetra_rt_candidates_degree_1, only: &
        evaluate_rt_candidates_degree_1
    use fortfem_generated_tetra_rt_candidates_degree_2, only: &
        evaluate_rt_candidates_degree_2
    use fortfem_generated_tetra_rt_candidates_degree_3, only: &
        evaluate_rt_candidates_degree_3
    use fortfem_generated_tetra_rt_candidates_degree_4, only: &
        evaluate_rt_candidates_degree_4
    use fortfem_generated_tetra_rt_coefficients, only: &
        load_tetra_rt_coefficients
    implicit none
    private

    type :: tetra_rt_t
        integer :: degree = -1
        integer :: dof_count = 0
        real(dp), allocatable :: coefficients(:, :)
    end type tetra_rt_t

    interface assignment(=)
        module procedure assign_tetra_rt
    end interface

    public :: assignment(=)
    public :: evaluate_tetra_rt
    public :: initialize_tetra_rt
    public :: tetra_rt_dof_count
    public :: tetra_rt_t

contains

    subroutine initialize_tetra_rt(degree, basis, status)
        integer, intent(in) :: degree
        type(tetra_rt_t), intent(out) :: basis
        integer, intent(out) :: status

        basis%degree = -1
        basis%dof_count = 0
        if (allocated(basis%coefficients)) deallocate(basis%coefficients)
        status = 1
        if (degree < 0 .or. degree > 4) return
        basis%dof_count = 3*(degree + 1)*(degree + 2)*(degree + 3)/6 + &
            (degree + 1)*(degree + 2)/2
        call load_tetra_rt_coefficients( &
            degree, basis%coefficients, status)
        if (status /= 0) return
        if (size(basis%coefficients, 1) /= basis%dof_count .or. &
            size(basis%coefficients, 2) /= basis%dof_count) then
            status = 2
            return
        end if
        basis%degree = degree
        status = 0
    end subroutine initialize_tetra_rt

    subroutine evaluate_tetra_rt( &
            basis, point, values, divergences, status)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(:, :), divergences(:)
        integer, intent(out) :: status

        real(dp), allocatable :: candidate_divergences(:)
        real(dp), allocatable :: candidate_values(:, :)
        real(dp) :: tolerance

        values = 0.0_dp
        divergences = 0.0_dp
        status = 1
        if (basis%degree < 0 .or. basis%dof_count < 1) return
        if (size(values, 1) /= 3 .or. &
            size(values, 2) /= basis%dof_count) return
        if (size(divergences) /= basis%dof_count) return
        tolerance = 64.0_dp*epsilon(1.0_dp)
        if (any(point < -tolerance) .or. &
            sum(point) > 1.0_dp + tolerance) return
        allocate(candidate_values(3, basis%dof_count))
        allocate(candidate_divergences(basis%dof_count))
        call evaluate_candidates( &
            basis, point, candidate_values, candidate_divergences)
        values = matmul(candidate_values, basis%coefficients)
        divergences = matmul( &
            candidate_divergences, basis%coefficients)
        status = 0
    end subroutine evaluate_tetra_rt

    pure integer function tetra_rt_dof_count(basis) result(dof_count)
        type(tetra_rt_t), intent(in) :: basis

        dof_count = basis%dof_count
    end function tetra_rt_dof_count

    subroutine evaluate_candidates(basis, point, values, divergences)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(:, :), divergences(:)

        select case (basis%degree)
        case (0)
            call evaluate_rt_candidates_degree_0( &
                point(1), point(2), point(3), values, divergences)
        case (1)
            call evaluate_rt_candidates_degree_1( &
                point(1), point(2), point(3), values, divergences)
        case (2)
            call evaluate_rt_candidates_degree_2( &
                point(1), point(2), point(3), values, divergences)
        case (3)
            call evaluate_rt_candidates_degree_3( &
                point(1), point(2), point(3), values, divergences)
        case (4)
            call evaluate_rt_candidates_degree_4( &
                point(1), point(2), point(3), values, divergences)
        end select
    end subroutine evaluate_candidates

    subroutine assign_tetra_rt(left, right)
        type(tetra_rt_t), intent(out) :: left
        type(tetra_rt_t), intent(in) :: right

        left%degree = right%degree
        left%dof_count = right%dof_count
        if (allocated(left%coefficients)) deallocate(left%coefficients)
        if (allocated(right%coefficients)) then
            allocate(left%coefficients, source=right%coefficients)
        end if
    end subroutine assign_tetra_rt
end module fortfem_tetra_rt_arbitrary_order
