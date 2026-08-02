module fortfem_solver_resource_budget
    !! Neutral wall-time, memory, and repetition budget contract.
    !!
    !! A caller owns the measurements and decides how to obtain them.  This
    !! module only validates fixed positive limits and reports whether one
    !! measured solver run fits those limits.  It deliberately has no timing,
    !! allocation, MPI, or solver implementation side effects.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    type, public :: solver_resource_budget_t
        private
        real(dp) :: max_wall_seconds = 0.0_dp
        integer(int64) :: max_peak_memory_bytes = 0_int64
        integer :: max_repetitions = 0
    end type solver_resource_budget_t

    public :: initialize_solver_resource_budget
    public :: validate_solver_resource_budget
    public :: evaluate_solver_resource_usage

contains

    subroutine initialize_solver_resource_budget( &
            budget, max_wall_seconds, max_peak_memory_bytes, max_repetitions, &
            status)
        type(solver_resource_budget_t), intent(out) :: budget
        real(dp), intent(in) :: max_wall_seconds
        integer(int64), intent(in) :: max_peak_memory_bytes
        integer, intent(in) :: max_repetitions
        type(fortsparse_status_t), intent(out) :: status

        budget = solver_resource_budget_t()
        budget%max_wall_seconds = max_wall_seconds
        budget%max_peak_memory_bytes = max_peak_memory_bytes
        budget%max_repetitions = max_repetitions
        call validate_solver_resource_budget(budget, status)
        if (status%code /= FORTSPARSE_OK) budget = solver_resource_budget_t()
    end subroutine initialize_solver_resource_budget

    subroutine validate_solver_resource_budget(budget, status)
        type(solver_resource_budget_t), intent(in) :: budget
        type(fortsparse_status_t), intent(out) :: status

        if (.not. ieee_is_finite(budget%max_wall_seconds) .or. &
            budget%max_wall_seconds <= 0.0_dp .or. &
            budget%max_peak_memory_bytes <= 0_int64 .or. &
            budget%max_repetitions < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "solver resource budget requires positive limits")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_solver_resource_budget

    subroutine evaluate_solver_resource_usage( &
            budget, wall_seconds, peak_memory_bytes, repetitions, &
            within_budget, status)
        type(solver_resource_budget_t), intent(in) :: budget
        real(dp), intent(in) :: wall_seconds
        integer(int64), intent(in) :: peak_memory_bytes
        integer, intent(in) :: repetitions
        logical, intent(out) :: within_budget
        type(fortsparse_status_t), intent(out) :: status

        within_budget = .false.
        call validate_solver_resource_budget(budget, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. ieee_is_finite(wall_seconds) .or. wall_seconds < 0.0_dp .or. &
            peak_memory_bytes < 0_int64 .or. repetitions < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "solver resource usage has invalid measurements")
            return
        end if
        within_budget = wall_seconds <= budget%max_wall_seconds .and. &
            peak_memory_bytes <= budget%max_peak_memory_bytes .and. &
            repetitions <= budget%max_repetitions
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_solver_resource_usage

end module fortfem_solver_resource_budget
