program test_paper_magnetic_box_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_magnetic_box_3d
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: reference = 0.07367135328_dp
    type(fortsparse_status_t) :: status
    real(dp) :: coarse_value, fine_value
    integer :: coarse_dofs, fine_dofs
    logical :: all_passed

    all_passed = .true.
    call solve_magnetic_box_3d( &
        6, 4, coarse_value, coarse_dofs, status)
    call record_condition(status%code == 0, &
        "Coarse magnetic-paper box solve succeeds")
    call solve_magnetic_box_3d( &
        12, 8, fine_value, fine_dofs, status)
    call record_condition(status%code == 0, &
        "Paper-resolution magnetic box solve succeeds")
    call record_condition(abs(fine_value - reference) / reference < 0.05_dp, &
        "Magnetic box centre value matches the membrane series")
    call record_condition(abs(fine_value - reference) < &
        abs(coarse_value - reference), &
        "Magnetic box centre value improves under refinement")
    call record_condition(fine_dofs > coarse_dofs, &
        "Magnetic box refinement increases edge degrees of freedom")

    call check_summary("Paper magnetic three-dimensional box")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_paper_magnetic_box_3d
