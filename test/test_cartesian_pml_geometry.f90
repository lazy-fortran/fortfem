program test_cartesian_pml_geometry
    use check, only: check_condition, check_summary
    use fortfem_api, only: build_cartesian_pml_element_stretch
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(2, 9) = reshape([ &
        0.2_dp, 0.2_dp, 0.4_dp, 0.2_dp, 0.2_dp, 0.4_dp, &
        1.2_dp, 0.2_dp, 1.4_dp, 0.2_dp, 1.2_dp, 0.4_dp, &
        1.2_dp, 1.2_dp, 1.4_dp, 1.2_dp, 1.2_dp, 1.4_dp], [2, 9])
    integer, parameter :: cells(3, 3) = reshape([ &
        1, 2, 3, 4, 5, 6, 7, 8, 9], [3, 3])
    real(dp), parameter :: physical_min(2) = [0.0_dp, 0.0_dp]
    real(dp), parameter :: physical_max(2) = [1.0_dp, 1.0_dp]
    real(dp), parameter :: outer_min(2) = [-0.5_dp, -0.5_dp]
    real(dp), parameter :: outer_max(2) = [1.5_dp, 1.5_dp]
    complex(dp), allocatable :: stretch(:, :)
    real(dp) :: expected_sigma
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call build_cartesian_pml_element_stretch( &
        vertices, cells, physical_min, physical_max, outer_min, outer_max, &
        2.0_dp, 4.0_dp, 2, stretch, status)
    call record_condition(status == 0, "Cartesian PML geometry is accepted")
    call record_condition(maxval(abs(stretch(:, 1) - &
        cmplx(1.0_dp, 0.0_dp, dp))) < 1.0e-14_dp, &
        "physical-domain element retains unit stretch")
    expected_sigma = 4.0_dp*((4.0_dp/15.0_dp)/0.5_dp)**2
    call record_condition(abs(stretch(1, 2) - &
        cmplx(1.0_dp, expected_sigma/2.0_dp, dp)) < 1.0e-14_dp .and. &
        abs(stretch(2, 2) - cmplx(1.0_dp, 0.0_dp, dp)) < 1.0e-14_dp, &
        "face-layer element receives the polynomial normal stretch")
    call record_condition(abs(stretch(1, 3) - &
        cmplx(1.0_dp, expected_sigma/2.0_dp, dp)) < 1.0e-14_dp .and. &
        abs(stretch(2, 3) - &
        cmplx(1.0_dp, expected_sigma/2.0_dp, dp)) < 1.0e-14_dp, &
        "corner-layer element receives independent coordinate stretches")
    call check_summary("Cartesian PML geometry")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_cartesian_pml_geometry
