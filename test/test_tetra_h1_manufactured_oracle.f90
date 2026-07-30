program test_tetra_h1_manufactured_oracle
    use check, only: check_condition, check_summary
    use fortfem_generated_tetra_h1_oracle, only: generated_tetra_h1_oracle
    use fortfem_kinds, only: dp
    implicit none

    real(dp) :: values(5)
    logical :: all_passed

    all_passed = .true.
    call generated_tetra_h1_oracle(0.1_dp, 0.2_dp, 0.3_dp, values)
    call record_condition( &
        abs(values(1) - 0.0024_dp) < 1.0e-15_dp, &
        "quartic tetrahedral bubble value")
    call record_condition( &
        maxval(abs(values(2:4) - [0.018_dp, 0.006_dp, 0.002_dp])) < &
        1.0e-15_dp, &
        "quartic tetrahedral bubble gradient")
    call record_condition( &
        abs(values(5) - 0.22_dp) < 1.0e-14_dp, &
        "source equals minus the analytical Laplacian")

    call check_boundary_trace()
    call check_summary("generated tetrahedral H1 manufactured oracle")
    if (.not. all_passed) error stop 1

contains

    subroutine check_boundary_trace()
        real(dp) :: boundary_values(5)

        call generated_tetra_h1_oracle(0.0_dp, 0.2_dp, 0.3_dp, boundary_values)
        call record_condition( &
            abs(boundary_values(1)) < 1.0e-15_dp, "zero trace on x=0")
        call generated_tetra_h1_oracle(0.1_dp, 0.0_dp, 0.3_dp, boundary_values)
        call record_condition( &
            abs(boundary_values(1)) < 1.0e-15_dp, "zero trace on y=0")
        call generated_tetra_h1_oracle(0.1_dp, 0.2_dp, 0.0_dp, boundary_values)
        call record_condition( &
            abs(boundary_values(1)) < 1.0e-15_dp, "zero trace on z=0")
        call generated_tetra_h1_oracle(0.2_dp, 0.3_dp, 0.5_dp, boundary_values)
        call record_condition( &
            abs(boundary_values(1)) < 1.0e-15_dp, &
            "zero trace on x+y+z=1")
    end subroutine check_boundary_trace

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_h1_manufactured_oracle
