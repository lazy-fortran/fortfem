program test_interoperability_records
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use check, only: check_condition, check_summary
    use fortfem_interoperability_records, only: interoperability_record_t, &
        read_interoperability_records
    implicit none

    character(len=*), parameter :: fixture = &
        "build/test_interoperability_records.csv"
    type(interoperability_record_t), allocatable :: records(:)
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call write_fixture(fixture)
    call read_interoperability_records(fixture, records, status)

    call record_condition(status == 0, "valid CSV is accepted")
    call record_condition(size(records) == 2, "all refinement records are read")
    call record_condition(records(1)%cells == 32, "cell count")
    call record_condition(records(2)%dofs == 81, "degree-of-freedom count")
    call record_condition( &
        abs(records(1)%h - 0.25_dp) < 1.0e-15_dp, "mesh spacing")
    call record_condition( &
        abs(records(2)%l2_error - 2.5e-2_dp) < 1.0e-15_dp, "L2 error")
    call record_condition( &
        abs(records(2)%graph_error - 1.25e-1_dp) < 1.0e-15_dp, &
        "graph error")
    call record_condition( &
        abs(records(1)%total_seconds - 3.0e-2_dp) < 1.0e-15_dp, &
        "total time")

    call write_bad_header(fixture)
    call read_interoperability_records(fixture, records, status)
    call record_condition(status == 2, "schema mismatch is rejected")

    call check_summary("interoperability benchmark records")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    subroutine write_fixture(path)
        character(len=*), intent(in) :: path
        integer :: unit

        open (newunit=unit, file=path, status="replace", action="write")
        write (unit, "(a)") &
            "mesh_id,cells,dofs,h,l2_error,h1_error,"// &
            "assembly_seconds,solve_seconds,total_seconds"
        write (unit, "(a)") &
            "unit-square-n4,32,25,0.25,0.1,0.5,0.01,0.02,0.03"
        write (unit, "(a)") &
            "unit-square-n8,128,81,0.125,0.025,0.125,0.02,0.03,0.05"
        close (unit)
    end subroutine write_fixture

    subroutine write_bad_header(path)
        character(len=*), intent(in) :: path
        integer :: unit

        open (newunit=unit, file=path, status="replace", action="write")
        write (unit, "(a)") "mesh_id,cells"
        write (unit, "(a)") "unit-square-n4,32"
        close (unit)
    end subroutine write_bad_header

end program test_interoperability_records
