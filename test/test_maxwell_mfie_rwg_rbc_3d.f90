program test_maxwell_mfie_rwg_rbc_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_maxwell_mfie_rwg_rbc_3d
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: matrix(:, :), principal_depth1(:, :)
    complex(dp), allocatable :: principal_depth2(:, :), principal_depth3(:, :)
    real(dp) :: error12, error23, vertices(3, 4)
    integer :: status, triangles(3, 4)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    triangles(:, 1) = [1, 3, 2]
    triangles(:, 2) = [1, 2, 4]
    triangles(:, 3) = [1, 4, 3]
    triangles(:, 4) = [2, 3, 4]

    call assemble_maxwell_mfie_rwg_rbc_3d( &
        vertices, triangles, 1.3_dp, 2, 1.0e-14_dp, 1, principal_depth1, &
        matrix, status)
    call assemble_maxwell_mfie_rwg_rbc_3d( &
        vertices, triangles, 1.3_dp, 2, 1.0e-14_dp, 2, principal_depth2, &
        matrix, status)
    call assemble_maxwell_mfie_rwg_rbc_3d( &
        vertices, triangles, 1.3_dp, 2, 1.0e-14_dp, 3, principal_depth3, &
        matrix, status)
    error12 = sqrt(sum(abs(principal_depth2 - principal_depth1)**2))
    error23 = sqrt(sum(abs(principal_depth3 - principal_depth2)**2))
    call record_condition(status == 0 .and. error23 < 0.65_dp*error12, &
        "adaptive RWG-RBC principal value converges under touching refinement")
    call record_condition(maxval(abs( &
        principal_depth3 - transpose(principal_depth3))) > 1.0e-5_dp, &
        "Maxwell magnetic principal value retains its nonsymmetric dual action")
    call record_condition(all(abs(matrix) < huge(1.0_dp)), &
        "MFIE jump-plus-principal-value matrix remains finite")

    call check_summary("Three-dimensional Maxwell RWG-RBC MFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_mfie_rwg_rbc_3d
