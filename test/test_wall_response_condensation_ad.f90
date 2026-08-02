program test_wall_response_condensation_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        condense_wall_response_blocks, condense_wall_response_blocks_jvp, &
        condense_wall_response_blocks_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: exterior_count = 2, wall_count = 2
    real(dp), parameter :: step = 1.0e-6_dp
    complex(dp) :: ee(exterior_count, exterior_count)
    complex(dp) :: ey(exterior_count, wall_count), ye(wall_count, exterior_count)
    complex(dp) :: yy(wall_count, wall_count)
    complex(dp) :: ee_dot(exterior_count, exterior_count)
    complex(dp) :: ey_dot(exterior_count, wall_count), ye_dot(wall_count, exterior_count)
    complex(dp) :: yy_dot(wall_count, wall_count)
    complex(dp) :: ee_bar(exterior_count, exterior_count)
    complex(dp) :: ey_bar(exterior_count, wall_count), ye_bar(wall_count, exterior_count)
    complex(dp) :: yy_bar(wall_count, wall_count)
    complex(dp) :: effective(exterior_count, exterior_count)
    complex(dp) :: effective_dot(exterior_count, exterior_count)
    complex(dp) :: effective_bar(exterior_count, exterior_count)
    complex(dp) :: effective_plus(exterior_count, exterior_count)
    complex(dp) :: effective_minus(exterior_count, exterior_count)
    complex(dp), allocatable :: ee_bar_alloc(:, :), ey_bar_alloc(:, :)
    complex(dp), allocatable :: ye_bar_alloc(:, :), yy_bar_alloc(:, :)
    real(dp) :: derivative_error, adjoint_error, lhs, rhs
    integer :: status
    logical :: all_passed

    all_passed = .true.
    ee = reshape([ &
        cmplx(2.0_dp, 0.1_dp, dp), cmplx(0.2_dp, -0.3_dp, dp), &
        cmplx(-0.1_dp, 0.4_dp, dp), cmplx(1.5_dp, -0.2_dp, dp)], [2, 2])
    ey = reshape([ &
        cmplx(0.3_dp, 0.2_dp, dp), cmplx(-0.1_dp, 0.4_dp, dp), &
        cmplx(0.2_dp, -0.3_dp, dp), cmplx(0.5_dp, 0.1_dp, dp)], [2, 2])
    ye = reshape([ &
        cmplx(0.1_dp, -0.2_dp, dp), cmplx(0.4_dp, 0.3_dp, dp), &
        cmplx(-0.2_dp, 0.1_dp, dp), cmplx(0.2_dp, -0.4_dp, dp)], [2, 2])
    yy = reshape([ &
        cmplx(1.8_dp, 0.2_dp, dp), cmplx(0.1_dp, -0.1_dp, dp), &
        cmplx(-0.2_dp, 0.3_dp, dp), cmplx(1.3_dp, -0.2_dp, dp)], [2, 2])
    ee_dot = cmplx(reshape([0.03_dp, -0.02_dp, 0.01_dp, 0.04_dp], [2, 2]), 0.0_dp, dp)
    ey_dot = cmplx(reshape([0.02_dp, 0.01_dp, -0.03_dp, 0.02_dp], [2, 2]), 0.0_dp, dp)
    ye_dot = cmplx(reshape([-0.01_dp, 0.02_dp, 0.03_dp, -0.02_dp], [2, 2]), 0.0_dp, dp)
    yy_dot = cmplx(reshape([0.01_dp, -0.02_dp, 0.02_dp, 0.03_dp], [2, 2]), 0.0_dp, dp)

    call condense_wall_response_blocks(ee, ey, ye, yy, effective, status)
    call record_condition(status == 0, &
        "wall response Schur complement accepts nonsingular blocks")
    call condense_wall_response_blocks_jvp( &
        ee, ey, ye, yy, ee_dot, ey_dot, ye_dot, yy_dot, effective_dot, status)
    call condense_wall_response_blocks( &
        ee + step*ee_dot, ey + step*ey_dot, ye + step*ye_dot, yy + step*yy_dot, &
        effective_plus, status)
    call condense_wall_response_blocks( &
        ee - step*ee_dot, ey - step*ey_dot, ye - step*ye_dot, yy - step*yy_dot, &
        effective_minus, status)
    derivative_error = maxval(abs(effective_dot - &
        (effective_plus - effective_minus)/(2.0_dp*step)))
    call record_condition(status == 0 .and. derivative_error < 2.0e-7_dp, &
        "wall response Schur-complement JVP matches central reassembly")

    effective_bar = reshape([ &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.1_dp, 0.6_dp, dp), cmplx(-0.2_dp, 0.5_dp, dp)], [2, 2])
    call condense_wall_response_blocks_vjp( &
        ee, ey, ye, yy, effective_bar, ee_bar_alloc, ey_bar_alloc, ye_bar_alloc, &
        yy_bar_alloc, status)
    lhs = real(sum(conjg(effective_bar)*effective_dot), dp)
    rhs = real(sum(conjg(ee_bar_alloc)*ee_dot) + &
        sum(conjg(ey_bar_alloc)*ey_dot) + sum(conjg(ye_bar_alloc)*ye_dot) + &
        sum(conjg(yy_bar_alloc)*yy_dot), dp)
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "wall response Schur-complement VJP satisfies the real complex adjoint")

    call check_summary("wall response condensation")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_wall_response_condensation_ad
