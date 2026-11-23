program test_parameter_continuation
    !! fortfem_parameter_continuation on a scalar family whose Newton basin
    !! shrinks with the parameter: atan(K(s) (x - c(s))) = 0, c(s) = 1 + s,
    !! K(s) = kmax^s. The solution is x = c(s); undamped Newton converges
    !! only from |K (x - c)| < 1.39, so a stage from the solution at s_ok
    !! converges only for K(s) (s - s_ok) < 1.39.
    !!
    !! 1. kmax = 10^1.5: the direct solve at s = 1 fails, the driver reaches
    !!    s = 1 with halved steps, and another largest step gives the same
    !!    solution.
    !! 2. kmax = 10^4: the steps needed fall below min_step; the driver
    !!    stalls and returns the last converged state.
    !! 3. No base: the solve at s = 0 fails and the start is returned.
    !! 4. try_target_first: one stage when the start is in the basin of the
    !!    target.
    !! 5. Invalid options are refused.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_parameter_continuation, only: continuation_options_t, &
        continuation_report_t, continuation_run, CONTINUATION_OK, &
        CONTINUATION_NO_BASE, CONTINUATION_STALLED, CONTINUATION_BAD_INPUT
    implicit none
    real(dp) :: kmax
    type(continuation_options_t) :: options
    type(continuation_report_t) :: report
    real(dp) :: x(1), first
    integer :: status, steps

    ! 1
    kmax = 10.0_dp**1.5_dp
    x = 1.9_dp
    call stage(1.0_dp, x, steps, status)
    call check_condition(status /= 0, "the direct solve at the target fails")
    x = 1.9_dp
    call continuation_run(stage, x, options, report, status)
    call check_condition(status == CONTINUATION_OK .and. abs(x(1) - 2) < 1.0e-12_dp, &
        "the continuation reaches the target")
    call check_condition(report%failed_stages > 0 .and. report%reached == 1.0_dp .and. &
        report%stages > report%failed_stages + 1, "with halved steps")
    print '(a,i0,a,i0,a,i0)', 'kmax 10^1.5: stages ', report%stages, ', failed ', &
        report%failed_stages, ', Newton steps ', report%total_steps
    first = x(1)
    x = 1.9_dp
    options%max_step = 0.02_dp
    call continuation_run(stage, x, options, report, status)
    call check_condition(status == CONTINUATION_OK .and. report%failed_stages == 0 .and. &
        abs(x(1) - first) < 1.0e-12_dp, "the end state does not depend on the path")
    options%max_step = 1.0_dp

    ! 2
    kmax = 1.0e4_dp
    x = 1.9_dp
    call continuation_run(stage, x, options, report, status)
    call check_condition(status == CONTINUATION_STALLED .and. report%reached < 1 .and. &
        report%reached > 0 .and. abs(x(1) - (1 + report%reached)) < 1.0e-12_dp, &
        "a basin below the smallest step: stalled at the last converged state")

    ! 3
    x = 5.0_dp
    call continuation_run(stage, x, options, report, status)
    call check_condition(status == CONTINUATION_NO_BASE .and. x(1) == 5.0_dp .and. &
        report%reached == -1.0_dp, "no base solution: the start is returned")

    ! 4
    options%try_target_first = .true.
    x = 2.00001_dp
    call continuation_run(stage, x, options, report, status)
    call check_condition(status == CONTINUATION_OK .and. report%stages == 1 .and. &
        abs(x(1) - 2) < 1.0e-12_dp, "the target first: one stage from a start in its basin")
    x = 1.9_dp
    kmax = 10.0_dp**1.5_dp
    call continuation_run(stage, x, options, report, status)
    call check_condition(status == CONTINUATION_OK .and. report%failed_stages > 0 .and. &
        abs(x(1) - 2) < 1.0e-12_dp, "the target first, then the base and the path")
    options%try_target_first = .false.

    ! 5
    options%min_step = 0.0_dp
    call continuation_run(stage, x, options, report, status)
    call check_condition(status == CONTINUATION_BAD_INPUT, "a zero smallest step is refused")
    call check_summary("Parameter continuation")

contains

    subroutine stage(s, y, newton_steps, stage_status)
        !! Undamped Newton on atan(K (y - c)) = 0, at most 30 steps.
        real(dp), intent(in) :: s
        real(dp), intent(inout) :: y(:)
        integer, intent(out) :: newton_steps, stage_status
        real(dp) :: k, c, t

        k = kmax**s
        c = 1 + s
        stage_status = 1
        do newton_steps = 0, 30
            t = k*(y(1) - c)
            if (.not. abs(t) < 1.0e6_dp) return
            if (abs(atan(t)) < 1.0e-14_dp) then
                stage_status = 0
                return
            end if
            y(1) = y(1) - atan(t)*(1 + t**2)/k
        end do
    end subroutine stage

end program test_parameter_continuation
