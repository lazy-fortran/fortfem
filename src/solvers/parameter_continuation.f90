module fortfem_parameter_continuation
    !! Natural-parameter continuation with step control
    !! (doc/CONSTRAINED_NEWTON.md): a problem family P(s), s in [0, 1], with
    !! P(0) easy and P(1) the target (a stiff parameter raised from order
    !! one, an anisotropy, a load). The caller owns the nonlinear solve of
    !! one stage; this driver owns only the sequence of s:
    !!
    !! - P(0) is solved from the start. With try_target_first P(1) is tried
    !!   from the start before, and P(0) only when that solve fails (a stiff
    !!   parameter usually makes the direct solve the most expensive one);
    !! - from the last converged state at s_ok the stage s_ok + ds is
    !!   solved (ds = max_step at first); a converged stage doubles ds up to
    !!   max_step, a failed one halves it and is retried from the same
    !!   state, and the continuation fails when a stage fails with ds below
    !!   min_step.
    !!
    !! The end state is a solution of P(1) and does not depend on the path
    !! unless P(1) has several solutions. No stage is set by hand.
    use fortfem_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: CONTINUATION_OK = 0, CONTINUATION_NO_BASE = 1, &
                                  CONTINUATION_STALLED = 2, CONTINUATION_BAD_INPUT = 3

    type, public :: continuation_options_t
        logical :: try_target_first = .false.
        real(dp) :: max_step = 1.0_dp, min_step = 1.0_dp/64
    end type continuation_options_t

    type, public :: continuation_report_t
        !! stages solved (failed ones included), of them failed, and the
        !! nonlinear steps of all of them
        integer :: stages = 0, failed_stages = 0, total_steps = 0
        !! the parameter of the last converged stage (-1: none)
        real(dp) :: reached = -1.0_dp
    end type continuation_report_t

    abstract interface
        subroutine continuation_stage(s, x, steps, status)
            !! Solves P(s) from x. status 0: x is the solution; otherwise x
            !! is discarded. steps: nonlinear steps spent (for the report).
            import :: dp
            real(dp), intent(in) :: s
            real(dp), intent(inout) :: x(:)
            integer, intent(out) :: steps, status
        end subroutine continuation_stage
    end interface

    public :: continuation_stage, continuation_run

contains

    subroutine continuation_run(stage, x, options, report, status)
        !! From x to the solution of P(1); on failure x is the last
        !! converged state (the start when none converged).
        procedure(continuation_stage) :: stage
        real(dp), intent(inout) :: x(:)
        type(continuation_options_t), intent(in) :: options
        type(continuation_report_t), intent(out) :: report
        integer, intent(out) :: status
        real(dp), allocatable :: xs(:)
        real(dp) :: s, s_ok, ds
        integer :: st

        status = CONTINUATION_BAD_INPUT
        if (.not. (options%max_step > 0 .and. options%min_step > 0)) return
        if (options%try_target_first) then
            call attempt(1.0_dp)
            if (st == 0) then
                status = CONTINUATION_OK
                return
            end if
        end if
        status = CONTINUATION_NO_BASE
        call attempt(0.0_dp)
        if (st /= 0) return
        status = CONTINUATION_STALLED
        s_ok = 0
        ds = min(options%max_step, 1.0_dp)
        do while (s_ok < 1)
            s = min(1.0_dp, s_ok + ds)
            call attempt(s)
            if (st == 0) then
                s_ok = s
                ds = min(options%max_step, 2*ds)
            else
                if (ds < options%min_step) return
                ds = ds/2
            end if
        end do
        status = CONTINUATION_OK

    contains

        subroutine attempt(at)
            real(dp), intent(in) :: at
            integer :: steps

            xs = x
            steps = 0
            call stage(at, xs, steps, st)
            report%stages = report%stages + 1
            report%total_steps = report%total_steps + steps
            if (st == 0) then
                x = xs
                report%reached = at
            else
                report%failed_stages = report%failed_stages + 1
            end if
        end subroutine attempt
    end subroutine continuation_run

end module fortfem_parameter_continuation
