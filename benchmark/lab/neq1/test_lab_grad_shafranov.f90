program test_lab_grad_shafranov
    ! NEQ1 Solov'ev oracle (Cerfon & Freidberg 2010) on the curved through-axis
    ! O-grid: the oracle satisfies Delta* psibar = (1 - A) x^2 + A and vanishes
    ! at the three shaping points; the free-axis Newton-Krylov solve converges
    ! and its flux, axis position and axis value approach the exact solution.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_newton_krylov, only: newton_options_t, newton_info_t, &
                                     newton_krylov_solve
    use fortfem_tensor_operator, only: tensor_l2_error
    use fortfem_lab_gs, only: gs_case_t, gs_system_t, read_gs_case, &
                              solovev_psibar, solovev_psi, build_gs_system, &
                              gs_initial_guess, gs_safety_factor
    implicit none
    type(gs_case_t) :: case
    character(512) :: path
    real(dp) :: errors(2), axis_errors(2)
    integer :: status, failures, k
    failures = 0
    path = 'benchmark/lab/neq1/solovev_iter.nml' ! default: run from the repository root
    if (command_argument_count() >= 1) call get_command_argument(1, path)
    call read_gs_case(trim(path), case, status)
    call expect(status == 0, 'read case namelist')
    if (status /= 0) error stop 'test_lab_grad_shafranov failed'
    call check_oracle()
    do k = 1, 2
        call solve(2*k + 2, errors(k), axis_errors(k))
    end do
    call expect(errors(2) < 0.1_dp*errors(1), 'p=6 flux error well below p=4')
    call expect(errors(2) < 5.e-5_dp, 'p=6 relative flux error < 5e-5')
    call expect(axis_errors(2) < 1.e-7_dp, 'p=6 axis position error < 1e-7 m')
    if (failures > 0) error stop 'test_lab_grad_shafranov failed'
    print '(a)', 'test_lab_grad_shafranov passed'
contains
    subroutine expect(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', label
        end if
    end subroutine

    subroutine check_oracle()
        real(dp) :: x, y, h, lap, worst
        integer :: i
        worst = 0
        h = 1.e-3_dp
        do i = 1, 7
            x = 1 + 0.2_dp*sin(real(i, dp))
            y = 0.3_dp*cos(real(2*i, dp))
            ! Delta* = x d/dx (1/x d/dx) + d2/dy2, fourth-order differences.
            lap = (-p(x + 2*h, y) + 16*p(x + h, y) - 30*p(x, y) + &
                   16*p(x - h, y) - p(x - 2*h, y))/(12*h**2) - &
                  (-p(x + 2*h, y) + 8*p(x + h, y) - 8*p(x - h, y) + &
                   p(x - 2*h, y))/(12*h*x) + &
                  (-p(x, y + 2*h) + 16*p(x, y + h) - 30*p(x, y) + &
                   16*p(x, y - h) - p(x, y - 2*h))/(12*h**2)
            worst = max(worst, abs(lap - ((1 - case%a)*x**2 + case%a)))
        end do
        call expect(worst < 1.e-7_dp, 'Solovev oracle satisfies Delta*')
        call expect(abs(p(1 + case%eps, 0.0_dp)) < 1.e-14_dp, 'outer point')
        call expect(abs(p(1 - case%eps, 0.0_dp)) < 1.e-14_dp, 'inner point')
        call expect(abs(p(1 - case%delta*case%eps, case%kappa*case%eps)) < &
                    1.e-14_dp, 'top point')
    end subroutine

    function p(x, y) result(value)
        real(dp), intent(in) :: x, y
        real(dp) :: value
        value = solovev_psibar(case, x, y)
    end function

    subroutine solve(degree, relative, axis_error)
        integer, intent(in) :: degree
        real(dp), intent(out) :: relative, axis_error
        type(gs_system_t) :: system
        type(newton_options_t) :: options
        type(newton_info_t) :: info
        real(dp), allocatable :: u(:), zero(:)
        real(dp) :: q
        integer :: status, its
        call build_gs_system(case, 2, degree, system, status)
        call expect(status == 0, 'build Grad-Shafranov system')
        allocate (u(system%op%mesh%ndof), zero(system%op%mesh%ndof))
        call gs_initial_guess(system, u, 1.e-3_dp, its, status)
        call expect(status == 0, 'flat-profile initial guess')
        options%rtol = 1.e-11_dp
        options%atol = 1.e-14_dp
        call newton_krylov_solve(system, u, options, info)
        call expect(info%status == 0, 'free-axis Newton-Krylov converged')
        call expect(info%newton <= 6, 'at most six Newton steps')
        zero = 0
        relative = tensor_l2_error(system%op, u, solovev_psi)/ &
                   tensor_l2_error(system%op, zero, solovev_psi)
        axis_error = norm2(system%axis_x - [case%r_axis, case%z_axis])
        call gs_safety_factor(system, u, 0.5_dp*system%psia, 128, q, status)
        call expect(status == 0, 'safety-factor contour')
        print '(a,i0,a,i0,a,i0,a,i0,a,es10.3,a,es10.3,a,f10.6)', 'p=', degree, &
            ' dofs=', system%op%mesh%ndof, ' newton=', info%newton, ' krylov=', &
            info%krylov, ' rel_l2=', relative, ' axis=', axis_error, &
            ' q(psin=0.5)=', q
    end subroutine
end program
