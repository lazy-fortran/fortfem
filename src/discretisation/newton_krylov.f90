module fortfem_newton_krylov
    ! Inexact Newton-Krylov driver for nonlinear systems F(u) = 0 that are only
    ! available through callbacks: residual, linearisation at an iterate,
    ! Jacobian action and preconditioner action. Right-preconditioned flexible
    ! GMRES(m), Eisenstat-Walker forcing (choice 2) and a backtracking line
    ! search on ||F||. Physics-free; operators are supplied by the extension.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: nonlinear_system_t, newton_options_t, newton_info_t
    public :: newton_krylov_solve, NEWTON_PHASES

    integer, parameter :: NEWTON_PHASES = 5
    integer, parameter :: PHASE_RESIDUAL = 1, PHASE_LINEARISE = 2, &
                          PHASE_JACOBIAN = 3, PHASE_PRECONDITION = 4, &
                          PHASE_KRYLOV = 5

    type, abstract :: nonlinear_system_t
    contains
        procedure(residual_interface), deferred :: residual
        procedure(linearise_interface), deferred :: linearise
        procedure(apply_interface), deferred :: jacobian
        procedure(apply_interface), deferred :: precondition
    end type

    abstract interface
        subroutine residual_interface(this, u, r, status)
            import nonlinear_system_t, dp
            class(nonlinear_system_t), intent(inout) :: this
            real(dp), intent(in) :: u(:)
            real(dp), intent(out) :: r(:)
            integer, intent(out) :: status
        end subroutine
        subroutine linearise_interface(this, u, status)
            import nonlinear_system_t, dp
            class(nonlinear_system_t), intent(inout) :: this
            real(dp), intent(in) :: u(:)
            integer, intent(out) :: status
        end subroutine
        subroutine apply_interface(this, v, w, status)
            import nonlinear_system_t, dp
            class(nonlinear_system_t), intent(inout) :: this
            real(dp), intent(in) :: v(:)
            real(dp), intent(out) :: w(:)
            integer, intent(out) :: status
        end subroutine
    end interface

    type :: newton_options_t
        real(dp) :: rtol = 1.e-10_dp ! stop when ||F|| <= max(atol, rtol*||F0||)
        real(dp) :: atol = 0.0_dp
        integer :: max_newton = 50
        integer :: restart = 30 ! GMRES restart length
        integer :: max_krylov = 300 ! Krylov iterations per Newton step
        real(dp) :: eta0 = 0.5_dp ! initial forcing term
        real(dp) :: eta_max = 0.9_dp
        real(dp) :: ew_gamma = 0.9_dp
        real(dp) :: ew_alpha = 2.0_dp
        logical :: constant_forcing = .false. ! use eta0 at every step
        real(dp) :: armijo = 1.e-4_dp
        real(dp) :: lambda_min = 1.e-4_dp
        integer :: max_backtracks = 20
    end type

    type :: newton_info_t
        ! status: 0 converged, 1 invalid input, 2 Newton limit, 3 line search
        ! failure, 4 callback error, 5 non-finite residual or Krylov breakdown
        integer :: status = 1, newton = 0, krylov = 0, residuals = 0
        integer :: jacobians = 0, preconditioners = 0, backtracks = 0
        real(dp) :: initial = huge(1.0_dp), residual = huge(1.0_dp)
        real(dp) :: time(NEWTON_PHASES) = 0 ! residual, linearise, J, M, Krylov rest
        real(dp) :: total = 0
        real(dp), allocatable :: history(:), forcing(:)
        integer, allocatable :: krylov_history(:)
    end type

contains

    function wall() result(t)
        real(dp) :: t
        integer(int64) :: count, rate
        call system_clock(count, rate)
        t = real(count, dp)/real(rate, dp)
    end function

    function pdot(a, b) result(s)
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: s
        integer :: i
        s = 0
        !$omp parallel do schedule(static) reduction(+:s)
        do i = 1, size(a)
            s = s + a(i)*b(i)
        end do
        !$omp end parallel do
    end function

    function pnorm(a) result(s)
        real(dp), intent(in) :: a(:)
        real(dp) :: s
        s = sqrt(pdot(a, a))
    end function

    subroutine paxpy(alpha, x, y)
        real(dp), intent(in) :: alpha, x(:)
        real(dp), intent(inout) :: y(:)
        integer :: i
        !$omp parallel do schedule(static)
        do i = 1, size(x)
            y(i) = y(i) + alpha*x(i)
        end do
        !$omp end parallel do
    end subroutine

    subroutine call_residual(system, u, r, info, status)
        class(nonlinear_system_t), intent(inout) :: system
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: r(:)
        type(newton_info_t), intent(inout) :: info
        integer, intent(out) :: status
        real(dp) :: t0
        t0 = wall()
        call system%residual(u, r, status)
        info%time(PHASE_RESIDUAL) = info%time(PHASE_RESIDUAL) + wall() - t0
        info%residuals = info%residuals + 1
    end subroutine

    subroutine fgmres(system, b, x, tol, options, info, iterations, status)
        ! Solve J x = b with right preconditioning, x0 = 0, to ||b - J x|| <= tol.
        class(nonlinear_system_t), intent(inout) :: system
        real(dp), intent(in) :: b(:), tol
        real(dp), intent(out) :: x(:)
        type(newton_options_t), intent(in) :: options
        type(newton_info_t), intent(inout) :: info
        integer, intent(out) :: iterations, status
        real(dp), allocatable :: v(:, :), z(:, :), h(:, :), g(:), cs(:), sn(:), &
                                 y(:), r(:), w(:)
        real(dp) :: beta, temp, t0, t1, resid
        integer :: m, n, i, j, k
        n = size(b)
        m = max(1, options%restart)
        allocate (v(n, m + 1), z(n, m), h(m + 1, m), g(m + 1), cs(m), sn(m), &
                  y(m), r(n), w(n))
        t1 = wall()
        x = 0
        r = b
        beta = pnorm(r)
        iterations = 0
        status = 0
        resid = beta
        do while (resid > tol .and. iterations < options%max_krylov)
            v(:, 1) = r/beta
            g = 0
            g(1) = beta
            h = 0
            k = 0
            do j = 1, m
                t0 = wall()
                call system%precondition(v(:, j), z(:, j), status)
                info%time(PHASE_PRECONDITION) = info%time(PHASE_PRECONDITION) + &
                                                wall() - t0
                info%preconditioners = info%preconditioners + 1
                if (status /= 0) then
                    status = 4
                    return
                end if
                t0 = wall()
                call system%jacobian(z(:, j), w, status)
                info%time(PHASE_JACOBIAN) = info%time(PHASE_JACOBIAN) + wall() - t0
                info%jacobians = info%jacobians + 1
                if (status /= 0) then
                    status = 4
                    return
                end if
                ! Modified Gram-Schmidt with one reorthogonalisation pass.
                do i = 1, j
                    h(i, j) = pdot(w, v(:, i))
                    call paxpy(-h(i, j), v(:, i), w)
                end do
                do i = 1, j
                    temp = pdot(w, v(:, i))
                    h(i, j) = h(i, j) + temp
                    call paxpy(-temp, v(:, i), w)
                end do
                h(j + 1, j) = pnorm(w)
                if (.not. ieee_is_finite(h(j + 1, j))) then
                    status = 5
                    return
                end if
                if (h(j + 1, j) > 0) v(:, j + 1) = w/h(j + 1, j)
                do i = 1, j - 1
                    temp = cs(i)*h(i, j) + sn(i)*h(i + 1, j)
                    h(i + 1, j) = -sn(i)*h(i, j) + cs(i)*h(i + 1, j)
                    h(i, j) = temp
                end do
                temp = hypot(h(j, j), h(j + 1, j))
                if (temp <= 0) then
                    status = 5
                    return
                end if
                cs(j) = h(j, j)/temp
                sn(j) = h(j + 1, j)/temp
                h(j, j) = temp
                h(j + 1, j) = 0
                g(j + 1) = -sn(j)*g(j)
                g(j) = cs(j)*g(j)
                iterations = iterations + 1
                k = j
                resid = abs(g(j + 1))
                if (resid <= tol .or. iterations >= options%max_krylov) exit
                if (h(j, j) == 0) exit
            end do
            do i = k, 1, -1
                y(i) = (g(i) - dot_product(h(i, i + 1:k), y(i + 1:k)))/h(i, i)
            end do
            do i = 1, k
                call paxpy(y(i), z(:, i), x)
            end do
            if (resid <= tol .or. iterations >= options%max_krylov) exit
            ! Restart from the true residual.
            t0 = wall()
            call system%jacobian(x, w, status)
            info%time(PHASE_JACOBIAN) = info%time(PHASE_JACOBIAN) + wall() - t0
            info%jacobians = info%jacobians + 1
            if (status /= 0) then
                status = 4
                return
            end if
            r = b - w
            beta = pnorm(r)
            resid = beta
            if (beta <= 0) exit
        end do
        info%time(PHASE_KRYLOV) = info%time(PHASE_KRYLOV) + wall() - t1
    end subroutine

    subroutine newton_krylov_solve(system, u, options, info)
        class(nonlinear_system_t), intent(inout) :: system
        real(dp), intent(inout) :: u(:)
        type(newton_options_t), intent(in) :: options
        type(newton_info_t), intent(out) :: info
        real(dp), allocatable :: f(:), d(:), trial(:), ftrial(:), rhs(:)
        real(dp) :: fnorm, target, eta, eta_previous, lambda, ftrial_norm, &
                    previous_norm, t_start, t0, safeguard, lambda_next, phase_sum
        integer :: n, status, its, k
        t_start = wall()
        n = size(u)
        allocate (info%history(0:options%max_newton), &
                  info%forcing(options%max_newton), &
                  info%krylov_history(options%max_newton))
        info%history = 0; info%forcing = 0; info%krylov_history = 0
        if (n < 1 .or. options%rtol < 0 .or. options%atol < 0) return
        allocate (f(n), d(n), trial(n), ftrial(n), rhs(n))
        call call_residual(system, u, f, info, status)
        if (status /= 0) then
            info%status = 4
            return
        end if
        fnorm = pnorm(f)
        if (.not. ieee_is_finite(fnorm)) then
            info%status = 5
            return
        end if
        info%initial = fnorm
        info%residual = fnorm
        info%history(0) = fnorm
        target = max(options%atol, options%rtol*fnorm)
        eta = min(options%eta0, options%eta_max)
        eta_previous = eta
        previous_norm = fnorm
        info%status = 2
        do k = 1, options%max_newton
            if (fnorm <= target) then
                info%status = 0
                exit
            end if
            t0 = wall()
            call system%linearise(u, status)
            info%time(PHASE_LINEARISE) = info%time(PHASE_LINEARISE) + wall() - t0
            if (status /= 0) then
                info%status = 4
                exit
            end if
            if (k > 1 .and. .not. options%constant_forcing) then
                ! Eisenstat-Walker choice 2 with the standard safeguards.
                eta = options%ew_gamma*(fnorm/previous_norm)**options%ew_alpha
                safeguard = options%ew_gamma*eta_previous**options%ew_alpha
                if (safeguard > 0.1_dp) eta = max(eta, safeguard)
                eta = min(eta, options%eta_max)
            end if
            ! Do not oversolve beyond the nonlinear target.
            eta = max(eta, 0.5_dp*target/fnorm)
            eta = min(eta, options%eta_max)
            info%forcing(k) = eta
            rhs = -f
            call fgmres(system, rhs, d, eta*fnorm, options, info, its, status)
            info%krylov = info%krylov + its
            info%krylov_history(k) = its
            if (status /= 0) then
                info%status = status
                exit
            end if
            lambda = 1
            do
                trial = u + lambda*d
                call call_residual(system, trial, ftrial, info, status)
                if (status /= 0) then
                    ftrial_norm = huge(1.0_dp)
                else
                    ftrial_norm = pnorm(ftrial)
                    if (.not. ieee_is_finite(ftrial_norm)) ftrial_norm = huge(1.0_dp)
                end if
                if (ftrial_norm <= (1 - options%armijo*lambda*(1 - eta))*fnorm) exit
                if (lambda <= options%lambda_min) exit
                info%backtracks = info%backtracks + 1
                ! Safeguarded parabolic model of ||F(u + lambda d)||^2.
                lambda_next = 0.5_dp*lambda
                if (ftrial_norm < huge(1.0_dp)) then
                    lambda_next = lambda**2*fnorm**2/(ftrial_norm**2 - &
                                                      fnorm**2 + 2*lambda*fnorm**2)
                end if
                lambda = min(0.5_dp*lambda, max(0.1_dp*lambda, lambda_next))
            end do
            if (ftrial_norm > (1 - options%armijo*lambda*(1 - eta))*fnorm) then
                info%status = 3
                exit
            end if
            u = trial
            f = ftrial
            previous_norm = fnorm
            fnorm = ftrial_norm
            eta_previous = eta
            info%newton = k
            info%history(k) = fnorm
            info%residual = fnorm
        end do
        if (info%status == 2 .and. fnorm <= target) info%status = 0
        info%total = wall() - t_start
        ! Krylov time excludes the separately recorded operator phases.
        phase_sum = info%time(PHASE_JACOBIAN) + info%time(PHASE_PRECONDITION)
        info%time(PHASE_KRYLOV) = max(0.0_dp, info%time(PHASE_KRYLOV) - phase_sum)
    end subroutine
end module fortfem_newton_krylov
