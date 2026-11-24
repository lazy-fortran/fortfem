module fortfem_constrained_newton
    !! Newton's method for a stationary point of an energy W(z) under linear
    !! equality constraints C z = d and fixed slots (doc/CONSTRAINED_NEWTON.md).
    !! A step solves the KKT system
    !!
    !!     [ H + shift D   C^T ] [ dz     ]   [ -g      ]
    !!     [ C             0   ] [ lambda ] = [ d - C z ]
    !!
    !! on the free slots by the symmetric LDL^T factorisation of
    !! fortfem_sparse_direct (H the Hessian given as its upper triangle, D a
    !! multiple of the identity). The factorisation gives the inertia, and
    !! with m independent constraints the Hessian on the tangent space of the
    !! constraints has (negative eigenvalues of the KKT matrix) - m negative
    !! eigenvalues, which is reported at convergence. A null direction of H
    !! that the constraints cross (a continuous symmetry of W with constraints
    !! transverse to its orbits) leaves the KKT matrix regular.
    !!
    !! Step acceptance (newton_options_t%acceptance):
    !!
    !! - NEWTON_ACCEPT_ENERGY: the full step is taken when W decreases by a
    !!   fraction of the decrease of the quadratic model; otherwise the shift
    !!   is raised. For minima.
    !! - NEWTON_ACCEPT_ENERGY_RELAXED: as NEWTON_ACCEPT_ENERGY, with relaxed
    !!   steps (the watchdog technique of Chamberlain et al. 1982). A Newton
    !!   step along a curved valley of a stiff term leaves the valley by the
    !!   curvature and pays the stiff energy for it, so the energy test alone
    !!   forces short steps although the next step returns to the valley. A
    !!   step that fails the energy test while the matrix is positive on the
    !!   constraints is therefore taken on trial, and kept when the step
    !!   after it brings W below the value required of the first; otherwise
    !!   the state returns to the start of the pair, the shift is raised, and
    !!   no further trial is made for a number of accepted steps that doubles
    !!   with every failed pair.
    !! - NEWTON_ACCEPT_GRADIENT: the step a dz, a = 1, 1/2, ..., is taken
    !!   when the scaled gradient of the Lagrangian falls below (1 - a/4)
    !!   times its value. For the stationary point next to the start,
    !!   whatever the signature of the Hessian.
    !!
    !! Converged when the gradient of the Lagrangian g + C^T lambda, scaled
    !! by the square roots of the Hessian diagonal, is below the tolerance.
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
        sparse_direct_factor_symmetric_csc, sparse_direct_refactor_csc, &
        sparse_direct_solve_factored, sparse_direct_inertia, sparse_direct_free, &
        SPARSE_DIRECT_LOWER
    use fortsparse_csr, only: csr_t
    implicit none
    private

    integer, parameter, public :: NEWTON_ACCEPT_ENERGY = 1, &
                                  NEWTON_ACCEPT_ENERGY_RELAXED = 2, &
                                  NEWTON_ACCEPT_GRADIENT = 3
    !! status of constrained_newton
    integer, parameter, public :: NEWTON_CONVERGED = 0, NEWTON_NO_STEP = 1, &
                                  NEWTON_STEP_LIMIT = 2, NEWTON_BAD_INPUT = 3

    type, public :: newton_constraints_t
        !! Rows of C by rows (CSR): sum_k val(k) z(col(k)) = target(row) for k
        !! in ptr(row):ptr(row + 1) - 1.
        integer, allocatable :: ptr(:), col(:)
        real(dp), allocatable :: val(:), target(:)
    end type newton_constraints_t

    type, public :: newton_options_t
        real(dp) :: tolerance = 1.0e-9_dp
        integer :: max_steps = 60
        integer :: acceptance = NEWTON_ACCEPT_ENERGY
        !! inertia of the unshifted KKT matrix at the converged state
        logical :: inertia = .true.
        logical :: verbose = .false.
        !! Retain damping after a poorly predicted energy step.
        logical :: model_agreement_damping = .false.
        !! Initial multiplier of the smallest absolute free Hessian diagonal
        !! for energy minimization. The gradient search starts unshifted.
        real(dp) :: initial_damping = 1.0e-9_dp
    end type newton_options_t

    type, public :: newton_report_t
        !! accepted steps, and those of them accepted by the natural
        !! monotonicity test after the energy test had failed
        integer :: steps = 0, relaxed_steps = 0, rejected_pairs = 0
        real(dp) :: energy = 0.0_dp
        !! scaled gradient of the Lagrangian and largest constraint defect
        real(dp) :: gradient_norm = 0.0_dp, constraint_norm = 0.0_dp
        !! largest multiplier (zero for constraints transverse to an exact
        !! symmetry)
        real(dp) :: multiplier_norm = 0.0_dp
        logical :: converged = .false.
        !! Eigenvalue counts of the Hessian on the tangent space of the
        !! constraints at the converged state (-1: not computed) and whether
        !! every pivot was clear of its rounding threshold.
        integer :: negative_eigenvalues = -1, zero_eigenvalues = -1
        logical :: inertia_reliable = .false.
        !! seconds in the system with Hessian, the system without, the
        !! factorisations and solves
        real(dp) :: time_hessian = 0.0_dp, time_energy = 0.0_dp, time_solve = 0.0_dp
    end type newton_report_t

    abstract interface
        subroutine newton_system(z, energy, gradient, hessian)
            !! W at z; with gradient also dW/dz; with hessian also the UPPER
            !! triangle (columns >= row, every diagonal entry present) of
            !! the Hessian by rows, on a pattern that does not depend on z.
            import :: dp, csr_t
            real(dp), intent(in) :: z(:)
            real(dp), intent(out) :: energy
            real(dp), intent(out), optional :: gradient(:)
            type(csr_t), intent(out), optional :: hessian
        end subroutine newton_system
    end interface

    type :: kkt_t
        type(sparse_direct_factor_t) :: factor
        logical :: factored = .false.
        integer :: n = 0, nf = 0, m = 0
        integer, allocatable :: free(:), local(:)
        !! C on the free slots by columns: rows and the positions of the
        !! values in newton_constraints_t%val
        integer, allocatable :: cptr(:), crow(:), csrc(:)
        !! constraints with a nonzero row on the free slots at the last
        !! factorisation
        integer :: active = 0
        logical, allocatable :: live(:)
    end type kkt_t

    abstract interface
        subroutine newton_constraint_update(z, constraints)
            !! Constraints that follow the state: new values and targets of
            !! the rows at the iterate z, on the same pattern (ptr and col
            !! are not changed).
            import :: dp, newton_constraints_t
            real(dp), intent(in) :: z(:)
            type(newton_constraints_t), intent(inout) :: constraints
        end subroutine newton_constraint_update
    end interface

    public :: newton_system, newton_constraint_update, newton_constraints_init, &
              newton_constraints_add, constrained_newton

contains

    subroutine newton_constraints_init(c)
        type(newton_constraints_t), intent(out) :: c

        allocate (c%ptr(1), source=1)
        allocate (c%col(0), c%val(0), c%target(0))
    end subroutine newton_constraints_init

    subroutine newton_constraints_add(c, cols, vals, target)
        !! Adds the row sum vals z(cols) = target.
        type(newton_constraints_t), intent(inout) :: c
        integer, intent(in) :: cols(:)
        real(dp), intent(in) :: vals(:), target

        if (.not. allocated(c%ptr)) call newton_constraints_init(c)
        c%col = [c%col, cols]
        c%val = [c%val, vals]
        c%target = [c%target, target]
        c%ptr = [c%ptr, size(c%col) + 1]
    end subroutine newton_constraints_add

    subroutine constrained_newton(system, z, fixed, constraints, options, report, &
                                  status, update)
        !! From z to a stationary point of W on {C z = d}, the slots fixed
        !! kept. A start off the constraints is moved onto them by the first
        !! full step. A row of C that vanishes on the free slots is inactive.
        !!
        !! With update the constraints follow the state: update is called at
        !! every iterate z_k and returns C(z_k) and d(z_k). With d = C(z_k)
        !! z_k every step is in the null space of C(z_k), e.g. orthogonal to
        !! the generators of a symmetry at the current state; the limit then
        !! satisfies g + C(z)^T lambda = 0 with the rows at the limit.
        procedure(newton_system) :: system
        real(dp), intent(inout) :: z(:)
        integer, intent(in) :: fixed(:)
        type(newton_constraints_t), intent(inout) :: constraints
        type(newton_options_t), intent(in) :: options
        type(newton_report_t), intent(out) :: report
        integer, intent(out) :: status
        procedure(newton_constraint_update), optional :: update
        type(kkt_t) :: kkt
        type(csr_t) :: hessian
        real(dp), allocatable :: grad(:), gt(:), diag(:), delta(:), hd(:), zt(:), &
                                 lambda(:), defect(:), lag(:)
        real(dp), allocatable :: zw(:)
        real(dp) :: w, wt, lam, dmin, pred, alpha, shift, ww, predw, lamw, reference, needed
        integer :: it, n, m, st, nneg, nzero, npos
        logical :: accepted, sure, tested, watch
        integer :: blocked, pause
        integer(int64) :: c0, rate

        n = size(z)
        m = 0
        if (allocated(constraints%target)) m = size(constraints%target)
        status = NEWTON_BAD_INPUT
        if (options%acceptance < NEWTON_ACCEPT_ENERGY .or. &
            options%acceptance > NEWTON_ACCEPT_GRADIENT) return
        if (options%acceptance /= NEWTON_ACCEPT_GRADIENT) then
            if (.not. ieee_is_finite(options%initial_damping)) return
            if (options%initial_damping <= 0) return
        end if
        if (any(fixed < 1) .or. any(fixed > n)) return
        if (m > 0) then
            if (any(constraints%col < 1) .or. any(constraints%col > n)) return
        end if
        call kkt_init(kkt, n, fixed, constraints, m, st)
        if (st /= 0) return
        allocate (grad(n), gt(n), diag(n), delta(n), hd(n), zt(n), lambda(m), defect(m), &
                  lag(n))
        call system_clock(c0, rate)
        lam = 0.0_dp
        if (options%acceptance /= NEWTON_ACCEPT_GRADIENT) lam = options%initial_damping
        status = NEWTON_STEP_LIMIT
        watch = .false.
        blocked = 0
        pause = 0
        do it = 0, options%max_steps
            call tick()
            call system(z, w, grad, hessian)
            call tock(report%time_hessian)
            call diagonal(hessian, diag)
            dmin = minval(abs(diag(kkt%free)))
            if (present(update)) call update(z, constraints)
            defect = constraint_defect(constraints, m, z)
            report%energy = w
            accepted = .false.
            tested = .false.
            do
                shift = lam*dmin
                call tick()
                call kkt_solve(kkt, hessian, shift, constraints, -grad, -defect, &
                               .true., delta, lambda, st)
                call tock(report%time_solve)
                if (st /= 0) then
                    ! singular: only a shift can help
                    if (options%acceptance == NEWTON_ACCEPT_GRADIENT .and. lam == 0) &
                        lam = 1.0e-10_dp
                    lam = 10*lam
                    if (lam > 1.0e12_dp) exit
                    cycle
                end if
                if (.not. tested) then
                    ! gradient of the Lagrangian with the multipliers of this
                    ! step: it vanishes with the step
                    tested = .true.
                    lag = grad + constraint_transpose(constraints, m, n, lambda)
                    report%gradient_norm = &
                        norm2(lag(kkt%free)/sqrt(abs(diag(kkt%free)) + tiny(1.0_dp)))
                    report%constraint_norm = maxval(abs(defect), 1)
                    if (m == 0) report%constraint_norm = 0
                    report%multiplier_norm = maxval(abs(lambda), 1)
                    if (m == 0) report%multiplier_norm = 0
                    if (options%verbose) then
                        call sparse_direct_inertia(kkt%factor, nneg, nzero, npos, st, sure)
                        print '(a,i4,a,es20.12,a,es10.2,a,es9.1,a,es9.1,a,i0,a,es9.2)', &
                            '    step', it, '  W', w, '  |g|', report%gradient_norm, &
                            '  shift', lam, '  multiplier', report%multiplier_norm, &
                            '  negative ', nneg - kkt%active, '  |dz|', norm2(delta)
                    end if
                    if (report%gradient_norm < options%tolerance .and. &
                        report%constraint_norm <= 1.0e-10_dp*(1 + maxval(abs(z)))) then
                        report%converged = .true.
                        status = NEWTON_CONVERGED
                        call inertia_at_convergence()
                        call kkt_free(kkt)
                        return
                    end if
                    if (it == options%max_steps) then
                        call kkt_free(kkt)
                        return
                    end if
                end if
                if (options%acceptance /= NEWTON_ACCEPT_GRADIENT) then
                    call upper_matvec(hessian, delta, hd)
                    pred = -(dot_product(grad, delta) + 0.5_dp*dot_product(delta, hd))
                    zt = z + delta
                    call tick()
                    call system(zt, wt)
                    call tock(report%time_energy)
                    reference = w
                    needed = pred
                    if (watch) then
                        reference = ww
                        needed = predw
                    end if
                    if (pred > 0 .and. ieee_is_finite(wt) .and. &
                        reference - wt > 1.0e-4_dp*needed - 1.0e-14_dp*abs(reference)) then
                        if (watch) then
                            report%relaxed_steps = report%relaxed_steps + 1
                            pause = pause/2
                        end if
                        watch = .false.
                        blocked = max(blocked - 1, 0)
                        if (options%model_agreement_damping) then
                            if (reference - wt > 0.75_dp*needed) lam = max(lam/10, 1.0e-12_dp)
                        else
                            lam = max(lam/100, 1.0e-12_dp)
                        end if
                        accepted = .true.
                        exit
                    end if
                    if (watch) then
                        ! the relaxed step did not pay off: back to its start
                        ! with a larger shift
                        zt = zw
                        lam = 10*lamw
                        watch = .false.
                        ! no further relaxed step for a number of accepted
                        ! steps that doubles with every failure
                        pause = 2*pause + 1
                        blocked = pause
                        report%steps = report%steps - 2
                        report%rejected_pairs = report%rejected_pairs + 1
                        accepted = .true.
                        exit
                    end if
                    if (options%acceptance == NEWTON_ACCEPT_ENERGY_RELAXED .and. &
                        blocked == 0 .and. pred > 0 .and. ieee_is_finite(wt)) then
                        ! W did not decrease as predicted although the matrix
                        ! is positive on the constraints: the step may have
                        ! left a curved valley of a stiff term, to which the
                        ! next step returns. Take it and judge the pair.
                        call sparse_direct_inertia(kkt%factor, nneg, nzero, npos, st, sure)
                        if (st == 0 .and. nneg == kkt%active .and. nzero == 0) then
                            watch = .true.
                            zw = z
                            ww = w
                            predw = pred
                            lamw = lam
                            accepted = .true.
                            exit
                        end if
                    end if
                    lam = 10*lam
                    if (lam > 1.0e12_dp) exit
                else
                    alpha = 1
                    do while (alpha > 1.0e-4_dp)
                        zt = z + alpha*delta
                        accepted = gradient_falls(alpha)
                        if (accepted) exit
                        alpha = alpha/2
                    end do
                    if (accepted) then
                        lam = lam/100
                        if (lam < 1.0e-12_dp) lam = 0
                        exit
                    end if
                    exit
                end if
            end do
            if (.not. accepted) then
                status = NEWTON_NO_STEP
                call kkt_free(kkt)
                return
            end if
            z = zt
            report%steps = report%steps + 1
        end do
        call kkt_free(kkt)

    contains

        subroutine tick()
            call system_clock(c0)
        end subroutine tick

        subroutine tock(acc)
            real(dp), intent(inout) :: acc
            integer(int64) :: now

            call system_clock(now)
            acc = acc + real(now - c0, dp)/rate
        end subroutine tock

        logical function gradient_falls(a)
            !! The scaled gradient of the Lagrangian at zt = z + a delta,
            !! with the multipliers of the step, is below (1 - a/4) times
            !! its value at z.
            real(dp), intent(in) :: a

            gradient_falls = .false.
            call tick()
            call system(zt, wt, gt)
            call tock(report%time_energy)
            if (.not. (ieee_is_finite(wt) .and. all(ieee_is_finite(gt)))) return
            lag = gt + constraint_transpose(constraints, m, n, lambda)
            gradient_falls = &
                norm2(lag(kkt%free)/sqrt(abs(diag(kkt%free)) + tiny(1.0_dp))) < &
                (1 - a/4)*report%gradient_norm
        end function gradient_falls

        subroutine inertia_at_convergence()
            !! Of the unshifted KKT matrix at the final state.
            if (.not. options%inertia) return
            call tick()
            if (shift /= 0) call kkt_solve(kkt, hessian, 0.0_dp, constraints, -grad, &
                                           -defect, .true., delta, lambda, st)
            ! an exactly singular matrix still has its counts
            call sparse_direct_inertia(kkt%factor, nneg, nzero, npos, st, sure)
            call tock(report%time_solve)
            if (st /= 0) return
            report%negative_eigenvalues = nneg - kkt%active
            report%zero_eigenvalues = nzero
            report%inertia_reliable = sure
            if (options%verbose) print '(a,i0,a,i0,a,i0,a,i0,a,l1)', &
                '    inertia of the KKT matrix (negative, zero, positive): ', nneg, &
                ', ', nzero, ', ', npos, '  constraints ', kkt%active, '  reliable ', sure
        end subroutine inertia_at_convergence
    end subroutine constrained_newton

    subroutine kkt_init(kkt, n, fixed, c, m, status)
        !! status /= 0: a row names a slot twice.
        type(kkt_t), intent(out) :: kkt
        integer, intent(in) :: n, fixed(:), m
        type(newton_constraints_t), intent(in) :: c
        integer, intent(out) :: status
        integer, allocatable :: count(:), seen(:)
        integer :: i, r, p, j

        status = 0
        kkt%n = n
        kkt%m = m
        allocate (kkt%local(n), seen(n))
        kkt%local = 1
        kkt%local(fixed) = 0
        kkt%free = pack([(i, i=1, n)], kkt%local == 1)
        kkt%nf = size(kkt%free)
        kkt%local = 0
        kkt%local(kkt%free) = [(i, i=1, kkt%nf)]
        allocate (kkt%cptr(kkt%nf + 1), count(kkt%nf))
        count = 0
        seen = 0
        do r = 1, m
            do p = c%ptr(r), c%ptr(r + 1) - 1
                if (seen(c%col(p)) == r) status = 1
                seen(c%col(p)) = r
                j = kkt%local(c%col(p))
                if (j > 0) count(j) = count(j) + 1
            end do
        end do
        kkt%cptr(1) = 1
        do j = 1, kkt%nf
            kkt%cptr(j + 1) = kkt%cptr(j) + count(j)
        end do
        allocate (kkt%crow(kkt%cptr(kkt%nf + 1) - 1), kkt%csrc(kkt%cptr(kkt%nf + 1) - 1))
        count = 0
        ! rows ascending within a column
        do r = 1, m
            do p = c%ptr(r), c%ptr(r + 1) - 1
                j = kkt%local(c%col(p))
                if (j == 0) cycle
                i = kkt%cptr(j) + count(j)
                kkt%crow(i) = r
                kkt%csrc(i) = p
                count(j) = count(j) + 1
            end do
        end do
    end subroutine kkt_init

    subroutine kkt_solve(kkt, h, shift, c, rhs, rhs_c, factorise, x, lambda, status)
        !! x on the free slots (zero on the fixed ones) and the multipliers
        !! of K (x, lambda) = (rhs, rhs_c); with factorise K is factorised
        !! first from the upper triangle h and the shift of its diagonal.
        type(kkt_t), intent(inout) :: kkt
        type(csr_t), intent(in) :: h
        real(dp), intent(in) :: shift, rhs(:), rhs_c(:)
        type(newton_constraints_t), intent(in) :: c
        logical, intent(in) :: factorise
        real(dp), intent(out) :: x(:), lambda(:)
        integer, intent(out) :: status
        integer, allocatable :: colptr(:), rowind(:)
        real(dp), allocatable :: val(:), b(:), y(:)
        integer :: nf, m, i, j, p, q, nz

        nf = kkt%nf
        m = kkt%m
        status = 0
        if (factorise) then
            nz = size(h%val) + size(kkt%csrc) + m
            allocate (colptr(nf + m + 1), rowind(nz), val(nz))
            if (.not. allocated(kkt%live)) allocate (kkt%live(m))
            kkt%live = .false.
            colptr(1) = 1
            nz = 0
            ! the rows of the upper triangle are the columns of the lower one
            do j = 1, nf
                i = kkt%free(j)
                do p = h%row_ptr(i), h%row_ptr(i + 1) - 1
                    q = kkt%local(h%col(p))
                    if (q == 0) cycle
                    nz = nz + 1
                    rowind(nz) = q
                    val(nz) = h%val(p)
                    if (q == j) val(nz) = val(nz) + shift
                end do
                do p = kkt%cptr(j), kkt%cptr(j + 1) - 1
                    nz = nz + 1
                    rowind(nz) = nf + kkt%crow(p)
                    val(nz) = c%val(kkt%csrc(p))
                    if (val(nz) /= 0) kkt%live(kkt%crow(p)) = .true.
                end do
                colptr(j + 1) = nz + 1
            end do
            do j = 1, m
                nz = nz + 1
                rowind(nz) = nf + j
                ! an inactive row keeps its multiplier at zero
                val(nz) = merge(0.0_dp, 1.0_dp, kkt%live(j))
                colptr(nf + j + 1) = nz + 1
            end do
            kkt%active = count(kkt%live)
            if (kkt%factored) then
                call sparse_direct_refactor_csc(kkt%factor, nf + m, colptr, &
                                                rowind(:nz), val(:nz), status)
            else
                call sparse_direct_factor_symmetric_csc(kkt%factor, nf + m, colptr, &
                    rowind(:nz), val(:nz), status, SPARSE_DIRECT_LOWER)
                kkt%factored = .true.
            end if
            if (status /= 0) return
        end if
        allocate (b(nf + m), y(nf + m))
        b(:nf) = rhs(kkt%free)
        b(nf + 1:) = rhs_c
        where (.not. kkt%live) b(nf + 1:) = 0
        call sparse_direct_solve_factored(kkt%factor, b, y, status)
        x = 0
        lambda = 0
        if (status /= 0) return
        if (.not. all(ieee_is_finite(y))) then
            status = -1
            return
        end if
        x(kkt%free) = y(:nf)
        lambda = y(nf + 1:)
    end subroutine kkt_solve

    subroutine kkt_free(kkt)
        type(kkt_t), intent(inout) :: kkt

        if (kkt%factored) call sparse_direct_free(kkt%factor)
        kkt%factored = .false.
    end subroutine kkt_free

    pure function constraint_defect(c, m, z) result(r)
        !! C z - d
        type(newton_constraints_t), intent(in) :: c
        integer, intent(in) :: m
        real(dp), intent(in) :: z(:)
        real(dp) :: r(m)
        integer :: row

        do row = 1, m
            r(row) = dot_product(c%val(c%ptr(row):c%ptr(row + 1) - 1), &
                                 z(c%col(c%ptr(row):c%ptr(row + 1) - 1))) - c%target(row)
        end do
    end function constraint_defect

    pure function constraint_transpose(c, m, n, lambda) result(v)
        !! C^T lambda
        type(newton_constraints_t), intent(in) :: c
        integer, intent(in) :: m, n
        real(dp), intent(in) :: lambda(:)
        real(dp) :: v(n)
        integer :: row, p

        v = 0
        do row = 1, m
            do p = c%ptr(row), c%ptr(row + 1) - 1
                v(c%col(p)) = v(c%col(p)) + c%val(p)*lambda(row)
            end do
        end do
    end function constraint_transpose

    pure subroutine upper_matvec(a, x, y)
        !! y = A x for the symmetric A given by its upper triangle by rows.
        type(csr_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: i, p, j

        y = 0
        do i = 1, a%nrow
            do p = a%row_ptr(i), a%row_ptr(i + 1) - 1
                j = a%col(p)
                y(i) = y(i) + a%val(p)*x(j)
                if (j /= i) y(j) = y(j) + a%val(p)*x(i)
            end do
        end do
    end subroutine upper_matvec

    pure subroutine diagonal(a, d)
        type(csr_t), intent(in) :: a
        real(dp), intent(out) :: d(:)
        integer :: i, p

        d = 0
        do i = 1, a%nrow
            do p = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (a%col(p) == i) d(i) = d(i) + a%val(p)
            end do
        end do
    end subroutine diagonal

end module fortfem_constrained_newton
