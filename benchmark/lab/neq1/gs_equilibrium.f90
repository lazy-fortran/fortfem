module fortfem_lab_gs
    ! NEQ1 consumer: fixed-boundary axisymmetric Grad-Shafranov equilibrium on
    ! the lab tensor operators with a free magnetic axis. Weak form
    !   integral (1/R) grad psi . grad v = integral g(R, psi, psia) v,
    ! psi = 0 on the boundary, psia = min psi located on the field (no axis
    ! DOFs). The Newton Jacobian is K - M[dg/dpsi] - c e_axis^T, where the
    ! rank-one term carries the free-axis normalisation. Physics lives here;
    ! the Newton-Krylov driver and tensor helpers are generic FortFEM code.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortfem_newton_krylov, only: nonlinear_system_t
    use fortfem_tensor_mesh, only: tensor_mesh_t, ogrid_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                                       apply_tensor_operator
    use fortfem_tensor_solver, only: tensor_multigrid_t, build_tensor_multigrid, &
                                     tensor_vcycle
    use fortfem_tensor_pointwise, only: tensor_quadrature_values, &
                                        tensor_weighted_mass, tensor_weighted_load
    use fortfem_tensor_points, only: tensor_locator_t, build_tensor_locator, &
                                     locate_tensor_point, evaluate_tensor_point, &
                                     tensor_point_basis, tensor_field_extremum
    use fortfem_lab_gs_source_p1f1, only: lab_gs_source_p1f1
    use fortfem_lab_gs_source_p2f3, only: lab_gs_source_p2f3
    implicit none
    private
    public :: gs_case_t, gs_system_t, read_gs_case, solovev_psibar, solovev_psi
    public :: build_gs_system, gs_initial_guess, gs_point, gs_safety_factor
    public :: set_active_case, gs_profile_variant

    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: gs_case_t
        real(dp) :: r0 = 1, f0 = 1, psi0 = 1, mu0 = 4.e-7_dp*pi, a = 0
        real(dp) :: eps = 0.32_dp, kappa = 1.7_dp, delta = 0.33_dp
        real(dp) :: c(7) = 0, psi_axis = 0, r_axis = 1, z_axis = 0
        real(dp) :: p_axis = 0, cf = 0
        integer :: alpha_p = 1, alpha_f = 1
    end type

    type, extends(nonlinear_system_t) :: gs_system_t
        type(gs_case_t) :: case
        type(tensor_operator_t) :: op
        type(tensor_multigrid_t) :: mg
        type(tensor_locator_t) :: locator
        integer :: variant = 1, axis_cell = 0
        real(dp) :: axis_xi(3) = 0, axis_x(2) = 0, psia = 0, psib = 0
        real(dp), allocatable :: psiq(:, :), g(:, :), gpsi(:, :), gpsia(:, :)
        real(dp), allocatable :: cvec(:), axis_basis(:), evaluated(:)
        integer, allocatable :: axis_dofs(:)
        real(dp) :: time_axis = 0, time_setup(4) = 0 ! mesh, operator, mg, locator
        integer :: axis_searches = 0
    contains
        procedure :: residual => gs_residual
        procedure :: linearise => gs_linearise
        procedure :: jacobian => gs_jacobian
        procedure :: precondition => gs_precondition
    end type

    type(gs_case_t), save :: active

contains

    function wall() result(t)
        real(dp) :: t
        integer(int64) :: count, rate
        call system_clock(count, rate)
        t = real(count, dp)/real(rate, dp)
    end function

    subroutine read_gs_case(path, case, status)
        character(*), intent(in) :: path
        type(gs_case_t), intent(out) :: case
        integer, intent(out) :: status
        real(dp) :: r0, f0, psi0, mu0, a, eps, kappa, delta, c(7), psi_axis, &
                    r_axis, z_axis, p_axis, cf, alpha_p, alpha_f
        integer :: unit
        namelist /solovev/ r0, f0, psi0, mu0, a, eps, kappa, delta, c, psi_axis, &
            r_axis, z_axis, p_axis, cf, alpha_p, alpha_f
        r0 = case%r0; f0 = case%f0; psi0 = case%psi0; mu0 = case%mu0; a = case%a
        eps = case%eps; kappa = case%kappa; delta = case%delta; c = case%c
        psi_axis = case%psi_axis; r_axis = case%r_axis; z_axis = case%z_axis
        p_axis = case%p_axis; cf = case%cf
        alpha_p = real(case%alpha_p, dp); alpha_f = real(case%alpha_f, dp)
        open (newunit=unit, file=path, status='old', action='read', iostat=status)
        if (status /= 0) return
        read (unit, nml=solovev, iostat=status)
        close (unit)
        if (status /= 0) return
        case = gs_case_t(r0, f0, psi0, mu0, a, eps, kappa, delta, c, psi_axis, &
                         r_axis, z_axis, p_axis, cf, nint(alpha_p), nint(alpha_f))
    end subroutine

    subroutine set_active_case(case)
        type(gs_case_t), intent(in) :: case
        active = case
    end subroutine

    integer function gs_profile_variant(case) result(variant)
        type(gs_case_t), intent(in) :: case
        variant = 0
        if (case%alpha_p == 1 .and. case%alpha_f == 1) variant = 1
        if (case%alpha_p == 2 .and. case%alpha_f == 3) variant = 2
    end function

    pure function solovev_psibar(case, x, y) result(psi)
        ! Cerfon & Freidberg 2010 (doi:10.1063/1.3328818), up-down symmetric
        ! homogeneous set plus particular solution, normalised x = R/R0.
        type(gs_case_t), intent(in) :: case
        real(dp), intent(in) :: x, y
        real(dp) :: psi, lx, h(7)
        lx = log(x)
        h(1) = 1
        h(2) = x**2
        h(3) = y**2 - x**2*lx
        h(4) = x**4 - 4*x**2*y**2
        h(5) = 2*y**4 - 9*y**2*x**2 + 3*x**4*lx - 12*x**2*y**2*lx
        h(6) = x**6 - 12*x**4*y**2 + 8*x**2*y**4
        h(7) = 8*y**6 - 140*y**4*x**2 + 75*y**2*x**4 - 15*x**6*lx + &
               180*x**4*y**2*lx - 120*x**2*y**4*lx
        psi = x**4/8 + case%a*(x**2*lx/2 - x**4/8) + dot_product(case%c, h)
    end function

    function solovev_psi(x) result(psi)
        real(dp), intent(in) :: x(:)
        real(dp) :: psi
        psi = active%psi0*solovev_psibar(active, x(1)/active%r0, x(2)/active%r0)
    end function

    subroutine solovev_boundary(angle, x)
        ! psibar = 0 along the elliptic-direction ray from the geometric centre.
        real(dp), intent(in) :: angle
        real(dp), intent(out) :: x(2)
        real(dp) :: d(2), lo, hi, mid
        integer :: it
        d = [cos(angle), active%kappa*sin(angle)]
        d = d/norm2(d)
        lo = 0
        hi = 0
        do it = 1, 400
            hi = hi + 0.005_dp
            if (solovev_psibar(active, 1 + hi*d(1), hi*d(2)) > 0) exit
            lo = hi
        end do
        do it = 1, 200
            mid = 0.5_dp*(lo + hi)
            if (solovev_psibar(active, 1 + mid*d(1), mid*d(2)) > 0) then
                hi = mid
            else
                lo = mid
            end if
            if (hi - lo <= 4*epsilon(1.0_dp)*hi) exit
        end do
        x = active%r0*[1 + 0.5_dp*(lo + hi)*d(1), 0.5_dp*(lo + hi)*d(2)]
    end subroutine

    subroutine inverse_r(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(k, 1)
            k(i, i) = 1/x(1)
        end do
        c = 0
    end subroutine

    subroutine build_gs_system(case, n, p, system, status)
        type(gs_case_t), intent(in) :: case
        integer, intent(in) :: n, p
        type(gs_system_t), intent(out) :: system
        integer, intent(out) :: status
        type(tensor_mesh_t) :: mesh
        real(dp) :: core, t0
        core = 0.32_dp
        call set_active_case(case)
        system%case = case
        system%variant = gs_profile_variant(case)
        status = 1
        if (system%variant == 0) return
        t0 = wall()
        call ogrid_tensor_mesh(2, n, p, mesh, status, boundary=solovev_boundary, &
                               center=[case%r0, 0.0_dp], &
                               core_half=core*case%eps*case%r0*[1.0_dp, case%kappa])
        if (status /= 0) return
        system%time_setup(1) = wall() - t0
        t0 = wall()
        call build_tensor_operator(mesh, inverse_r, system%op, status)
        if (status /= 0) return
        system%time_setup(2) = wall() - t0
        t0 = wall()
        call build_tensor_multigrid(mesh, inverse_r, system%mg, status)
        if (status /= 0) return
        system%time_setup(3) = wall() - t0
        t0 = wall()
        call build_tensor_locator(mesh, system%locator, status)
        if (status /= 0) return
        system%time_setup(4) = wall() - t0
        associate (mq => product(system%op%no), ne => mesh%nelem)
            allocate (system%psiq(mq, ne), system%g(mq, ne), system%gpsi(mq, ne), &
                      system%gpsia(mq, ne))
        end associate
        allocate (system%cvec(mesh%ndof), system%axis_basis(size(mesh%dof, 1)), &
                  system%axis_dofs(size(mesh%dof, 1)))
    end subroutine

    subroutine evaluate_source(this, u, status)
        ! Quadrature-point flux, free axis and generated source kernel at u.
        ! The driver linearises at the last accepted residual point, so the
        ! state is reused when u is unchanged.
        class(gs_system_t), intent(inout) :: this
        real(dp), intent(in) :: u(:)
        integer, intent(out) :: status
        real(dp) :: t0, x(2)
        integer :: e, q
        status = 0
        if (allocated(this%evaluated)) then
            if (all(this%evaluated == u)) return
        end if
        call tensor_quadrature_values(this%op, u, this%psiq)
        t0 = wall()
        call tensor_field_extremum(this%op%mesh, this%locator, u, this%axis_cell, &
                                   this%axis_xi, x, this%psia, status)
        this%time_axis = this%time_axis + wall() - t0
        this%axis_searches = this%axis_searches + 1
        if (status /= 0) return
        this%axis_x = x
        if (this%psia >= this%psib) then
            status = 6
            return
        end if
        !$omp parallel do schedule(static) private(q)
        do e = 1, this%op%mesh%nelem
            do q = 1, size(this%psiq, 1)
                if (this%variant == 1) then
                    call lab_gs_source_p1f1(this%op%points(1, q, e), this%psiq(q, e), &
                                            this%psia, this%psib, this%case%mu0, &
                                            this%case%p_axis, this%case%cf, &
                                            this%g(q, e), this%gpsi(q, e), &
                                            this%gpsia(q, e))
                else
                    call lab_gs_source_p2f3(this%op%points(1, q, e), this%psiq(q, e), &
                                            this%psia, this%psib, this%case%mu0, &
                                            this%case%p_axis, this%case%cf, &
                                            this%g(q, e), this%gpsi(q, e), &
                                            this%gpsia(q, e))
                end if
            end do
        end do
        !$omp end parallel do
        this%evaluated = u
    end subroutine

    subroutine gs_residual(this, u, r, status)
        class(gs_system_t), intent(inout) :: this
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: r(:)
        integer, intent(out) :: status
        real(dp), allocatable :: load(:)
        call evaluate_source(this, u, status)
        if (status /= 0) return
        allocate (load(size(u)))
        call tensor_weighted_load(this%op, this%g, load)
        call apply_tensor_operator(this%op, u, r)
        ! Boundary rows of r are u - psib from the eliminated operator.
        r = r - load
        where (this%op%mesh%boundary) r = r - this%psib
    end subroutine

    subroutine gs_linearise(this, u, status)
        class(gs_system_t), intent(inout) :: this
        real(dp), intent(in) :: u(:)
        integer, intent(out) :: status
        call evaluate_source(this, u, status)
        if (status /= 0) return
        call tensor_weighted_load(this%op, this%gpsia, this%cvec)
        call tensor_point_basis(this%op%mesh, this%locator, this%axis_xi, &
                                this%axis_basis)
        this%axis_dofs = this%op%mesh%dof(:, this%axis_cell)
    end subroutine

    subroutine gs_jacobian(this, v, w, status)
        class(gs_system_t), intent(inout) :: this
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: w(:)
        integer, intent(out) :: status
        real(dp), allocatable :: mv(:)
        real(dp) :: va
        allocate (mv(size(v)))
        call apply_tensor_operator(this%op, v, w)
        if (this%variant /= 1) then
            call tensor_weighted_mass(this%op, this%gpsi, v, mv)
            w = w - mv
        end if
        va = dot_product(this%axis_basis, v(this%axis_dofs))
        w = w - va*this%cvec
        status = 0
    end subroutine

    subroutine gs_precondition(this, v, w, status)
        class(gs_system_t), intent(inout) :: this
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: w(:)
        integer, intent(out) :: status
        call tensor_vcycle(this%mg, 1, v, w, status)
    end subroutine

    subroutine gs_initial_guess(this, u, tol, iterations, status)
        ! Flat-profile start: with Delta = psib - psia = 1 the source at
        ! psin = 0 is g0(R); U = K^-1 L[g0]. Since g scales as 1/Delta for the
        ! normalised-profile family, u0 = U/Delta with Delta^2 = -min U.
        class(gs_system_t), intent(inout) :: this
        real(dp), intent(out) :: u(:)
        real(dp), intent(in) :: tol
        integer, intent(out) :: iterations, status
        real(dp), allocatable :: rhs(:), r(:), z(:), p(:), ap(:)
        real(dp) :: rho, previous, alpha, initial, delta
        integer :: e, q, it
        if (allocated(this%evaluated)) deallocate (this%evaluated)
        !$omp parallel do schedule(static) private(q)
        do e = 1, this%op%mesh%nelem
            do q = 1, size(this%g, 1)
                if (this%variant == 1) then
                    call lab_gs_source_p1f1(this%op%points(1, q, e), -1.0_dp, &
                                            -1.0_dp, 0.0_dp, this%case%mu0, &
                                            this%case%p_axis, this%case%cf, &
                                            this%g(q, e), this%gpsi(q, e), &
                                            this%gpsia(q, e))
                else
                    call lab_gs_source_p2f3(this%op%points(1, q, e), -1.0_dp, &
                                            -1.0_dp, 0.0_dp, this%case%mu0, &
                                            this%case%p_axis, this%case%cf, &
                                            this%g(q, e), this%gpsi(q, e), &
                                            this%gpsia(q, e))
                end if
            end do
        end do
        !$omp end parallel do
        allocate (rhs(size(u)), r(size(u)), z(size(u)), p(size(u)), ap(size(u)))
        call tensor_weighted_load(this%op, this%g, rhs)
        u = 0
        r = rhs
        initial = norm2(rhs)
        previous = 0
        iterations = 0
        status = 1
        do it = 1, 200
            call tensor_vcycle(this%mg, 1, r, z, status)
            if (status /= 0) return
            rho = dot_product(r, z)
            if (it == 1) then
                p = z
            else
                p = z + (rho/previous)*p
            end if
            call apply_tensor_operator(this%op, p, ap)
            alpha = rho/dot_product(p, ap)
            u = u + alpha*p
            r = r - alpha*ap
            previous = rho
            iterations = it
            if (norm2(r) <= tol*initial) exit
        end do
        status = 1
        if (minval(u) >= 0) return
        delta = sqrt(-minval(u))
        u = u/delta
        status = 0
    end subroutine

    subroutine gs_point(this, u, x, cell, psi, gradient, status)
        class(gs_system_t), intent(in) :: this
        real(dp), intent(in) :: u(:), x(2)
        integer, intent(inout) :: cell
        real(dp), intent(out) :: psi, gradient(2)
        integer, intent(out) :: status
        real(dp) :: xi(2)
        call locate_tensor_point(this%op%mesh, this%locator, x, cell, xi, status)
        if (status /= 0) return
        call evaluate_tensor_point(this%op%mesh, this%locator, u, cell, xi, psi, &
                                   gradient, status)
    end subroutine

    subroutine gs_safety_factor(this, u, level, rays, q, status, start)
        ! q = F/(2 pi) contour integral dl/(R |grad psi|) on the level set
        ! psi = level, parametrised by the ray angle about the located axis:
        ! dl/|grad psi| = rho dtheta / |d psi/d rho| (trapezoid, periodic).
        ! Each ray root is bracketed (points outside the mesh count as above
        ! the level) and refined by Newton steps safeguarded with bisection.
        class(gs_system_t), intent(in) :: this
        real(dp), intent(in) :: u(:), level
        integer, intent(in) :: rays
        real(dp), intent(out) :: q
        integer, intent(out) :: status
        real(dp), intent(in), optional :: start
        real(dp) :: rho, theta, e(2), x(2), psi, gradient(2), slope, total, &
                    psin, fsq, lo, hi, next
        integer :: k, it, cell
        rho = 0.5_dp*this%case%eps*this%case%r0
        if (present(start)) rho = start
        cell = this%axis_cell
        total = 0
        status = 1
        do k = 0, rays - 1
            theta = 2*pi*real(k, dp)/real(rays, dp)
            e = [cos(theta), sin(theta)]
            lo = 0
            hi = -1
            do it = 1, 80
                x = this%axis_x + rho*e
                call gs_point(this, u, x, cell, psi, gradient, status)
                if (status /= 0 .or. psi >= level) then
                    hi = rho
                    exit
                end if
                lo = rho
                rho = 1.25_dp*rho
            end do
            if (hi < 0) then
                status = 1
                return
            end if
            rho = 0.5_dp*(lo + hi)
            do it = 1, 200
                x = this%axis_x + rho*e
                call gs_point(this, u, x, cell, psi, gradient, status)
                if (status /= 0) then
                    hi = rho
                    rho = 0.5_dp*(lo + hi)
                    cycle
                end if
                if (psi < level) then
                    lo = rho
                else
                    hi = rho
                end if
                slope = dot_product(gradient, e)
                next = -1
                if (slope > 0) next = rho + (level - psi)/slope
                if (next <= lo .or. next >= hi) next = 0.5_dp*(lo + hi)
                if (abs(next - rho) <= 1.e-14_dp*rho .or. &
                    hi - lo <= 1.e-14_dp*hi) exit
                rho = next
            end do
            x = this%axis_x + rho*e
            call gs_point(this, u, x, cell, psi, gradient, status)
            if (status /= 0) return
            slope = dot_product(gradient, e)
            if (slope <= 0) then
                status = 1
                return
            end if
            total = total + rho/(x(1)*slope)
        end do
        psin = (level - this%psia)/(this%psib - this%psia)
        fsq = this%case%f0**2 + this%case%cf*(1 - psin)**this%case%alpha_f
        q = sqrt(fsq)/(2*pi)*total*(2*pi/real(rays, dp))
        status = 0
    end subroutine
end module fortfem_lab_gs
