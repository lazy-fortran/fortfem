module uns1_island_problem
    ! Single-helicity island chain of M. Hoelzl, S. Guenter, Q. Yu, K. Lackner,
    ! Phys. Plasmas 14, 052501 (2007), doi:10.1063/1.2722980, Sec. II.A:
    ! periodic cylinder, a = 1, B_phi = 1, q(r) = 0.9 sqrt(1 + (r/0.48)^4),
    ! Psi_pert = Psi_0 r^2 (1 - r^2)^2 cos(m theta - n phi), m/n = 4/3.
    ! R0 = 3.3 is inferred: the paper's w/w_c = 1.1, 1.9, 3.4, 6.0 at
    ! chi ratio 1e6..1e9 with w_c = r_s ratio^(-1/4) sqrt(8/(eps_s s_s n)) give
    ! R0/a = 3.295. Psi_0 is fixed by the paper's island width w = 0.068 a.
    ! Single helicity reduces the problem exactly to the helical-frame plane
    ! (X, Y) with u = theta - (n/m) phi: grad T = M (T_X, T_Y), M^T M =
    ! I + alpha^2 [[Y^2, -XY], [-XY, X^2]], alpha = n/(m R0), and the field
    ! direction enters as c = M^T b = (z x grad psi*)/|B| with the helical flux
    ! psi* = Psi - alpha r^2/2. Then K2 = chi_perp M^T M + (chi_par - chi_perp) c c^T.
    ! Manufactured T (constant on field lines plus an O(1/ratio) perturbation):
    ! a C-infinity core-to-wall profile that is flat at the separatrix, axis and
    ! wall, plus a bump peaked at the four O-points; T = 0 on the wall.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: problem_setup, coefficient, perp_coefficient, exact, flux, psistar
    public :: grad_exact, grad_psistar, opoint, xpoint, width, chain_level, direction
    real(dp), parameter, public :: c48 = 0.48_dp, r0 = 3.3_dp, wtarget = 0.068_dp
    integer, parameter, public :: mpol = 4, ntor = 3
    real(dp), parameter :: alpha = real(ntor, dp)/(mpol*r0), pi = acos(-1.0_dp)
    real(dp), parameter :: tsep = 0.4_dp, tcore = 0.6_dp, bump = 0.3_dp, g1 = 0.7_dp
    real(dp), public :: amplitude = 0, ratio = 1, psi_sep = 0, psi_o = 0, &
                        psi_wall = 0, psi_c = 0, r_x = 0, r_o = 0
contains
    subroutine problem_setup(chi_ratio)
        real(dp), intent(in) :: chi_ratio
        real(dp) :: lo, hi
        integer :: it
        ratio = chi_ratio
        lo = 1.e-7_dp; hi = 1.e-2_dp
        do it = 1, 200
            amplitude = (lo + hi)/2
            if (width() < wtarget) then
                lo = amplitude
            else
                hi = amplitude
            end if
            if (hi - lo < 1.e-17_dp) exit
        end do
        amplitude = (lo + hi)/2
        r_x = ridge(pi/mpol)
        r_o = ridge(0.0_dp)
        psi_sep = psistar(r_x*[cos(pi/mpol), sin(pi/mpol)])
        psi_o = psistar([r_o, 0.0_dp])
        psi_wall = psistar([1.0_dp, 0.0_dp])
        psi_c = psi_sep - 3*(psi_o - psi_sep)
    end subroutine

    real(dp) function ridge(angle)
        real(dp), intent(in) :: angle
        real(dp), parameter :: g = (sqrt(5.0_dp) - 1)/2
        real(dp) :: lo, hi, x1, x2, f1, f2
        integer :: it
        lo = 0.3_dp; hi = 0.7_dp
        x1 = hi - g*(hi - lo); x2 = lo + g*(hi - lo)
        f1 = psistar(x1*[cos(angle), sin(angle)])
        f2 = psistar(x2*[cos(angle), sin(angle)])
        do it = 1, 200
            if (hi - lo < 1.e-15_dp) exit
            if (f1 < f2) then
                lo = x1; x1 = x2; f1 = f2; x2 = lo + g*(hi - lo)
                f2 = psistar(x2*[cos(angle), sin(angle)])
            else
                hi = x2; x2 = x1; f2 = f1; x1 = hi - g*(hi - lo)
                f1 = psistar(x1*[cos(angle), sin(angle)])
            end if
        end do
        ridge = (lo + hi)/2
    end function

    real(dp) function width()
        ! Largest radial extent of the separatrix: along the O-point ray.
        real(dp) :: rx, ro, sep, lo, hi, inner, outer
        integer :: it
        rx = ridge(pi/mpol)
        ro = ridge(0.0_dp)
        sep = psistar(rx*[cos(pi/mpol), sin(pi/mpol)])
        lo = 0.2_dp; hi = ro
        do it = 1, 200
            inner = (lo + hi)/2
            if (psistar([inner, 0.0_dp]) < sep) then
                lo = inner
            else
                hi = inner
            end if
        end do
        lo = ro; hi = 0.95_dp
        do it = 1, 200
            outer = (lo + hi)/2
            if (psistar([outer, 0.0_dp]) < sep) then
                hi = outer
            else
                lo = outer
            end if
        end do
        width = outer - inner
    end function

    subroutine opoint(k, x)
        integer, intent(in) :: k
        real(dp), intent(out) :: x(2)
        x = r_o*[cos((k - 1)*2*pi/mpol), sin((k - 1)*2*pi/mpol)]
    end subroutine

    subroutine xpoint(k, x)
        integer, intent(in) :: k
        real(dp), intent(out) :: x(2)
        x = r_x*[cos((k - 0.5_dp)*2*pi/mpol), sin((k - 0.5_dp)*2*pi/mpol)]
    end subroutine

    subroutine flux_function(x, psi, grad)
        ! Psi and grad Psi (poloidal flux), r^2 cos(4u) = (X^4-6X^2Y^2+Y^4)/r^2.
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: psi, grad(2)
        real(dp) :: r2, q, poly, gpoly(2), h, gh(2), env
        r2 = x(1)**2 + x(2)**2
        q = 0.9_dp*sqrt(1 + (r2/c48**2)**2)
        psi = c48**2/(2*0.9_dp*r0)*asinh(r2/c48**2)
        grad = x/(r0*q)
        poly = x(1)**4 - 6*x(1)**2*x(2)**2 + x(2)**4
        gpoly = [4*x(1)**3 - 12*x(1)*x(2)**2, -12*x(1)**2*x(2) + 4*x(2)**3]
        if (r2 > 0) then
            h = poly/r2
            gh = gpoly/r2 - 2*poly*x/r2**2
        else
            h = 0
            gh = 0
        end if
        env = (1 - r2)**2
        psi = psi + amplitude*env*h
        grad = grad + amplitude*(env*gh - 4*(1 - r2)*h*x)
    end subroutine

    function psistar(x) result(value)
        real(dp), intent(in) :: x(2)
        real(dp) :: value, g(2)
        call flux_function(x, value, g)
        value = value - alpha*(x(1)**2 + x(2)**2)/2
    end function

    function chain_level(x) result(value)
        real(dp), intent(in) :: x(2)
        real(dp) :: value
        value = psistar(x)
    end function

    subroutine grad_psistar(x, g, babs)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: g(2), babs
        real(dp) :: psi
        call flux_function(x, psi, g)
        babs = sqrt(1 + g(1)**2 + g(2)**2)
        g = g - alpha*x
    end subroutine

    subroutine metric(x, m)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: m(2, 2)
        m(1, 1) = 1 + alpha**2*x(2)**2
        m(2, 2) = 1 + alpha**2*x(1)**2
        m(1, 2) = -alpha**2*x(1)*x(2)
        m(2, 1) = m(1, 2)
    end subroutine

    subroutine direction(x, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: c(2)
        real(dp) :: g(2), babs
        call grad_psistar(x(1:2), g, babs)
        c = [-g(2), g(1)]/babs
    end subroutine

    subroutine coefficient(x, k, reaction)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), reaction
        real(dp) :: m(2, 2), c(2)
        integer :: i, j
        call metric(x(1:2), m)
        call direction(x, c)
        do j = 1, 2
            do i = 1, 2
                k(i, j) = m(i, j) + (ratio - 1)*c(i)*c(j)
            end do
        end do
        reaction = 0
    end subroutine

    subroutine perp_coefficient(x, k, reaction)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), reaction
        real(dp) :: m(2, 2)
        call metric(x(1:2), m)
        k(1:2, 1:2) = m
        reaction = 0
    end subroutine

    pure real(dp) function phi(t)
        real(dp), intent(in) :: t
        phi = 0
        if (t > 0) phi = exp(-1/t)
    end function

    pure real(dp) function dphi(t)
        real(dp), intent(in) :: t
        dphi = 0
        if (t > 0) dphi = exp(-1/t)/t**2
    end function

    pure subroutine step(t, value, slope)
        ! C-infinity step: 0 for t <= 0, 1 for t >= 1.
        real(dp), intent(in) :: t
        real(dp), intent(out) :: value, slope
        real(dp) :: a, b
        a = phi(t); b = phi(1 - t)
        value = a/(a + b)
        slope = (dphi(t)*b + a*dphi(1 - t))/(a + b)**2
    end subroutine

    subroutine profile(x, value, dvalue)
        ! T without the perturbation, and dT/dpsi*.
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: value, dvalue
        real(dp) :: psi, s, v, dv, gb
        psi = psistar(x)
        value = tsep
        dvalue = 0
        if (norm2(x) < r_x) then
            call step((psi_sep - psi)/psi_sep/g1, v, dv)
            value = value + tcore*v
            dvalue = dvalue - tcore*dv/(psi_sep*g1)
        else
            call step((psi_sep - psi)/(psi_sep - psi_wall), v, dv)
            value = value - tsep*v
            dvalue = dvalue + tsep*dv/(psi_sep - psi_wall)
        end if
        s = (psi - psi_c)/(psi_o - psi_c)
        if (s > 0) then
            gb = exp(1 - 1/s)
            value = value + bump*gb
            dvalue = dvalue + bump*gb/(s**2*(psi_o - psi_c))
        end if
    end subroutine

    function exact(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value, d, r2
        call profile(x(1:2), value, d)
        r2 = x(1)**2 + x(2)**2
        value = value + (1 - r2)**2*x(1)*x(2)/ratio
    end function

    subroutine perturbation(x, g)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: g(2)
        real(dp) :: r2
        r2 = x(1)**2 + x(2)**2
        g = (1 - r2)**2*[x(2), x(1)] - 4*(1 - r2)*x(1)*x(2)*x
    end subroutine

    subroutine grad_exact(x, g)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: g(2)
        real(dp) :: value, d, gp(2), babs, gg(2)
        call profile(x, value, d)
        call grad_psistar(x, gp, babs)
        call perturbation(x, gg)
        g = d*gp + gg/ratio
    end subroutine

    subroutine flux(x, f)
        ! F = K2 grad T with c . grad(profile) = 0 cancelled analytically.
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:)
        real(dp) :: m(2, 2), c(2), g(2), gg(2)
        call metric(x(1:2), m)
        call direction(x, c)
        call grad_exact(x(1:2), g)
        call perturbation(x(1:2), gg)
        f(1:2) = matmul(m, g) + (1 - 1/ratio)*c*dot_product(c, gg)
    end subroutine
end module

module uns1_amg_coarse
    ! FortSparse SA-AMG (lazy-fortran/fortsparse#29) as the p-multigrid coarse
    ! solve: one fixed V-cycle, set up on the first call for a given matrix.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortsparse_amg, only: amg_params_t, amg_handle_t, amg_setup_csr, amg_apply, &
                              amg_destroy
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK
    implicit none
    private
    public :: amg_coarse, amg_reset, amg_cycles
    type(amg_handle_t), save :: handle
    logical, save :: ready = .false.
    integer, save :: amg_cycles = 1
contains
    subroutine amg_reset()
        if (ready) call amg_destroy(handle)
        ready = .false.
    end subroutine

    subroutine amg_coarse(matrix, rhs, x, status)
        type(csr_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: status
        type(amg_params_t) :: params
        type(fortsparse_status_t) :: info
        real(dp), allocatable :: r(:), z(:)
        integer :: i, j, k
        status = 1
        if (.not. ready) then
            call amg_setup_csr(matrix%n, matrix%rowptr, matrix%colidx, matrix%val, &
                               params, handle, info)
            if (info%code /= FORTSPARSE_OK) return
            ready = .true.
        end if
        x = 0
        allocate (r(size(rhs)), z(size(rhs)))
        r = rhs
        do k = 1, amg_cycles
            call amg_apply(handle, r, z, info)
            if (info%code /= FORTSPARSE_OK) return
            x = x + z
            if (k == amg_cycles) exit
            r = rhs
            do i = 1, matrix%n
                do j = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                    r(i) = r(i) - matrix%val(j)*x(matrix%colidx(j))
                end do
            end do
        end do
        status = 0
    end subroutine
end module

program island_diffusion
    ! Usage: island_diffusion KIND LEVEL P RATIO [COARSE [TOL]]
    !   KIND ogrid|island, COARSE ldlt|amg|direct. Writes one CSV row to stdout.
    !   direct factors the assembled fine operator (the system PCG solves) by
    !   sparse LDLT, so the error columns are free of iteration error.
    !   island_diffusion check   verifies the analytic gradients by central FD.
    use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    use uns1_island_problem
    use uns1_amg_coarse, only: amg_coarse, amg_reset
    use fortfem_tensor_mesh, only: tensor_mesh_t, ogrid_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                                       apply_tensor_operator, tensor_load_flux, &
                                       tensor_l2_error
    use fortfem_tensor_contractions, only: tensor_contract
    use fortfem_tensor_solver, only: tensor_solve_info_t, solve_tensor_diffusion
    use fortfem_tensor_operator, only: assemble_tensor_operator
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
                                     sparse_direct_factor_csc, &
                                     sparse_direct_solve_factored, sparse_direct_free
    use fortsparse, only: sparse_options_t, FORTSPARSE_SYM_INDEFINITE
    use fortfem_tensor_mesh_tools, only: tensor_mesh_quality_t, tensor_mesh_check
    use fortfem_island_quad_mesh, only: island_chain_t, island_mesh_counts_t, &
                                        island_chain_quad_mesh
    implicit none
    character(len=32) :: kind, coarse, word
    type(tensor_mesh_t) :: mesh
    type(tensor_operator_t) :: op, perp
    type(tensor_solve_info_t) :: info
    type(tensor_mesh_quality_t) :: quality
    type(island_chain_t) :: chain
    type(island_mesh_counts_t) :: counts
    real(dp), allocatable :: rhs(:), x(:), y(:)
    real(dp) :: chi, tol, t0, t1, t2, rel, norm_t, inner, amp, etot, eperp, nodal, &
                interp
    integer :: level, p, status, k, i, extra
    call get_command_argument(1, kind)
    if (trim(kind) == 'check') then
        call check_gradients()
        stop
    end if
    if (command_argument_count() < 4) error stop 'KIND LEVEL P RATIO [COARSE [TOL]]'
    call get_command_argument(2, word); read (word, *) level
    call get_command_argument(3, word); read (word, *) p
    call get_command_argument(4, word); read (word, *) chi
    coarse = 'ldlt'
    if (command_argument_count() >= 5) call get_command_argument(5, coarse)
    tol = 1.e-10_dp
    if (command_argument_count() >= 6) then
        call get_command_argument(6, word); read (word, *) tol
    end if
    call problem_setup(chi)
    t0 = omp_get_wtime()
    select case (trim(kind))
    case ('ogrid')
        ! MF1 five-block O-grid of the unit disk: unaligned, through the axis.
        call ogrid_tensor_mesh(2, level, p, mesh, status, q=p)
    case ('island')
        allocate (chain%opoint(2, mpol), chain%xpoint(2, mpol))
        do k = 1, mpol
            call opoint(k, chain%opoint(:, k))
            call xpoint(k, chain%xpoint(:, k))
        end do
        chain%r_core = 0.15_dp
        chain%r_inner = 0.3_dp
        chain%r_outer = 0.75_dp
        chain%r_wall = 1
        counts%half = level
        counts%core_rings = level
        counts%inner_rows = level + 1
        counts%island_rings = level
        counts%outer_rows = level + 1
        counts%wall_rows = level
        call island_chain_quad_mesh(chain_level, chain, counts, p, p, mesh, status)
    case default
        error stop 'KIND must be ogrid or island'
    end select
    if (status /= 0) error stop 'mesh'
    call get_environment_variable('UNS1_MESH_OUTLINE', word, status=status)
    if (status == 0) call dump_outline(mesh, trim(word))
    call tensor_mesh_check(mesh, quality)
    if (quality%inverted_cells /= 0) error stop 'inverted cells'
    ! Load and diagnostics on a rule with two extra points per direction
    ! (UNS1_LOAD_EXTRA): the collocated p+1 rule under-integrates the steep
    ! manufactured flux on distorted cells; the stiffness keeps p+1 points.
    extra = 2
    call get_environment_variable('UNS1_LOAD_EXTRA', word, status=status)
    if (status == 0) read (word, *) extra
    call build_tensor_operator(mesh, coefficient, op, status, nq=p + 1 + extra)
    if (status /= 0) error stop 'operator'
    call tensor_load_flux(op, flux, rhs)
    allocate (x(mesh%ndof), y(mesh%ndof))
    x = 0
    t1 = omp_get_wtime()
    if (trim(coarse) == 'direct') then
        call direct_solve(mesh, rhs, x, info)
    else if (trim(coarse) == 'amg') then
        call solve_tensor_diffusion(mesh, coefficient, rhs, tol, x, info, &
                                    maxiter=500, coarse_hook=amg_coarse)
        call amg_reset()
    else
        call solve_tensor_diffusion(mesh, coefficient, rhs, tol, x, info, maxiter=500)
    end if
    t2 = omp_get_wtime()
    norm_t = tensor_l2_error(op, 0*x, exact)
    rel = tensor_l2_error(op, x, exact)/norm_t
    call weighted_inner(op, x, inner)
    amp = inner/norm_t**2
    nodal = 0
    do i = 1, mesh%ndof
        y(i) = exact(mesh%coordinates(:, i))
        nodal = max(nodal, abs(x(i) - y(i)))
    end do
    ! Nodal interpolant: a best-approximation proxy independent of the solve.
    interp = tensor_l2_error(op, y, exact)/norm_t
    call get_environment_variable('UNS1_CELL_ERRORS', word, status=status)
    if (status == 0) call dump_cells(op, x, trim(word))
    ! Spurious parallel power: total minus perpendicular energy of T_h.
    call build_tensor_operator(mesh, perp_coefficient, perp, status)
    call apply_tensor_operator(op, x, y, eliminate=.false.)
    etot = dot_product(x, y)
    call apply_tensor_operator(perp, x, y, eliminate=.false.)
    eperp = dot_product(x, y)
    write (output_unit, '(a,",",a,",",i0,",",i0,",",i0,",",es9.2,",",i0,",",i0,'// &
           '",",es12.5,",",es12.5,",",es12.5,",",es12.5,",",i0,",",i0,",",es10.3,'// &
           '",",f10.4,",",f10.4,",",i0,",",i0,",",es10.3,",",f8.5,",",es12.5)') &
        trim(kind), trim(coarse), level, p, p, chi, mesh%ndof, mesh%nelem, rel, &
        nodal, 1/amp - 1, (etot - eperp)/eperp, info%iterations, info%status, &
        info%residual, t1 - t0, t2 - t1, omp_get_max_threads(), peak_rss_kb(), &
        quality%min_scaled_jacobian, width(), interp
contains
    subroutine direct_solve(mesh, rhs, x, info)
        ! Same discrete operator and quadrature as solve_tensor_diffusion.
        type(tensor_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(inout) :: x(:)
        type(tensor_solve_info_t), intent(out) :: info
        type(tensor_operator_t) :: stiff
        type(csr_matrix_t) :: a
        type(sparse_direct_factor_t) :: factor
        type(sparse_options_t) :: options
        integer, allocatable :: ptr(:), rows(:), fill(:)
        real(dp), allocatable :: values(:), r(:)
        integer :: n, i, j, k, pos, status
        info%status = 1
        call build_tensor_operator(mesh, coefficient, stiff, status)
        if (status /= 0) return
        call assemble_tensor_operator(stiff, a)
        n = a%n
        allocate (ptr(n + 1), fill(n), rows(size(a%val)), values(size(a%val)), r(n))
        ptr = 0
        do k = 1, size(a%colidx)
            ptr(a%colidx(k) + 1) = ptr(a%colidx(k) + 1) + 1
        end do
        ptr(1) = 1
        do i = 1, n
            ptr(i + 1) = ptr(i + 1) + ptr(i)
        end do
        fill = ptr(1:n)
        do i = 1, n
            do k = a%rowptr(i), a%rowptr(i + 1) - 1
                j = a%colidx(k)
                pos = fill(j)
                rows(pos) = i
                values(pos) = a%val(k)
                fill(j) = pos + 1
            end do
        end do
        options%symmetry = FORTSPARSE_SYM_INDEFINITE
        call sparse_direct_factor_csc(factor, n, ptr, rows, values, status, options)
        if (status /= 0) return
        call sparse_direct_solve_factored(factor, rhs, x, status)
        call sparse_direct_free(factor)
        if (status /= 0) return
        call apply_tensor_operator(stiff, x, r)
        info%residual = norm2(r - rhs)/max(norm2(rhs), tiny(1.0_dp))
        info%iterations = 0
        info%status = 0
    end subroutine

    subroutine weighted_inner(op, x, value)
        ! integral of T_h * T on the quadrature of the operator.
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value
        real(dp), allocatable :: v(:), local(:), t1(:), t2(:), b(:, :)
        real(dp) :: identity(1, 1)
        integer :: e, q, ws
        ws = maxval([op%ni, op%no])**3
        allocate (v(product(op%no)), local(product(op%ni)), t1(ws), t2(ws))
        b = op%line%b(:, :, 1)
        identity = 1
        value = 0
        do e = 1, op%mesh%nelem
            local = x(op%mesh%dof(:, e))
            call tensor_contract(op%ni, op%no, b, b, identity, local, v, t1, t2)
            do q = 1, product(op%no)
                value = value + op%volume(q, e)*v(q)*exact(op%points(:, q, e))
            end do
        end do
    end subroutine

    subroutine dump_cells(op, x, path)
        ! Per cell: centroid, area and squared L2 error, for error maps.
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: x(:)
        character(len=*), intent(in) :: path
        real(dp), allocatable :: v(:), local(:), t1(:), t2(:), b(:, :)
        real(dp) :: identity(1, 1), err, area, centre(2)
        integer :: e, q, ws, u
        ws = maxval([op%ni, op%no])**3
        allocate (v(product(op%no)), local(product(op%ni)), t1(ws), t2(ws))
        b = op%line%b(:, :, 1)
        identity = 1
        open (newunit=u, file=path, action='write')
        do e = 1, op%mesh%nelem
            local = x(op%mesh%dof(:, e))
            call tensor_contract(op%ni, op%no, b, b, identity, local, v, t1, t2)
            err = 0; area = 0; centre = 0
            do q = 1, product(op%no)
                err = err + op%volume(q, e)*(v(q) - exact(op%points(:, q, e)))**2
                area = area + op%volume(q, e)
                centre = centre + op%volume(q, e)*op%points(:, q, e)
            end do
            write (u, '(i0,4es14.6)') e, centre/area, area, err
        end do
        close (u)
    end subroutine

    subroutine dump_outline(mesh, path)
        ! Cell boundaries as closed polylines of geometry nodes, one cell a line.
        type(tensor_mesh_t), intent(in) :: mesh
        character(len=*), intent(in) :: path
        integer :: e, u, q, i, a
        integer, allocatable :: loop(:)
        q = mesh%geometry_degree
        allocate (loop(4*q + 1))
        do i = 0, q - 1
            loop(1 + i) = 1 + i
            loop(1 + q + i) = (q + 1)*i + q + 1
            loop(1 + 2*q + i) = (q + 1)*q + q + 1 - i
            loop(1 + 3*q + i) = (q + 1)*(q - i) + 1
        end do
        loop(4*q + 1) = 1
        open (newunit=u, file=path, action='write')
        do e = 1, mesh%nelem
            write (u, '(*(es13.5))') (mesh%geometry(:, loop(a), e), a=1, size(loop))
        end do
        close (u)
    end subroutine

    integer function peak_rss_kb()
        character(len=256) :: line
        integer :: u, ios
        peak_rss_kb = -1
        open (newunit=u, file='/proc/self/status', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:6) == 'VmHWM:') read (line(7:), *) peak_rss_kb
        end do
        close (u)
    end function

    subroutine check_gradients()
        ! Central differences of psi* and T against the analytic gradients,
        ! plus the field-line constancy c . grad(profile) = 0.
        real(dp) :: xp(2), g(2), fd(2), babs, h, worst, c(2), e1(2), e2(2), &
                    worstc, f(2), kmat(2, 2), react, fref(2), worstf
        integer :: s
        call problem_setup(1.e6_dp)
        print '(a,es14.7,a,es12.5,a,f9.6,a,f9.6)', 'Psi_0=', amplitude, &
            ' psi_sep=', psi_sep, ' r_X=', r_x, ' w=', width()
        h = 1.e-6_dp
        worst = 0; worstc = 0; worstf = 0
        do s = 1, 4000
            xp = 0.97_dp*[cos(2.399963_dp*s), sin(2.399963_dp*s)]*sqrt(s/4000.0_dp)
            if (norm2(xp) < 1.e-3_dp) cycle
            call grad_psistar(xp, g, babs)
            e1 = [h, 0.0_dp]; e2 = [0.0_dp, h]
            fd = [psistar(xp + e1) - psistar(xp - e1), &
                  psistar(xp + e2) - psistar(xp - e2)]/(2*h)
            worst = max(worst, norm2(fd - g)/max(norm2(g), 1.e-3_dp))
            call grad_exact(xp, g)
            fd = [exact(xp + e1) - exact(xp - e1), &
                  exact(xp + e2) - exact(xp - e2)]/(2*h)
            worst = max(worst, norm2(fd - g)/max(norm2(g), 1.e-1_dp))
            call direction(xp, c)
            ratio = huge(1.0_dp)
            call grad_exact(xp, g)
            worstc = max(worstc, abs(dot_product(c, g))/max(norm2(g), 1.e-12_dp))
            ratio = 1.e6_dp
            ! The cancelled flux equals K2 grad T evaluated directly.
            call flux(xp, f)
            call coefficient(xp, kmat, react)
            call grad_exact(xp, g)
            fref = matmul(kmat, g)
            worstf = max(worstf, norm2(f - fref)/max(norm2(fref), 1.e-3_dp))
        end do
        print '(a,es9.2,a,es9.2,a,es9.2)', 'max FD gradient error ', worst, &
            '  max |c.grad profile|/|grad| ', worstc, '  flux identity ', worstf
        if (worst > 1.e-6_dp .or. worstc > 1.e-12_dp .or. worstf > 1.e-6_dp) &
            error stop 'gradient check'
        if (abs(width() - wtarget) > 1.e-12_dp) error stop 'island width'
    end subroutine
end program
