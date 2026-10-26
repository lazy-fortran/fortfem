module fourier_equilibration_fixture
    use fortfem_kinds, only: dp
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t, &
        meridian_fourier_source_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), chi = 0.75_dp
    integer, parameter :: modes = 2, nb = 5, nphi = 8
    real(dp), parameter :: amplitude(nb) = [1.0_dp, 0.3_dp, 0.0_dp, 0.0_dp, 0.2_dp]
    type, extends(meridian_fourier_tensor_t) :: conductivity_t
    contains
        procedure :: tensor
    end type
    type, extends(meridian_fourier_source_t) :: source_t
    contains
        procedure :: source
    end type
contains
    subroutine tensor(self, t, x, k, status)
        class(conductivity_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: status
        integer :: j, c
        k = 0
        do j = 1, size(x, 2)
            do c = 1, 3
                k(c, c, j) = 1 + 0.2_dp*cos(x(3, j)) + 0.1_dp*sin(x(3, j))
            end do
        end do
        status = 0
        if (t < 1) status = 1
        associate (unused => self)
        end associate
    end subroutine

    pure subroutine bubble(x, b, g, lap)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: b, g(2), lap
        real(dp) :: r, z
        r = x(1) - 3; z = x(2)
        b = (0.25_dp - r*r)*(0.25_dp - z*z)
        g = [-2*r*(0.25_dp - z*z), -2*z*(0.25_dp - r*r)]
        lap = -1 + 2*(r*r + z*z)
    end subroutine

    subroutine source(self, t, x, f0, f, status)
        class(source_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        real(dp) :: b, g(2), lap, phi, a, ap, app, k, kp
        integer :: j
        ! Strong source of u = bubble * (1 + .3 cos(phi) + .2 sin(2 phi))
        ! on a straight cylinder of unit metric radius, with zero boundary.
        f = 0
        do j = 1, size(x, 2)
            call bubble(x(:2, j), b, g, lap)
            phi = x(3, j)
            a = 1 + 0.3_dp*cos(phi) + 0.2_dp*sin(2*phi)
            ap = -0.3_dp*sin(phi) + 0.4_dp*cos(2*phi)
            app = -0.3_dp*cos(phi) - 0.8_dp*sin(2*phi)
            k = 1 + 0.2_dp*cos(phi) + 0.1_dp*sin(phi)
            kp = -0.2_dp*sin(phi) + 0.1_dp*cos(phi)
            f0(j) = -k*(lap*a + b*app) - kp*b*ap
        end do
        status = 0
        if (t < 1) status = 1
        associate (unused => self)
        end associate
    end subroutine
end module

program test_meridian_fourier_equilibration_slow
    use fourier_equilibration_fixture, only: dp, pi, chi, modes, nb, nphi, &
        amplitude, conductivity_t, source_t, bubble
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_lagrange_nodes, meridian_qgrad
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_diffusion, &
        meridian_fourier_load
    use fortfem_meridian_fourier_equilibration, only: fourier_equilibrated_t, &
        meridian_fourier_equilibrated_estimate
    use fortfem_sparse_direct, only: sparse_direct_solve_once
    use fortsparse, only: csc_t, fortsparse_status_t, status_ok
    use fortfem_equilibrated_flux_2d, only: EQUILIBRATION_BAD_INPUT
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    type(mesh_2d_t) :: mesh
    type(meridian_mesh_t) :: mm
    type(meridian_feec_t) :: fs
    type(conductivity_t) :: coef
    type(source_t) :: src
    type(fourier_equilibrated_t) :: est
    type(csc_t) :: a
    type(fortsparse_status_t) :: ss
    real(dp), allocatable :: u(:), rhs(:), xy(:, :), lower(:), xi(:), eta(:), w(:)
    logical, allocatable :: mask(:)
    real(dp) :: b, g(2), lap, err, prev, rate, solve_time, bound_time
    integer :: n, lev, st, i, j, p, block, repeat, repeats, study_n
    integer(int64) :: c0, c1, hz
    character(32) :: arg

    repeats = 1; study_n = 0
    call get_command_argument(1, arg)
    if (len_trim(arg) > 0) read(arg, *) study_n
    call get_command_argument(2, arg)
    if (len_trim(arg) > 0) read(arg, *) repeats
    call system_clock(count_rate=hz)
    ! Exact P4 field, nonzero strong source, both cosine and sine modes.
    call build(2, 4)
    call meridian_lagrange_nodes(fs, xy)
    allocate (u(nb*fs%n0))
    do i = 1, fs%n0
        call bubble(xy(:, i), b, g, lap)
        do block = 1, nb
            u((block - 1)*fs%n0 + i) = amplitude(block)*b
        end do
    end do
    call estimate()
    print '(a,4es12.4)', 'exact bound/flux/osc/div ', est%bound, est%flux_defect, &
        est%oscillation, est%divergence_defect
    if (est%bound > 2.0e-9_dp) error stop 'exact recovery'
    if (maxval(est%mode_divergence) > 2.0e-9_dp) error stop 'mode balance'
    call estimate(9)
    if (est%bound > 2.0e-9_dp) error stop 'odd-grid exact recovery'
    ! The sampled K drops below 1: a declared chi_lower = 1 is rejected.
    call meridian_fourier_equilibrated_estimate(fs, modes, nphi, coef, src, u, &
        1.0_dp, mesh%vertices, mesh%triangles, xi, eta, w, lower, &
        1/(pi*sqrt(2*chi)), est, st)
    if (st /= EQUILIBRATION_BAD_INPUT) error stop 'ellipticity contract'
    deallocate (u)

    prev = huge(1.0_dp)
    do lev = 1, 3
        n = 2**(lev + 1)
        if (study_n > 0) n = study_n
        call build(n, 2)
        call system_clock(c0)
        call meridian_fourier_diffusion(fs, modes, nphi, coef, a, st)
        if (st /= 0) error stop 'assembly'
        call meridian_fourier_load(fs, modes, nphi, src, rhs, st)
        if (st /= 0) error stop 'load'
        allocate (mask(nb*fs%n0), u(nb*fs%n0))
        do block = 1, nb
            mask((block - 1)*fs%n0 + 1:block*fs%n0) = fs%sp(0)%bnd
        end do
        ! Homogeneous Dirichlet elimination; retain the sparse pattern.
        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (.not. mask(j) .and. .not. mask(a%row_idx(p))) cycle
                a%val(p) = 0
                if (j == a%row_idx(p)) a%val(p) = 1
            end do
        end do
        where (mask) rhs = 0
        call sparse_direct_solve_once(a, rhs, u, ss)
        if (.not. status_ok(ss)) error stop 'solve'
        call system_clock(c1)
        solve_time = real(c1 - c0, dp)/hz
        call system_clock(c0)
        do repeat = 1, repeats
            call estimate()
        end do
        call system_clock(c1)
        bound_time = real(c1 - c0, dp)/(hz*repeats)
        call energy_error(err)
        rate = 0
        if (lev > 1) rate = log(prev/err)/log(2.0_dp)
        print '(a,i4,a,i7,7(a,es12.4))', 'n ', n, ' dof ', size(u), &
            ' error ', err, ' bound ', est%bound, ' index ', est%bound/err, &
            ' rate ', rate, ' solve_s ', solve_time, ' bound_s ', bound_time, &
            ' cost ', bound_time/solve_time
        if (est%bound < err) error stop 'containment'
        if (est%bound > 1.6_dp*err) error stop 'efficiency'
        if (maxval(est%mode_divergence) > 1.0e-8_dp) error stop 'mode divergence'
        if (study_n == 0) then
            if (lev > 1) then
                if (rate < 1.8_dp) error stop 'rate'
            end if
        end if
        prev = err
        deallocate (u, mask)
        if (study_n > 0) exit
    end do
    print '(a)', 'PASS Fourier-coupled equilibrated majorant'
contains
    subroutine build(n, order)
        integer, intent(in) :: n, order
        call mesh%create_rectangular(n + 1, n + 1, 2.5_dp, 3.5_dp, -0.5_dp, 0.5_dp)
        call meridian_mesh_straight(mesh, 1, mm, st)
        if (st /= 0) error stop 'mesh'
        call meridian_feec_build(mm, order, fs, st, cylinder_r0=1.0_dp)
        if (st /= 0) error stop 'complex'
        call triangle_duffy_quadrature(2*order + 4, xi, eta, w, st)
        if (st /= 0) error stop 'rule'
        if (allocated(lower)) deallocate (lower)
        allocate (lower(mesh%n_triangles))
        lower = chi
    end subroutine

    subroutine estimate(phi_points)
        integer, intent(in), optional :: phi_points
        integer :: np
        np = nphi
        if (present(phi_points)) np = phi_points
        call meridian_fourier_equilibrated_estimate(fs, modes, np, coef, src, u, &
            chi, mesh%vertices, mesh%triangles, xi, eta, w, lower, &
            1/(pi*sqrt(2*chi)), est, st)
        if (st /= 0) error stop 'estimate'
    end subroutine

    subroutine energy_error(err)
        real(dp), intent(out) :: err
        real(dp) :: x(2), wq, val(fs%nl0), grad(2, fs%nl0), gh(3), &
            b, g(2), lap, phi, factor, deriv, weight, e(3), sumsq, c, s, ucos, usin
        integer :: t, q, j, n, off
        ! Independent direct trigonometric sum on twice the assembly phi grid.
        sumsq = 0
        do t = 1, mesh%n_triangles
            do q = 1, fs%nq
                call meridian_qgrad(fs, t, q, x, wq, val, grad, st)
                call bubble(x, b, g, lap)
                do j = 1, 2*nphi
                    phi = 2*pi*(j - 1)/(2*nphi)
                    factor = 1 + 0.3_dp*cos(phi) + 0.2_dp*sin(2*phi)
                    deriv = -0.3_dp*sin(phi) + 0.4_dp*cos(2*phi)
                    gh = 0
                    gh([1, 3]) = matmul(grad, u(fs%sp(0)%dof(:, t)))
                    do n = 1, modes
                        c = cos(n*phi); s = sin(n*phi)
                        off = (2*n - 1)*fs%n0
                        ucos = dot_product(val, u(off + fs%sp(0)%dof(:, t)))
                        gh([1, 3]) = gh([1, 3]) &
                            + c*matmul(grad, u(off + fs%sp(0)%dof(:, t)))
                        off = 2*n*fs%n0
                        usin = dot_product(val, u(off + fs%sp(0)%dof(:, t)))
                        gh([1, 3]) = gh([1, 3]) &
                            + s*matmul(grad, u(off + fs%sp(0)%dof(:, t)))
                        gh(2) = gh(2) + n*(-ucos*s + usin*c)
                    end do
                    e = [g(1)*factor, b*deriv, g(2)*factor] - gh
                    weight = 1 + 0.2_dp*cos(phi) + 0.1_dp*sin(phi)
                    sumsq = sumsq + wq*weight*dot_product(e, e)/(2*nphi)
                end do
            end do
        end do
        err = sqrt(sumsq)
    end subroutine
end program
