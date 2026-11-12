module meridian_fourier_test_coefficients
    use fortfem_kinds, only: dp
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t, &
        meridian_fourier_source_t, meridian_fourier_anisotropic_t
    implicit none
    private
    public :: identity_t, helical_t, polynomial_source_t, trig_source_t
    public :: helical_direction_t, ellipse_t, ellipse_source_t, ellipse_psi
    public :: exact_u, exact_grad, helical_tensor

    type, extends(meridian_fourier_tensor_t) :: identity_t
    contains
        procedure :: tensor => identity_tensor
    end type

    type, extends(meridian_fourier_tensor_t) :: helical_t
    contains
        procedure :: tensor => helical_eval
    end type

    !> helical_tensor in factored form: k_perp = 1e-3, k_par = 1 + 1e-3.
    type, extends(meridian_fourier_anisotropic_t) :: helical_direction_t
    contains
        procedure :: direction => helical_direction
    end type

    !> b along the level sets of ellipse_psi, k_par/k_perp = 1e12.
    type, extends(meridian_fourier_anisotropic_t) :: ellipse_t
    contains
        procedure :: direction => ellipse_direction
    end type

    !> Weak flux f = k_perp grad ellipse_psi (b . grad psi = 0).
    type, extends(meridian_fourier_source_t) :: ellipse_source_t
    contains
        procedure :: source => ellipse_source
    end type

    !> Weak flux f = K grad u of the trigonometric polynomial exact_u.
    type, extends(meridian_fourier_source_t) :: polynomial_source_t
    contains
        procedure :: source => polynomial_source
    end type

    !> f0 = g(R, Z) cos(2 phi) + h(R, Z) sin(phi), no flux.
    type, extends(meridian_fourier_source_t) :: trig_source_t
    contains
        procedure :: source => trig_source
    end type

contains

    subroutine identity_tensor(self, t, x, k, status)
        class(identity_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: status
        integer :: p
        k = 0.0_dp
        do p = 1, size(x, 2)
            k(1, 1, p) = 1.0_dp; k(2, 2, p) = 1.0_dp; k(3, 3, p) = 1.0_dp
        end do
        status = 0
    end subroutine identity_tensor

    pure subroutine helical_tensor(x, k)
        !! K = 1e-3 I + b b for a smooth phi-dependent unit vector b.
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: k(3, 3)
        real(dp) :: b(3)
        integer :: i
        b = [0.3_dp*x(2)*cos(2*x(3)), 3.0_dp/x(1), &
            0.2_dp*sin(x(3)) + (x(1) - 3.0_dp)]
        b = b/norm2(b)
        do i = 1, 3
            k(:, i) = b*b(i)
            k(i, i) = k(i, i) + 1.0e-3_dp
        end do
    end subroutine helical_tensor

    subroutine helical_eval(self, t, x, k, status)
        class(helical_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: status
        integer :: p
        do p = 1, size(x, 2)
            call helical_tensor(x(:, p), k(:, :, p))
        end do
        status = 0
    end subroutine helical_eval

    subroutine helical_direction(self, t, x, b, k_par, k_perp, status)
        class(helical_direction_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: b(:, :), k_par(:), k_perp(:)
        integer, intent(out) :: status
        integer :: p
        do p = 1, size(x, 2)
            b(:, p) = [0.3_dp*x(2, p)*cos(2*x(3, p)), 3.0_dp/x(1, p), &
                0.2_dp*sin(x(3, p)) + (x(1, p) - 3.0_dp)]
            b(:, p) = b(:, p)/norm2(b(:, p))
        end do
        k_perp = 1.0e-3_dp
        k_par = 1.0_dp + 1.0e-3_dp
        status = 0
    end subroutine helical_direction

    pure real(dp) function ellipse_psi(x) result(psi)
        real(dp), intent(in) :: x(2)
        psi = (x(1) - 3.0_dp)**2 + 2*x(2)**2
    end function ellipse_psi

    subroutine ellipse_direction(self, t, x, b, k_par, k_perp, status)
        class(ellipse_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: b(:, :), k_par(:), k_perp(:)
        integer, intent(out) :: status
        integer :: p
        do p = 1, size(x, 2)
            b(:, p) = [-4*x(2, p), 0.0_dp, 2*(x(1, p) - 3.0_dp)]
            b(:, p) = b(:, p)/norm2(b(:, p))
        end do
        k_perp = 1.0_dp
        k_par = 1.0e12_dp
        status = 0
    end subroutine ellipse_direction

    subroutine ellipse_source(self, t, x, f0, f, status)
        class(ellipse_source_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        integer :: p
        f0 = 0.0_dp
        do p = 1, size(x, 2)
            f(:, p) = [2*(x(1, p) - 3.0_dp), 0.0_dp, 4*x(2, p)]
        end do
        status = 0
    end subroutine ellipse_source

    pure real(dp) function exact_u(x) result(u)
        !! Quadratic in (R, Z) per Fourier component, modes 0, 1 and 2.
        real(dp), intent(in) :: x(3)
        u = x(1)**2 - x(1)*x(2) + 0.5_dp &
            + (x(2)**2 + 0.3_dp*x(1))*cos(x(3)) &
            + (x(1)*x(2) - 1.0_dp)*sin(2*x(3))
    end function exact_u

    pure function exact_grad(x) result(g)
        !! Physical gradient (d_R, (1/R) d_phi, d_Z).
        real(dp), intent(in) :: x(3)
        real(dp) :: g(3)
        g(1) = 2*x(1) - x(2) + 0.3_dp*cos(x(3)) + x(2)*sin(2*x(3))
        g(2) = (-(x(2)**2 + 0.3_dp*x(1))*sin(x(3)) &
            + 2*(x(1)*x(2) - 1.0_dp)*cos(2*x(3)))/x(1)
        g(3) = -x(1) + 2*x(2)*cos(x(3)) + x(1)*sin(2*x(3))
    end function exact_grad

    subroutine polynomial_source(self, t, x, f0, f, status)
        class(polynomial_source_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        real(dp) :: k(3, 3)
        integer :: p
        f0 = 0.0_dp
        do p = 1, size(x, 2)
            call helical_tensor(x(:, p), k)
            f(:, p) = matmul(k, exact_grad(x(:, p)))
        end do
        status = 0
    end subroutine polynomial_source

    subroutine trig_source(self, t, x, f0, f, status)
        class(trig_source_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        integer :: p
        f = 0.0_dp
        do p = 1, size(x, 2)
            f0(p) = exp(x(1, p) - 3.0_dp)*cos(2*x(3, p)) + x(2, p)**3*sin(x(3, p))
        end do
        status = 0
    end subroutine trig_source
end module meridian_fourier_test_coefficients

program test_meridian_fourier_scalar
    !! Real Fourier extension of the meridional V0 space:
    !!  - with K = I the mode blocks equal the per-n Dirichlet forms
    !!    grad_n^H M_1 grad_n (n = 0) and half of them (cos, sin; n > 0),
    !!    with zero cross-mode blocks;
    !!  - the scalar load equals half of the per-n load of each component;
    !!  - a Galerkin solve with a phi-dependent anisotropic tensor
    !!    (ratio 1e3) reproduces a trigonometric polynomial with quadratic
    !!    (R, Z) components, which lies in the P2 space, to roundoff;
    !!  - the matrix-free factored apply equals the assembled form;
    !!  - at k_par/k_perp = 1e12 with b along the level sets of a quadratic
    !!    psi, refinement with the factored residual returns psi to
    !!    roundoff.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_stiffness, meridian_load, meridian_lagrange_nodes, &
        meridian_real_to_complex
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_diffusion, &
        meridian_fourier_load, meridian_fourier_index, &
        meridian_fourier_anisotropic_apply
    use fortsparse, only: csc_t, csc_z_t, csc_matvec
    use meridian_test_fields, only: solve_z, zero_rows
    use meridian_test_meshes, only: make_annulus_mesh
    use meridian_fourier_test_coefficients, only: identity_t, helical_t, &
        polynomial_source_t, trig_source_t, exact_u, helical_direction_t, &
        ellipse_t, ellipse_source_t, ellipse_psi
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: modes = 2, nphi = 9
    type(meridian_mesh_t) :: mm
    type(meridian_feec_t) :: fs
    type(identity_t) :: iso
    type(helical_t) :: hel
    type(polynomial_source_t) :: psrc
    type(trig_source_t) :: tsrc
    type(helical_direction_t) :: hdir
    type(ellipse_t) :: ell
    type(ellipse_source_t) :: esrc
    type(csc_t) :: a, az_check_a
    type(csc_z_t) :: kn, az
    real(dp), allocatable :: v(:), full(:), y(:), ref(:), b(:), xy(:, :)
    complex(dp), allocatable :: zb(:), zx(:), zr(:)
    logical, allocatable :: mask(:)
    real(dp) :: err, scale
    integer :: st, n, part, n0, nb, i, j, n_other
    logical :: ok, all_passed

    all_passed = .true.
    call make_annulus_mesh(3, 12, 0.3_dp, 1.0_dp, 3.0_dp, 1, mm)
    call meridian_feec_build(mm, 2, fs, st)
    call record(st == 0, "build P2 complex")
    n0 = fs%n0
    nb = 2*modes + 1

    ! Isotropic blocks against the per-n Dirichlet forms.
    call meridian_fourier_diffusion(fs, modes, nphi, iso, a, st)
    call record(st == 0, "isotropic Fourier diffusion assembly")
    allocate (v(n0), full(n0*nb), y(n0*nb), ref(n0))
    do i = 1, n0
        v(i) = sin(0.7_dp*i) + 0.1_dp*cos(1.3_dp*i*i)
    end do
    err = 0.0_dp; scale = 0.0_dp
    do n = 0, modes
        call meridian_stiffness(fs, 0, n, kn, st)
        ref = real(csc_matvec(kn, cmplx(v, 0.0_dp, dp)), dp)
        if (n > 0) ref = 0.5_dp*ref
        scale = max(scale, maxval(abs(ref)))
        do part = 0, merge(0, 1, n == 0)
            full = 0.0_dp
            full(meridian_fourier_index(n0, n, part, 1): &
                meridian_fourier_index(n0, n, part, n0)) = v
            y = csc_matvec(a, full)
            do j = 0, modes
                do i = 0, merge(0, 1, j == 0)
                    associate (blk => y(meridian_fourier_index(n0, j, i, 1): &
                            meridian_fourier_index(n0, j, i, n0)))
                        if (j == n .and. i == part) then
                            err = max(err, maxval(abs(blk - ref)))
                        else
                            err = max(err, maxval(abs(blk)))
                        end if
                    end associate
                end do
            end do
        end do
    end do
    print '(a, es10.3)', "  isotropic block defect / scale: ", err/scale
    call record(err <= 1.0e-12_dp*scale, "K = I blocks equal per-n forms")

    ! Scalar load against the per-n loads.
    call meridian_fourier_load(fs, modes, nphi, tsrc, b, st)
    call record(st == 0, "Fourier load assembly")
    call meridian_load(fs, 0, g_cos, zb, st)
    err = maxval(abs(b(meridian_fourier_index(n0, 2, 0, 1): &
        meridian_fourier_index(n0, 2, 0, n0)) - 0.5_dp*real(zb, dp)))
    scale = maxval(abs(real(zb, dp)))
    call meridian_load(fs, 0, h_sin, zb, st)
    err = max(err, maxval(abs(b(meridian_fourier_index(n0, 1, 1, 1): &
        meridian_fourier_index(n0, 1, 1, n0)) - 0.5_dp*real(zb, dp))))
    n_other = 0
    do n = 0, modes
        do part = 0, merge(0, 1, n == 0)
            if (n == 2 .and. part == 0) cycle
            if (n == 1 .and. part == 1) cycle
            err = max(err, maxval(abs(b(meridian_fourier_index(n0, n, part, 1): &
                meridian_fourier_index(n0, n, part, n0)))))
            n_other = n_other + 1
        end do
    end do
    print '(a, es10.3)', "  load defect / scale: ", err/scale
    call record(err <= 1.0e-12_dp*scale .and. n_other == 3, &
        "scalar load equals half the per-n loads")

    ! Anisotropic, phi-dependent Galerkin reproduction.
    call meridian_fourier_diffusion(fs, modes, nphi, hel, a, st)
    call record(st == 0, "anisotropic Fourier diffusion assembly")
    call meridian_fourier_load(fs, modes, nphi, psrc, b, st)
    call record(st == 0, "weak flux load assembly")
    call meridian_lagrange_nodes(fs, xy)
    allocate (mask(n0*nb), zx(n0*nb), zr(n0*nb))
    deallocate (full); allocate (full(n0*nb))
    full = 0.0_dp
    do i = 1, n0
        full(meridian_fourier_index(n0, 0, 0, i)) = component(xy(:, i), 0, 0)
        do n = 1, modes
            do part = 0, 1
                full(meridian_fourier_index(n0, n, part, i)) = &
                    component(xy(:, i), n, part)
            end do
        end do
    end do
    do part = 0, nb - 1
        mask(part*n0 + 1:(part + 1)*n0) = fs%sp(0)%bnd
    end do
    y = full
    where (.not. mask) y = 0.0_dp
    zr = cmplx(b - csc_matvec(a, y), 0.0_dp, dp)
    where (mask) zr = 0.0_dp
    call meridian_real_to_complex(a, az)
    call zero_rows(az, mask)
    call solve_z(az, zr, zx, ok)
    call record(ok, "anisotropic solve")
    y = y + real(zx, dp)
    err = maxval(abs(y - full))/maxval(abs(full))
    print '(a, es10.3)', "  trigonometric-polynomial reproduction error: ", err
    call record(err <= 1.0e-9_dp, "anisotropic Galerkin reproduces P2 x trig")

    ! Factored apply against the assembled form of the same tensor.
    call meridian_fourier_diffusion(fs, modes, nphi, hdir, a, st)
    call record(st == 0, "factored tensor assembly")
    y = csc_matvec(a, full)
    call meridian_fourier_anisotropic_apply(fs, modes, nphi, hdir, full, ref, st)
    call record(st == 0, "factored apply")
    err = maxval(abs(ref - y))/maxval(abs(y))
    print '(a, es10.3)', "  factored apply vs assembled: ", err
    call record(err <= 1.0e-13_dp, "factored apply equals the assembled form")
    call meridian_fourier_diffusion(fs, modes, nphi, hel, az_check_a, st)
    err = maxval(abs(csc_matvec(az_check_a, full) - y))/maxval(abs(y))
    call record(err <= 1.0e-13_dp, "default tensor equals k_perp I + (k_par - k_perp) b b")

    ! Ratio 1e12, exact discrete solution psi.
    call ellipse_refinement(err)
    print '(a, es10.3)', "  1e12 error with the factored residual: ", err
    call record(err <= 1.0e-12_dp, "factored refinement reaches roundoff at ratio 1e12")

    call check_summary("Meridian real Fourier scalar diffusion")
    if (.not. all_passed) error stop 1

contains

    subroutine ellipse_refinement(err)
        !! Dirichlet solve, axisymmetric block, u = psi on the boundary; the
        !! assembled factor only preconditions the refinement.
        real(dp), intent(out) :: err
        type(csc_t) :: ae
        type(csc_z_t) :: aez
        real(dp), allocatable :: u(:), rhs(:), ex(:), r(:)
        complex(dp), allocatable :: cr(:), cx(:)
        logical, allocatable :: fixed(:)
        integer :: it
        logical :: solved
        call meridian_fourier_diffusion(fs, 0, 1, ell, ae, st)
        call meridian_fourier_load(fs, 0, 1, esrc, rhs, st)
        allocate (ex(n0), fixed(n0), cr(n0), cx(n0), r(n0))
        do i = 1, n0
            ex(i) = ellipse_psi(xy(:, i))
        end do
        fixed = fs%sp(0)%bnd
        call meridian_real_to_complex(ae, aez)
        call zero_rows(aez, fixed)
        u = merge(ex, 0.0_dp, fixed)
        do it = 0, 8
            call meridian_fourier_anisotropic_apply(fs, 0, 1, ell, u, r, st)
            r = rhs - r
            cr = cmplx(merge(0.0_dp, r, fixed), 0.0_dp, dp)
            call solve_z(aez, cr, cx, solved)
            u = u + real(cx, dp)
        end do
        err = maxval(abs(u - ex))/maxval(abs(ex))
    end subroutine ellipse_refinement

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    subroutine g_cos(x, w)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: w(3)
        w = 0.0_dp
        w(1) = exp(x(1) - 3.0_dp)
    end subroutine g_cos

    subroutine h_sin(x, w)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: w(3)
        w = 0.0_dp
        w(1) = x(2)**3
    end subroutine h_sin

    real(dp) function component(x, n, part) result(c)
        !! Fourier component of exact_u at x by an exact DFT.
        real(dp), intent(in) :: x(2)
        integer, intent(in) :: n, part
        integer :: j
        real(dp) :: phi
        c = 0.0_dp
        do j = 0, 15
            phi = 2*pi*j/16
            if (n == 0) then
                c = c + exact_u([x, phi])/16
            else if (part == 0) then
                c = c + 2*exact_u([x, phi])*cos(n*phi)/16
            else
                c = c + 2*exact_u([x, phi])*sin(n*phi)/16
            end if
        end do
    end function component

end program test_meridian_fourier_scalar
