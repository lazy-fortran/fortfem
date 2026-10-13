module fortfem_meridian_fourier_scalar
    !! Real scalar fields on the meridional V0 space times a truncated real
    !! Fourier series in the toroidal (or axial) angle phi,
    !!   u(R, Z, phi) = u_0 + sum_{n=1}^N (u_n^c cos(n phi) + u_n^s sin(n phi)),
    !! with variable, phi-dependent symmetric tensor coefficients.  This is
    !! the mode-coupled counterpart of the per-n complex: a coefficient that
    !! depends on phi couples all retained modes.
    !!
    !! Layout: block b = 0 holds u_0, b = 2n - 1 the cosine and b = 2n the
    !! sine amplitudes of mode n; the global index of V0 dof i in block b is
    !! b*n0 + i.  The measure is R dR dZ (the 2 pi of the phi integral is
    !! dropped, as in fortfem_meridian_feec_n) times the phi mean
    !! (1/2 pi) int dphi, evaluated with nphi uniform points
    !! phi_j = 2 pi j / nphi.  That rule is exact for trigonometric
    !! polynomials of degree < nphi; nphi >= 2N + 1 is required so that the
    !! isotropic forms are exact, and a phi-dependent coefficient carries a
    !! declared phi-quadrature (aliasing) error.  Physical components are
    !! ordered (R, phi, Z) as in the per-n complex; the phi derivative
    !! contributes (1/R) d/dphi with the metric radius R (R0 on the
    !! periodic-cylinder chart).
    !!
    !! Forms: a(v, u) = < grad v . K grad u > and l(v) = < v f0 + grad v . f >.
    !! Coefficients are evaluated by batched callbacks per triangle, so a
    !! consumer can evaluate a field defined on the same mesh without point
    !! location.  Element loops are OpenMP parallel.
    use fortfem_kinds, only: dp
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_qgrad
    use fortsparse, only: csc_t, csc_from_triplet, fortsparse_status_t, status_ok
    implicit none
    private

    type, abstract, public :: meridian_fourier_tensor_t
    contains
        procedure(fourier_tensor_i), deferred :: tensor
    end type meridian_fourier_tensor_t

    type, abstract, public :: meridian_fourier_source_t
    contains
        procedure(fourier_source_i), deferred :: source
    end type meridian_fourier_source_t

    abstract interface
        subroutine fourier_tensor_i(self, t, x, k, status)
            !! Symmetric physical tensor k(3, 3, m) in (R, phi, Z) at the
            !! points x(3, m) = (R, Z, phi) of triangle t.
            import :: dp, meridian_fourier_tensor_t
            class(meridian_fourier_tensor_t), intent(in) :: self
            integer, intent(in) :: t
            real(dp), intent(in) :: x(:, :)
            real(dp), intent(out) :: k(:, :, :)
            integer, intent(out) :: status
        end subroutine fourier_tensor_i

        subroutine fourier_source_i(self, t, x, f0, f, status)
            !! Scalar density f0(m) and physical flux f(3, m) in (R, phi, Z)
            !! at the points x(3, m) = (R, Z, phi) of triangle t.
            import :: dp, meridian_fourier_source_t
            class(meridian_fourier_source_t), intent(in) :: self
            integer, intent(in) :: t
            real(dp), intent(in) :: x(:, :)
            real(dp), intent(out) :: f0(:), f(:, :)
            integer, intent(out) :: status
        end subroutine fourier_source_i
    end interface

    public :: meridian_fourier_diffusion, meridian_fourier_load
    public :: meridian_fourier_index, meridian_fourier_value

    integer, parameter :: par_min = 32
    real(dp), parameter :: pi = acos(-1.0_dp)

contains

    pure integer function meridian_fourier_index(n0, n, part, i) result(g)
        !! Global index of V0 dof i of mode n; part 0 cosine (and n = 0),
        !! part 1 sine.
        integer, intent(in) :: n0, n, part, i
        if (n == 0) then
            g = i
        else
            g = (2*n - 1 + part)*n0 + i
        end if
    end function meridian_fourier_index

    pure subroutine basis_trig(modes, phi, rr, val, grad, gv)
        !! Physical gradients gv(3, nl*(2N+1)) (R, phi, Z) of the local
        !! basis functions val(i) * trig_b(phi) at one point.
        integer, intent(in) :: modes
        real(dp), intent(in) :: phi, rr, val(:), grad(:, :)
        real(dp), intent(out) :: gv(:, :)
        real(dp) :: c, s
        integer :: n, nl, off

        nl = size(val)
        gv(1, :nl) = grad(1, :)
        gv(2, :nl) = 0.0_dp
        gv(3, :nl) = grad(2, :)
        do n = 1, modes
            c = cos(n*phi)
            s = sin(n*phi)
            off = (2*n - 1)*nl
            gv(1, off + 1:off + nl) = grad(1, :)*c
            gv(2, off + 1:off + nl) = -n*val*s/rr
            gv(3, off + 1:off + nl) = grad(2, :)*c
            off = off + nl
            gv(1, off + 1:off + nl) = grad(1, :)*s
            gv(2, off + 1:off + nl) = n*val*c/rr
            gv(3, off + 1:off + nl) = grad(2, :)*s
        end do
    end subroutine basis_trig

    pure subroutine basis_values(modes, phi, val, v)
        integer, intent(in) :: modes
        real(dp), intent(in) :: phi, val(:)
        real(dp), intent(out) :: v(:)
        integer :: n, nl, off

        nl = size(val)
        v(:nl) = val
        do n = 1, modes
            off = (2*n - 1)*nl
            v(off + 1:off + nl) = val*cos(n*phi)
            v(off + nl + 1:off + 2*nl) = val*sin(n*phi)
        end do
    end subroutine basis_values

    pure real(dp) function radius(fs, x) result(r)
        type(meridian_feec_t), intent(in) :: fs
        real(dp), intent(in) :: x(2)
        r = x(1)
        if (fs%cyl_r0 > 0.0_dp) r = fs%cyl_r0
    end function radius

    subroutine meridian_fourier_diffusion(fs, modes, nphi, coef, a, status)
        !! a_(gi, gj) = < grad phi_gi . K grad phi_gj > on the Fourier-V0
        !! space; symmetric, and positive definite on the complement of the
        !! constants for positive definite K.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes, nphi
        class(meridian_fourier_tensor_t), intent(in) :: coef
        type(csc_t), intent(out) :: a
        integer, intent(out) :: status

        integer, allocatable :: rows(:), cols(:), gdof(:)
        real(dp), allocatable :: vals(:), ml(:, :), gv(:, :), kg(:, :)
        real(dp), allocatable :: val(:, :), grad(:, :, :), xq(:, :), w(:)
        real(dp), allocatable :: kt(:, :, :)
        real(dp) :: x(2)
        type(fortsparse_status_t) :: fst
        integer :: nt, nl, nb, m, t, q, j, p, i, base, bad, st, np, ntot

        status = 1
        if (modes < 0 .or. nphi < 2*modes + 1) return
        nt = fs%mm%topology%n_triangles
        nl = fs%nl0
        nb = 2*modes + 1
        m = nl*nb
        np = fs%nq*nphi
        ntot = fs%n0*nb
        allocate (rows(m*m*nt), cols(m*m*nt), vals(m*m*nt))
        bad = 0
        !$omp parallel if(nt >= par_min) default(shared) &
        !$omp private(t, q, j, p, i, base, st, ml, gv, kg, val, grad, xq, w, &
        !$omp kt, gdof, x) reduction(max:bad)
        allocate (ml(m, m), gv(3, m), kg(3, m), val(nl, fs%nq), &
            grad(2, nl, fs%nq), xq(3, np), w(fs%nq), kt(3, 3, np), gdof(m))
        !$omp do schedule(dynamic, 4)
        do t = 1, nt
            ml = 0.0_dp
            st = 0
            do q = 1, fs%nq
                call meridian_qgrad(fs, t, q, x, w(q), val(:, q), grad(:, :, q), &
                    st)
                if (st /= 0) exit
                do j = 1, nphi
                    xq(:, (q - 1)*nphi + j) = [x, 2*pi*(j - 1)/nphi]
                end do
            end do
            if (st == 0) call coef%tensor(t, xq, kt, st)
            if (st /= 0) then
                bad = 1
                cycle
            end if
            do q = 1, fs%nq
                do j = 1, nphi
                    p = (q - 1)*nphi + j
                    call basis_trig(modes, xq(3, p), radius(fs, xq(:2, p)), &
                        val(:, q), grad(:, :, q), gv)
                    kg = matmul(kt(:, :, p), gv)
                    do i = 1, m
                        ml(:, i) = ml(:, i) + (w(q)/nphi)*(gv(1, :)*kg(1, i) &
                            + gv(2, :)*kg(2, i) + gv(3, :)*kg(3, i))
                    end do
                end do
            end do
            call local_dofs(t, gdof)
            base = (t - 1)*m*m
            do j = 1, m
                do i = 1, m
                    base = base + 1
                    rows(base) = gdof(i)
                    cols(base) = gdof(j)
                    vals(base) = ml(i, j)
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        if (bad /= 0) return
        call csc_from_triplet(ntot, ntot, rows, cols, vals, a, fst)
        if (status_ok(fst)) status = 0

    contains

        subroutine local_dofs(t, gdof)
            integer, intent(in) :: t
            integer, intent(out) :: gdof(:)
            integer :: b
            do b = 0, nb - 1
                gdof(b*nl + 1:(b + 1)*nl) = b*fs%n0 + fs%sp(0)%dof(:, t)
            end do
        end subroutine local_dofs
    end subroutine meridian_fourier_diffusion

    subroutine meridian_fourier_load(fs, modes, nphi, src, b, status)
        !! b_gi = < phi_gi f0 + grad phi_gi . f > on the Fourier-V0 space.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes, nphi
        class(meridian_fourier_source_t), intent(in) :: src
        real(dp), allocatable, intent(out) :: b(:)
        integer, intent(out) :: status

        real(dp), allocatable :: bl(:, :), gv(:, :), v(:), val(:), grad(:, :)
        real(dp), allocatable :: xq(:, :), w(:), f0(:), f(:, :)
        real(dp) :: x(2)
        integer :: nt, nl, nb, m, t, q, j, p, i, bad, st, np, blk

        status = 1
        if (modes < 0 .or. nphi < 2*modes + 1) return
        nt = fs%mm%topology%n_triangles
        nl = fs%nl0
        nb = 2*modes + 1
        m = nl*nb
        np = fs%nq*nphi
        allocate (bl(m, nt), b(fs%n0*nb))
        bad = 0
        !$omp parallel if(nt >= par_min) default(shared) &
        !$omp private(t, q, j, p, i, st, gv, v, val, grad, xq, w, f0, f, x) &
        !$omp reduction(max:bad)
        allocate (gv(3, m), v(m), val(nl), grad(2, nl), xq(3, np), w(np), &
            f0(np), f(3, np))
        !$omp do schedule(dynamic, 4)
        do t = 1, nt
            bl(:, t) = 0.0_dp
            st = 0
            do q = 1, fs%nq
                call meridian_qgrad(fs, t, q, x, w(q), val, grad, st)
                if (st /= 0) exit
                do j = 1, nphi
                    xq(:, (q - 1)*nphi + j) = [x, 2*pi*(j - 1)/nphi]
                end do
            end do
            if (st == 0) call src%source(t, xq, f0, f, st)
            if (st /= 0) then
                bad = 1
                cycle
            end if
            do q = 1, fs%nq
                call meridian_qgrad(fs, t, q, x, w(q), val, grad, st)
                do j = 1, nphi
                    p = (q - 1)*nphi + j
                    call basis_trig(modes, xq(3, p), radius(fs, x), val, grad, gv)
                    call basis_values(modes, xq(3, p), val, v)
                    do i = 1, m
                        bl(i, t) = bl(i, t) + (w(q)/nphi)*(v(i)*f0(p) &
                            + gv(1, i)*f(1, p) + gv(2, i)*f(2, p) + gv(3, i)*f(3, p))
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        if (bad /= 0) return
        b = 0.0_dp
        do t = 1, nt
            do blk = 0, nb - 1
                do i = 1, nl
                    p = blk*fs%n0 + fs%sp(0)%dof(i, t)
                    b(p) = b(p) + bl(blk*nl + i, t)
                end do
            end do
        end do
        status = 0
    end subroutine meridian_fourier_load

    pure real(dp) function meridian_fourier_value(n0, modes, u, dofs, val, phi) &
            result(s)
        !! u at a point with local V0 dofs(:) and basis values val(:).
        integer, intent(in) :: n0, modes, dofs(:)
        real(dp), intent(in) :: u(:), val(:), phi
        integer :: n

        s = dot_product(val, u(dofs))
        do n = 1, modes
            s = s + cos(n*phi)*dot_product(val, u((2*n - 1)*n0 + dofs)) &
                + sin(n*phi)*dot_product(val, u(2*n*n0 + dofs))
        end do
    end function meridian_fourier_value

end module fortfem_meridian_fourier_scalar
