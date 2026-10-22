module fortfem_meridian_fourier_curl
    !! Mode-coupled real-Fourier forms on the meridional curl complex
    !! V0 -> V1 -> V2 (Lagrange, Nedelec x Lagrange, Raviart-Thomas x DG) times
    !! a truncated real Fourier series in the toroidal (or axial) angle phi,
    !!   F(R, Z, phi) = F_0 + sum_{n=1}^N (F_n^c cos(n phi) + F_n^s sin(n phi)),
    !! with symmetric tensor coefficients that depend on (R, Z, phi).  It is
    !! the curl-complex counterpart of fortfem_meridian_fourier_scalar and
    !! uses the same layout: block b = 0 holds F_0, b = 2n - 1 the cosine and
    !! b = 2n the sine amplitudes of mode n; the global index of dof i of the
    !! per-n space in block b is b*ndof + i (meridian_fourier_index with
    !! n0 = ndof of the space: n0 for V0, nn + n0 for V1, nn + nd for V2).
    !!
    !! Measure and quadrature are those of the scalar module: R dR dZ (the
    !! 2 pi of the phi integral dropped) times the phi mean (1/2 pi) int dphi
    !! on nphi uniform points phi_j = 2 pi j / nphi.  The rule is exact for
    !! trigonometric polynomials of degree < nphi, so for a coefficient of
    !! Fourier degree d the forms are exact (in phi) for nphi >= 2N + d + 1;
    !! nphi >= 2N + 1 is required.  Physical components are (R, phi, Z).
    !!
    !! Operators.  With the per-n real operators grad_n = a + i n b and
    !! curl_n = a + i n b of fortfem_meridian_feec_n, the real-block
    !! operator acting on (cos, sin) amplitudes is
    !!   n = 0:   a,
    !!   n > 0:   [[a, n b], [-n b, a]]     (rows (cos, sin), columns (cos, sin)),
    !! the real form of multiplication by a + i n b on Re[(c - i s) e^{i n phi}].
    !! The operators do not couple modes, and since grad_n V0 is in V1 and
    !! curl_n V1 is in V2 exactly, curl_F grad_F = 0 holds as a matrix
    !! identity for every mode count.
    !!
    !! Forms (K is a symmetric tensor of physical components):
    !!   meridian_fourier_mass       M(v, u)   = < v . K u >, u, v in V1 or V2;
    !!   meridian_fourier_curl_curl  S(v, u)   = < curl v . K curl u >, V1 x V1
    !!                               (K acts on the V2 field B = curl A, e.g.
    !!                               a reluctivity tensor);
    !!   meridian_fourier_gradient_coupling
    !!                               C(v, lambda) = < v . K grad lambda >,
    !!                               v in V1, lambda in V0.
    !! The mass form is assembled element by element; the other two are
    !! the exact products C_F^T M_2(K) C_F and M_1(K) G_F because the discrete
    !! curl and gradient map the spaces into each other exactly.  A
    !! phi-dependent K couples the modes through its Fourier content only:
    !! block (b, b') is nonzero only if K has a harmonic of degree
    !! |n_b - n_b'| or n_b + n_b'.  Element loops are OpenMP parallel.
    use fortfem_kinds, only: dp
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_qpoint
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t
    use fortsparse, only: csc_t, csc_from_triplet, csc_matmul, csc_transpose, &
        fortsparse_status_t, status_ok
    implicit none
    private

    public :: meridian_fourier_mass, meridian_fourier_curl_curl
    public :: meridian_fourier_gradient_coupling
    public :: meridian_fourier_gradient, meridian_fourier_curl
    public :: meridian_fourier_identity_forms

    integer, parameter :: par_min = 32
    real(dp), parameter :: pi = acos(-1.0_dp)

contains

    subroutine meridian_fourier_identity_forms(fs, modes, c, s, mg, status)
        !! Exact phi integration for K=I. The already assembled meridian
        !! mass is copied with weights 1, 1/2, 1/2, ... . Storage is linear
        !! in the number of real Fourier blocks, including construction.
        !! This avoids the general tensor builder's all-mode-pair triplets.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes
        type(csc_t), intent(out) :: c, s, mg
        integer, intent(out) :: status
        type(csc_t) :: m, g, ct, mc
        type(fortsparse_status_t) :: st
        call meridian_fourier_curl(fs, modes, c, status)
        if (status /= 0) return
        call identity_mass(2, m)
        if (status /= 0) return
        call csc_matmul(m, c, mc, st)
        status = 1
        if (.not. status_ok(st)) return
        call csc_transpose(c, ct, st)
        if (.not. status_ok(st)) return
        call csc_matmul(ct, mc, s, st)
        if (.not. status_ok(st)) return
        call meridian_fourier_gradient(fs, modes, g, status)
        if (status /= 0) return
        call identity_mass(1, m)
        if (status /= 0) return
        status = 1
        call csc_matmul(m, g, mg, st)
        if (status_ok(st)) status = 0
    contains
        subroutine identity_mass(space, matrix)
            integer, intent(in) :: space
            type(csc_t), intent(out) :: matrix
            integer :: b, j, p, offset, pos, nb, nd
            real(dp) :: weight
            nb = 2*modes + 1; nd = fs%sp(space)%ndof
            matrix%nrow = nb*nd; matrix%ncol = nb*nd
            matrix%nnz = nb*fs%mass(space)%nnz
            allocate (matrix%col_ptr(matrix%ncol + 1), &
                matrix%row_idx(matrix%nnz), matrix%val(matrix%nnz))
            !$omp parallel do default(shared) private(b, j, p, offset, pos, weight) &
            !$omp& schedule(static)
            do b = 0, nb - 1
                offset = b*fs%mass(space)%nnz
                weight = merge(1.0_dp, 0.5_dp, b == 0)
                do j = 1, nd
                    matrix%col_ptr(b*nd + j) = fs%mass(space)%col_ptr(j) + offset
                    do p = fs%mass(space)%col_ptr(j), fs%mass(space)%col_ptr(j + 1) - 1
                        pos = offset + p
                        matrix%row_idx(pos) = b*nd + fs%mass(space)%row_idx(p)
                        matrix%val(pos) = weight*fs%mass(space)%val(p)
                    end do
                end do
            end do
            !$omp end parallel do
            matrix%col_ptr(matrix%ncol + 1) = matrix%nnz + 1
            status = 0
        end subroutine
    end subroutine

    subroutine block_operator(a, b, modes, op, status)
        !! Real-block operator of a + i n b on the (cos, sin) amplitudes.
        type(csc_t), intent(in) :: a, b
        integer, intent(in) :: modes
        type(csc_t), intent(out) :: op
        integer, intent(out) :: status

        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: fst
        integer :: nr, nc, cap, pos, n, j, p

        nr = a%nrow
        nc = a%ncol
        cap = a%nnz + modes*(2*a%nnz + 2*b%nnz)
        allocate (rows(cap), cols(cap), vals(cap))
        pos = 0
        call put(a, 0, 0, 1.0_dp)
        do n = 1, modes
            call put(a, (2*n - 1)*nr, (2*n - 1)*nc, 1.0_dp)
            call put(a, 2*n*nr, 2*n*nc, 1.0_dp)
            call put(b, (2*n - 1)*nr, 2*n*nc, real(n, dp))
            call put(b, 2*n*nr, (2*n - 1)*nc, -real(n, dp))
        end do
        call csc_from_triplet((2*modes + 1)*nr, (2*modes + 1)*nc, &
            rows(1:pos), cols(1:pos), vals(1:pos), op, fst)
        status = merge(0, 1, status_ok(fst))

    contains

        subroutine put(m, ro, co, s)
            type(csc_t), intent(in) :: m
            integer, intent(in) :: ro, co
            real(dp), intent(in) :: s
            do j = 1, m%ncol
                do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                    pos = pos + 1
                    rows(pos) = m%row_idx(p) + ro
                    cols(pos) = j + co
                    vals(pos) = s*m%val(p)
                end do
            end do
        end subroutine put
    end subroutine block_operator

    subroutine meridian_fourier_gradient(fs, modes, g, status)
        !! Real-block gradient V0 -> V1 on the Fourier spaces.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes
        type(csc_t), intent(out) :: g
        integer, intent(out) :: status

        status = 1
        if (modes < 0) return
        call block_operator(fs%grad%a, fs%grad%b, modes, g, status)
    end subroutine meridian_fourier_gradient

    subroutine meridian_fourier_curl(fs, modes, c, status)
        !! Real-block curl V1 -> V2 on the Fourier spaces.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes
        type(csc_t), intent(out) :: c
        integer, intent(out) :: status

        status = 1
        if (modes < 0) return
        call block_operator(fs%curl%a, fs%curl%b, modes, c, status)
    end subroutine meridian_fourier_curl

    subroutine meridian_fourier_mass(fs, space, modes, nphi, coef, m, status)
        !! M_(gi, gj) = < phi_gi . K phi_gj > on the Fourier-V1 (space = 1) or
        !! Fourier-V2 (space = 2) space; symmetric, and positive definite for
        !! positive definite K.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: space, modes, nphi
        class(meridian_fourier_tensor_t), intent(in) :: coef
        type(csc_t), intent(out) :: m
        integer, intent(out) :: status

        integer, allocatable :: rows(:), cols(:), gdof(:)
        real(dp), allocatable :: vals(:), ml(:, :), gl(:, :), phi(:, :)
        real(dp), allocatable :: kp(:, :), xq(:, :), w(:), kt(:, :, :)
        real(dp), allocatable :: phis(:, :, :), tr(:), tt(:, :), h(:, :)
        real(dp), allocatable :: hb(:, :)
        real(dp) :: x(2)
        type(fortsparse_status_t) :: fst
        integer :: nt, nl, nb, ml_n, t, q, j, p, i, b, b2, base, bad, st, np
        integer :: ndof

        status = 1
        if (space < 1 .or. space > 2) return
        if (modes < 0 .or. nphi < 2*modes + 1) return
        nt = fs%mm%topology%n_triangles
        nl = fs%sp(space)%nloc
        ndof = fs%sp(space)%ndof
        nb = 2*modes + 1
        ml_n = nl*nb
        np = fs%nq*nphi
        allocate (rows(ml_n*ml_n*nt), cols(ml_n*ml_n*nt), vals(ml_n*ml_n*nt))
        ! The quadrature angles phi_j are the same at every meridian point,
        ! so the phi sum factorises: H_j = sum_q w_q/nphi phi_q^T K phi_q,
        ! then block (b, b2) = sum_j tr_b(phi_j) tr_b2(phi_j) H_j, one GEMM
        ! with the angle table tt(j, (b, b2)).
        allocate (tt(nphi, nb*nb), tr(nb))
        do j = 1, nphi
            call trig_values(modes, 2*pi*(j - 1)/nphi, tr)
            do b2 = 0, nb - 1
                do b = 0, nb - 1
                    tt(j, b2*nb + b + 1) = tr(b + 1)*tr(b2 + 1)
                end do
            end do
        end do
        bad = 0
        !$omp parallel if(nt >= par_min) default(shared) &
        !$omp private(t, q, j, p, i, b, b2, base, st, ml, gl, phi, kp, xq, &
        !$omp w, kt, phis, gdof, x, h, hb) reduction(max:bad)
        allocate (ml(ml_n, ml_n), gl(nl, nl), phi(3, nl), kp(3, nl), &
            xq(3, np), w(fs%nq), kt(3, 3, np), phis(3, nl, fs%nq), &
            gdof(ml_n), h(nl*nl, nphi), hb(nl*nl, nb*nb))
        !$omp do schedule(dynamic, 4)
        do t = 1, nt
            ml = 0.0_dp
            st = 0
            do q = 1, fs%nq
                call meridian_qpoint(fs, space, t, q, x, w(q), phi, st)
                if (st /= 0) exit
                do i = 1, nl
                    phis(:, i, q) = real(fs%sp(space)%sgn(i, t), dp)*phi(:, i)
                end do
                do j = 1, nphi
                    xq(:, (q - 1)*nphi + j) = [x, 2*pi*(j - 1)/nphi]
                end do
            end do
            if (st == 0) call coef%tensor(t, xq, kt, st)
            if (st /= 0) then
                bad = 1
                cycle
            end if
            h = 0.0_dp
            do q = 1, fs%nq
                do j = 1, nphi
                    p = (q - 1)*nphi + j
                    kp = matmul(kt(:, :, p), phis(:, :, q))
                    gl = matmul(transpose(phis(:, :, q)), kp)
                    h(:, j) = h(:, j) + (w(q)/nphi)*reshape(gl, [nl*nl])
                end do
            end do
            hb = matmul(h, tt)
            do b2 = 0, nb - 1
                do b = 0, nb - 1
                    ml(b*nl + 1:(b + 1)*nl, b2*nl + 1:(b2 + 1)*nl) = &
                        reshape(hb(:, b2*nb + b + 1), [nl, nl])
                end do
            end do
            do b = 0, nb - 1
                gdof(b*nl + 1:(b + 1)*nl) = b*ndof + fs%sp(space)%dof(:, t)
            end do
            base = (t - 1)*ml_n*ml_n
            do j = 1, ml_n
                do i = 1, ml_n
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
        call csc_from_triplet(ndof*nb, ndof*nb, rows, cols, vals, m, fst)
        if (status_ok(fst)) status = 0
    end subroutine meridian_fourier_mass

    pure subroutine trig_values(modes, phi, tr)
        !! tr(b + 1) = 1, cos phi, sin phi, cos 2 phi, ...
        integer, intent(in) :: modes
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: tr(:)
        integer :: n

        tr(1) = 1.0_dp
        do n = 1, modes
            tr(2*n) = cos(n*phi)
            tr(2*n + 1) = sin(n*phi)
        end do
    end subroutine trig_values

    subroutine meridian_fourier_curl_curl(fs, modes, nphi, coef, s, status)
        !! S = C_F^T M_2(K) C_F: < curl v . K curl u > on Fourier-V1.
        !! Symmetric positive semi-definite with grad V0 in the null space.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes, nphi
        class(meridian_fourier_tensor_t), intent(in) :: coef
        type(csc_t), intent(out) :: s
        integer, intent(out) :: status

        type(csc_t) :: c, ct, m2, mc
        type(fortsparse_status_t) :: fst

        call meridian_fourier_curl(fs, modes, c, status)
        if (status /= 0) return
        call meridian_fourier_mass(fs, 2, modes, nphi, coef, m2, status)
        if (status /= 0) return
        status = 2
        call csc_matmul(m2, c, mc, fst)
        if (.not. status_ok(fst)) return
        call csc_transpose(c, ct, fst)
        if (.not. status_ok(fst)) return
        call csc_matmul(ct, mc, s, fst)
        if (status_ok(fst)) status = 0
    end subroutine meridian_fourier_curl_curl

    subroutine meridian_fourier_gradient_coupling(fs, modes, nphi, coef, cpl, &
            status)
        !! cpl = M_1(K) G_F: < v . K grad lambda >, rows Fourier-V1, columns
        !! Fourier-V0.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes, nphi
        class(meridian_fourier_tensor_t), intent(in) :: coef
        type(csc_t), intent(out) :: cpl
        integer, intent(out) :: status

        type(csc_t) :: g, m1
        type(fortsparse_status_t) :: fst

        call meridian_fourier_gradient(fs, modes, g, status)
        if (status /= 0) return
        call meridian_fourier_mass(fs, 1, modes, nphi, coef, m1, status)
        if (status /= 0) return
        status = 2
        call csc_matmul(m1, g, cpl, fst)
        if (status_ok(fst)) status = 0
    end subroutine meridian_fourier_gradient_coupling

end module fortfem_meridian_fourier_curl
