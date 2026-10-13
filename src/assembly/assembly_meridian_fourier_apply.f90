module fortfem_meridian_fourier_apply
    !! Matrix-free apply of the mode-coupled real-Fourier mass form of
    !! fortfem_meridian_fourier_curl, and of the curl-curl form built on it,
    !! without assembling the (2N+1)^2 mode blocks per element.
    !!
    !! M(v, u) = < v . K u > is applied element by element: the field of u at
    !! each meridional quadrature point is synthesised on the uniform phi grid
    !! phi_j = 2 pi (j - 1)/nphi by an inverse FFT of its 2N + 1 real Fourier
    !! amplitudes, K(R, Z, phi_j) is applied pointwise, and the result is
    !! projected back to the amplitudes by a forward FFT (fortnum_fft), then
    !! contracted with the basis.  This is the same quadrature as the
    !! assembled form, so the apply equals the matrix product to roundoff.
    !! K is evaluated once at setup (6 values per phi point, quadrature point
    !! and element); storage is O(nt nq nphi) plus O(nt nl N) per apply,
    !! against O(nt (nl N)^2) for the assembled matrix, and the apply costs
    !! O(nt nq nphi log nphi) instead of O(nt (nl N)^2).  Element loops are
    !! OpenMP parallel; the scatter is race free (dof-ordered gather).
    use fortfem_kinds, only: dp
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_qpoint
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t
    use fortnum_fft, only: fft_r2c, fft_c2c, fortnum_fft_plan_t, fft_plan_init, &
        fft_c2c_plan_init
    use fortsparse, only: csc_t, csc_matvec
    !$ use omp_lib, only: omp_get_max_threads
    implicit none
    private

    type, public :: meridian_fourier_op_t
        !! Setup data of one mass operator (space 1 or 2, modes, nphi, K).
        integer :: space = 0, modes = 0, nphi = 0
        integer :: nt = 0, nq = 0, nl = 0, ndof = 0
        integer, allocatable :: dof(:, :)       ! (nl, nt) global dof of the space
        real(dp), allocatable :: phis(:, :, :, :) ! (3, nl, nq, nt) signed basis
        real(dp), allocatable :: wq(:, :)       ! (nq, nt) weight / nphi
        real(dp), allocatable :: kk(:, :, :, :) ! (6, nphi, nq, nt): xx yy zz xy xz yz
        integer, allocatable :: dptr(:), dlist(:) ! dof -> flat element-local slots
        type(fortnum_fft_plan_t) :: rplan, cplan
    end type meridian_fourier_op_t

    public :: meridian_fourier_op_init, meridian_fourier_mass_apply
    public :: meridian_fourier_curl_curl_apply, meridian_fourier_op_bytes

    integer, parameter :: par_min = 32
    real(dp), parameter :: pi = acos(-1.0_dp)

contains

    subroutine meridian_fourier_op_init(fs, space, modes, nphi, coef, op, status)
        !! Build the apply data for M(K) on Fourier-V1 (space = 1) or V2
        !! (space = 2); same arguments and quadrature as meridian_fourier_mass.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: space, modes, nphi
        class(meridian_fourier_tensor_t), intent(in) :: coef
        type(meridian_fourier_op_t), intent(out) :: op
        integer, intent(out) :: status

        real(dp), allocatable :: phi(:, :), xq(:, :), kt(:, :, :), w(:)
        real(dp) :: x(2)
        integer, allocatable :: cnt(:)
        integer :: nt, nl, nq, t, q, j, i, st, bad, np, d, k

        status = 1
        if (space < 1 .or. space > 2) return
        if (modes < 0 .or. nphi < 2*modes + 1) return
        nt = fs%mm%topology%n_triangles
        nl = fs%sp(space)%nloc
        nq = fs%nq
        np = nq*nphi
        op%space = space; op%modes = modes; op%nphi = nphi
        op%nt = nt; op%nq = nq; op%nl = nl; op%ndof = fs%sp(space)%ndof
        allocate (op%dof(nl, nt), op%phis(3, nl, nq, nt), op%wq(nq, nt), &
            op%kk(6, nphi, nq, nt))
        op%dof = fs%sp(space)%dof(:, 1:nt)
        bad = 0
        !$omp parallel if(nt >= par_min) default(shared) &
        !$omp private(t, q, j, i, st, phi, xq, kt, w, x) reduction(max:bad)
        allocate (phi(3, nl), xq(3, np), kt(3, 3, np), w(nq))
        !$omp do schedule(dynamic, 4)
        do t = 1, nt
            st = 0
            do q = 1, nq
                call meridian_qpoint(fs, space, t, q, x, w(q), phi, st)
                if (st /= 0) exit
                do i = 1, nl
                    op%phis(:, i, q, t) = real(fs%sp(space)%sgn(i, t), dp)*phi(:, i)
                end do
                op%wq(q, t) = w(q)/nphi
                do j = 1, nphi
                    xq(:, (q - 1)*nphi + j) = [x, 2*pi*(j - 1)/nphi]
                end do
            end do
            if (st == 0) call coef%tensor(t, xq, kt, st)
            if (st /= 0) then
                bad = 1
                cycle
            end if
            do q = 1, nq
                do j = 1, nphi
                    op%kk(:, j, q, t) = [kt(1, 1, (q - 1)*nphi + j), &
                        kt(2, 2, (q - 1)*nphi + j), kt(3, 3, (q - 1)*nphi + j), &
                        kt(1, 2, (q - 1)*nphi + j), kt(1, 3, (q - 1)*nphi + j), &
                        kt(2, 3, (q - 1)*nphi + j)]
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        if (bad /= 0) return
        ! dof -> element-local slots (CSR), for a race-free scatter
        allocate (cnt(op%ndof + 1), op%dptr(op%ndof + 1), op%dlist(nl*nt))
        cnt = 0
        do t = 1, nt
            do i = 1, nl
                cnt(op%dof(i, t)) = cnt(op%dof(i, t)) + 1
            end do
        end do
        op%dptr(1) = 1
        do d = 1, op%ndof
            op%dptr(d + 1) = op%dptr(d) + cnt(d)
        end do
        cnt = 0
        do t = 1, nt
            do i = 1, nl
                d = op%dof(i, t)
                k = op%dptr(d) + cnt(d)
                cnt(d) = cnt(d) + 1
                op%dlist(k) = (t - 1)*nl + i
            end do
        end do
        call fft_plan_init(op%rplan, nphi)
        call fft_c2c_plan_init(op%cplan, nphi)
        status = 0
    end subroutine meridian_fourier_op_init

    function meridian_fourier_op_bytes(op) result(nbytes)
        !! Persistent storage of the setup data (excluding FFT plans).
        type(meridian_fourier_op_t), intent(in) :: op
        integer(8) :: nbytes

        nbytes = 8_8*(size(op%phis, kind=8) + size(op%wq, kind=8) &
            + size(op%kk, kind=8)) + 4_8*(size(op%dof, kind=8) &
            + size(op%dptr, kind=8) + size(op%dlist, kind=8))
    end function meridian_fourier_op_bytes

    subroutine meridian_fourier_mass_apply(op, u, y)
        !! y = M(K) u for u of size ndof*(2 modes + 1); equals the product with
        !! the matrix of meridian_fourier_mass to roundoff.
        type(meridian_fourier_op_t), intent(in) :: op
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: y(:)

        real(dp), allocatable :: ye(:, :)
        real(dp), allocatable :: ul(:, :), yl(:, :), pm(:, :), r(:, :), v(:, :)
        real(dp), allocatable :: rr(:)
        complex(dp), allocatable :: z(:), c(:)
        integer :: nb, nl, nq, nphi, modes, t, q, i, b, j, n, comp, d, k, ndof
        real(dp) :: kx(6), vv(3), ww(3)

        modes = op%modes; nphi = op%nphi; nb = 2*modes + 1
        nl = op%nl; nq = op%nq; ndof = op%ndof
        allocate (ye(nl*nb, op%nt))
        !$omp parallel if(op%nt >= par_min) default(shared) &
        !$omp private(t, q, i, b, j, n, comp, ul, yl, pm, r, v, rr, z, c, kx, vv, ww)
        allocate (ul(nl, nb), yl(nl, nb), pm(3, nb), r(3, nb), v(3, nphi), &
            z(nphi), c(nphi/2 + 1), rr(nphi))
        !$omp do schedule(dynamic, 8)
        do t = 1, op%nt
            do b = 0, nb - 1
                do i = 1, nl
                    ul(i, b + 1) = u(b*ndof + op%dof(i, t))
                end do
            end do
            yl = 0.0_dp
            do q = 1, nq
                pm = matmul(op%phis(:, :, q, t), ul)
                ! synthesis: amplitudes -> values on the phi grid
                do comp = 1, 3
                    z = (0.0_dp, 0.0_dp)
                    z(1) = cmplx(pm(comp, 1), 0.0_dp, dp)
                    do n = 1, modes
                        z(n + 1) = cmplx(pm(comp, 2*n), -pm(comp, 2*n + 1), dp)
                    end do
                    call fft_c2c(z, 1, op%cplan)
                    v(comp, :) = real(z, dp)
                end do
                ! pointwise tensor
                do j = 1, nphi
                    kx = op%kk(:, j, q, t)
                    vv = v(:, j)
                    ww(1) = kx(1)*vv(1) + kx(4)*vv(2) + kx(5)*vv(3)
                    ww(2) = kx(4)*vv(1) + kx(2)*vv(2) + kx(6)*vv(3)
                    ww(3) = kx(5)*vv(1) + kx(6)*vv(2) + kx(3)*vv(3)
                    v(:, j) = ww
                end do
                ! analysis: grid values -> amplitudes (cos, sin parts)
                do comp = 1, 3
                    rr = v(comp, :)
                    call fft_r2c(rr, c, op%rplan)
                    r(comp, 1) = real(c(1), dp)
                    do n = 1, modes
                        r(comp, 2*n) = real(c(n + 1), dp)
                        r(comp, 2*n + 1) = -aimag(c(n + 1))
                    end do
                end do
                yl = yl + op%wq(q, t)*matmul(transpose(op%phis(:, :, q, t)), r)
            end do
            do b = 0, nb - 1
                ye(b*nl + 1:(b + 1)*nl, t) = yl(:, b + 1)
            end do
        end do
        !$omp end do
        !$omp end parallel
        ! race-free scatter by dof
        !$omp parallel do if(ndof >= 4*par_min) default(shared) &
        !$omp private(b, d, k, i, t) schedule(static)
        do d = 1, ndof
            do b = 0, nb - 1
                y(b*ndof + d) = 0.0_dp
                do k = op%dptr(d), op%dptr(d + 1) - 1
                    i = mod(op%dlist(k) - 1, nl) + 1
                    t = (op%dlist(k) - 1)/nl + 1
                    y(b*ndof + d) = y(b*ndof + d) + ye(b*nl + i, t)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine meridian_fourier_mass_apply

    subroutine meridian_fourier_curl_curl_apply(op, c, ct, u, y)
        !! y = C^T M_2(K) C u with C, C^T the Fourier curl and its transpose
        !! (meridian_fourier_curl and csc_transpose of it); op is a space-2
        !! operator.  Equals the matrix of meridian_fourier_curl_curl applied
        !! to u to roundoff.
        type(meridian_fourier_op_t), intent(in) :: op
        type(csc_t), intent(in) :: c, ct
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: y(:)
        real(dp), allocatable :: cu(:), mcu(:)

        allocate (cu(c%nrow), mcu(c%nrow))
        cu = csc_matvec(c, u)
        call meridian_fourier_mass_apply(op, cu, mcu)
        y = csc_matvec(ct, mcu)
    end subroutine meridian_fourier_curl_curl_apply

end module fortfem_meridian_fourier_apply
