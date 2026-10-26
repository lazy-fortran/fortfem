program fortfem_hcurl_hx_ladder
    !! Iteration and wall-time ladder of the Hiptmair-Xu preconditioner for
    !! the per-mode meridional curl-curl system K_n = rot^T M rot + n^2
    !! M_(1/R) on N_k (gauge R A_phi = 0), against the native sparse
    !! Cholesky.  Geometry: the NN1 vacuum annulus of KIN6D, minor radii
    !! 0.3 - 1.2 m about R0 = 6 m, P2 isoparametric, nr radial and 4 nr
    !! angular cells.
    !!
    !! usage: fortfem_hcurl_hx_ladder K NR N SMOOTHER AUX [nodirect]
    !!   SMOOTHER: patch | cheb, AUX: amg | direct; nodirect skips the
    !!   Cholesky reference (its columns are then 0); optional 7th argument:
    !!   CG tolerance (default max(1e-10, 10 x direct floor)).
    !! prints one CSV row:
    !!   k,nr,n,smoother,aux,threads,ndof,setup_s,solve_s,its,factor_s,
    !!   direct_solve_s,floor,cg_res,rel_diff,patch_mb
    !$ use omp_lib, only: omp_get_max_threads
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_from_map
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build
    use fortfem_meridian_hcurl_auxiliary, only: meridian_hcurl_aux_t, &
        meridian_hcurl_aux_build, meridian_hcurl_mode_matrix
    use fortfem_hcurl_auxiliary_space, only: hx_preconditioner_t, hx_options_t, &
        hx_space_t, hx_setup, hx_destroy, hx_pcg, hx_masked_csr, &
        HX_SMOOTHER_PATCH, HX_SMOOTHER_CHEBYSHEV, HX_AUX_AMG, HX_AUX_DIRECT
    use fortsparse, only: csc_t, csr_t, csr_matvec, csc_normalize_symmetric_raw, &
        sparse_solver_t, sparse_factor, sparse_solve, sparse_destroy, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_NATIVE, &
        FORTSPARSE_SYM_SPD
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 6.0_dp
    character(32) :: arg, smoother_name, aux_name, extra
    integer :: k, nr, n, st, info, its, i, threads
    type(meridian_mesh_t) :: mm
    type(meridian_feec_t) :: fs
    type(meridian_hcurl_aux_t) :: aux
    type(csc_t) :: kmat, m
    type(csr_t) :: op
    type(hx_preconditioner_t) :: pc
    type(hx_options_t) :: opt
    type(hx_space_t) :: spaces(3)
    type(sparse_solver_t) :: direct
    type(fortsparse_status_t) :: fst
    logical, allocatable :: fixed(:)
    integer, allocatable :: cp(:), ri(:)
    real(dp), allocatable :: b(:), x(:), xd(:), v(:), rd(:)
    real(dp) :: t0, t_setup, t_solve, t_factor, t_dsolve, floor, res, corr, tol

    call get_command_argument(1, arg); read (arg, *) k
    call get_command_argument(2, arg); read (arg, *) nr
    call get_command_argument(3, arg); read (arg, *) n
    call get_command_argument(4, smoother_name)
    call get_command_argument(5, aux_name)
    call get_command_argument(6, extra)
    threads = 1
    !$ threads = omp_get_max_threads()

    call annulus(nr, mm)
    call meridian_feec_build(mm, k, fs, st)
    if (st /= 0) error stop "feec build"
    call meridian_hcurl_mode_matrix(fs, n, kmat, st)
    if (st /= 0) error stop "mode matrix"
    fixed = fs%sp(1)%bnd(:fs%nn)
    call hx_masked_csr(kmat, fixed, op, st)
    allocate (b(fs%nn), x(fs%nn), xd(fs%nn), rd(fs%nn))
    do i = 1, fs%nn
        b(i) = sin(0.37_dp*i) + 0.2_dp*cos(2.1_dp*i)
    end do
    where (fixed) b = 0.0_dp
    b = b/norm2(b)

    t_factor = 0.0_dp; t_dsolve = 0.0_dp; floor = 0.0_dp; xd = 0.0_dp
    if (trim(extra) /= 'nodirect') call direct_reference()
    tol = max(1.0e-10_dp, 10*floor)
    if (command_argument_count() >= 7) then
        call get_command_argument(7, arg)
        read (arg, *) tol
    end if

    t0 = wtime()
    call meridian_hcurl_aux_build(fs, aux, st)
    if (st /= 0) error stop "aux build"
    spaces(1)%p = aux%gradient
    spaces(2)%p = aux%interp(1)
    spaces(3)%p = aux%interp(2)
    opt%smoother = HX_SMOOTHER_PATCH
    if (trim(smoother_name) == 'cheb') opt%smoother = HX_SMOOTHER_CHEBYSHEV
    opt%aux_solver = HX_AUX_AMG
    if (trim(aux_name) == 'direct') opt%aux_solver = HX_AUX_DIRECT
    if (k > 1) then
        call hx_setup(pc, kmat, fixed, spaces, st, opt, aux%transfer, &
            aux%edge_fixed, aux%patch_ptr, aux%patch_dof)
    else
        call hx_setup(pc, kmat, fixed, spaces, st, opt, &
            patch_ptr=aux%patch_ptr, patch_dof=aux%patch_dof)
    end if
    if (st /= 0) error stop "hx setup"
    t_setup = wtime() - t0
    x = 0.0_dp
    t0 = wtime()
    call hx_pcg(op, pc, b, x, tol, 2000, its, res, info)
    t_solve = wtime() - t0
    if (info /= 0) its = -its
    write (*, '(i0,",",i0,",",i0,",",a,",",a,",",i0,",",i0,9(",",es11.4))') &
        k, nr, n, trim(smoother_name), trim(aux_name), threads, fs%nn, &
        t_setup, t_solve, real(its, dp), t_factor, t_dsolve, floor, res, &
        norm2(x - xd)/max(norm2(xd), tiny(1.0_dp)), pc%patch_bytes/1.0e6_dp
    call hx_destroy(pc)

contains

    subroutine direct_reference()
        !! Native supernodal Cholesky of the same masked operator.
        cp = op%row_ptr
        ri = op%col
        v = op%val
        call csc_normalize_symmetric_raw(fs%nn, cp, ri, v, 1.0e-12_dp, corr, fst)
        m%nrow = fs%nn; m%ncol = fs%nn; m%nnz = size(v)
        m%col_ptr = cp; m%row_idx = ri; m%val = v
        direct%backend_id = FORTSPARSE_BACKEND_NATIVE
        direct%options%symmetry = FORTSPARSE_SYM_SPD
        t0 = wtime()
        call sparse_factor(direct, m, fst)
        if (.not. status_ok(fst)) error stop "direct factor"
        t_factor = wtime() - t0
        t0 = wtime()
        call sparse_solve(direct, b, xd, fst)
        t_dsolve = wtime() - t0
        call sparse_destroy(direct)
        call csr_matvec(op, xd, rd)
        floor = norm2(b - rd)/norm2(b)
    end subroutine direct_reference

    real(dp) function wtime()
        integer(8) :: c, r
        call system_clock(c, r)
        wtime = real(c, dp)/real(r, dp)
    end function wtime

    subroutine annulus_map(p, xx)
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: xx(2)
        xx = [r0 + p(1)*cos(p(2)), p(1)*sin(p(2))]
    end subroutine annulus_map

    subroutine annulus(nr, mm)
        integer, intent(in) :: nr
        type(meridian_mesh_t), intent(out) :: mm
        type(mesh_2d_t) :: topo
        real(dp), allocatable :: pv(:, :, :)
        real(dp) :: dr, dth, xx(2)
        real(dp), parameter :: r1 = 0.3_dp, r2 = 1.2_dp
        integer :: i, j, jn, kk, nth, status
        nth = 4*nr
        dr = (r2 - r1)/real(nr, dp)
        dth = 2.0_dp*pi/real(nth, dp)
        topo%n_vertices = (nr + 1)*nth
        topo%n_triangles = 2*nr*nth
        topo%has_triangles = .true.
        allocate (topo%vertices(2, topo%n_vertices))
        allocate (topo%triangles(3, topo%n_triangles))
        allocate (pv(2, 3, topo%n_triangles))
        do i = 0, nr
            do j = 0, nth - 1
                call annulus_map([r1 + i*dr, j*dth], xx)
                topo%vertices(:, i*nth + j + 1) = xx
            end do
        end do
        kk = 0
        do i = 0, nr - 1
            do j = 0, nth - 1
                jn = mod(j + 1, nth)
                kk = kk + 1
                topo%triangles(:, kk) = [i*nth + j + 1, (i + 1)*nth + j + 1, &
                    (i + 1)*nth + jn + 1]
                pv(:, :, kk) = reshape([r1 + i*dr, j*dth, r1 + (i + 1)*dr, &
                    j*dth, r1 + (i + 1)*dr, (j + 1)*dth], [2, 3])
                kk = kk + 1
                topo%triangles(:, kk) = [i*nth + j + 1, (i + 1)*nth + jn + 1, &
                    i*nth + jn + 1]
                pv(:, :, kk) = reshape([r1 + i*dr, j*dth, r1 + (i + 1)*dr, &
                    (j + 1)*dth, r1 + i*dr, (j + 1)*dth], [2, 3])
            end do
        end do
        call meridian_mesh_from_map(topo, pv, 2, mm, status, annulus_map)
        if (status /= 0) error stop "annulus mesh failed"
    end subroutine annulus

end program fortfem_hcurl_hx_ladder
