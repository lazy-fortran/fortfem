program test_meridian_hcurl_auxiliary
    !! Hiptmair-Xu auxiliary space on the meridional Nedelec space N_k:
    !!  - rot(E G x) = 0: the embedded Whitney gradients are discrete
    !!    gradients of N_k (checks E, its signs and G on a curved mesh),
    !!  - E Pi reproduces a constant vector field on a straight mesh,
    !!  - HX-preconditioned CG for K_n = rot^T M rot + n^2 M_(1/R) converges
    !!    in a bounded number of iterations on h and p ladders and returns
    !!    the direct solution.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_qpoint
    use fortfem_meridian_hcurl_auxiliary, only: meridian_hcurl_aux_t, &
        meridian_hcurl_aux_build, meridian_hcurl_mode_matrix
    use fortfem_hcurl_auxiliary_space, only: hx_preconditioner_t, hx_options_t, &
        hx_space_t, hx_setup, hx_destroy, hx_pcg, hx_masked_csr, &
        HX_SMOOTHER_PATCH, HX_AUX_AMG
    use fortsparse, only: csc_t, csr_t, csr_matvec, csc_matvec, &
        csc_normalize_symmetric_raw, &
        sparse_solver_t, sparse_factor, sparse_solve, sparse_destroy, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_NATIVE, &
        FORTSPARSE_SYM_SPD
    use meridian_test_meshes, only: make_annulus_mesh
    implicit none

    logical :: all_passed

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    call test_gradient_embedding()
    call test_constant_reproduction()
    call test_pcg_ladder()
    call check_summary("Meridian H(curl) Hiptmair-Xu auxiliary space")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    subroutine build(nr, k, geo, fs, aux)
        integer, intent(in) :: nr, k, geo
        type(meridian_feec_t), intent(out) :: fs
        type(meridian_hcurl_aux_t), intent(out) :: aux
        type(meridian_mesh_t) :: mm
        integer :: st
        call make_annulus_mesh(nr, 4*nr, 0.3_dp, 1.2_dp, 6.0_dp, geo, mm)
        call meridian_feec_build(mm, k, fs, st)
        if (st /= 0) error stop "feec build"
        call meridian_hcurl_aux_build(fs, aux, st)
        if (st /= 0) error stop "aux build"
    end subroutine build

    subroutine test_gradient_embedding()
        type(meridian_feec_t) :: fs
        type(meridian_hcurl_aux_t) :: aux
        real(dp), allocatable :: x(:), a(:), full(:), b(:)
        integer :: k, i
        real(dp) :: rot, scale, floor
        do k = 1, 5
            call build(3, k, 2, fs, aux)
            allocate (x(aux%nv))
            do i = 1, aux%nv
                x(i) = sin(0.7_dp*i) + cos(1.3_dp*i**2)
            end do
            a = csc_matvec(aux%transfer, csc_matvec(aux%gradient, x))
            allocate (full(fs%sp(1)%ndof))
            full = 0.0_dp
            full(:fs%nn) = a
            b = csc_matvec(fs%curl%a, full)
            rot = maxval(abs(b(fs%nn + 1:)))
            scale = maxval(abs(a))
            floor = complex_floor(fs)
            print '(a,i0,a,2es10.2)', "  k = ", k, &
                " rot(E G x)/|E G x|, complex floor rot(G_k y)/|G_k y| ", &
                rot/scale, floor
            call record(rot <= max(1.0e-12_dp, 10*floor)*scale .and. &
                scale > 0.1_dp, "rot(E G x) = 0 to the complex rounding floor")
            deallocate (x, full)
        end do
    end subroutine test_gradient_embedding

    real(dp) function complex_floor(fs) result(floor)
        !! Rounding floor of the N_k complex itself: rot of FortFEM's own
        !! discrete gradient of a V0 field, relative to its size.  It grows
        !! with k (monomial-moment reference basis).
        type(meridian_feec_t), intent(in) :: fs
        real(dp), allocatable :: y(:), ga(:), b(:)
        integer :: i
        allocate (y(fs%n0))
        do i = 1, fs%n0
            y(i) = sin(0.7_dp*i) + cos(1.3_dp*i**2)
        end do
        ga = csc_matvec(fs%grad%a, y)
        b = csc_matvec(fs%curl%a, ga)
        floor = maxval(abs(b(fs%nn + 1:)))/maxval(abs(ga(:fs%nn)))
    end function complex_floor

    subroutine test_constant_reproduction()
        type(meridian_feec_t) :: fs
        type(meridian_hcurl_aux_t) :: aux
        real(dp), allocatable :: u(:), phi(:, :), cv(:)
        real(dp) :: x(2), w, val(2), err
        real(dp), parameter :: c(2) = [0.7_dp, -1.9_dp]
        integer :: k, t, q, l, st
        do k = 1, 5
            call build(3, k, 1, fs, aux)
            allocate (cv(aux%nv))
            cv = c(1)
            u = csc_matvec(aux%interp(1), cv)
            cv = c(2)
            u = u + csc_matvec(aux%interp(2), cv)
            u = csc_matvec(aux%transfer, u)
            allocate (phi(3, fs%sp(1)%nloc))
            err = 0.0_dp
            do t = 1, fs%mm%topology%n_triangles
                do q = 1, fs%nq
                    call meridian_qpoint(fs, 1, t, q, x, w, phi, st)
                    val = 0.0_dp
                    do l = 1, fs%nln
                        val = val + fs%sp(1)%sgn(l, t)*u(fs%sp(1)%dof(l, t)) &
                            *[phi(1, l), phi(3, l)]
                    end do
                    err = max(err, maxval(abs(val - c)))
                end do
            end do
            print '(a,i0,a,2es10.2)', "  k = ", k, &
                " max|E Pi c - c|, complex floor ", err, complex_floor(fs)
            call record(err < max(1.0e-12_dp, 10*complex_floor(fs))*norm2(c), &
                "E Pi reproduces constant fields to the complex rounding floor")
            deallocate (cv, phi)
        end do
    end subroutine test_constant_reproduction

    subroutine test_pcg_ladder()
        !! Stop at max(1e-10, 10 r_d), r_d the relative residual of the
        !! direct solve (the conditioning floor of K_n grows with k).
        integer :: k, nr, its(3), j
        real(dp) :: dif, floor
        do k = 1, 4
            j = 0
            do nr = 2, 8, 2
                if (nr == 6) cycle
                j = j + 1
                call solve_case(nr, k, 1, its(j), dif, floor)
                call record(dif < max(1.0e-8_dp, 1.0e3_dp*floor), &
                    "HX-PCG returns the direct solution")
            end do
            print '(a,i0,a,3(1x,i0))', "  k = ", k, " HX-PCG iterations nr 2/4/8:", &
                its
            call record(maxval(its) <= 40, "bounded HX-PCG iterations")
            call record(its(3) <= its(1) + 8, "iterations nearly h-independent")
        end do
    end subroutine test_pcg_ladder

    subroutine solve_case(nr, k, n, its, dif, floor)
        integer, intent(in) :: nr, k, n
        integer, intent(out) :: its
        real(dp), intent(out) :: dif, floor
        type(meridian_feec_t) :: fs
        type(meridian_hcurl_aux_t) :: aux
        type(csc_t) :: kmat, m
        type(hx_preconditioner_t) :: pc
        type(hx_options_t) :: opt
        type(hx_space_t) :: spaces(3)
        type(csr_t) :: op
        type(sparse_solver_t) :: direct
        type(fortsparse_status_t) :: fst
        logical, allocatable :: fixed(:)
        integer, allocatable :: cp(:), ri(:)
        real(dp), allocatable :: b(:), x(:), xd(:), v(:), rd(:)
        real(dp) :: res, corr, tol
        integer :: st, info, i

        call build(nr, k, 2, fs, aux)
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
        ! direct oracle on the same masked operator (CSR of a symmetric
        ! matrix is its CSC)
        cp = op%row_ptr
        ri = op%col
        v = op%val
        call csc_normalize_symmetric_raw(fs%nn, cp, ri, v, 1.0e-12_dp, corr, fst)
        m%nrow = fs%nn; m%ncol = fs%nn; m%nnz = size(v)
        m%col_ptr = cp; m%row_idx = ri; m%val = v
        direct%backend_id = FORTSPARSE_BACKEND_NATIVE
        direct%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(direct, m, fst)
        if (.not. status_ok(fst)) error stop "direct factor"
        call sparse_solve(direct, b, xd, fst)
        call sparse_destroy(direct)
        call csr_matvec(op, xd, rd)
        floor = norm2(b - rd)/norm2(b)
        tol = max(1.0e-10_dp, 10*floor)

        spaces(1)%p = aux%gradient
        spaces(2)%p = aux%interp(1)
        spaces(3)%p = aux%interp(2)
        opt%smoother = HX_SMOOTHER_PATCH
        opt%aux_solver = HX_AUX_AMG
        if (k > 1) then
            call hx_setup(pc, kmat, fixed, spaces, st, opt, aux%transfer, &
                aux%edge_fixed, aux%patch_ptr, aux%patch_dof)
        else
            call hx_setup(pc, kmat, fixed, spaces, st, opt, &
                patch_ptr=aux%patch_ptr, patch_dof=aux%patch_dof)
        end if
        if (st /= 0) error stop "hx setup"
        x = 0.0_dp
        call hx_pcg(op, pc, b, x, tol, 500, its, res, info)
        if (info /= 0) its = 999
        dif = norm2(x - xd)/norm2(xd)
        print '(a,2i3,a,i4,3es10.2)', "  k nr", k, nr, &
            " its, |x - x_d|/|x_d|, direct floor, CG residual", its, dif, floor, res
        call hx_destroy(pc)
    end subroutine solve_case

end program test_meridian_hcurl_auxiliary
