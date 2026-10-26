module fortfem_hcurl_auxiliary_space
    !! Hiptmair-Xu auxiliary-space preconditioner for symmetric positive
    !! definite H(curl) systems A = curl^T M_a curl + M_b (Hiptmair and Xu,
    !! SIAM J. Numer. Anal. 45 (2007) 2483, doi:10.1137/060660588).
    !!
    !! Purely algebraic: the caller supplies the fine operator, its fixed
    !! (Dirichlet) dofs and the prolongations of the auxiliary spaces into
    !! the lowest-order edge space N_1, e.g. the discrete gradient G of the
    !! nodal P1 space and the componentwise nodal interpolation Pi_d of
    !! (P1)^dim.  For a high-order fine space N_k an optional transfer E:
    !! N_1 -> N_k (the exact embedding) adds a p-level on top:
    !!
    !!   level k : smoother S_k (vertex-patch Schwarz or Chebyshev)
    !!   level 1 : A_1 = E^T A E, smoother S_1 and the auxiliary corrections
    !!             P_j (P_j^T A_1 P_j)^{-1} P_j^T, j = 1..m, visited in the
    !!             palindromic order 1..m..1 (multiplicative), or an exact
    !!             factorisation of A_1.
    !!
    !! Every stage is symmetric, so the cycle is an SPD preconditioner for
    !! CG.  Auxiliary problems are scalar nodal problems solved by FortSparse
    !! smoothed aggregation (one V-cycle) or by a native sparse factor.
    !! The vertex-patch smoother (overlapping additive Schwarz over the dofs
    !! supported in each vertex star) is the robust choice in the polynomial
    !! degree (Arnold, Falk and Winther, Numer. Math. 85 (2000) 197; Schoeberl,
    !! Melenk, Pechstein and Zaglmayr, IMA J. Numer. Anal. 28 (2008) 1).
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_t, csr_t, csr_matvec, csc_matmul, csc_transpose, &
        csc_from_triplet, csr_smoother_t, smoother_setup, smoother_apply, &
        SMOOTHER_CHEBYSHEV, SMOOTHER_L1_JACOBI, amg_params_t, amg_handle_t, &
        amg_setup_csr, amg_apply, amg_destroy, sparse_solver_t, sparse_factor, &
        sparse_solve, sparse_destroy, fortsparse_status_t, status_ok, &
        FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_SYM_SPD
    implicit none
    private

    integer, parameter, public :: HX_SMOOTHER_PATCH = 1
    integer, parameter, public :: HX_SMOOTHER_CHEBYSHEV = 2
    integer, parameter, public :: HX_SMOOTHER_L1_JACOBI = 3
    integer, parameter, public :: HX_AUX_AMG = 1, HX_AUX_DIRECT = 2
    integer, parameter, public :: HX_COARSE_HX = 1, HX_COARSE_DIRECT = 2

    type, public :: hx_options_t
        integer :: smoother = HX_SMOOTHER_PATCH      !! fine-level smoother
        integer :: smoothing_steps = 1               !! per pre/post smoothing
        integer :: chebyshev_degree = 2
        integer :: aux_solver = HX_AUX_AMG
        integer :: coarse = HX_COARSE_HX
    end type hx_options_t

    !> Prolongation of one auxiliary space into the lowest-order space.
    type, public :: hx_space_t
        type(csc_t) :: p
    end type hx_space_t

    type :: patch_smoother_t
        integer :: npatch = 0
        integer, allocatable :: ptr(:), dof(:)
        integer(8), allocatable :: fptr(:)
        real(dp), allocatable :: factor(:), buf(:)
        integer, allocatable :: gptr(:), gidx(:)
        real(dp) :: omega = 1.0_dp
    end type patch_smoother_t

    type :: aux_solver_t
        integer :: kind = 0
        type(csr_t) :: a
        type(amg_handle_t) :: amg
        type(sparse_solver_t) :: direct
    end type aux_solver_t

    type :: aux_level_t
        type(csr_t) :: p, pt
        type(aux_solver_t) :: solver
        real(dp), allocatable :: r(:), y(:)
    end type aux_level_t

    type, public :: hx_preconditioner_t
        logical :: ready = .false.
        logical :: two_level = .false.
        type(hx_options_t) :: opt
        integer :: n = 0, n1 = 0
        logical, allocatable :: fixed(:)
        type(csr_t) :: a, e, et
        type(csr_smoother_t) :: cheb
        type(patch_smoother_t) :: patch
        type(csr_t) :: a1
        type(csr_smoother_t) :: cheb1
        type(aux_level_t), allocatable :: aux(:)
        type(aux_solver_t) :: coarse
        real(dp), allocatable :: res(:), tmp(:), rc(:), yc(:), res1(:), tmp1(:)
        integer :: applications = 0
        real(dp) :: patch_bytes = 0.0_dp
    end type hx_preconditioner_t

    public :: hx_setup, hx_destroy, hx_apply, hx_masked_csr, hx_csc_to_csr, hx_pcg

    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf
        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dpotrs
    end interface

contains

    subroutine hx_setup(pc, a, fixed, spaces, status, options, transfer, &
            coarse_fixed, patch_ptr, patch_dof)
        !! a: symmetric positive definite on the free dofs (CSC, n x n);
        !! fixed(n): Dirichlet dofs (identity rows and columns are imposed).
        !! spaces(j)%p: auxiliary prolongations into N_1 (n1 x m_j).  Without
        !! transfer the fine space is N_1 (n1 = n).  transfer: E (n x n1) with
        !! coarse_fixed(n1).  patch_ptr/patch_dof: overlapping fine-level
        !! patches (required for HX_SMOOTHER_PATCH on the fine level).
        type(hx_preconditioner_t), intent(out) :: pc
        type(csc_t), intent(in) :: a
        logical, intent(in) :: fixed(:)
        type(hx_space_t), intent(in) :: spaces(:)
        integer, intent(out) :: status
        type(hx_options_t), intent(in), optional :: options
        type(csc_t), intent(in), optional :: transfer
        logical, intent(in), optional :: coarse_fixed(:)
        integer, intent(in), optional :: patch_ptr(:), patch_dof(:)

        type(csc_t) :: am, em, ae, a1c, pm, ap, apc
        logical, allocatable :: cfix(:), nofix(:)
        type(fortsparse_status_t) :: fst
        integer :: j, n

        status = 1
        if (present(options)) pc%opt = options
        n = size(fixed)
        if (a%nrow /= n .or. a%ncol /= n) return
        if (size(spaces) < 1) return
        pc%n = n
        pc%fixed = fixed
        pc%two_level = present(transfer)
        call mask_csc(a, fixed, fixed, .true., am, status)
        if (status /= 0) return
        call hx_csc_to_csr(am, pc%a, status)
        if (status /= 0) return

        if (pc%two_level) then
            if (.not. present(coarse_fixed)) return
            if (transfer%nrow /= n) return
            pc%n1 = transfer%ncol
            if (size(coarse_fixed) /= pc%n1) return
            cfix = coarse_fixed
            ! E with fixed fine rows and fixed coarse columns removed.
            call mask_csc(transfer, fixed, cfix, .false., em, status)
            if (status /= 0) return
            call hx_csc_to_csr(em, pc%e, status)
            if (status /= 0) return
            call transpose_csr(em, pc%et, status)
            if (status /= 0) return
            call csc_matmul(am, em, ae, fst)
            if (.not. status_ok(fst)) return
            call csc_transpose(em, a1c, fst)
            if (.not. status_ok(fst)) return
            call csc_matmul(a1c, ae, ap, fst)
            if (.not. status_ok(fst)) return
            call symmetric_regularised(ap, a1c, status)
            if (status /= 0) return
            select case (pc%opt%smoother)
            case (HX_SMOOTHER_PATCH)
                if (.not. (present(patch_ptr) .and. present(patch_dof))) then
                    status = 1
                    return
                end if
                call patch_setup(pc%a, fixed, patch_ptr, patch_dof, &
                    pc%patch, pc%patch_bytes, status)
            case (HX_SMOOTHER_CHEBYSHEV)
                call smoother_setup(pc%a, SMOOTHER_CHEBYSHEV, pc%cheb, fst, &
                    degree=pc%opt%chebyshev_degree)
                status = merge(0, 1, status_ok(fst))
            case (HX_SMOOTHER_L1_JACOBI)
                call smoother_setup(pc%a, SMOOTHER_L1_JACOBI, pc%cheb, fst, &
                    degree=pc%opt%chebyshev_degree)
                status = merge(0, 1, status_ok(fst))
            case default
                status = 1
            end select
            if (status /= 0) return
        else
            pc%n1 = n
            cfix = fixed
            a1c = am
            if (pc%opt%smoother == HX_SMOOTHER_PATCH .and. present(patch_ptr) &
                .and. present(patch_dof)) then
                ! Lowest order: the patch smoother acts on N_1 directly.
                call patch_setup(pc%a, fixed, patch_ptr, patch_dof, &
                    pc%patch, pc%patch_bytes, status)
                if (status /= 0) return
            end if
        end if
        call hx_csc_to_csr(a1c, pc%a1, status)
        if (status /= 0) return

        if (pc%opt%coarse == HX_COARSE_DIRECT) then
            call aux_solver_setup(a1c, HX_AUX_DIRECT, pc%coarse, status)
            if (status /= 0) return
        else
            if (.not. (pc%opt%smoother == HX_SMOOTHER_PATCH .and. &
                .not. pc%two_level .and. pc%patch%npatch > 0)) then
                call smoother_setup(pc%a1, SMOOTHER_CHEBYSHEV, pc%cheb1, fst, &
                    degree=pc%opt%chebyshev_degree)
                if (.not. status_ok(fst)) then
                    status = 1
                    return
                end if
            end if
            allocate (pc%aux(size(spaces)))
            do j = 1, size(spaces)
                if (spaces(j)%p%nrow /= pc%n1) then
                    status = 1
                    return
                end if
                if (allocated(nofix)) deallocate (nofix)
                allocate (nofix(spaces(j)%p%ncol))
                nofix = .false.
                call mask_csc(spaces(j)%p, cfix, nofix, .false., pm, status)
                if (status /= 0) return
                call hx_csc_to_csr(pm, pc%aux(j)%p, status)
                if (status /= 0) return
                call transpose_csr(pm, pc%aux(j)%pt, status)
                if (status /= 0) return
                call csc_matmul(a1c, pm, ae, fst)
                if (.not. status_ok(fst)) return
                call csc_transpose(pm, em, fst)
                if (.not. status_ok(fst)) return
                call csc_matmul(em, ae, ap, fst)
                if (.not. status_ok(fst)) return
                call symmetric_regularised(ap, apc, status)
                if (status /= 0) return
                call aux_solver_setup(apc, pc%opt%aux_solver, pc%aux(j)%solver, &
                    status)
                if (status /= 0) return
                allocate (pc%aux(j)%r(pm%ncol), pc%aux(j)%y(pm%ncol))
            end do
        end if
        allocate (pc%res(n), pc%tmp(n), pc%rc(pc%n1), pc%yc(pc%n1), &
            pc%res1(pc%n1), pc%tmp1(pc%n1))
        pc%ready = .true.
        status = 0
    end subroutine hx_setup

    subroutine hx_destroy(pc)
        type(hx_preconditioner_t), intent(inout) :: pc
        integer :: j
        if (allocated(pc%aux)) then
            do j = 1, size(pc%aux)
                call aux_solver_free(pc%aux(j)%solver)
            end do
        end if
        call aux_solver_free(pc%coarse)
        pc%ready = .false.
    end subroutine hx_destroy

    subroutine hx_apply(self, input, output)
        !! output = B input, one symmetric cycle.
        type(hx_preconditioner_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        integer :: s, i

        output = 0.0_dp
        if (.not. self%ready) return
        self%applications = self%applications + 1
        if (self%two_level) then
            do s = 1, self%opt%smoothing_steps
                call fine_smooth(self, input, output)
            end do
            call residual(self%a, input, output, self%res)
            call csr_matvec(self%et, self%res, self%rc)
            call level1_cycle(self, self%rc, self%yc)
            call csr_matvec(self%e, self%yc, self%tmp)
            !$omp parallel do if(size(output) >= 4096)
            do i = 1, size(output)
                output(i) = output(i) + self%tmp(i)
            end do
            do s = 1, self%opt%smoothing_steps
                call fine_smooth(self, input, output)
            end do
        else
            call level1_cycle(self, input, output)
        end if
        where (self%fixed) output = input
    end subroutine hx_apply

    subroutine fine_smooth(self, b, x)
        type(hx_preconditioner_t), intent(inout) :: self
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        type(fortsparse_status_t) :: fst
        integer :: i

        if (self%opt%smoother == HX_SMOOTHER_PATCH) then
            call residual(self%a, b, x, self%res)
            call patch_apply(self%patch, self%res, self%tmp)
            !$omp parallel do if(size(x) >= 4096)
            do i = 1, size(x)
                x(i) = x(i) + self%patch%omega*self%tmp(i)
            end do
        else
            call smoother_apply(self%a, self%cheb, b, x, fst)
        end if
    end subroutine fine_smooth

    subroutine level1_cycle(self, b, y)
        !! S_1, P_1 .. P_m .. P_1, S_1 on A_1 (or an exact coarse solve).
        type(hx_preconditioner_t), intent(inout) :: self
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: y(:)
        integer :: j, m, v, i

        y = 0.0_dp
        if (self%opt%coarse == HX_COARSE_DIRECT) then
            call aux_solve(self%coarse, b, y)
            return
        end if
        call level1_smooth(self, b, y)
        m = size(self%aux)
        do v = 1, 2*m - 1
            j = v
            if (v > m) j = 2*m - v
            call residual(self%a1, b, y, self%res1)
            call csr_matvec(self%aux(j)%pt, self%res1, self%aux(j)%r)
            call aux_solve(self%aux(j)%solver, self%aux(j)%r, self%aux(j)%y)
            call csr_matvec(self%aux(j)%p, self%aux(j)%y, self%tmp1)
            !$omp parallel do if(size(y) >= 4096)
            do i = 1, size(y)
                y(i) = y(i) + self%tmp1(i)
            end do
        end do
        call level1_smooth(self, b, y)
    end subroutine level1_cycle

    subroutine level1_smooth(self, b, y)
        type(hx_preconditioner_t), intent(inout) :: self
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: y(:)
        type(fortsparse_status_t) :: fst
        integer :: s, i

        do s = 1, self%opt%smoothing_steps
            if (.not. self%two_level .and. self%patch%npatch > 0) then
                call residual(self%a1, b, y, self%res1)
                call patch_apply(self%patch, self%res1, self%tmp1)
                !$omp parallel do if(size(y) >= 4096)
                do i = 1, size(y)
                    y(i) = y(i) + self%patch%omega*self%tmp1(i)
                end do
            else
                call smoother_apply(self%a1, self%cheb1, b, y, fst)
            end if
        end do
    end subroutine level1_smooth

    subroutine residual(a, b, x, r)
        type(csr_t), intent(in) :: a
        real(dp), intent(in) :: b(:), x(:)
        real(dp), intent(out) :: r(:)
        integer :: i
        call csr_matvec(a, x, r)
        !$omp parallel do if(size(r) >= 4096)
        do i = 1, size(r)
            r(i) = b(i) - r(i)
        end do
    end subroutine residual

    ! ------------------------------------------------------------------
    ! Auxiliary (and exact coarse) solvers
    ! ------------------------------------------------------------------

    subroutine aux_solver_setup(a, kind, s, status)
        type(csc_t), intent(in) :: a
        integer, intent(in) :: kind
        type(aux_solver_t), intent(inout) :: s
        integer, intent(out) :: status
        type(amg_params_t) :: params
        type(fortsparse_status_t) :: fst

        status = 1
        s%kind = kind
        select case (kind)
        case (HX_AUX_AMG)
            call hx_csc_to_csr(a, s%a, status)
            if (status /= 0) return
            call amg_setup_csr(s%a%nrow, s%a%row_ptr, s%a%col, s%a%val, params, &
                s%amg, fst)
            status = merge(0, 1, status_ok(fst))
        case (HX_AUX_DIRECT)
            s%direct%backend_id = FORTSPARSE_BACKEND_NATIVE
            s%direct%options%symmetry = FORTSPARSE_SYM_SPD
            call sparse_factor(s%direct, a, fst)
            status = merge(0, 1, status_ok(fst))
        end select
    end subroutine aux_solver_setup

    subroutine aux_solve(s, r, y)
        type(aux_solver_t), intent(inout) :: s
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: y(:)
        type(fortsparse_status_t) :: fst
        select case (s%kind)
        case (HX_AUX_AMG)
            call amg_apply(s%amg, r, y, fst)
        case (HX_AUX_DIRECT)
            call sparse_solve(s%direct, r, y, fst)
        case default
            y = 0.0_dp
        end select
    end subroutine aux_solve

    subroutine aux_solver_free(s)
        type(aux_solver_t), intent(inout) :: s
        select case (s%kind)
        case (HX_AUX_AMG)
            call amg_destroy(s%amg)
        case (HX_AUX_DIRECT)
            call sparse_destroy(s%direct)
        end select
        s%kind = 0
    end subroutine aux_solver_free

    ! ------------------------------------------------------------------
    ! Overlapping patch (additive Schwarz) smoother
    ! ------------------------------------------------------------------

    subroutine patch_setup(a, fixed, ptr, dof, ps, bytes, status)
        !! Dense Cholesky factors of the free-dof blocks A_P of every patch;
        !! application z = sum_P R_P^T A_P^{-1} R_P r, gathered per dof
        !! (no atomics, thread-count independent).  omega = 1/lambda_max(BA)
        !! from a power iteration.
        type(csr_t), intent(in) :: a
        logical, intent(in) :: fixed(:)
        integer, intent(in) :: ptr(:), dof(:)
        type(patch_smoother_t), intent(inout) :: ps
        real(dp), intent(out) :: bytes
        integer, intent(out) :: status
        integer, allocatable :: marker(:), cnt(:)
        real(dp), allocatable :: x(:), y(:), z(:)
        integer :: np, p, i, j, k, m, n, q, info, fail, it
        real(dp) :: lam, nrm

        status = 1
        bytes = 0.0_dp
        np = size(ptr) - 1
        n = a%nrow
        if (np < 1) return
        ! free dofs of every patch
        allocate (ps%ptr(np + 1))
        ps%ptr(1) = 1
        do p = 1, np
            m = 0
            do k = ptr(p), ptr(p + 1) - 1
                if (dof(k) < 1 .or. dof(k) > n) return
                if (.not. fixed(dof(k))) m = m + 1
            end do
            ps%ptr(p + 1) = ps%ptr(p) + m
        end do
        allocate (ps%dof(ps%ptr(np + 1) - 1))
        do p = 1, np
            m = ps%ptr(p) - 1
            do k = ptr(p), ptr(p + 1) - 1
                if (fixed(dof(k))) cycle
                m = m + 1
                ps%dof(m) = dof(k)
            end do
        end do
        ps%npatch = np
        allocate (ps%fptr(np + 1))
        ps%fptr(1) = 1
        do p = 1, np
            m = ps%ptr(p + 1) - ps%ptr(p)
            ps%fptr(p + 1) = ps%fptr(p) + int(m, 8)*int(m, 8)
        end do
        allocate (ps%factor(ps%fptr(np + 1) - 1), ps%buf(size(ps%dof)))
        bytes = 8.0_dp*real(size(ps%factor), dp)
        fail = 0
        !$omp parallel default(shared) private(marker, p, m, i, j, k, q, info)
        allocate (marker(n))
        marker = 0
        !$omp do schedule(dynamic, 16) reduction(max:fail)
        do p = 1, np
            m = ps%ptr(p + 1) - ps%ptr(p)
            if (m == 0) cycle
            do i = 1, m
                marker(ps%dof(ps%ptr(p) + i - 1)) = i
            end do
            ps%factor(ps%fptr(p):ps%fptr(p + 1) - 1) = 0.0_dp
            do j = 1, m
                q = ps%dof(ps%ptr(p) + j - 1)
                do k = a%row_ptr(q), a%row_ptr(q + 1) - 1
                    i = marker(a%col(k))
                    if (i > 0) ps%factor(ps%fptr(p) + int(j - 1, 8)*m + i - 1) &
                        = a%val(k)
                end do
            end do
            call dpotrf('L', m, ps%factor(ps%fptr(p)), m, info)
            if (info /= 0) fail = max(fail, 1)
            do i = 1, m
                marker(ps%dof(ps%ptr(p) + i - 1)) = 0
            end do
        end do
        !$omp end do
        deallocate (marker)
        !$omp end parallel
        if (fail /= 0) return
        ! dof -> buffer positions
        allocate (cnt(n + 1))
        cnt = 0
        do k = 1, size(ps%dof)
            cnt(ps%dof(k) + 1) = cnt(ps%dof(k) + 1) + 1
        end do
        cnt(1) = 1
        do i = 2, n + 1
            cnt(i) = cnt(i) + cnt(i - 1)
        end do
        ps%gptr = cnt
        allocate (ps%gidx(size(ps%dof)))
        do k = 1, size(ps%dof)
            ps%gidx(cnt(ps%dof(k))) = k
            cnt(ps%dof(k)) = cnt(ps%dof(k)) + 1
        end do
        ! omega from lambda_max(B A)
        allocate (x(n), y(n), z(n))
        do i = 1, n
            x(i) = sin(1.2345_dp*real(i, dp)) + 0.5_dp
            if (fixed(i)) x(i) = 0.0_dp
        end do
        lam = 1.0_dp
        do it = 1, 20
            nrm = norm2(x)
            if (nrm <= tiny(1.0_dp)) exit
            x = x/nrm
            call csr_matvec(a, x, y)
            call patch_apply(ps, y, z)
            lam = norm2(z)
            x = z
        end do
        ps%omega = 1.0_dp/max(lam, tiny(1.0_dp))
        status = 0
    end subroutine patch_setup

    subroutine patch_apply(ps, r, z)
        type(patch_smoother_t), intent(inout) :: ps
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer :: p, m, k, i, info
        real(dp) :: s

        !$omp parallel do schedule(dynamic, 16) private(m, k, info)
        do p = 1, ps%npatch
            m = ps%ptr(p + 1) - ps%ptr(p)
            if (m == 0) cycle
            do k = ps%ptr(p), ps%ptr(p + 1) - 1
                ps%buf(k) = r(ps%dof(k))
            end do
            call dpotrs('L', m, 1, ps%factor(ps%fptr(p)), m, ps%buf(ps%ptr(p)), &
                m, info)
        end do
        !$omp end parallel do
        !$omp parallel do private(s, k) if(size(z) >= 4096)
        do i = 1, size(z)
            s = 0.0_dp
            do k = ps%gptr(i), ps%gptr(i + 1) - 1
                s = s + ps%buf(ps%gidx(k))
            end do
            z(i) = s
        end do
        !$omp end parallel do
    end subroutine patch_apply

    ! ------------------------------------------------------------------
    ! Sparse helpers
    ! ------------------------------------------------------------------

    subroutine hx_masked_csr(a, fixed, r, status)
        !! CSR of a with identity rows and columns at the fixed dofs.
        type(csc_t), intent(in) :: a
        logical, intent(in) :: fixed(:)
        type(csr_t), intent(out) :: r
        integer, intent(out) :: status
        type(csc_t) :: m
        call mask_csc(a, fixed, fixed, .true., m, status)
        if (status /= 0) return
        call hx_csc_to_csr(m, r, status)
    end subroutine hx_masked_csr

    subroutine mask_csc(a, rowfix, colfix, identity, m, status)
        !! Drop rows rowfix and columns colfix; identity: put 1 on the
        !! diagonal of every dropped (square) row.
        type(csc_t), intent(in) :: a
        logical, intent(in) :: rowfix(:), colfix(:)
        logical, intent(in) :: identity
        type(csc_t), intent(out) :: m
        integer, intent(out) :: status
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:)
        type(fortsparse_status_t) :: fst
        integer :: j, p, c

        status = 1
        if (size(rowfix) /= a%nrow .or. size(colfix) /= a%ncol) return
        allocate (row(a%nnz + a%ncol), col(a%nnz + a%ncol), val(a%nnz + a%ncol))
        c = 0
        do j = 1, a%ncol
            if (identity .and. colfix(j)) then
                c = c + 1
                row(c) = j; col(c) = j; val(c) = 1.0_dp
                cycle
            end if
            if (colfix(j)) cycle
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (rowfix(a%row_idx(p))) cycle
                c = c + 1
                row(c) = a%row_idx(p); col(c) = j; val(c) = a%val(p)
            end do
        end do
        call csc_from_triplet(a%nrow, a%ncol, row(:c), col(:c), val(:c), m, fst)
        if (status_ok(fst)) status = 0
    end subroutine mask_csc

    subroutine symmetric_regularised(a, s, status)
        !! s = (a + a^T)/2 + shift diag(a), with a unit diagonal on empty rows
        !! (dofs that the auxiliary map does not reach); exact symmetry for
        !! the SPD factor.  P^T A P is only semidefinite when P has a kernel
        !! (possible for the componentwise interpolation on coarse meshes);
        !! the relative shift makes it definite and leaves the correction
        !! P (P^T A P)^{-1} P^T unchanged on ker P.
        type(csc_t), intent(in) :: a
        type(csc_t), intent(out) :: s
        integer, intent(out) :: status
        integer, allocatable :: row(:), col(:)
        real(dp), allocatable :: val(:)
        logical, allocatable :: live(:)
        type(fortsparse_status_t) :: fst
        integer :: j, p, c
        real(dp), parameter :: shift = 1.0e-8_dp

        status = 1
        if (a%nrow /= a%ncol) return
        allocate (row(2*a%nnz + 2*a%ncol), col(2*a%nnz + 2*a%ncol), &
            val(2*a%nnz + 2*a%ncol), live(a%ncol))
        live = .false.
        c = 0
        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (a%row_idx(p) == j .and. a%val(p) > 0.0_dp) then
                    live(j) = .true.
                    c = c + 1
                    row(c) = j; col(c) = j; val(c) = shift*a%val(p)
                end if
                c = c + 1
                row(c) = a%row_idx(p); col(c) = j; val(c) = 0.5_dp*a%val(p)
                c = c + 1
                row(c) = j; col(c) = a%row_idx(p); val(c) = 0.5_dp*a%val(p)
            end do
        end do
        do j = 1, a%ncol
            if (live(j)) cycle
            c = c + 1
            row(c) = j; col(c) = j; val(c) = 1.0_dp
        end do
        call csc_from_triplet(a%nrow, a%ncol, row(:c), col(:c), val(:c), s, fst)
        if (status_ok(fst)) status = 0
    end subroutine symmetric_regularised

    subroutine hx_csc_to_csr(a, r, status)
        !! CSR of a (sorted columns) via the transpose.
        type(csc_t), intent(in) :: a
        type(csr_t), intent(out) :: r
        integer, intent(out) :: status
        type(csc_t) :: t
        type(fortsparse_status_t) :: fst
        status = 1
        call csc_transpose(a, t, fst)
        if (.not. status_ok(fst)) return
        r%nrow = a%nrow
        r%ncol = a%ncol
        r%row_ptr = t%col_ptr
        r%col = t%row_idx(:t%nnz)
        r%val = t%val(:t%nnz)
        status = 0
    end subroutine hx_csc_to_csr

    subroutine transpose_csr(a, r, status)
        !! CSR of a^T (the CSC arrays of a reinterpreted).
        type(csc_t), intent(in) :: a
        type(csr_t), intent(out) :: r
        integer, intent(out) :: status
        type(csc_t) :: t
        type(fortsparse_status_t) :: fst
        status = 1
        ! canonical (sorted) row indices: transpose twice
        call csc_transpose(a, t, fst)
        if (.not. status_ok(fst)) return
        call hx_csc_to_csr(t, r, status)
    end subroutine transpose_csr

    subroutine hx_pcg(a, pc, b, x, tolerance, max_iterations, iterations, &
            relative_residual, info)
        !! Preconditioned CG for the masked operator a (hx_masked_csr) with
        !! the HX cycle; stops at |b - A x| <= tolerance |b| (true residual
        !! checked on exit).  x is the initial guess.  info = 0 on
        !! convergence, 1 at max_iterations, 2 on breakdown.
        !! The object-based CG of FortNum (real_conjugate_gradient with
        !! real_linear_operator_t) replaces this once the FortNum pin has it.
        type(csr_t), intent(in) :: a
        type(hx_preconditioner_t), intent(inout) :: pc
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations
        integer, intent(out) :: iterations, info
        real(dp), intent(out) :: relative_residual
        real(dp), allocatable :: r(:), z(:), p(:), q(:)
        real(dp) :: bnorm, rho, rho_old, alpha, pq
        integer :: n

        n = size(b)
        allocate (r(n), z(n), p(n), q(n))
        iterations = 0
        info = 1
        bnorm = norm2(b)
        relative_residual = 0.0_dp
        if (bnorm <= 0.0_dp) then
            x = 0.0_dp
            info = 0
            return
        end if
        call residual(a, b, x, r)
        relative_residual = pnorm(r)/bnorm
        if (relative_residual <= tolerance) then
            info = 0
            return
        end if
        call hx_apply(pc, r, z)
        rho = pdot(r, z)
        p = z
        do while (iterations < max_iterations)
            call csr_matvec(a, p, q)
            pq = pdot(p, q)
            if (.not. (pq > 0.0_dp)) then
                info = 2
                exit
            end if
            alpha = rho/pq
            call paxpy(alpha, p, x)
            call paxpy(-alpha, q, r)
            iterations = iterations + 1
            relative_residual = pnorm(r)/bnorm
            if (relative_residual <= tolerance) then
                info = 0
                exit
            end if
            call hx_apply(pc, r, z)
            rho_old = rho
            rho = pdot(r, z)
            call pxpby(z, rho/rho_old, p)
        end do
        call residual(a, b, x, r)
        relative_residual = pnorm(r)/bnorm
        if (info == 0 .and. relative_residual > 10.0_dp*tolerance) info = 1
    end subroutine hx_pcg

    real(dp) function pdot(u, v) result(s)
        real(dp), intent(in) :: u(:), v(:)
        integer :: i
        s = 0.0_dp
        !$omp parallel do reduction(+:s) if(size(u) >= 4096)
        do i = 1, size(u)
            s = s + u(i)*v(i)
        end do
    end function pdot

    real(dp) function pnorm(u) result(s)
        real(dp), intent(in) :: u(:)
        s = sqrt(pdot(u, u))
    end function pnorm

    subroutine paxpy(alpha, u, v)
        real(dp), intent(in) :: alpha, u(:)
        real(dp), intent(inout) :: v(:)
        integer :: i
        !$omp parallel do if(size(u) >= 4096)
        do i = 1, size(u)
            v(i) = v(i) + alpha*u(i)
        end do
    end subroutine paxpy

    subroutine pxpby(u, beta, v)
        real(dp), intent(in) :: u(:), beta
        real(dp), intent(inout) :: v(:)
        integer :: i
        !$omp parallel do if(size(u) >= 4096)
        do i = 1, size(u)
            v(i) = u(i) + beta*v(i)
        end do
    end subroutine pxpby

end module fortfem_hcurl_auxiliary_space
