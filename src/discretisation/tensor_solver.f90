module fortfem_tensor_solver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t, build_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                     apply_tensor_operator, tensor_diagonal, assemble_tensor_operator, &
                                       diffusion_coefficient
    use fortfem_tensor_contractions, only: tensor_contract
    use fortfem_tensor_seam, only: tensor_seam_t, identity_tensor_seam, &
                                   coarsen_tensor_seam, seam_expand, seam_fold, &
                                   apply_seam_operator, seam_diagonal, &
                                   assemble_seam_operator
    use fortfem_gll_line, only: gll_rule, lagrange_table
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
              sparse_direct_factor_csc, sparse_direct_solve_factored, sparse_direct_free
    use fortsparse, only: sparse_options_t, FORTSPARSE_SYM_SPD
    implicit none
    private
    public :: tensor_solve_info_t, tensor_multigrid_t, build_tensor_multigrid
    public :: solve_tensor_diffusion, tensor_vcycle, free_tensor_multigrid
    public :: tensor_transfer, coarse_solver_hook, tensor_smooth
    public :: build_line_multigrid, solve_tensor_pcg
    type :: tensor_solve_info_t
        integer :: status = 1, iterations = 0
        real(dp) :: residual = huge(1.0_dp)
    end type
    type :: tensor_level_t
        type(tensor_operator_t) :: op
        !! Vectors of this level live on the seam's reduced (primary) DOFs.
        type(tensor_seam_t) :: seam
        real(dp), allocatable :: diagonal(:), interpolation(:, :)
        integer, allocatable :: owner(:)
        real(dp) :: lambda_max = 1
        ! Nested line-table hierarchies: prolongation to this level from the
        ! next coarser one and its transpose (both CSR).
        logical :: algebraic = .false.
        type(csr_matrix_t) :: prolong, restrict
    end type
    type :: tensor_multigrid_t
        type(tensor_level_t), allocatable :: level(:)
        type(csr_matrix_t) :: coarse
        type(sparse_direct_factor_t) :: factor
        integer :: nlevels = 0, smoothing = 3
    end type
    abstract interface
        subroutine coarse_solver_hook(matrix, rhs, x, status)
            import csr_matrix_t, dp
            type(csr_matrix_t), intent(in) :: matrix
            real(dp), intent(in) :: rhs(:)
            real(dp), intent(out) :: x(:)
            integer, intent(out) :: status
        end subroutine
    end interface
contains
    subroutine build_tensor_multigrid(mesh, coefficient, mg, status, external_coarse, &
                                      seam)
        type(tensor_mesh_t), intent(in) :: mesh
        procedure(diffusion_coefficient) :: coefficient
        type(tensor_multigrid_t), intent(out) :: mg
        integer, intent(out) :: status
        logical, intent(in), optional :: external_coarse
        !! Nonconforming periodic seam; coarse levels constrain consistently.
        type(tensor_seam_t), intent(in), optional :: seam
        type(tensor_mesh_t) :: coarse_mesh
        real(dp), allocatable :: xf(:), wf(:), xc(:), wc(:), d(:, :)
        integer :: p, l, e, a, n, i
        logical :: external
        external = .false.
        if (present(external_coarse)) external = external_coarse
        p = mesh%degree
        mg%nlevels = 1
        do while (p > 1)
            p = (p + 1)/2
            mg%nlevels = mg%nlevels + 1
        end do
        allocate (mg%level(mg%nlevels))
        call build_tensor_operator(mesh, coefficient, mg%level(1)%op, status)
        if (status /= 0) return
        if (present(seam)) then
            mg%level(1)%seam = seam
        else
            call identity_tensor_seam(mesh, mg%level(1)%seam)
        end if
        do l = 1, mg%nlevels
            call seam_diagonal(mg%level(l)%op, mg%level(l)%seam, mg%level(l)%diagonal)
            if (any(mg%level(l)%diagonal <= 0)) then
                status = 3
                return
            end if
            call estimate_spectrum(mg%level(l))
            if (l == mg%nlevels) exit
            p = (mg%level(l)%op%mesh%degree + 1)/2
            call build_tensor_mesh(mesh%dim, p, mesh%geometry_degree, &
                                   mesh%vertices, mesh%geometry, coarse_mesh, status)
            if (status /= 0) return
            call coarse_boundary(mg%level(l)%op%mesh, coarse_mesh)
            call build_tensor_operator(coarse_mesh, coefficient, &
                                       mg%level(l + 1)%op, status)
            if (status /= 0) return
            call coarsen_tensor_seam(mg%level(l)%seam, mg%level(l + 1)%op%mesh, &
                                     mg%level(l + 1)%seam, status)
            if (status /= 0) return
            call gll_rule(mg%level(l)%op%mesh%degree, xf, wf, status)
            call gll_rule(p, xc, wc, status)
            call lagrange_table(xc, xf, mg%level(l)%interpolation, d)
            n = mg%level(l)%op%mesh%ndof
            allocate (mg%level(l)%owner(n))
            mg%level(l)%owner = 0
            do e = 1, mesh%nelem
                do a = 1, size(mg%level(l)%op%mesh%dof, 1)
                    i = mg%level(l)%op%mesh%dof(a, e)
                    if (mg%level(l)%owner(i) == 0) mg%level(l)%owner(i) = e
                end do
            end do
        end do
        call factor_coarse(mg, external, status)
    end subroutine

    subroutine build_line_multigrid(meshes, transfers, coefficient, mg, status, &
                                    external_coarse, smoothing)
        ! Nested hierarchy for line-table meshes (h- or mode-coarsening): meshes
        ! from fine to coarse, transfers(l) = prolongation level l+1 -> l (CSR).
        type(tensor_mesh_t), intent(in) :: meshes(:)
        type(csr_matrix_t), intent(in) :: transfers(:)
        procedure(diffusion_coefficient) :: coefficient
        type(tensor_multigrid_t), intent(out) :: mg
        integer, intent(out) :: status
        logical, intent(in), optional :: external_coarse
        integer, intent(in), optional :: smoothing
        integer :: l, i, k, j, n
        integer, allocatable :: count(:)
        logical :: external
        status = 1
        external = .false.
        if (present(external_coarse)) external = external_coarse
        if (present(smoothing)) mg%smoothing = smoothing
        mg%nlevels = size(meshes)
        if (mg%nlevels < 1 .or. size(transfers) < mg%nlevels - 1) return
        allocate (mg%level(mg%nlevels))
        do l = 1, mg%nlevels
            call build_tensor_operator(meshes(l), coefficient, mg%level(l)%op, status)
            if (status /= 0) return
            call identity_tensor_seam(meshes(l), mg%level(l)%seam)
            call tensor_diagonal(mg%level(l)%op, mg%level(l)%diagonal)
            if (any(mg%level(l)%diagonal <= 0)) then
                status = 3
                return
            end if
            call estimate_spectrum(mg%level(l))
            if (l == mg%nlevels) exit
            status = 1
            if (transfers(l)%n /= meshes(l)%ndof) return
            if (maxval(transfers(l)%colidx) > meshes(l + 1)%ndof) return
            mg%level(l)%algebraic = .true.
            mg%level(l)%prolong = transfers(l)
            ! Explicit transpose for a race-free row-parallel restriction.
            n = meshes(l + 1)%ndof
            associate (r => mg%level(l)%restrict, pp => transfers(l))
                r%n = n
                allocate (r%rowptr(n + 1), r%colidx(size(pp%colidx)), &
                          r%val(size(pp%val)), count(n + 1))
                count = 0
                do k = 1, size(pp%colidx)
                    count(pp%colidx(k) + 1) = count(pp%colidx(k) + 1) + 1
                end do
                count(1) = 1
                do i = 1, n
                    count(i + 1) = count(i + 1) + count(i)
                end do
                r%rowptr = count
                do i = 1, pp%n
                    do k = pp%rowptr(i), pp%rowptr(i + 1) - 1
                        j = pp%colidx(k)
                        r%colidx(count(j)) = i
                        r%val(count(j)) = pp%val(k)
                        count(j) = count(j) + 1
                    end do
                end do
                deallocate (count)
            end associate
        end do
        call factor_coarse(mg, external, status)
    end subroutine

    subroutine factor_coarse(mg, external, status)
        type(tensor_multigrid_t), intent(inout) :: mg
        logical, intent(in) :: external
        integer, intent(out) :: status
        type(sparse_options_t) :: options
        integer, allocatable :: ptr(:), rows(:), fill(:)
        real(dp), allocatable :: values(:)
        integer :: n, i, j, a, pos
        call assemble_seam_operator(mg%level(mg%nlevels)%op, &
                                    mg%level(mg%nlevels)%seam, mg%coarse)
        if (external) then
            status = 0
            return
        end if
        ! CSR -> CSC without dropping small entries or introducing a dense matrix.
        n = mg%coarse%n
        allocate (ptr(n + 1), fill(n), rows(size(mg%coarse%val)), &
                  values(size(mg%coarse%val)))
        ptr = 0
        do i = 1, size(mg%coarse%colidx)
            j = mg%coarse%colidx(i)
            ptr(j + 1) = ptr(j + 1) + 1
        end do
        ptr(1) = 1
        do i = 1, n
            ptr(i + 1) = ptr(i + 1) + ptr(i)
        end do
        fill = ptr(1:n)
        do i = 1, n
            do a = mg%coarse%rowptr(i), mg%coarse%rowptr(i + 1) - 1
                j = mg%coarse%colidx(a)
                pos = fill(j)
                rows(pos) = i
                values(pos) = mg%coarse%val(a)
                fill(j) = pos + 1
            end do
        end do
        options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_direct_factor_csc(mg%factor, n, ptr, rows, values, status, options)
    end subroutine

    subroutine coarse_boundary(fine, coarse)
        type(tensor_mesh_t), intent(in) :: fine
        type(tensor_mesh_t), intent(inout) :: coarse
        integer :: e, dir, side, a, b, n, p, ia(3), ib(3)
        logical :: boundary
        n = fine%degree + 1; p = coarse%degree + 1
        do e = 1, fine%nelem
            do dir = 1, fine%dim
                do side = 0, 1
                    boundary = .true.
                    do a = 1, n**fine%dim
                        ia = [modulo(a - 1, n), modulo((a - 1)/n, n), (a - 1)/n**2]
                        if (ia(dir) /= side*(n - 1)) cycle
                        if (.not. fine%boundary(fine%dof(a, e))) boundary = .false.
                    end do
                    if (.not. boundary) cycle
                    do b = 1, p**fine%dim
                        ib = [modulo(b - 1, p), modulo((b - 1)/p, p), (b - 1)/p**2]
                        if (ib(dir) == side*(p - 1)) coarse%boundary(coarse%dof(b, e)) &
                            = .true.
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine estimate_spectrum(level)
        type(tensor_level_t), intent(inout) :: level
        real(dp), allocatable :: v(:), previous(:), work(:), product(:), scale(:)
        real(dp) :: alpha, beta, last, upper
        integer :: i, it, n
        n = level%seam%nreduced
        allocate (v(n), previous(n), work(n), product(n))
        scale = sqrt(level%diagonal)
        do i = 1, n
            v(i) = sin(real(137*i, dp))
        end do
        where (level%seam%boundary) v = 0
        beta = norm2(v)
        if (beta <= tiny(1.0_dp)) return
        v = v/beta
        previous = 0; last = 0; upper = 0
        do it = 1, 12
            work = v/scale
            call apply_seam_operator(level%op, level%seam, work, product)
            product = product/scale - last*previous
            alpha = dot_product(v, product)
            product = product - alpha*v
            beta = norm2(product)
            upper = max(upper, alpha + last + beta)
            if (beta < 1.e-13_dp) exit
            previous = v; v = product/beta; last = beta
        end do
        level%lambda_max = 1.15_dp*upper
    end subroutine

    subroutine tensor_transfer(mg, l, coarse, fine, restrict)
        !! p-transfer between the reduced spaces of levels l and l+1:
        !! prolongation R_f I P_c and its transpose P_c^T I^T R_f^T.
        type(tensor_multigrid_t), intent(in) :: mg
        integer, intent(in) :: l
        real(dp), intent(inout) :: coarse(:), fine(:)
        logical, intent(in) :: restrict
        real(dp), allocatable :: cf(:), ff(:)
        associate (sf => mg%level(l)%seam, sc => mg%level(l + 1)%seam)
            if (sf%nslave == 0 .and. sc%nslave == 0) then
                call transfer_full(mg, l, coarse, fine, restrict)
                return
            end if
            allocate (cf(sc%nfull), ff(sf%nfull))
            if (restrict) then
                ff = 0
                ff(sf%primary) = fine
                call transfer_full(mg, l, cf, ff, .true.)
                call seam_fold(sc, cf, coarse)
                where (sc%boundary) coarse = 0
            else
                call seam_expand(sc, coarse, cf)
                call transfer_full(mg, l, cf, ff, .false.)
                fine = ff(sf%primary)
            end if
        end associate
    end subroutine

    subroutine transfer_full(mg, l, coarse, fine, restrict)
        type(tensor_multigrid_t), intent(in) :: mg
        integer, intent(in) :: l
        real(dp), intent(inout) :: coarse(:), fine(:)
        logical, intent(in) :: restrict
        integer :: nf(3), nc(3), ws, e, a, id, colour, ce
        if (mg%level(l)%algebraic) then
            call csr_transfer(mg%level(l), restrict, coarse, fine)
            if (restrict) then
                where (mg%level(l + 1)%op%mesh%boundary) coarse = 0
            else
                where (mg%level(l)%op%mesh%boundary) fine = 0
            end if
            return
        end if
        nf = mg%level(l)%op%ni; nc = mg%level(l + 1)%op%ni
        ws = maxval(nf)**3
        if (restrict) then
            coarse = 0
        else
            fine = 0
        end if
        !$omp parallel default(shared) private(e,a,id,colour,ce)
        block
            real(dp), allocatable :: u(:), v(:), t1(:), t2(:), b(:, :)
            real(dp) :: identity(1, 1)
            allocate (u(ws), v(ws), t1(ws), t2(ws))
            identity = 1
            if (restrict) then
                b = transpose(mg%level(l)%interpolation)
            else
                b = mg%level(l)%interpolation
            end if
            do colour = 1, mg%level(l + 1)%op%ncolours
                !$omp do schedule(static)
                do ce = mg%level(l + 1)%op%colour_ptr(colour), &
                    mg%level(l + 1)%op%colour_ptr(colour + 1) - 1
                    e = mg%level(l + 1)%op%colour_cells(ce)
                    if (restrict) then
                        u = 0
                        do a = 1, product(nf)
                            id = mg%level(l)%op%mesh%dof(a, e)
                            if (mg%level(l)%owner(id) /= e) cycle
                            if (mg%level(l)%op%mesh%boundary(id)) cycle
                            u(a) = fine(id)
                        end do
                        if (mg%level(l)%op%mesh%dim == 2) then
                            call tensor_contract(nf, nc, b, b, identity, u, v, t1, t2)
                        else
                            call tensor_contract(nf, nc, b, b, b, u, v, t1, t2)
                        end if
                        do a = 1, product(nc)
                            id = mg%level(l + 1)%op%mesh%dof(a, e)
                            coarse(id) = coarse(id) + v(a)
                        end do
                    else
                        u(1:product(nc)) = coarse(mg%level(l + 1)%op%mesh%dof(:, e))
                        if (mg%level(l)%op%mesh%dim == 2) then
                            call tensor_contract(nc, nf, b, b, identity, u, v, t1, t2)
                        else
                            call tensor_contract(nc, nf, b, b, b, u, v, t1, t2)
                        end if
                        do a = 1, product(nf)
                            id = mg%level(l)%op%mesh%dof(a, e)
                            if (mg%level(l)%owner(id) == e) fine(id) = v(a)
                        end do
                    end if
                end do
                !$omp end do
            end do
        end block
        !$omp end parallel
        if (restrict) then
            where (mg%level(l + 1)%op%mesh%boundary) coarse = 0
        else
            where (mg%level(l)%op%mesh%boundary) fine = 0
        end if
    end subroutine

    subroutine csr_transfer(level, restrict, coarse, fine)
        type(tensor_level_t), intent(in) :: level
        logical, intent(in) :: restrict
        real(dp), intent(inout) :: coarse(:), fine(:)
        real(dp), allocatable :: masked(:)
        integer :: i, k
        if (restrict) then
            masked = fine
            where (level%op%mesh%boundary) masked = 0
            !$omp parallel do schedule(static) private(k)
            do i = 1, level%restrict%n
                coarse(i) = 0
                do k = level%restrict%rowptr(i), level%restrict%rowptr(i + 1) - 1
                    coarse(i) = coarse(i) + level%restrict%val(k)* &
                                masked(level%restrict%colidx(k))
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) private(k)
            do i = 1, level%prolong%n
                fine(i) = 0
                do k = level%prolong%rowptr(i), level%prolong%rowptr(i + 1) - 1
                fine(i) = fine(i) + level%prolong%val(k)*coarse(level%prolong%colidx(k))
                end do
            end do
            !$omp end parallel do
        end if
    end subroutine

    subroutine tensor_smooth(op, diagonal, lambda_max, rhs, x, steps, chebyshev, seam)
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: diagonal(:), lambda_max, rhs(:)
        real(dp), intent(inout) :: x(:)
        integer, intent(in) :: steps
        logical, intent(in) :: chebyshev
        type(tensor_seam_t), intent(in), optional :: seam
        real(dp) :: theta, delta, omega
        real(dp), allocatable :: ax(:)
        integer :: j
        allocate (ax(size(x)))
        theta = 0.55_dp*lambda_max
        delta = 0.45_dp*lambda_max
        do j = 1, steps
            omega = 1/theta
            if (chebyshev) omega = 1/(theta - delta* &
                              cos(acos(-1.0_dp)*(real(j, dp) - 0.5_dp)/real(steps, dp)))
            if (present(seam)) then
                call apply_seam_operator(op, seam, x, ax)
            else
                call apply_tensor_operator(op, x, ax)
            end if
            x = x + omega*(rhs - ax)/diagonal
        end do
    end subroutine

    recursive subroutine tensor_vcycle(mg, l, rhs, x, status, coarse_hook)
        type(tensor_multigrid_t), intent(inout) :: mg
        integer, intent(in) :: l
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: status
        procedure(coarse_solver_hook), optional :: coarse_hook
        real(dp), allocatable :: r(:), ax(:), rc(:), xc(:), correction(:)
        status = 0
        x = 0
        if (l == mg%nlevels) then
            if (present(coarse_hook)) then
                call coarse_hook(mg%coarse, rhs, x, status)
            else
                call sparse_direct_solve_factored(mg%factor, rhs, x, status)
            end if
            return
        end if
        allocate (r(size(rhs)), ax(size(rhs)), correction(size(rhs)), &
                  rc(mg%level(l + 1)%seam%nreduced), xc(mg%level(l + 1)%seam%nreduced))
        call tensor_smooth(mg%level(l)%op, mg%level(l)%diagonal, &
                           mg%level(l)%lambda_max, rhs, x, mg%smoothing, .true., &
                           seam=mg%level(l)%seam)
        call apply_seam_operator(mg%level(l)%op, mg%level(l)%seam, x, ax)
        r = rhs - ax
        call tensor_transfer(mg, l, rc, r, .true.)
        call tensor_vcycle(mg, l + 1, rc, xc, status, coarse_hook)
        if (status /= 0) return
        call tensor_transfer(mg, l, xc, correction, .false.)
        x = x + correction
        call tensor_smooth(mg%level(l)%op, mg%level(l)%diagonal, &
                           mg%level(l)%lambda_max, rhs, x, mg%smoothing, .true., &
                           seam=mg%level(l)%seam)
    end subroutine

    subroutine solve_tensor_diffusion(mesh, coefficient, rhs, tol, x, info, &
                                      bc, maxiter, coarse_hook, seam)
        !! rhs, bc and x are full-mesh vectors also with a seam; the solve runs
        !! on the seam's primary DOFs and returns x = P x_reduced.
        type(tensor_mesh_t), intent(in) :: mesh
        procedure(diffusion_coefficient) :: coefficient
        real(dp), intent(in) :: rhs(:), tol
        real(dp), intent(inout) :: x(:)
        type(tensor_solve_info_t), intent(out) :: info
        real(dp), intent(in), optional :: bc(:)
        integer, intent(in), optional :: maxiter
        procedure(coarse_solver_hook), optional :: coarse_hook
        type(tensor_seam_t), intent(in), optional :: seam
        type(tensor_multigrid_t) :: mg
        integer :: status, n
        info%status = 1
        n = mesh%ndof
        if (size(x) /= n .or. size(rhs) /= n .or. tol <= 0) return
        if (present(seam)) then
            if (seam%nfull /= n) return
        end if
        call build_tensor_multigrid(mesh, coefficient, mg, status, &
                                    external_coarse=present(coarse_hook), seam=seam)
        if (status /= 0) then
            info%status = status
            call free_tensor_multigrid(mg)
            return
        end if
        call solve_tensor_pcg(mg, rhs, tol, x, info, bc, maxiter, coarse_hook)
        call free_tensor_multigrid(mg)
    end subroutine

    subroutine solve_tensor_pcg(mg, rhs, tol, x, info, bc, maxiter, coarse_hook)
        ! Multigrid-preconditioned CG on level 1 of a built hierarchy. rhs, bc
        ! and x are full-mesh vectors; CG runs on the level's reduced DOFs.
        type(tensor_multigrid_t), intent(inout) :: mg
        real(dp), intent(in) :: rhs(:), tol
        real(dp), intent(inout) :: x(:)
        type(tensor_solve_info_t), intent(out) :: info
        real(dp), intent(in), optional :: bc(:)
        integer, intent(in), optional :: maxiter
        procedure(coarse_solver_hook), optional :: coarse_hook
        real(dp), allocatable :: b(:), r(:), z(:), p(:), ap(:), lift(:), u(:)
        real(dp) :: rho, previous, alpha, denom, initial
        integer :: limit, it, status, n, nfull
        info%status = 1
        if (mg%nlevels < 1) return
        nfull = mg%level(1)%op%mesh%ndof
        if (size(x) /= nfull .or. size(rhs) /= nfull .or. tol <= 0) return
        limit = 300
        if (present(maxiter)) limit = maxiter
        associate (level => mg%level(1))
            n = level%seam%nreduced
            allocate (b(n), r(n), z(n), p(n), ap(n), lift(n), u(n))
            call reduce(rhs, b, .true.)
            call reduce(x, u, .false.)
            lift = 0
            if (present(bc)) then
                if (size(bc) /= nfull) return
                call reduce(bc, z, .false.)
                where (level%seam%boundary) lift = z
                call apply_seam_operator(level%op, level%seam, lift, ap, &
                                         eliminate=.false.)
                b = b - ap
            end if
            where (level%seam%boundary)
                b = 0
                u = 0
            end where
            call apply_seam_operator(level%op, level%seam, u, ap)
            r = b - ap
            initial = max(norm2(b), tiny(1.0_dp))
            info%residual = norm2(r)/initial
            if (info%residual <= tol) then
                info%status = 0
            else
                previous = 0
                do it = 1, limit
                    call tensor_vcycle(mg, 1, r, z, status, coarse_hook)
                    if (status /= 0) exit
                    rho = dot_product(r, z)
                    if (rho <= 0) exit
                    if (it == 1) then
                        p = z
                    else
                        p = z + (rho/previous)*p
                    end if
                    call apply_seam_operator(level%op, level%seam, p, ap)
                    denom = dot_product(p, ap)
                    if (denom <= 0) exit
                    alpha = rho/denom
                    u = u + alpha*p
                    r = r - alpha*ap
                    info%iterations = it
                    info%residual = norm2(r)/initial
                    if (info%residual <= tol) then
                        ! Confirm using the actual operator, not only the recurrence.
                        call apply_seam_operator(level%op, level%seam, u, ap)
                        r = b - ap
                        info%residual = norm2(r)/initial
                        if (info%residual <= tol) then
                            info%status = 0
                            exit
                        end if
                    end if
                    previous = rho
                end do
            end if
            u = u + lift
            if (level%seam%nslave == 0) then
                x = u
            else
                call seam_expand(level%seam, u, x)
            end if
        end associate
    contains
        subroutine reduce(full, reduced, fold)
            !! Loads fold with P^T; primary values are restricted by injection.
            real(dp), intent(in) :: full(:)
            real(dp), intent(out) :: reduced(:)
            logical, intent(in) :: fold
            if (mg%level(1)%seam%nslave == 0) then
                reduced = full
            else if (fold) then
                call seam_fold(mg%level(1)%seam, full, reduced)
            else
                reduced = full(mg%level(1)%seam%primary)
            end if
        end subroutine
    end subroutine

    subroutine free_tensor_multigrid(mg)
        type(tensor_multigrid_t), intent(inout) :: mg
        call sparse_direct_free(mg%factor)
        if (allocated(mg%level)) deallocate (mg%level)
        mg%nlevels = 0
    end subroutine
end module
