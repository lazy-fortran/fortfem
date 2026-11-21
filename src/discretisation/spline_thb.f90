module fortfem_spline_thb
    !! Truncated hierarchical (THB) refinement of C0 multipatch spline spaces
    !! (doc/design/multipatch_spline_spaces.md, THB refinement).
    !!
    !! Level spaces V_0 c V_1 c ... are the multipatch spaces on dyadically
    !! refined breakpoints, all blocks (THB_ISOTROPIC) or only the block
    !! directions flagged in layout%radial_refine (THB_RADIAL: the annulus
    !! radius, or the radial directions of an X-point layout); the level is
    !! chosen when it is created. Omega_m are unions of level-m cells, Omega_0 the whole
    !! domain. A level-m function is active if its support lies in Omega_m
    !! but not in Omega_(m+1), and enters truncated: its V_(m+1) coefficients
    !! of functions with support in Omega_(m+1) are dropped, level by level
    !! (Giannelli, Juttler, Speleers, CAGD 29, 2012). Refinement follows the
    !! class-2 admissible recursion of Buffa, Gantner, Giannelli,
    !! Praetorius, Vazquez (arXiv:2107.02023), so every truncated function
    !! overlaps at most two levels of active cells.
    !!
    !! The hierarchical space is an ordinary spline_space_t whose patches
    !! are the (block, level) pairs and whose cells are the active cells; its
    !! extraction rows hold each truncated function on the active cells of
    !! every level it reaches. The geometry of every level is the exact knot
    !! insertion of the level-0 map.
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_interval_bernstein_jets, only: interval_bernstein_jet
    use fortfem_gll_line, only: gauss_rule
    use fortfem_spline_layout, only: block_layout_t
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_cell_broken, spline_block_t, SPLINE_C0
    use fortfem_spline_eval, only: spline_quadrature_t, spline_evaluate_at
    use fortfem_spline_energy, only: pointwise_energy_t, spline_energy, &
        spline_energy_hessian_csr
    use fortsparse, only: pcg_params_t, pcg_handle_t, fortsparse_status_t, &
        pcg_setup_csr, pcg_solve, pcg_destroy, PCG_JACOBI, FORTSPARSE_OK
    implicit none
    private

    integer, parameter, public :: THB_ISOTROPIC = 0
    integer, parameter, public :: THB_RADIAL = 1

    type :: rows_t
        integer, allocatable :: ptr(:), col(:)
        real(dp), allocatable :: val(:)
    end type rows_t

    type :: thb_level_t
        type(spline_space_t) :: space
        type(spline_space_options_t) :: options
        integer :: direction = THB_ISOTROPIC
        logical, allocatable :: omega(:), refined(:) ! (cell)
        integer, allocatable :: parent(:) ! (cell): cell of the coarser level
        type(rows_t) :: children ! cell -> cells of the finer level
        type(rows_t) :: cell_dofs ! cell -> functions nonzero on it
        type(rows_t) :: support ! function -> cells
        type(rows_t) :: dof_rows ! function -> (broken, weight)
        type(rows_t) :: prolong ! function -> finer functions (coefficients)
    end type thb_level_t

    type, public :: spline_thb_t
        type(block_layout_t) :: layout
        type(thb_level_t), allocatable :: level(:) ! level(1) is V_0
        type(spline_space_t) :: space
        integer, allocatable :: cell_level(:), cell_index(:) ! (cell of space)
    end type spline_thb_t

    type, extends(pointwise_energy_t) :: projection_kernel_t
        real(dp), allocatable :: g(:)
    contains
        procedure :: evaluate => projection_evaluate
    end type projection_kernel_t

    public :: spline_thb_init
    public :: spline_thb_refine
    public :: spline_thb_mark
    public :: spline_project

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    subroutine spline_thb_init(layout, options, thb, status)
        !! Hierarchy with the single level V_0. Status 1: not a C0 space or
        !! explicit multiplicities, otherwise the spline_space_build status.
        type(block_layout_t), intent(in) :: layout
        type(spline_space_options_t), intent(in) :: options
        type(spline_thb_t), intent(out) :: thb
        integer, intent(out) :: status

        integer :: k, d

        status = 1
        if (options%regularity /= SPLINE_C0) return
        if (allocated(options%knots)) then
            do k = 1, size(options%knots, 2)
                do d = 1, 2
                    if (allocated(options%knots(d, k)%multiplicity)) return
                end do
            end do
        end if
        thb%layout = layout
        allocate(thb%level(1))
        call spline_space_build(layout, options, thb%level(1)%space, status)
        if (status /= 0) return
        associate (lev => thb%level(1))
            lev%options = options
            do k = 1, layout%block_count
                do d = 1, 2
                    lev%options%knots(d, k)%reversed = .false.
                    if (d == 1) then
                        lev%options%knots(d, k)%breakpoints = lev%space%blocks(k)%breaks1
                    else
                        lev%options%knots(d, k)%breakpoints = lev%space%blocks(k)%breaks2
                    end if
                end do
            end do
            allocate(lev%omega(lev%space%cell_count), lev%refined(lev%space%cell_count))
            lev%omega = .true.
            lev%refined = .false.
            call level_tables(lev)
        end associate
        call assemble(thb, status)
    end subroutine spline_thb_init

    subroutine spline_thb_refine(thb, marked, status, direction)
        !! Refine the marked cells of thb%space and the cells the class-2
        !! admissibility requires. direction (THB_ISOTROPIC default or
        !! THB_RADIAL) applies to a level created by this call. Status 1:
        !! invalid cells or direction, 2: a level could not be built.
        type(spline_thb_t), intent(inout) :: thb
        integer, intent(in) :: marked(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: direction

        logical, allocatable :: flag(:, :)
        integer, allocatable :: stack(:, :)
        integer :: dir, h, m, c, top, nlev, maxcell, q, r, i, pcell, nstack

        status = 1
        dir = THB_ISOTROPIC
        if (present(direction)) dir = direction
        if (dir /= THB_ISOTROPIC .and. dir /= THB_RADIAL) return
        if (dir == THB_RADIAL .and. .not. any(thb%layout%radial_refine)) return
        if (any(marked < 1) .or. any(marked > thb%space%cell_count)) return
        nlev = size(thb%level)
        maxcell = maxval([(thb%level(m)%space%cell_count, m = 1, nlev)])
        allocate(flag(maxcell, nlev), stack(2, max(16, 4*size(marked))))
        flag = .false.
        nstack = 0
        do i = 1, size(marked)
            h = marked(i)
            call push_cell(thb%cell_level(h), thb%cell_index(h))
        end do
        do while (nstack > 0)
            m = stack(1, nstack)
            c = stack(2, nstack)
            nstack = nstack - 1
            if (flag(c, m)) cycle
            flag(c, m) = .true.
            if (m == 1) cycle
            associate (coarse => thb%level(m - 1))
                pcell = thb%level(m)%parent(c)
                do q = coarse%cell_dofs%ptr(pcell), coarse%cell_dofs%ptr(pcell + 1) - 1
                    i = coarse%cell_dofs%col(q)
                    do r = coarse%support%ptr(i), coarse%support%ptr(i + 1) - 1
                        top = coarse%support%col(r)
                        if (.not. coarse%omega(top) .or. coarse%refined(top)) cycle
                        if (flag(top, m - 1)) cycle
                        call push_cell(m - 1, top)
                    end do
                end do
            end associate
        end do
        if (any(flag(:, nlev))) then
            call add_level(thb, dir, status)
            if (status /= 0) then
                status = 2
                return
            end if
        end if
        do m = 1, min(nlev, size(thb%level) - 1)
            associate (lev => thb%level(m), fine => thb%level(m + 1))
                do c = 1, lev%space%cell_count
                    if (.not. flag(c, m)) cycle
                    if (.not. lev%omega(c) .or. lev%refined(c)) cycle
                    lev%refined(c) = .true.
                    do q = lev%children%ptr(c), lev%children%ptr(c + 1) - 1
                        fine%omega(lev%children%col(q)) = .true.
                    end do
                end do
            end associate
        end do
        call assemble(thb, status)

    contains

        subroutine push_cell(level, cell)
            integer, intent(in) :: level, cell
            integer, allocatable :: grown(:, :)

            if (nstack == size(stack, 2)) then
                allocate(grown(2, 2*nstack))
                grown(:, :nstack) = stack
                call move_alloc(grown, stack)
            end if
            nstack = nstack + 1
            stack(:, nstack) = [level, cell]
        end subroutine push_cell

    end subroutine spline_thb_refine

    subroutine spline_thb_mark(indicator, theta, marked)
        !! Dorfler marking: the fewest cells with the largest indicators eta
        !! such that their sum of eta^2 reaches theta times the total.
        real(dp), intent(in) :: indicator(:), theta
        integer, allocatable, intent(out) :: marked(:)

        integer, allocatable :: order(:)
        real(dp) :: total, partial
        integer :: i, j, x, n

        n = size(indicator)
        order = [(i, i = 1, n)]
        do i = 2, n
            x = order(i)
            j = i - 1
            do while (j >= 1)
                if (indicator(order(j)) >= indicator(x)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = x
        end do
        total = sum(indicator**2)
        partial = 0.0_dp
        do i = 1, n
            partial = partial + indicator(order(i))**2
            if (partial >= theta*total) exit
        end do
        marked = order(:min(i, n))
    end subroutine spline_thb_mark

    subroutine spline_project(source, source_coef, target, quad, coef, status)
        !! L2 projection of the fields source_coef(nfield, source dofs) onto
        !! target, both on the same layout, with the target quadrature. Exact
        !! when the source space lies in the target space (refinement).
        !! Status 1: invalid input, 2: the source cannot be evaluated at a
        !! point, 3: the mass solve failed.
        type(spline_space_t), intent(in) :: source, target
        real(dp), intent(in) :: source_coef(:, :)
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(out) :: coef(:, :)
        integer, intent(out) :: status

        type(projection_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        type(pcg_params_t) :: params
        type(pcg_handle_t) :: handle
        type(fortsparse_status_t) :: info
        real(dp), allocatable :: values(:, :), zero(:, :), g(:, :), jet(:), x(:)
        real(dp) :: energy, residual
        integer :: nf, q, f, ierr, iterations

        status = 1
        nf = size(source_coef, 1)
        if (size(coef, 1) /= nf .or. size(coef, 2) /= target%dof_count) return
        if (size(source_coef, 2) /= source%dof_count) return
        allocate(values(nf, quad%point_count), jet(nf))
        status = 2
        do q = 1, quad%point_count
            call spline_evaluate_at(source, source_coef, quad%block(q), quad%xi(:, q), &
                0, jet, ierr)
            if (ierr /= 0) return
            values(:, q) = jet
        end do
        kernel%order = 0
        kernel%g = values(1, :)
        allocate(zero(1, target%dof_count), g(1, target%dof_count))
        zero = 0.0_dp
        status = 3
        call spline_energy_hessian_csr(target, quad, kernel, zero, matrix, ierr)
        if (ierr /= 0) return
        params%preconditioner = PCG_JACOBI
        params%rtol = 1.0e-14_dp
        params%max_iterations = 20*matrix%n + 100
        call pcg_setup_csr(matrix%n, matrix%rowptr, matrix%colidx, matrix%val, params, &
            handle, info)
        if (info%code /= FORTSPARSE_OK) return
        allocate(x(matrix%n))
        do f = 1, nf
            kernel%g = values(f, :)
            call spline_energy(target, quad, kernel, zero, energy, ierr, g)
            if (ierr /= 0) exit
            x = 0.0_dp
            call pcg_solve(handle, params, -g(1, :), x, info, iterations, residual)
            if (info%code /= FORTSPARSE_OK) exit
            coef(f, :) = x
        end do
        call pcg_destroy(handle)
        if (ierr /= 0 .or. info%code /= FORTSPARSE_OK) return
        status = 0
    end subroutine spline_project

    subroutine projection_evaluate(self, cell, first_point, npoint, jet, weight, &
            need, w, gradient, hessian, status)
        class(projection_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status

        integer :: q

        associate (unused => cell, unused_weight => weight)
        end associate
        do q = 1, npoint
            w(q) = 0.5_dp*jet(1, q)**2 - self%g(first_point + q - 1)*jet(1, q)
            if (need >= 1) gradient(1, q) = jet(1, q) - self%g(first_point + q - 1)
            if (need >= 2) hessian(1, q) = 1.0_dp
        end do
        status = 0
    end subroutine projection_evaluate

    subroutine add_level(thb, direction, status)
        !! V_(m+1) from V_m: refined breakpoints, geometry by knot insertion,
        !! parent and child maps, supports and the prolongation V_m -> V_(m+1).
        type(spline_thb_t), intent(inout) :: thb
        integer, intent(in) :: direction
        integer, intent(out) :: status

        type(thb_level_t), allocatable :: levels(:)
        integer :: m, k, d

        m = size(thb%level)
        allocate(levels(m + 1))
        levels(:m) = thb%level
        call move_alloc(levels, thb%level)
        associate (coarse => thb%level(m), fine => thb%level(m + 1))
            fine%direction = direction
            fine%options = coarse%options
            do k = 1, thb%layout%block_count
                do d = 1, 2
                    if (thb%layout%periodic(k) .and. d == 1) cycle
                    if (direction == THB_RADIAL .and. &
                        .not. thb%layout%radial_refine(d, k)) cycle
                    call refine_breakpoints(fine%options%knots(d, k)%breakpoints)
                end do
            end do
            call spline_space_build(thb%layout, fine%options, fine%space, status)
            if (status /= 0) return
            call insert_geometry(coarse%space, fine%space, thb%layout, status)
            if (status /= 0) return
            allocate(fine%omega(fine%space%cell_count), &
                fine%refined(fine%space%cell_count))
            fine%omega = .false.
            fine%refined = .false.
            call parent_map(coarse%space, fine%space, fine%parent, coarse%children)
            call level_tables(fine)
            call build_prolongation(coarse, fine, thb%layout, status)
        end associate
    end subroutine add_level

    subroutine refine_breakpoints(breaks)
        real(dp), allocatable, intent(inout) :: breaks(:)
        real(dp), allocatable :: fine(:)
        integer :: i, n

        n = size(breaks)
        allocate(fine(2*n - 1))
        do i = 1, n - 1
            fine(2*i - 1) = breaks(i)
            fine(2*i) = 0.5_dp*(breaks(i) + breaks(i + 1))
        end do
        fine(2*n - 1) = breaks(n)
        call move_alloc(fine, breaks)
    end subroutine refine_breakpoints

    subroutine level_tables(lev)
        !! Rows of E by function, functions per cell, and supports.
        type(thb_level_t), intent(inout) :: lev

        integer, allocatable :: broken(:), counts(:), seen(:)
        integer :: p, cell, a, q, i, n, nloc, pass

        associate (space => lev%space)
            call rows_from_broken(space, lev%dof_rows)
            p = space%degree
            nloc = (p + 1)**2
            allocate(broken(nloc), seen(space%dof_count), counts(space%cell_count + 1))
            allocate(lev%cell_dofs%ptr(space%cell_count + 1))
            do pass = 1, 2
                seen = 0
                lev%cell_dofs%ptr(1) = 1
                do cell = 1, space%cell_count
                    call spline_cell_broken(space, cell, broken)
                    n = 0
                    do a = 1, nloc
                        do q = space%broken_ptr(broken(a)), &
                                space%broken_ptr(broken(a) + 1) - 1
                            i = space%broken_dof(q)
                            if (space%broken_weight(q) == 0.0_dp .or. seen(i) == cell) &
                                cycle
                            seen(i) = cell
                            if (pass == 2) lev%cell_dofs%col(lev%cell_dofs%ptr(cell) + n) = i
                            n = n + 1
                        end do
                    end do
                    lev%cell_dofs%ptr(cell + 1) = lev%cell_dofs%ptr(cell) + n
                end do
                if (pass == 1) allocate(lev%cell_dofs%col(lev%cell_dofs%ptr( &
                    space%cell_count + 1) - 1))
            end do
            call transpose_rows(lev%cell_dofs, space%dof_count, lev%support)
        end associate
    end subroutine level_tables

    subroutine rows_from_broken(space, rows)
        type(spline_space_t), intent(in) :: space
        type(rows_t), intent(out) :: rows

        integer, allocatable :: next(:)
        integer :: j, q, i

        allocate(rows%ptr(space%dof_count + 1), next(space%dof_count + 1))
        rows%ptr = 0
        do q = 1, size(space%broken_dof)
            rows%ptr(space%broken_dof(q) + 1) = rows%ptr(space%broken_dof(q) + 1) + 1
        end do
        rows%ptr(1) = 1
        do i = 2, space%dof_count + 1
            rows%ptr(i) = rows%ptr(i) + rows%ptr(i - 1)
        end do
        allocate(rows%col(size(space%broken_dof)), rows%val(size(space%broken_dof)))
        next = rows%ptr
        do j = 1, space%broken_count
            do q = space%broken_ptr(j), space%broken_ptr(j + 1) - 1
                i = space%broken_dof(q)
                rows%col(next(i)) = j
                rows%val(next(i)) = space%broken_weight(q)
                next(i) = next(i) + 1
            end do
        end do
    end subroutine rows_from_broken

    subroutine transpose_rows(a, ncol, at)
        type(rows_t), intent(in) :: a
        integer, intent(in) :: ncol
        type(rows_t), intent(out) :: at

        integer, allocatable :: next(:)
        integer :: i, q, j
        logical :: values

        values = allocated(a%val)
        allocate(at%ptr(ncol + 1), next(ncol + 1))
        at%ptr = 0
        do q = 1, size(a%col)
            at%ptr(a%col(q) + 1) = at%ptr(a%col(q) + 1) + 1
        end do
        at%ptr(1) = 1
        do j = 2, ncol + 1
            at%ptr(j) = at%ptr(j) + at%ptr(j - 1)
        end do
        allocate(at%col(size(a%col)))
        if (values) allocate(at%val(size(a%col)))
        next = at%ptr
        do i = 1, size(a%ptr) - 1
            do q = a%ptr(i), a%ptr(i + 1) - 1
                j = a%col(q)
                at%col(next(j)) = i
                if (values) at%val(next(j)) = a%val(q)
                next(j) = next(j) + 1
            end do
        end do
    end subroutine transpose_rows

    subroutine parent_map(coarse, fine, parent, children)
        !! Coarser cell of each fine cell (by the element mid point) and the
        !! inverse map.
        type(spline_space_t), intent(in) :: coarse, fine
        integer, allocatable, intent(out) :: parent(:)
        type(rows_t), intent(out) :: children

        type(rows_t) :: up
        real(dp) :: x(2)
        integer :: cell, k, e(2), ec(2)

        allocate(parent(fine%cell_count), up%ptr(fine%cell_count + 1), &
            up%col(fine%cell_count))
        do cell = 1, fine%cell_count
            k = fine%cell_block(cell)
            e = fine%cell_element(:, cell)
            x(1) = 0.5_dp*sum(fine%blocks(k)%breaks1(e(1):e(1) + 1))
            x(2) = 0.5_dp*sum(fine%blocks(k)%breaks2(e(2):e(2) + 1))
            ec(1) = locate(coarse%blocks(k)%breaks1, x(1))
            ec(2) = locate(coarse%blocks(k)%breaks2, x(2))
            parent(cell) = coarse%element_cell(coarse%patch_cell_offset(k) + ec(1) + &
                coarse%blocks(k)%element_count(1)*(ec(2) - 1))
            up%ptr(cell) = cell
            up%col(cell) = parent(cell)
        end do
        up%ptr(fine%cell_count + 1) = fine%cell_count + 1
        call transpose_rows(up, coarse%cell_count, children)
    end subroutine parent_map

    subroutine insert_geometry(coarse, fine, layout, status)
        !! Fine broken control points of the coarse block maps.
        type(spline_space_t), intent(in) :: coarse
        type(spline_space_t), intent(inout) :: fine
        type(block_layout_t), intent(in) :: layout
        integer, intent(out) :: status

        real(dp), allocatable :: a1(:, :), a2(:, :), g(:, :, :), h(:, :, :)
        integer :: k, nc(2), nf(2), i, j, d

        do k = 1, layout%block_count
            call insertion_1d(coarse%blocks(k), fine%blocks(k), 1, layout%periodic(k), &
                coarse%degree, a1, status)
            if (status /= 0) return
            call insertion_1d(coarse%blocks(k), fine%blocks(k), 2, .false., &
                coarse%degree, a2, status)
            if (status /= 0) return
            nc = [size(a1, 2), size(a2, 2)]
            nf = [size(a1, 1), size(a2, 1)]
            allocate(g(2, nc(1), nc(2)), h(2, nf(1), nf(2)))
            do j = 1, nc(2)
                do i = 1, nc(1)
                    g(:, i, j) = coarse%geometry(:, coarse%blocks(k)%offset + i + &
                        coarse%blocks(k)%n(1)*(j - 1))
                end do
            end do
            do d = 1, 2
                h(d, :, :) = matmul(a1, matmul(g(d, :, :), transpose(a2)))
            end do
            do j = 1, fine%blocks(k)%n(2)
                do i = 1, fine%blocks(k)%n(1)
                    fine%geometry(:, fine%blocks(k)%offset + i + &
                        fine%blocks(k)%n(1)*(j - 1)) = h(:, modulo(i - 1, nf(1)) + 1, j)
                end do
            end do
            deallocate(g, h)
        end do
    end subroutine insert_geometry

    subroutine insertion_1d(coarse, fine, direction, periodic, degree, a, status)
        !! a(f, i): coefficient of the fine function f in the coarse function
        !! i along one block direction (periodic functions folded), by least
        !! squares at the fine Gauss points (exact, the spaces are nested).
        type(spline_block_t), intent(in) :: coarse, fine
        integer, intent(in) :: direction, degree
        logical, intent(in) :: periodic
        real(dp), allocatable, intent(out) :: a(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: x(:), w(:), vf(:, :), vc(:, :), gram(:, :), breaks(:)
        integer, allocatable :: pivot(:)
        integer :: nfun_f, nfun_c, e, q, row, nq

        nq = degree + 1
        call gauss_rule(nq, x, w, status)
        if (status /= 0) return
        nfun_f = folded_count(fine, direction, periodic)
        nfun_c = folded_count(coarse, direction, periodic)
        if (direction == 1) then
            breaks = fine%breaks1
        else
            breaks = fine%breaks2
        end if
        allocate(vf(nq*(size(breaks) - 1), nfun_f), vc(nq*(size(breaks) - 1), nfun_c))
        row = 0
        do e = 1, size(breaks) - 1
            do q = 1, nq
                row = row + 1
                call direction_values(fine, direction, periodic, degree, breaks(e) + &
                    0.5_dp*(x(q) + 1.0_dp)*(breaks(e + 1) - breaks(e)), vf(row, :))
                call direction_values(coarse, direction, periodic, degree, breaks(e) + &
                    0.5_dp*(x(q) + 1.0_dp)*(breaks(e + 1) - breaks(e)), vc(row, :))
            end do
        end do
        gram = matmul(transpose(vf), vf)
        a = matmul(transpose(vf), vc)
        allocate(pivot(nfun_f))
        call dgesv(nfun_f, nfun_c, gram, nfun_f, pivot, a, nfun_f, status)
        where (abs(a) < 1.0e-14_dp) a = 0.0_dp
    end subroutine insertion_1d

    pure integer function folded_count(blk, direction, periodic) result(n)
        type(spline_block_t), intent(in) :: blk
        integer, intent(in) :: direction
        logical, intent(in) :: periodic

        if (periodic .and. direction == 1) then
            n = blk%element_count(1)
        else
            n = blk%n(direction)
        end if
    end function folded_count

    subroutine direction_values(blk, direction, periodic, degree, x, values)
        !! All (folded) functions of one block direction at x.
        type(spline_block_t), intent(in) :: blk
        integer, intent(in) :: direction, degree
        logical, intent(in) :: periodic
        real(dp), intent(in) :: x
        real(dp), intent(out) :: values(:)

        real(dp) :: b(degree + 1), d1(degree + 1), d2(degree + 1), local(degree + 1), t
        integer :: e, a, i, n, ierr

        values = 0.0_dp
        n = size(values)
        if (direction == 1) then
            e = locate(blk%breaks1, x)
            t = (x - blk%breaks1(e))/(blk%breaks1(e + 1) - blk%breaks1(e))
            call interval_bernstein_jet(degree, t, b, d1, d2, ierr)
            local = matmul(blk%extraction1(:, :, e), b)
            do a = 1, degree + 1
                i = blk%first1(e) + a - 1
                if (periodic) i = modulo(i - 1, n) + 1
                values(i) = values(i) + local(a)
            end do
        else
            e = locate(blk%breaks2, x)
            t = (x - blk%breaks2(e))/(blk%breaks2(e + 1) - blk%breaks2(e))
            call interval_bernstein_jet(degree, t, b, d1, d2, ierr)
            local = matmul(blk%extraction2(:, :, e), b)
            do a = 1, degree + 1
                i = blk%first2(e) + a - 1
                values(i) = values(i) + local(a)
            end do
        end if
    end subroutine direction_values

    subroutine build_prolongation(coarse, fine, layout, status)
        !! coarse%prolong(I, K): coefficient of the fine function K in the
        !! coarse function I. Each fine function is read off one of its
        !! broken functions whose extraction column holds it alone.
        type(thb_level_t), intent(inout) :: coarse
        type(thb_level_t), intent(in) :: fine
        type(block_layout_t), intent(in) :: layout
        integer, intent(out) :: status

        type :: matrix_t
            real(dp), allocatable :: a(:, :)
        end type matrix_t
        type(rows_t) :: by_fine
        type(matrix_t), allocatable :: a1(:), a2(:)
        integer, allocatable :: rep(:), touched(:), owner(:)
        real(dp), allocatable :: rep_w(:), acc(:)
        integer :: j, kf, k, i1, i2, a, b, jc, q, ntouch, nloc1, i, block, nnz

        status = 2
        associate (fs => fine%space, cs => coarse%space)
            allocate(rep(fs%dof_count), rep_w(fs%dof_count), owner(fs%broken_count))
            rep = 0
            do block = 1, fs%block_count
                owner(fs%blocks(block)%offset + 1:fs%blocks(block)%offset + &
                    product(fs%blocks(block)%n)) = block
            end do
            do j = 1, fs%broken_count
                if (fs%broken_ptr(j + 1) - fs%broken_ptr(j) /= 1) cycle
                kf = fs%broken_dof(fs%broken_ptr(j))
                if (rep(kf) /= 0) cycle
                rep(kf) = j
                rep_w(kf) = fs%broken_weight(fs%broken_ptr(j))
            end do
            if (any(rep == 0)) return
            allocate(a1(layout%block_count), a2(layout%block_count))
            do k = 1, layout%block_count
                call insertion_1d(cs%blocks(k), fs%blocks(k), 1, layout%periodic(k), &
                    cs%degree, a1(k)%a, status)
                if (status /= 0) return
                call insertion_1d(cs%blocks(k), fs%blocks(k), 2, .false., cs%degree, &
                    a2(k)%a, status)
                if (status /= 0) return
            end do
            allocate(acc(cs%dof_count), touched(cs%dof_count))
            acc = 0.0_dp
            nnz = 0
            allocate(by_fine%ptr(fs%dof_count + 1), by_fine%col(1024), &
                by_fine%val(1024))
            by_fine%ptr(1) = 1
            do kf = 1, fs%dof_count
                j = rep(kf)
                k = owner(j)
                nloc1 = fs%blocks(k)%n(1)
                i1 = modulo(j - fs%blocks(k)%offset - 1, nloc1) + 1
                i2 = (j - fs%blocks(k)%offset - 1)/nloc1 + 1
                if (layout%periodic(k)) i1 = modulo(i1 - 1, size(a1(k)%a, 1)) + 1
                ntouch = 0
                do b = 1, size(a2(k)%a, 2)
                    if (a2(k)%a(i2, b) == 0.0_dp) cycle
                    do a = 1, size(a1(k)%a, 2)
                        if (a1(k)%a(i1, a) == 0.0_dp) cycle
                        jc = cs%blocks(k)%offset + a + cs%blocks(k)%n(1)*(b - 1)
                        do q = cs%broken_ptr(jc), cs%broken_ptr(jc + 1) - 1
                            i = cs%broken_dof(q)
                            if (acc(i) == 0.0_dp) then
                                ntouch = ntouch + 1
                                touched(ntouch) = i
                            end if
                            acc(i) = acc(i) + a1(k)%a(i1, a)*a2(k)%a(i2, b)* &
                                cs%broken_weight(q)/rep_w(kf)
                            if (acc(i) == 0.0_dp) acc(i) = tiny(1.0_dp)
                        end do
                    end do
                end do
                do q = 1, ntouch
                    i = touched(q)
                    if (abs(acc(i)) > 1.0e-14_dp) then
                        if (nnz == size(by_fine%col)) call grow(by_fine)
                        nnz = nnz + 1
                        by_fine%col(nnz) = i
                        by_fine%val(nnz) = acc(i)
                    end if
                    acc(i) = 0.0_dp
                end do
                by_fine%ptr(kf + 1) = nnz + 1
            end do
            by_fine%col = by_fine%col(:nnz)
            by_fine%val = by_fine%val(:nnz)
            call transpose_rows(by_fine, cs%dof_count, coarse%prolong)
        end associate
        status = 0
    end subroutine build_prolongation

    subroutine grow(rows)
        type(rows_t), intent(inout) :: rows
        integer, allocatable :: col(:)
        real(dp), allocatable :: val(:)

        allocate(col(2*size(rows%col)), val(2*size(rows%col)))
        col(:size(rows%col)) = rows%col
        val(:size(rows%col)) = rows%val
        call move_alloc(col, rows%col)
        call move_alloc(val, rows%val)
    end subroutine grow

    subroutine assemble(thb, status)
        !! Hierarchical space from the level spaces and the domains Omega_m.
        type(spline_thb_t), intent(inout) :: thb
        integer, intent(out) :: status

        type(spline_space_t) :: s
        integer, allocatable :: level_offset(:), cell_offset(:), hier(:)
        integer, allocatable :: tri_dof(:), tri_broken(:), vec(:), nvec(:)
        integer, allocatable :: btouch(:), counts(:)
        real(dp), allocatable :: tri_w(:), val(:), nval(:), bacc(:), work_acc(:)
        integer, allocatable :: work_touched(:)
        logical, allocatable :: bactive(:)
        integer :: nlev, nb, m, n, k, c, i, q, r, f, cell, ntri, nbt, kk, j
        integer :: total_broken, total_cells, broken(64), nloc, pos

        nlev = size(thb%level)
        nb = thb%layout%block_count
        allocate(level_offset(nlev + 1), cell_offset(nlev + 1))
        level_offset(1) = 0
        cell_offset(1) = 0
        do m = 1, nlev
            level_offset(m + 1) = level_offset(m) + thb%level(m)%space%broken_count
            cell_offset(m + 1) = cell_offset(m) + thb%level(m)%space%cell_count
        end do
        total_broken = level_offset(nlev + 1)
        total_cells = cell_offset(nlev + 1)
        nloc = (thb%level(1)%space%degree + 1)**2
        allocate(hier(total_cells), bactive(total_broken))
        hier = 0
        bactive = .false.
        n = 0
        do m = 1, nlev
            associate (lev => thb%level(m))
                do c = 1, lev%space%cell_count
                    if (.not. lev%omega(c) .or. lev%refined(c)) cycle
                    n = n + 1
                    hier(cell_offset(m) + c) = n
                    call spline_cell_broken(lev%space, c, broken(:nloc))
                    bactive(level_offset(m) + broken(:nloc)) = .true.
                end do
            end associate
        end do
        s%cell_count = n
        if (allocated(thb%cell_level)) deallocate(thb%cell_level, thb%cell_index)
        allocate(thb%cell_level(n), thb%cell_index(n))
        do m = 1, nlev
            do c = 1, thb%level(m)%space%cell_count
                if (hier(cell_offset(m) + c) == 0) cycle
                thb%cell_level(hier(cell_offset(m) + c)) = m
                thb%cell_index(hier(cell_offset(m) + c)) = c
            end do
        end do
        allocate(tri_dof(1024), tri_broken(1024), tri_w(1024), bacc(total_broken), &
            btouch(total_broken))
        allocate(work_acc(maxval([(thb%level(m)%space%dof_count, m = 1, nlev)])), &
            work_touched(maxval([(thb%level(m)%space%dof_count, m = 1, nlev)])))
        work_acc = 0.0_dp
        bacc = 0.0_dp
        ntri = 0
        f = 0
        do m = 1, nlev
            associate (lev => thb%level(m))
                do i = 1, lev%space%dof_count
                    if (.not. function_active(lev, i)) cycle
                    f = f + 1
                    nbt = 0
                    vec = [i]
                    val = [1.0_dp]
                    do kk = m, nlev
                        associate (cur => thb%level(kk))
                            do q = 1, size(vec)
                                do r = cur%dof_rows%ptr(vec(q)), &
                                        cur%dof_rows%ptr(vec(q) + 1) - 1
                                    j = level_offset(kk) + cur%dof_rows%col(r)
                                    if (.not. bactive(j)) cycle
                                    if (bacc(j) == 0.0_dp) then
                                        nbt = nbt + 1
                                        btouch(nbt) = j
                                    end if
                                    bacc(j) = bacc(j) + val(q)*cur%dof_rows%val(r)
                                    if (bacc(j) == 0.0_dp) bacc(j) = tiny(1.0_dp)
                                end do
                            end do
                            if (kk == nlev) exit
                            call truncate_step(cur, thb%level(kk + 1), vec, val, nvec, &
                                nval, work_acc, work_touched)
                            call move_alloc(nvec, vec)
                            call move_alloc(nval, val)
                            if (size(vec) == 0) exit
                        end associate
                    end do
                    do q = 1, nbt
                        j = btouch(q)
                        if (abs(bacc(j)) > 1.0e-14_dp) call push(f, j, bacc(j))
                        bacc(j) = 0.0_dp
                    end do
                end do
            end associate
        end do
        s%degree = thb%level(1)%space%degree
        s%regularity = SPLINE_C0
        s%layout = thb%layout
        s%layout_block_count = nb
        s%block_count = nb*nlev
        s%broken_count = total_broken
        s%dof_count = f
        allocate(s%blocks(nb*nlev), s%patch_block(nb*nlev), s%patch_level(nb*nlev), &
            s%patch_cell_offset(nb*nlev), s%element_cell(total_cells), &
            s%geometry(2, total_broken), s%cell_block(n), s%cell_element(2, n))
        do m = 1, nlev
            associate (lev => thb%level(m))
                do k = 1, nb
                    pos = k + nb*(m - 1)
                    s%blocks(pos) = lev%space%blocks(k)
                    s%blocks(pos)%offset = level_offset(m) + lev%space%blocks(k)%offset
                    s%patch_block(pos) = k
                    s%patch_level(pos) = m - 1
                    s%patch_cell_offset(pos) = cell_offset(m) + &
                        lev%space%patch_cell_offset(k)
                end do
                s%geometry(:, level_offset(m) + 1:level_offset(m + 1)) = lev%space%geometry
                do c = 1, lev%space%cell_count
                    cell = hier(cell_offset(m) + c)
                    s%element_cell(cell_offset(m) + c) = cell
                    if (cell == 0) cycle
                    s%cell_block(cell) = lev%space%cell_block(c) + nb*(m - 1)
                    s%cell_element(:, cell) = lev%space%cell_element(:, c)
                end do
            end associate
        end do
        allocate(counts(total_broken + 1))
        counts = 0
        do q = 1, ntri
            counts(tri_broken(q) + 1) = counts(tri_broken(q) + 1) + 1
        end do
        counts(1) = 1
        do j = 2, total_broken + 1
            counts(j) = counts(j) + counts(j - 1)
        end do
        s%broken_ptr = counts
        allocate(s%broken_dof(ntri), s%broken_weight(ntri))
        do q = 1, ntri
            pos = counts(tri_broken(q))
            counts(tri_broken(q)) = counts(tri_broken(q)) + 1
            s%broken_dof(pos) = tri_dof(q)
            s%broken_weight(pos) = tri_w(q)
        end do
        s%c0_interfaces = thb%level(1)%space%c0_interfaces
        thb%space = s
        status = 0

    contains

        subroutine push(dof, b, w)
            integer, intent(in) :: dof, b
            real(dp), intent(in) :: w
            integer, allocatable :: gi(:), gb(:)
            real(dp), allocatable :: gw(:)

            if (ntri == size(tri_dof)) then
                allocate(gi(2*ntri), gb(2*ntri), gw(2*ntri))
                gi(:ntri) = tri_dof
                gb(:ntri) = tri_broken
                gw(:ntri) = tri_w
                call move_alloc(gi, tri_dof)
                call move_alloc(gb, tri_broken)
                call move_alloc(gw, tri_w)
            end if
            ntri = ntri + 1
            tri_dof(ntri) = dof
            tri_broken(ntri) = b
            tri_w(ntri) = w
        end subroutine push

    end subroutine assemble

    logical function function_active(lev, i) result(active)
        !! Support in Omega_m and not in Omega_(m+1).
        type(thb_level_t), intent(in) :: lev
        integer, intent(in) :: i

        integer :: q, c
        logical :: all_refined

        active = .false.
        all_refined = .true.
        do q = lev%support%ptr(i), lev%support%ptr(i + 1) - 1
            c = lev%support%col(q)
            if (.not. lev%omega(c)) return
            if (.not. lev%refined(c)) all_refined = .false.
        end do
        active = .not. all_refined
    end function function_active

    subroutine truncate_step(coarse, fine, vec, val, nvec, nval, acc, touched)
        !! Fine coefficients of sum val(q) phi_vec(q) without the functions
        !! whose support lies in Omega of the fine level. acc (zero on entry
        !! and exit) and touched are work arrays of at least the fine size.
        type(thb_level_t), intent(in) :: coarse, fine
        integer, intent(in) :: vec(:)
        real(dp), intent(in) :: val(:)
        integer, allocatable, intent(out) :: nvec(:)
        real(dp), allocatable, intent(out) :: nval(:)
        real(dp), intent(inout) :: acc(:)
        integer, intent(inout) :: touched(:)

        integer :: q, r, k, nt, c, n
        logical :: inside

        nt = 0
        do q = 1, size(vec)
            do r = coarse%prolong%ptr(vec(q)), coarse%prolong%ptr(vec(q) + 1) - 1
                k = coarse%prolong%col(r)
                if (acc(k) == 0.0_dp) then
                    nt = nt + 1
                    touched(nt) = k
                end if
                acc(k) = acc(k) + val(q)*coarse%prolong%val(r)
                if (acc(k) == 0.0_dp) acc(k) = tiny(1.0_dp)
            end do
        end do
        allocate(nvec(nt), nval(nt))
        n = 0
        do q = 1, nt
            k = touched(q)
            inside = .true.
            do r = fine%support%ptr(k), fine%support%ptr(k + 1) - 1
                c = fine%support%col(r)
                if (.not. fine%omega(c)) then
                    inside = .false.
                    exit
                end if
            end do
            if (.not. inside .and. abs(acc(k)) > 1.0e-14_dp) then
                n = n + 1
                nvec(n) = k
                nval(n) = acc(k)
            end if
            acc(k) = 0.0_dp
        end do
        nvec = nvec(:n)
        nval = nval(:n)
    end subroutine truncate_step

    pure integer function locate(breaks, x) result(e)
        real(dp), intent(in) :: breaks(:), x

        do e = 1, size(breaks) - 2
            if (x < breaks(e + 1)) return
        end do
        e = size(breaks) - 1
    end function locate

end module fortfem_spline_thb
