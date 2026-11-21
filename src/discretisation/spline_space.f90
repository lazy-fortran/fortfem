module fortfem_spline_space
    !! Multipatch spline spaces on a block layout
    !! (doc/design/multipatch_spline_spaces.md).
    !!
    !! Every block carries tensor B-splines on open knot vectors. The broken
    !! space is the direct sum of these patch spaces, numbered block by block
    !! with xi1 fastest. The global space is fixed by the sparse extraction
    !! E, stored by broken function: broken_dof(broken_ptr(j):broken_ptr(j+1)-1)
    !! are the global functions I with weights E(I, j), so that the broken
    !! coefficients are E^T c. The block maps (1) of the design note are kept
    !! as broken control points.
    !!
    !! A periodic block (the annulus) carries, along xi1, the n + p unclamped
    !! B-splines of the periodic extension, and E identifies the functions i
    !! and i + n. Its xi1 breakpoints are the images of the breakpoints of
    !! the core sides tagged LAYOUT_ANNULUS_INTERFACE, and the core traces
    !! there are the restrictions of the annulus traces (C0 coupling).
    use fortfem_kinds, only: dp
    use fortfem_spline_extraction_1d, only: build_graded_breakpoints, &
        build_open_knot_vector, build_bezier_extraction_1d, &
        build_periodic_extraction_1d, SPLINE_GRID_UNIFORM
    use fortfem_spline_layout, only: block_layout_t, layout_side_point, &
        LAYOUT_ANNULUS_INTERFACE
    use fortfem_multipatch_dof_graph, only: build_multipatch_signed_dof_map
    use fortfem_spline_c1, only: c1_block_knots_t, build_c1_extraction
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    implicit none
    private

    integer, parameter, public :: SPLINE_C0 = 0
    integer, parameter, public :: SPLINE_C1 = 1

    type, public :: knot_spec_t
        !! Breakpoints on [0, 1]: either element_count graded cells or explicit
        !! breakpoints. multiplicity (interior entries) lowers the regularity
        !! at chosen breakpoints. reversed mirrors a grading, s -> 1 - s.
        integer :: element_count = 1
        integer :: grading = SPLINE_GRID_UNIFORM
        logical :: reversed = .false.
        real(dp), allocatable :: breakpoints(:)
        integer, allocatable :: multiplicity(:)
    end type knot_spec_t

    type, public :: spline_space_options_t
        integer :: degree = 3
        integer :: regularity = SPLINE_C0
        type(knot_spec_t), allocatable :: knots(:, :) ! (2, block)
        integer, allocatable :: block_degree(:, :) ! (2, block), optional
    end type spline_space_options_t

    type, public :: spline_block_t
        integer :: n(2) = 0
        integer :: element_count(2) = 0
        integer :: offset = 0 ! broken index of the first function - 1
        real(dp), allocatable :: knots1(:), knots2(:)
        real(dp), allocatable :: breaks1(:), breaks2(:)
        real(dp), allocatable :: extraction1(:, :, :), extraction2(:, :, :)
        integer, allocatable :: first1(:), first2(:)
    end type spline_block_t

    type, public :: spline_space_t
        integer :: degree = 0
        integer :: regularity = SPLINE_C0
        integer :: block_count = 0
        integer :: broken_count = 0
        integer :: dof_count = 0
        integer :: cell_count = 0
        type(block_layout_t) :: layout
        type(spline_block_t), allocatable :: blocks(:)
        integer, allocatable :: cell_block(:) ! (cell)
        integer, allocatable :: cell_element(:, :) ! (2, cell)
        integer, allocatable :: broken_ptr(:) ! (broken + 1)
        integer, allocatable :: broken_dof(:)
        real(dp), allocatable :: broken_weight(:)
        real(dp), allocatable :: geometry(:, :) ! (2, broken)
        !> C0 interfaces (block, side, neighbour block, neighbour side): all
        !> interior sides of a C0 space, the side_c0 sides of a C1 space, and
        !> the core sides on an annulus interface (neighbour side 3).
        integer, allocatable :: c0_interfaces(:, :) ! (4, interface)
    end type spline_space_t

    public :: spline_space_build
    public :: spline_boundary_dofs
    public :: spline_cell_broken
    public :: spline_greville
    public :: spline_space_update_geometry

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

contains

    subroutine spline_space_build(layout, options, space, status)
        !! Status 1: invalid options, 2: knot construction failed, 3: knot
        !! vectors of two blocks disagree along a shared side, 4: unknown
        !! regularity, 5: geometry construction failed, 6: the C1 space needs
        !! bilinear blocks, degree >= 3 and interior multiplicity >= 2,
        !! 7: an interface is not analysis-suitable G1.
        !! C1 spaces default to interior multiplicity 2 (C^(p-2) in blocks).
        !! The xi1 knots of a periodic block follow from its interface, or,
        !! without one, from options%knots(1, k) without multiplicity.
        type(block_layout_t), intent(in) :: layout
        type(spline_space_options_t), intent(in) :: options
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status

        integer :: k, d, offset, cell, e1, e2, nb, pass
        logical :: periodic

        status = 1
        nb = layout%block_count
        if (nb < 1 .or. options%degree < 1) return
        if (.not. allocated(options%knots)) return
        if (any(shape(options%knots) /= [2, nb])) return
        if (allocated(options%block_degree)) then
            if (any(options%block_degree /= options%degree)) return
        end if
        status = 4
        if (options%regularity /= SPLINE_C0 .and. options%regularity /= SPLINE_C1) &
            return
        space%degree = options%degree
        space%regularity = options%regularity
        space%block_count = nb
        space%layout = layout
        allocate(space%blocks(nb))
        do pass = 1, 2
            do k = 1, nb
                periodic = layout%periodic(k)
                if (periodic .neqv. pass == 2) cycle
                do d = 1, 2
                    if (periodic .and. d == 1) then
                        call build_periodic_direction(space, k, options%knots(1, k), &
                            status)
                    else
                        call build_block_direction(options%knots(d, k), &
                            options%degree, merge(2, 1, options%regularity == &
                            SPLINE_C1), space%blocks(k), d, status)
                    end if
                    if (status /= 0) then
                        status = 2
                        return
                    end if
                end do
            end do
        end do
        offset = 0
        do k = 1, nb
            space%blocks(k)%offset = offset
            offset = offset + product(space%blocks(k)%n)
        end do
        space%broken_count = offset
        status = 3
        if (.not. sides_conform(space)) return
        if (space%regularity == SPLINE_C1) then
            call build_c1_space(space, status)
        else
            call build_c0_extraction(space, status)
            if (status == 0 .and. any(layout%periodic)) &
                call constrain_periodic(space, status)
        end if
        if (status /= 0) return
        space%cell_count = sum([(product(space%blocks(k)%element_count), &
            k = 1, nb)])
        allocate(space%cell_block(space%cell_count), &
            space%cell_element(2, space%cell_count))
        cell = 0
        do k = 1, nb
            do e2 = 1, space%blocks(k)%element_count(2)
                do e1 = 1, space%blocks(k)%element_count(1)
                    cell = cell + 1
                    space%cell_block(cell) = k
                    space%cell_element(:, cell) = [e1, e2]
                end do
            end do
        end do
        call build_geometry(space, status)
        if (status /= 0) then
            status = 5
            return
        end if
        call list_c0_interfaces(space)
        status = 0
    end subroutine spline_space_build

    subroutine list_c0_interfaces(space)
        type(spline_space_t), intent(inout) :: space

        integer :: k, s, k2, s2, n, pass, ka

        ka = 0
        if (any(space%layout%periodic)) ka = findloc(space%layout%periodic, .true., 1)
        do pass = 1, 2
            n = 0
            do k = 1, space%block_count
                do s = 1, 4
                    k2 = space%layout%neighbour(s, k)
                    s2 = space%layout%neighbour_side(s, k)
                    if (k2 /= 0) then
                        if (k2 < k .or. (k2 == k .and. s2 < s)) cycle
                        if (space%regularity == SPLINE_C1 .and. &
                            .not. space%layout%side_c0(s, k)) cycle
                    else if (space%layout%side_tag(s, k) == &
                            LAYOUT_ANNULUS_INTERFACE .and. k /= ka .and. ka /= 0) then
                        k2 = ka
                        s2 = 3
                    else
                        cycle
                    end if
                    n = n + 1
                    if (pass == 2) space%c0_interfaces(:, n) = [k, s, k2, s2]
                end do
            end do
            if (pass == 1) allocate(space%c0_interfaces(4, n))
        end do
    end subroutine list_c0_interfaces

    subroutine spline_space_update_geometry(space, layout, status)
        !! Replace the block maps after the layout geometry changed (moved
        !! vertices or curves, e.g. by layout_set_curve) on the same topology.
        !! C0 spaces keep their extraction, so coefficients keep their meaning.
        !! C1 spaces rebuild it (the gluing data follow the vertices) with the
        !! same dimension. Quadrature data must be rebuilt.
        !! Status 1: the topology differs, 5: geometry construction failed.
        type(spline_space_t), intent(inout) :: space
        type(block_layout_t), intent(in) :: layout
        integer, intent(out) :: status

        status = 1
        if (layout%block_count /= space%layout%block_count) return
        if (any(layout%corners /= space%layout%corners)) return
        if (any(layout%side_tag /= space%layout%side_tag)) return
        if (any(layout%neighbour /= space%layout%neighbour)) return
        if (any(layout%periodic .neqv. space%layout%periodic)) return
        if (any(abs(layout%side_u0 - space%layout%side_u0) > 0.0_dp)) return
        if (any(abs(layout%side_du - space%layout%side_du) > 0.0_dp)) return
        space%layout = layout
        if (space%regularity == SPLINE_C1) then
            deallocate(space%broken_ptr, space%broken_dof, space%broken_weight)
            call build_c1_space(space, status)
            if (status /= 0) return
        end if
        deallocate(space%geometry)
        call build_geometry(space, status)
        if (status /= 0) status = 5
    end subroutine spline_space_update_geometry

    subroutine spline_cell_broken(space, cell, broken)
        !! Broken indices of the (p + 1)^2 patch functions on a cell, local
        !! index a + (p + 1)(b - 1) with a along xi1.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: cell
        integer, intent(out) :: broken(:)

        integer :: k, a, b, i, j, p

        p = space%degree
        k = space%cell_block(cell)
        associate (blk => space%blocks(k))
            do b = 1, p + 1
                j = blk%first2(space%cell_element(2, cell)) + b - 1
                do a = 1, p + 1
                    i = blk%first1(space%cell_element(1, cell)) + a - 1
                    broken(a + (p + 1)*(b - 1)) = blk%offset + i + blk%n(1)*(j - 1)
                end do
            end do
        end associate
    end subroutine spline_cell_broken

    subroutine spline_boundary_dofs(space, tag, dofs, field, field_count)
        !! Global functions that do not vanish on the sides with the given tag.
        !! With field and field_count the indices refer to coefficients stored
        !! as coef(field_count, dof_count).
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: tag
        integer, allocatable, intent(out) :: dofs(:)
        integer, intent(in), optional :: field, field_count

        logical, allocatable :: mark(:)
        integer, allocatable :: trace(:)
        integer :: k, s, t, j, q

        allocate(mark(space%dof_count))
        mark = .false.
        do k = 1, space%block_count
            do s = 1, 4
                if (space%layout%side_tag(s, k) /= tag) cycle
                call side_trace(space%blocks(k), s, trace)
                do t = 1, size(trace)
                    j = space%blocks(k)%offset + trace(t)
                    do q = space%broken_ptr(j), space%broken_ptr(j + 1) - 1
                        if (space%broken_weight(q) /= 0.0_dp) &
                            mark(space%broken_dof(q)) = .true.
                    end do
                end do
            end do
        end do
        dofs = pack([(j, j = 1, space%dof_count)], mark)
        if (present(field) .and. present(field_count)) &
            dofs = field + field_count*(dofs - 1)
    end subroutine spline_boundary_dofs

    pure subroutine spline_greville(knots, degree, greville)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree
        real(dp), intent(out) :: greville(:)

        integer :: i

        do i = 1, size(knots) - degree - 1
            greville(i) = sum(knots(i + 1:i + degree))/real(degree, dp)
        end do
    end subroutine spline_greville

    subroutine build_block_direction(spec, degree, default_multiplicity, blk, &
            direction, status)
        type(knot_spec_t), intent(in) :: spec
        integer, intent(in) :: degree, default_multiplicity, direction
        type(spline_block_t), intent(inout) :: blk
        integer, intent(out) :: status

        real(dp), allocatable :: breaks(:), knots(:), extraction(:, :, :)
        real(dp), allocatable :: element_breaks(:)
        integer, allocatable :: multiplicity(:), first(:)
        integer :: n

        if (allocated(spec%breakpoints)) then
            breaks = spec%breakpoints
            status = 0
        else
            call build_graded_breakpoints(spec%element_count, spec%grading, &
                breaks, status)
            if (status /= 0) return
        end if
        if (spec%reversed) breaks = 1.0_dp - breaks(size(breaks):1:-1)
        n = size(breaks)
        if (allocated(spec%multiplicity)) then
            status = 1
            if (size(spec%multiplicity) /= n) return
            multiplicity = spec%multiplicity
            if (spec%reversed) multiplicity = multiplicity(n:1:-1)
        else
            allocate(multiplicity(n))
            multiplicity = min(default_multiplicity, degree)
        end if
        call build_open_knot_vector(degree, breaks, multiplicity, knots, status)
        if (status /= 0) return
        call build_bezier_extraction_1d(knots, degree, extraction, first, &
            element_breaks, status)
        if (status /= 0) return
        blk%n(direction) = size(knots) - degree - 1
        blk%element_count(direction) = size(extraction, 3)
        if (direction == 1) then
            blk%knots1 = knots
            blk%breaks1 = element_breaks
            blk%extraction1 = extraction
            blk%first1 = first
        else
            blk%knots2 = knots
            blk%breaks2 = element_breaks
            blk%extraction2 = extraction
            blk%first2 = first
        end if
    end subroutine build_block_direction

    logical function sides_conform(space) result(ok)
        !! Matching sides need identical knot vectors, mirrored when the
        !! neighbour runs the side in the opposite direction.
        type(spline_space_t), intent(in) :: space

        real(dp), allocatable :: mine(:), theirs(:)
        integer :: k, s, k2, s2

        ok = .false.
        do k = 1, space%block_count
            do s = 1, 4
                k2 = space%layout%neighbour(s, k)
                if (k2 == 0) cycle
                s2 = space%layout%neighbour_side(s, k)
                call side_knots(space%blocks(k), s, mine)
                call side_knots(space%blocks(k2), s2, theirs)
                if (space%layout%reversed(s, k)) &
                    theirs = 1.0_dp - theirs(size(theirs):1:-1)
                if (size(mine) /= size(theirs)) return
                if (maxval(abs(mine - theirs)) > 1.0e-12_dp) return
            end do
        end do
        ok = .true.
    end function sides_conform

    subroutine side_knots(blk, side, knots)
        type(spline_block_t), intent(in) :: blk
        integer, intent(in) :: side
        real(dp), allocatable, intent(out) :: knots(:)

        if (side == 1 .or. side == 3) then
            knots = blk%knots1
        else
            knots = blk%knots2
        end if
    end subroutine side_knots

    subroutine side_trace(blk, side, trace)
        !! Patch-local indices of the functions on a side, in side order.
        type(spline_block_t), intent(in) :: blk
        integer, intent(in) :: side
        integer, allocatable, intent(out) :: trace(:)

        integer :: i, j, n1, n2

        n1 = blk%n(1)
        n2 = blk%n(2)
        select case (side)
        case (1)
            trace = [(i, i = 1, n1)]
        case (2)
            trace = [(n1 + n1*(j - 1), j = 1, n2)]
        case (3)
            trace = [(i + n1*(n2 - 1), i = 1, n1)]
        case default
            trace = [(1 + n1*(j - 1), j = 1, n2)]
        end select
    end subroutine side_trace

    subroutine build_c0_extraction(space, status)
        type(spline_space_t), intent(inout) :: space
        integer, intent(out) :: status

        integer, allocatable :: offsets(:), left_patch(:), left_local(:)
        integer, allocatable :: right_patch(:), right_local(:), signs(:)
        integer, allocatable :: mine(:), theirs(:), local_to_global(:)
        integer :: k, s, k2, s2, nb, global_count

        nb = space%block_count
        offsets = [(space%blocks(k)%offset + 1, k = 1, nb), space%broken_count + 1]
        allocate(left_patch(0), left_local(0), right_patch(0), right_local(0))
        do k = 1, nb
            do s = 1, 4
                k2 = space%layout%neighbour(s, k)
                if (k2 == 0) cycle
                s2 = space%layout%neighbour_side(s, k)
                if (k2 < k .or. (k2 == k .and. s2 < s)) cycle
                call side_trace(space%blocks(k), s, mine)
                call side_trace(space%blocks(k2), s2, theirs)
                if (space%layout%reversed(s, k)) theirs = theirs(size(theirs):1:-1)
                left_patch = [left_patch, spread(k, 1, size(mine))]
                right_patch = [right_patch, spread(k2, 1, size(mine))]
                left_local = [left_local, mine]
                right_local = [right_local, theirs]
            end do
        end do
        allocate(signs(size(left_patch)), local_to_global(space%broken_count))
        signs = 1
        call build_multipatch_signed_dof_map(offsets, left_patch, left_local, &
            right_patch, right_local, signs, local_to_global, global_count, status)
        if (status /= 0) then
            status = 3
            return
        end if
        space%dof_count = global_count
        space%broken_ptr = [(k, k = 1, space%broken_count + 1)]
        space%broken_dof = abs(local_to_global)
        space%broken_weight = real(sign(1, local_to_global), dp)
        status = 0
    end subroutine build_c0_extraction

    subroutine build_c1_space(space, status)
        !! C1 extraction on the non-periodic blocks; a periodic annulus (the
        !! last block) is coupled to them in C0 afterwards.
        type(spline_space_t), intent(inout) :: space
        integer, intent(out) :: status

        type(c1_block_knots_t), allocatable :: knots(:)
        type(block_layout_t) :: core
        integer, allocatable :: dof(:), broken(:), counts(:), dof_vertex(:)
        real(dp), allocatable :: weight(:)
        integer :: k, q, j, position, nc, ka

        status = 6
        ka = 0
        if (count(space%layout%periodic) > 1) return
        if (any(space%layout%periodic)) then
            ka = findloc(space%layout%periodic, .true., 1)
            if (ka /= space%block_count) return
        end if
        nc = space%block_count - merge(1, 0, ka /= 0)
        call truncate_layout(space%layout, nc, core)
        allocate(knots(nc))
        do k = 1, nc
            knots(k)%knots1 = space%blocks(k)%knots1
            knots(k)%knots2 = space%blocks(k)%knots2
        end do
        call build_c1_extraction(core, space%degree, knots, &
            [(space%blocks(k)%offset, k = 1, nc)], dof, broken, weight, &
            space%dof_count, status, dof_vertex)
        if (status /= 0) then
            status = 5 + status
            return
        end if
        if (ka /= 0) then
            call couple_annulus_c1(space, ka, dof, broken, weight, space%dof_count, &
                dof_vertex, status)
            if (status /= 0) then
                status = 7
                return
            end if
        end if
        allocate(counts(space%broken_count + 1))
        counts = 0
        do q = 1, size(dof)
            counts(broken(q) + 1) = counts(broken(q) + 1) + 1
        end do
        counts(1) = 1
        do j = 2, space%broken_count + 1
            counts(j) = counts(j) + counts(j - 1)
        end do
        space%broken_ptr = counts
        allocate(space%broken_dof(size(dof)), space%broken_weight(size(dof)))
        do q = 1, size(dof)
            position = counts(broken(q))
            counts(broken(q)) = counts(broken(q)) + 1
            space%broken_dof(position) = dof(q)
            space%broken_weight(position) = weight(q)
        end do
        status = 0
    end subroutine build_c1_space

    subroutine truncate_layout(layout, nc, core)
        !! The first nc blocks of a layout (all vertices kept).
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: nc
        type(block_layout_t), intent(out) :: core

        core = layout
        core%block_count = nc
        core%corners = layout%corners(:, :nc)
        core%side_tag = layout%side_tag(:, :nc)
        core%neighbour = layout%neighbour(:, :nc)
        core%neighbour_side = layout%neighbour_side(:, :nc)
        core%reversed = layout%reversed(:, :nc)
        core%arc_side = layout%arc_side(:nc)
        core%arc_center = layout%arc_center(:, :nc)
        core%periodic = layout%periodic(:nc)
        core%side_curve = layout%side_curve(:, :nc)
        core%side_u0 = layout%side_u0(:, :nc)
        core%side_du = layout%side_du(:, :nc)
        core%side_c0 = layout%side_c0(:, :nc)
    end subroutine truncate_layout

    subroutine couple_annulus_c1(space, ka, dof, broken, weight, core_count, &
            dof_vertex, status)
        !! Append the periodic annulus ka to the C1 core triplets: the core
        !! functions with a trace on the interface are replaced by
        !! zero-trace combinations of the vertex functions at each interface
        !! vertex and by extensions of the annulus interface functions, whose
        !! core traces are the restrictions (16). Status 3: an annulus trace
        !! has no extension into the core.
        type(spline_space_t), intent(inout) :: space
        integer, intent(in) :: ka, core_count, dof_vertex(:)
        integer, allocatable, intent(inout) :: dof(:), broken(:)
        real(dp), allocatable, intent(inout) :: weight(:)
        integer, intent(out) :: status

        integer, allocatable :: row_ptr(:), row_broken(:), trace_row(:), cols(:)
        integer, allocatable :: trace(:), column(:), new_dof(:), sel_rows(:)
        integer, allocatable :: sel_cols(:), out_dof(:), out_broken(:)
        real(dp), allocatable :: row_weight(:), t(:, :), g(:, :), restriction(:, :)
        real(dp), allocatable :: a(:, :), x(:), y(:, :), out_weight(:)
        logical, allocatable :: has_trace(:), mark(:)
        integer :: q, k, s, m, r, ntrace, n, nr, n1, i, j, c, v, count, info, nvec
        integer :: nout
        real(dp) :: residual

        status = 3
        call rows_by_dof(dof, broken, weight, core_count, row_ptr, row_broken, &
            row_weight)
        allocate(trace_row(space%broken_count))
        trace_row = 0
        ntrace = 0
        do k = 1, space%block_count
            if (k == ka) cycle
            do s = 1, 4
                if (space%layout%side_tag(s, k) /= LAYOUT_ANNULUS_INTERFACE) cycle
                call side_trace(space%blocks(k), s, trace)
                do m = 1, size(trace)
                    ntrace = ntrace + 1
                    trace_row(space%blocks(k)%offset + trace(m)) = ntrace
                end do
            end do
        end do
        allocate(has_trace(core_count), column(core_count))
        has_trace = .false.
        do i = 1, core_count
            do q = row_ptr(i), row_ptr(i + 1) - 1
                if (trace_row(row_broken(q)) /= 0 .and. abs(row_weight(q)) > 0.0_dp) &
                    has_trace(i) = .true.
            end do
        end do
        cols = pack([(i, i = 1, core_count)], has_trace)
        column = 0
        column(cols) = [(c, c = 1, size(cols))]
        allocate(t(ntrace, size(cols)))
        t = 0.0_dp
        do c = 1, size(cols)
            i = cols(c)
            do q = row_ptr(i), row_ptr(i + 1) - 1
                r = trace_row(row_broken(q))
                if (r /= 0) t(r, c) = t(r, c) + row_weight(q)
            end do
        end do
        n = space%blocks(ka)%element_count(1)
        nr = space%blocks(ka)%n(2)
        n1 = space%blocks(ka)%n(1)
        allocate(g(ntrace, n))
        g = 0.0_dp
        r = 0
        do k = 1, space%block_count
            if (k == ka) cycle
            do s = 1, 4
                if (space%layout%side_tag(s, k) /= LAYOUT_ANNULUS_INTERFACE) cycle
                call interface_restriction(space, k, s, ka, restriction, info)
                if (info /= 0) return
                g(r + 1:r + size(restriction, 1), :) = restriction
                r = r + size(restriction, 1)
            end do
        end do
        nout = 0
        allocate(out_dof(1024), out_broken(1024), out_weight(1024), &
            new_dof(core_count))
        count = 0
        new_dof = 0
        do i = 1, core_count
            if (has_trace(i)) cycle
            count = count + 1
            new_dof(i) = count
            call emit(count, [i], [1.0_dp])
        end do
        do v = 1, space%layout%vertex_count
            sel_cols = pack([(c, c = 1, size(cols))], dof_vertex(cols) == v)
            if (size(sel_cols) == 0) cycle
            call null_space(t(:, sel_cols), y, nvec)
            do j = 1, nvec
                count = count + 1
                call emit(count, cols(sel_cols), y(:, j))
            end do
        end do
        associate (blk => space%blocks(ka))
            allocate(mark(ntrace))
            do i = 1, n
                do j = 1, nr
                    count = count + 1
                    call push_annulus(count, i, j)
                    if (j /= nr) cycle
                    mark = abs(g(:, i)) > 1.0e-14_dp
                    sel_cols = pack([(c, c = 1, size(cols))], &
                        any(abs(t(:, :)) > 0.0_dp .and. spread(mark, 2, size(cols)), 1))
                    mark = mark .or. any(abs(t(:, sel_cols)) > 0.0_dp, 2)
                    sel_rows = pack([(r, r = 1, ntrace)], mark)
                    a = t(sel_rows, sel_cols)
                    call min_norm_solve(a, g(sel_rows, i), x, residual)
                    if (residual > 1.0e-10_dp*max(1.0_dp, maxval(abs(g(:, i))))) &
                        return
                    call emit(count, cols(sel_cols), x)
                end do
            end do
        end associate
        dof = out_dof(:nout)
        broken = out_broken(:nout)
        weight = out_weight(:nout)
        space%dof_count = count
        status = 0

    contains

        subroutine emit(target, sources, factors)
            !! Row target += sum factors(l) * E(sources(l), :).
            integer, intent(in) :: target, sources(:)
            real(dp), intent(in) :: factors(:)
            integer :: l, p

            do l = 1, size(sources)
                if (abs(factors(l)) < 1.0e-14_dp) cycle
                do p = row_ptr(sources(l)), row_ptr(sources(l) + 1) - 1
                    call append(target, row_broken(p), factors(l)*row_weight(p))
                end do
            end do
        end subroutine emit

        subroutine push_annulus(target, i, j)
            integer, intent(in) :: target, i, j
            integer :: b

            b = space%blocks(ka)%offset + i + n1*(j - 1)
            call append(target, b, 1.0_dp)
            if (i <= space%degree) call append(target, b + n, 1.0_dp)
        end subroutine push_annulus

        subroutine append(target, b, w)
            integer, intent(in) :: target, b
            real(dp), intent(in) :: w
            integer, allocatable :: gi(:), gb(:)
            real(dp), allocatable :: gw(:)

            if (nout == size(out_dof)) then
                allocate(gi(2*nout), gb(2*nout), gw(2*nout))
                gi(:nout) = out_dof
                gb(:nout) = out_broken
                gw(:nout) = out_weight
                call move_alloc(gi, out_dof)
                call move_alloc(gb, out_broken)
                call move_alloc(gw, out_weight)
            end if
            nout = nout + 1
            out_dof(nout) = target
            out_broken(nout) = b
            out_weight(nout) = w
        end subroutine append

    end subroutine couple_annulus_c1

    subroutine rows_by_dof(dof, broken, weight, dof_count, ptr, row_broken, &
            row_weight)
        integer, intent(in) :: dof(:), broken(:), dof_count
        real(dp), intent(in) :: weight(:)
        integer, allocatable, intent(out) :: ptr(:), row_broken(:)
        real(dp), allocatable, intent(out) :: row_weight(:)

        integer, allocatable :: next(:)
        integer :: q, i

        allocate(ptr(dof_count + 1), next(dof_count + 1), row_broken(size(dof)), &
            row_weight(size(dof)))
        ptr = 0
        do q = 1, size(dof)
            ptr(dof(q) + 1) = ptr(dof(q) + 1) + 1
        end do
        ptr(1) = 1
        do i = 2, dof_count + 1
            ptr(i) = ptr(i) + ptr(i - 1)
        end do
        next = ptr
        do q = 1, size(dof)
            row_broken(next(dof(q))) = broken(q)
            row_weight(next(dof(q))) = weight(q)
            next(dof(q)) = next(dof(q)) + 1
        end do
    end subroutine rows_by_dof

    subroutine null_space(a, y, kept)
        !! Orthonormal null space y(:, 1:kept) of a (small, dense).
        real(dp), intent(in) :: a(:, :)
        real(dp), allocatable, intent(out) :: y(:, :)
        integer, intent(out) :: kept

        real(dp), allocatable :: gram(:, :), lambda(:), work(:)
        integer :: n, info

        n = size(a, 2)
        allocate(gram(n, n), lambda(n), work(10*n + 10))
        gram = matmul(transpose(a), a)
        call dsyev('V', 'U', n, gram, n, lambda, work, size(work), info)
        kept = 0
        if (info /= 0) then
            allocate(y(n, 0))
            return
        end if
        kept = count(lambda <= 1.0e-12_dp*max(maxval(lambda), 1.0e-300_dp))
        y = gram(:, :kept)
    end subroutine null_space

    subroutine min_norm_solve(a, b, x, residual)
        !! Minimal-norm least-squares solution of a x = b (small, dense).
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp), allocatable, intent(out) :: x(:)
        real(dp), intent(out) :: residual

        real(dp), allocatable :: gram(:, :), lambda(:), work(:), c(:)
        integer :: n, info, i

        n = size(a, 2)
        allocate(gram(n, n), lambda(n), work(10*n + 10), x(n))
        gram = matmul(transpose(a), a)
        call dsyev('V', 'U', n, gram, n, lambda, work, size(work), info)
        x = 0.0_dp
        residual = huge(1.0_dp)
        if (info /= 0) return
        c = matmul(transpose(gram), matmul(transpose(a), b))
        do i = 1, n
            if (lambda(i) > 1.0e-12_dp*maxval(lambda)) &
                x = x + gram(:, i)*c(i)/lambda(i)
        end do
        residual = maxval(abs(matmul(a, x) - b))
    end subroutine min_norm_solve

    subroutine build_geometry(space, status)
        !! Control points of the block maps (1): bilinear interpolation at the
        !! Greville points, or, with a curved side C, collocation of C at the
        !! Greville points blended linearly to the opposite side. C1 spaces
        !! use G = B + (1 - tau)^2 (C - chord) instead, B bilinear and tau
        !! the distance parameter from C, so that the transversal derivative
        !! on the opposite side is bilinear and that interface stays AS-G1.
        type(spline_space_t), intent(inout) :: space
        integer, intent(out) :: status

        real(dp), allocatable :: g1(:), g2(:), arc(:, :), ring(:, :)
        real(dp) :: corner(2, 4), x(2), tangent(2), opposite(2), tau
        integer :: k, i, j, n1, n2, side, ka
        logical :: blend

        blend = space%regularity == SPLINE_C1
        allocate(space%geometry(2, space%broken_count))
        status = 0
        ka = 0
        do k = 1, space%block_count
            if (.not. space%layout%periodic(k)) cycle
            ka = k
            call periodic_geometry(space, k, ring, status)
            if (status /= 0) return
        end do
        do k = 1, space%block_count
            if (space%layout%periodic(k)) cycle
            associate (blk => space%blocks(k), layout => space%layout)
                n1 = blk%n(1)
                n2 = blk%n(2)
                allocate(g1(n1), g2(n2))
                call spline_greville(blk%knots1, space%degree, g1)
                call spline_greville(blk%knots2, space%degree, g2)
                corner = layout%vertices(:, layout%corners(:, k))
                side = layout%arc_side(k)
                if (side /= 0 .and. ka /= 0) then
                    if (layout%side_tag(side, k) == LAYOUT_ANNULUS_INTERFACE) then
                        call collocate_interface(space, k, side, ka, ring, arc, status)
                        if (status /= 0) return
                    end if
                end if
                if (side /= 0 .and. .not. allocated(arc)) then
                    if (side == 1 .or. side == 3) then
                        call collocate_side(layout, k, side, blk%knots1, &
                            space%degree, g1, arc, status)
                    else
                        call collocate_side(layout, k, side, blk%knots2, &
                            space%degree, g2, arc, status)
                    end if
                    if (status /= 0) return
                end if
                do j = 1, n2
                    do i = 1, n1
                        select case (side)
                        case (0)
                            x = bilinear(corner, g1(i), g2(j))
                        case (1, 3)
                            if (blend) then
                                x = bilinear(corner, g1(i), g2(j)) + &
                                    quadratic_blossom(blk%knots2, space%degree, j, &
                                    side == 3)*(arc(:, i) - bilinear(corner, g1(i), &
                                    merge(0.0_dp, 1.0_dp, side == 1)))
                            else
                                call layout_side_point(layout, k, 4 - side, g1(i), &
                                    opposite, tangent)
                                tau = merge(g2(j), 1.0_dp - g2(j), side == 1)
                                x = (1.0_dp - tau)*arc(:, i) + tau*opposite
                            end if
                        case default
                            if (blend) then
                                x = bilinear(corner, g1(i), g2(j)) + &
                                    quadratic_blossom(blk%knots1, space%degree, i, &
                                    side == 2)*(arc(:, j) - bilinear(corner, &
                                    merge(0.0_dp, 1.0_dp, side == 4), g2(j)))
                            else
                                call layout_side_point(layout, k, 6 - side, g2(j), &
                                    opposite, tangent)
                                tau = merge(g1(i), 1.0_dp - g1(i), side == 4)
                                x = (1.0_dp - tau)*arc(:, j) + tau*opposite
                            end if
                        end select
                        space%geometry(:, blk%offset + i + n1*(j - 1)) = x
                    end do
                end do
                deallocate(g1, g2)
                if (allocated(arc)) deallocate(arc)
            end associate
        end do
    end subroutine build_geometry

    subroutine collocate_side(layout, block, side, knots, degree, greville, &
            control, status)
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: block, side, degree
        real(dp), intent(in) :: knots(:), greville(:)
        real(dp), allocatable, intent(out) :: control(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: matrix(:, :), rhs(:, :), values(:), derivatives(:)
        real(dp) :: tangent(2)
        integer, allocatable :: pivot(:)
        integer :: n, q

        n = size(greville)
        allocate(matrix(n, n), rhs(n, 2), pivot(n))
        do q = 1, n
            call evaluate_bspline_basis(knots, degree, greville(q), values, &
                derivatives, status)
            if (status /= 0) return
            matrix(q, :) = values
            call layout_side_point(layout, block, side, greville(q), rhs(q, :), &
                tangent)
        end do
        call dgesv(n, 2, matrix, n, pivot, rhs, n, status)
        if (status /= 0) return
        control = transpose(rhs)
    end subroutine collocate_side

    subroutine build_periodic_direction(space, k, spec, status)
        !! xi1 of the periodic block k: breakpoints from the interface sides
        !! of the blocks built so far, else from spec.
        type(spline_space_t), intent(inout) :: space
        integer, intent(in) :: k
        type(knot_spec_t), intent(in) :: spec
        integer, intent(out) :: status

        real(dp), allocatable :: breaks(:), side_breaks(:), merged(:)
        real(dp), allocatable :: knots(:), extraction(:, :, :)
        integer, allocatable :: first(:)
        real(dp) :: start, delta
        integer :: k2, s, i, n

        allocate(breaks(0))
        do k2 = 1, space%block_count
            if (space%layout%periodic(k2)) cycle
            do s = 1, 4
                if (space%layout%side_tag(s, k2) /= LAYOUT_ANNULUS_INTERFACE) cycle
                call interface_map(space%layout, k2, s, k, start, delta)
                if (s == 1 .or. s == 3) then
                    side_breaks = space%blocks(k2)%breaks1
                else
                    side_breaks = space%blocks(k2)%breaks2
                end if
                breaks = [breaks, modulo(start + delta*side_breaks, 1.0_dp)]
            end do
        end do
        if (size(breaks) == 0) then
            call build_graded_breakpoints(spec%element_count, spec%grading, &
                breaks, status)
            if (status /= 0) return
        end if
        breaks = [breaks, 0.0_dp, 1.0_dp]
        where (breaks > 1.0_dp - 1.0e-12_dp) breaks = 1.0_dp
        where (breaks < 1.0e-12_dp) breaks = 0.0_dp
        call sort_reals(breaks)
        allocate(merged(size(breaks)))
        n = 1
        merged(1) = breaks(1)
        do i = 2, size(breaks)
            if (breaks(i) - merged(n) > 1.0e-12_dp) then
                n = n + 1
                merged(n) = breaks(i)
            end if
        end do
        merged(n) = 1.0_dp
        call build_periodic_extraction_1d(space%degree, merged(:n), knots, &
            extraction, first, status)
        if (status /= 0) return
        associate (blk => space%blocks(k))
            blk%n(1) = size(knots) - space%degree - 1
            blk%element_count(1) = size(extraction, 3)
            blk%knots1 = knots
            blk%breaks1 = merged(:n)
            blk%extraction1 = extraction
            blk%first1 = first
        end associate
    end subroutine build_periodic_direction

    subroutine interface_map(layout, k, s, ka, start, delta)
        !! xi1 of the periodic block ka along side s of block k is
        !! start + delta*t, t the side parameter, taken modulo 1. The map is
        !! part of the topology and does not change with the geometry.
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: k, s, ka
        real(dp), intent(out) :: start, delta

        start = modulo(layout%side_u0(s, k) - layout%side_u0(3, ka), 1.0_dp)
        if (start > 1.0_dp - 1.0e-12_dp) start = 0.0_dp
        delta = layout%side_du(s, k)
    end subroutine interface_map

    subroutine periodic_values(blk, degree, x, values, status)
        !! The n periodic functions of a periodic block's xi1 at x (mod 1).
        type(spline_block_t), intent(in) :: blk
        integer, intent(in) :: degree
        real(dp), intent(in) :: x
        real(dp), allocatable, intent(out) :: values(:)
        integer, intent(out) :: status

        real(dp), allocatable :: full(:), derivatives(:)
        integer :: n

        n = blk%element_count(1)
        call evaluate_bspline_basis(blk%knots1, degree, modulo(x, 1.0_dp), full, &
            derivatives, status)
        if (status /= 0) return
        values = full(:n)
        values(:degree) = values(:degree) + full(n + 1:n + degree)
    end subroutine periodic_values

    subroutine constrain_periodic(space, status)
        !! Rewrite the one-to-one C0 extraction: the copies i + n of the
        !! periodic functions i, and the core functions on interface sides,
        !! become combinations of the remaining global functions.
        type(spline_space_t), intent(inout) :: space
        integer, intent(out) :: status

        type :: combination_t
            integer, allocatable :: dof(:)
            real(dp), allocatable :: weight(:)
        end type combination_t
        type(combination_t), allocatable :: map(:)
        integer, allocatable :: trace(:), new_index(:), ptr(:), dofs(:)
        real(dp), allocatable :: restriction(:, :), weights(:)
        integer :: k, ka, n, nr, n1, i, j, g, s, m, old_count, kept, q, nnz

        status = 3
        if (count(space%layout%periodic) /= 1) return
        ka = findloc(space%layout%periodic, .true., 1)
        old_count = space%dof_count
        allocate(map(old_count), new_index(old_count))
        new_index = 0
        n = space%blocks(ka)%element_count(1)
        nr = space%blocks(ka)%n(2)
        n1 = space%blocks(ka)%n(1)
        do j = 1, nr
            do i = n + 1, n1
                new_index(annulus_dof(i, j)) = -1
            end do
        end do
        do k = 1, space%block_count
            if (k == ka) cycle
            do s = 1, 4
                if (space%layout%side_tag(s, k) /= LAYOUT_ANNULUS_INTERFACE) cycle
                call side_trace(space%blocks(k), s, trace)
                do m = 1, size(trace)
                    new_index(space%broken_dof(space%blocks(k)%offset + trace(m))) = -1
                end do
            end do
        end do
        kept = 0
        do g = 1, old_count
            if (new_index(g) /= 0) cycle
            kept = kept + 1
            new_index(g) = kept
            map(g)%dof = [kept]
            map(g)%weight = [1.0_dp]
        end do
        do j = 1, nr
            do i = n + 1, n1
                map(annulus_dof(i, j))%dof = [new_index(annulus_dof(i - n, j))]
                map(annulus_dof(i, j))%weight = [1.0_dp]
            end do
        end do
        do k = 1, space%block_count
            if (k == ka) cycle
            do s = 1, 4
                if (space%layout%side_tag(s, k) /= LAYOUT_ANNULUS_INTERFACE) cycle
                call interface_restriction(space, k, s, ka, restriction, status)
                if (status /= 0) then
                    status = 3
                    return
                end if
                call side_trace(space%blocks(k), s, trace)
                do m = 1, size(trace)
                    q = space%blocks(k)%offset + trace(m)
                    g = space%broken_dof(q)
                    if (allocated(map(g)%dof)) cycle
                    dofs = [(new_index(annulus_dof(i, nr)), i = 1, n)]
                    weights = space%broken_weight(q)*restriction(m, :)
                    map(g)%dof = pack(dofs, abs(weights) > 1.0e-14_dp)
                    map(g)%weight = pack(weights, abs(weights) > 1.0e-14_dp)
                end do
            end do
        end do
        allocate(ptr(space%broken_count + 1))
        ptr(1) = 1
        do q = 1, space%broken_count
            ptr(q + 1) = ptr(q) + size(map(space%broken_dof(q))%dof)
        end do
        nnz = ptr(space%broken_count + 1) - 1
        if (allocated(dofs)) deallocate(dofs)
        if (allocated(weights)) deallocate(weights)
        allocate(dofs(nnz), weights(nnz))
        do q = 1, space%broken_count
            g = space%broken_dof(q)
            dofs(ptr(q):ptr(q + 1) - 1) = map(g)%dof
            weights(ptr(q):ptr(q + 1) - 1) = space%broken_weight(q)*map(g)%weight
        end do
        space%broken_ptr = ptr
        space%broken_dof = dofs
        space%broken_weight = weights
        space%dof_count = kept
        status = 0

    contains

        integer function annulus_dof(i, j) result(dof)
            integer, intent(in) :: i, j

            dof = space%broken_dof(space%blocks(ka)%offset + i + n1*(j - 1))
        end function annulus_dof

    end subroutine constrain_periodic

    subroutine interface_restriction(space, k, s, ka, restriction, status)
        !! restriction(m, i): coefficient of the side function m of block k
        !! in the restriction of the periodic function i of block ka, by
        !! interpolation at the side Greville points (exact, the knots match).
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: k, s, ka
        real(dp), allocatable, intent(out) :: restriction(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: knots(:), greville(:), matrix(:, :), values(:)
        real(dp), allocatable :: derivatives(:)
        integer, allocatable :: pivot(:)
        real(dp) :: start, delta
        integer :: nside, q

        call side_knots(space%blocks(k), s, knots)
        nside = size(knots) - space%degree - 1
        allocate(greville(nside), matrix(nside, nside), pivot(nside), &
            restriction(nside, space%blocks(ka)%element_count(1)))
        call spline_greville(knots, space%degree, greville)
        call interface_map(space%layout, k, s, ka, start, delta)
        do q = 1, nside
            call evaluate_bspline_basis(knots, space%degree, greville(q), values, &
                derivatives, status)
            if (status /= 0) return
            matrix(q, :) = values
            call periodic_values(space%blocks(ka), space%degree, &
                start + delta*greville(q), values, status)
            if (status /= 0) return
            restriction(q, :) = values
        end do
        call dgesv(nside, size(restriction, 2), matrix, nside, pivot, restriction, &
            nside, status)
    end subroutine interface_restriction

    subroutine periodic_geometry(space, k, ring, status)
        !! Annulus map (1 - xi2) w(xi1) + xi2 r(xi1), ruled between the
        !! periodic Greville interpolants w of the wall side (side 1) and r of
        !! the interface side (side 3); ring returns the control points of r.
        type(spline_space_t), intent(inout) :: space
        integer, intent(in) :: k
        real(dp), allocatable, intent(out) :: ring(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: g1(:), g2(:), matrix(:, :), rhs(:, :), values(:)
        real(dp), allocatable :: wall(:, :)
        integer, allocatable :: pivot(:)
        real(dp) :: tangent(2)
        integer :: n, i, j, q

        associate (blk => space%blocks(k), layout => space%layout)
            n = blk%element_count(1)
            allocate(g1(blk%n(1)), g2(blk%n(2)), matrix(n, n), rhs(n, 4), pivot(n))
            call spline_greville(blk%knots1, space%degree, g1)
            call spline_greville(blk%knots2, space%degree, g2)
            do q = 1, n
                call periodic_values(blk, space%degree, g1(q), values, status)
                if (status /= 0) return
                matrix(q, :) = values
                call layout_side_point(layout, k, 1, modulo(g1(q), 1.0_dp), &
                    rhs(q, 1:2), tangent)
                call layout_side_point(layout, k, 3, modulo(g1(q), 1.0_dp), &
                    rhs(q, 3:4), tangent)
            end do
            call dgesv(n, 4, matrix, n, pivot, rhs, n, status)
            if (status /= 0) return
            wall = transpose(rhs(:, 1:2))
            ring = transpose(rhs(:, 3:4))
            do j = 1, blk%n(2)
                do i = 1, blk%n(1)
                    q = modulo(i - 1, n) + 1
                    space%geometry(:, blk%offset + i + blk%n(1)*(j - 1)) = &
                        (1.0_dp - g2(j))*wall(:, q) + g2(j)*ring(:, q)
                end do
            end do
        end associate
    end subroutine periodic_geometry

    subroutine collocate_interface(space, k, side, ka, ring, control, status)
        !! Side control points of block k on the interface: the restriction
        !! of the annulus interface curve, so both maps agree exactly.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: k, side, ka
        real(dp), intent(in) :: ring(:, :)
        real(dp), allocatable, intent(out) :: control(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: restriction(:, :)

        call interface_restriction(space, k, side, ka, restriction, status)
        if (status /= 0) return
        control = matmul(ring, transpose(restriction))
    end subroutine collocate_interface

    pure subroutine sort_reals(a)
        real(dp), intent(inout) :: a(:)
        real(dp) :: x
        integer :: i, j

        do i = 2, size(a)
            x = a(i)
            j = i - 1
            do while (j >= 1)
                if (a(j) <= x) exit
                a(j + 1) = a(j)
                j = j - 1
            end do
            a(j + 1) = x
        end do
    end subroutine sort_reals

    pure real(dp) function quadratic_blossom(knots, degree, i, mirrored) result(c)
        !! B-spline coefficient i of (1 - x)^2, or of x^2 when mirrored (blossom).
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree, i
        logical, intent(in) :: mirrored

        real(dp) :: f(degree)
        integer :: a, b

        f = knots(i + 1:i + degree)
        if (.not. mirrored) f = 1.0_dp - f
        c = 0.0_dp
        do b = 2, degree
            do a = 1, b - 1
                c = c + f(a)*f(b)
            end do
        end do
        c = c/real(degree*(degree - 1)/2, dp)
    end function quadratic_blossom

    pure function bilinear(corner, u, v) result(x)
        real(dp), intent(in) :: corner(2, 4), u, v
        real(dp) :: x(2)

        x = (1.0_dp - u)*(1.0_dp - v)*corner(:, 1) + u*(1.0_dp - v)*corner(:, 2) + &
            u*v*corner(:, 3) + (1.0_dp - u)*v*corner(:, 4)
    end function bilinear

end module fortfem_spline_space
