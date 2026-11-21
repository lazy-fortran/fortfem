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
    use fortfem_kinds, only: dp
    use fortfem_spline_extraction_1d, only: build_graded_breakpoints, &
        build_open_knot_vector, build_bezier_extraction_1d, SPLINE_GRID_UNIFORM
    use fortfem_spline_layout, only: block_layout_t, layout_side_point
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
    end type spline_space_t

    public :: spline_space_build
    public :: spline_boundary_dofs
    public :: spline_cell_broken
    public :: spline_greville

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    subroutine spline_space_build(layout, options, space, status)
        !! Status 1: invalid options, 2: knot construction failed, 3: knot
        !! vectors of two blocks disagree along a shared side, 4: unknown
        !! regularity, 5: geometry construction failed, 6: the C1 space needs
        !! bilinear blocks, degree >= 3 and interior multiplicity >= 2,
        !! 7: an interface is not analysis-suitable G1.
        !! C1 spaces default to interior multiplicity 2 (C^(p-2) in blocks).
        type(block_layout_t), intent(in) :: layout
        type(spline_space_options_t), intent(in) :: options
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status

        integer :: k, d, offset, cell, e1, e2, nb

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
        offset = 0
        do k = 1, nb
            do d = 1, 2
                call build_block_direction(options%knots(d, k), options%degree, &
                    merge(2, 1, options%regularity == SPLINE_C1), space%blocks(k), d, &
                    status)
                if (status /= 0) then
                    status = 2
                    return
                end if
            end do
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
        status = 0
    end subroutine spline_space_build

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
        type(spline_space_t), intent(inout) :: space
        integer, intent(out) :: status

        type(c1_block_knots_t), allocatable :: knots(:)
        integer, allocatable :: dof(:), broken(:), counts(:)
        real(dp), allocatable :: weight(:)
        integer :: k, q, j, position

        allocate(knots(space%block_count))
        do k = 1, space%block_count
            knots(k)%knots1 = space%blocks(k)%knots1
            knots(k)%knots2 = space%blocks(k)%knots2
        end do
        call build_c1_extraction(space%layout, space%degree, knots, &
            [(space%blocks(k)%offset, k = 1, space%block_count)], dof, broken, &
            weight, space%dof_count, status)
        if (status /= 0) then
            status = 5 + status
            return
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

    subroutine build_geometry(space, status)
        !! Control points of the ruled block maps (1): bilinear interpolation
        !! at the Greville points, or, with an arc side, collocation of the arc
        !! at the Greville points blended linearly to the opposite side.
        type(spline_space_t), intent(inout) :: space
        integer, intent(out) :: status

        real(dp), allocatable :: g1(:), g2(:), arc(:, :)
        real(dp) :: corner(2, 4), x(2), tangent(2), opposite(2), tau
        integer :: k, i, j, n1, n2, side

        allocate(space%geometry(2, space%broken_count))
        status = 0
        do k = 1, space%block_count
            associate (blk => space%blocks(k), layout => space%layout)
                n1 = blk%n(1)
                n2 = blk%n(2)
                allocate(g1(n1), g2(n2))
                call spline_greville(blk%knots1, space%degree, g1)
                call spline_greville(blk%knots2, space%degree, g2)
                corner = layout%vertices(:, layout%corners(:, k))
                side = layout%arc_side(k)
                if (side /= 0) then
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
                            call layout_side_point(layout, k, 4 - side, g1(i), &
                                opposite, tangent)
                            tau = merge(g2(j), 1.0_dp - g2(j), side == 1)
                            x = (1.0_dp - tau)*arc(:, i) + tau*opposite
                        case default
                            call layout_side_point(layout, k, 6 - side, g2(j), &
                                opposite, tangent)
                            tau = merge(g1(i), 1.0_dp - g1(i), side == 4)
                            x = (1.0_dp - tau)*arc(:, j) + tau*opposite
                        end select
                        space%geometry(:, blk%offset + i + n1*(j - 1)) = x
                    end do
                end do
                deallocate(g1, g2)
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

    pure function bilinear(corner, u, v) result(x)
        real(dp), intent(in) :: corner(2, 4), u, v
        real(dp) :: x(2)

        x = (1.0_dp - u)*(1.0_dp - v)*corner(:, 1) + u*(1.0_dp - v)*corner(:, 2) + &
            u*v*corner(:, 3) + (1.0_dp - u)*v*corner(:, 4)
    end function bilinear

end module fortfem_spline_space
