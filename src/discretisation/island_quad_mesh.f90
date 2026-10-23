module fortfem_island_quad_mesh
    ! Contour-following curved quadrilateral mesh of a disk that contains one
    ! chain of m islands of a scalar level function (pure geometry, no physics):
    ! O-points are maxima, X-points saddles, and the separatrix is the level set
    ! through the X-points. Both separatrix branches and two bounding level sets
    ! (inside and outside the chain) are exact mesh lines; rows between them
    ! blend radially, so the separatrix kinks stay at X-point vertices instead
    ! of bending level sets running through cells. Every island is an O-grid
    ! of rings homothetic to its lens about the O-point, every X-point is a mesh
    ! vertex (each of its four sectors split by one edge), and the outer layer
    ! blends to a circular wall. All edges are shared curves sampled once at GLL
    ! points; cell geometry is the transfinite (Coons) interpolant of its edges.
    ! Assumption: along every ray from the centre the level increases from
    ! r_core to a single ridge and decreases from it to r_wall.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortfem_tensor_mesh, only: tensor_mesh_t, build_tensor_mesh
    use fortfem_gll_line, only: gll_rule
    use fortfem_tensor_mesh_tools, only: orient_tensor_cells, mark_boundary_faces
    implicit none
    private
    public :: island_chain_t, island_mesh_counts_t, island_chain_quad_mesh
    public :: level_function
    abstract interface
        function level_function(x) result(value)
            import dp
            real(dp), intent(in) :: x(2)
            real(dp) :: value
        end function
    end interface
    type :: island_chain_t
        real(dp) :: center(2) = 0
        real(dp), allocatable :: opoint(:, :), xpoint(:, :) ! (2, m)
        ! Core-square corner radius, bounding contours crossing the ray through
        ! the first O-point at r_inner and r_outer, and the circular wall.
        real(dp) :: r_core = 0.2_dp, r_inner = 0.35_dp, r_outer = 0.7_dp
        real(dp) :: r_wall = 1
    end type
    type :: island_mesh_counts_t
        integer :: half = 2 ! cells along half an island (X-point to midpoint)
        integer :: core_rings = 2, inner_rows = 3, island_rings = 2
        integer :: outer_rows = 3, wall_rows = 2
    end type
    integer, parameter :: ROW_SQUARE = 1, ROW_BLEND = 2, ROW_INNER = 3, &
                          ROW_OUTER = 4, ROW_CIRCLE = 5
contains
    subroutine island_chain_quad_mesh(level, chain, counts, p, q, mesh, status)
        procedure(level_function) :: level
        type(island_chain_t), intent(in) :: chain
        type(island_mesh_counts_t), intent(in) :: counts
        integer, intent(in) :: p, q
        type(tensor_mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp), allocatable :: z(:), w(:), vx(:, :), samples(:, :, :), geom(:, :, :), &
                                 ucol(:), uxs(:), rowlevel(:), rowfrac(:), &
                                 lens(:, :), ring(:, :, :), diamond(:, :, :)
        integer, allocatable :: rowkind(:), rowa(:), rowb(:), vid(:, :), &
                                edgev(:, :), table(:), cells(:, :), corners(:, :), &
                                order(:), lensid(:), ringid(:, :), diamondid(:, :)
        real(dp) :: corner(2, 0:3), sep, opolar(2), t, c0(2)
        integer :: m, h, n, nc, nrow, rin, rsep, rout, nv, ne, ncell, hsize, j, k, i, &
                   row, flipped, nisl, kk, j0, ib, jb, nmax, emax, cmax, so
        status = 1
        m = 0
        if (allocated(chain%xpoint)) m = size(chain%xpoint, 2)
        if (m < 1 .or. .not. allocated(chain%opoint)) return
        if (size(chain%opoint, 2) /= m) return
        h = counts%half
        if (h < 1 .or. modulo(m*h, 2) /= 0 .or. q < 1 .or. p < 1) return
        if (min(counts%core_rings, counts%inner_rows, counts%outer_rows, &
                counts%wall_rows) < 1 .or. counts%island_rings < 0) return
        if (.not. (0 < chain%r_core .and. chain%r_core < chain%r_inner .and. &
                   chain%r_outer < chain%r_wall)) return
        n = 2*m*h
        nc = n/4
        nisl = counts%island_rings
        call gll_rule(q, z, w, status)
        if (status /= 0) return
        status = 1
        ! X-points sorted by angle; columns are uniform in angle between them.
        allocate (uxs(m + 1), order(m), ucol(0:n))
        do k = 1, m
            uxs(k) = atan2(chain%xpoint(2, k) - chain%center(2), &
                           chain%xpoint(1, k) - chain%center(1))
            order(k) = k
        end do
        do k = 2, m
            j = k
            do while (j > 1)
                if (uxs(order(j)) >= uxs(order(j - 1))) exit
                i = order(j); order(j) = order(j - 1); order(j - 1) = i
                j = j - 1
            end do
        end do
        uxs(1:m) = uxs(order)
        uxs(m + 1) = uxs(1) + 2*pi
        do k = 1, m
            do i = 0, 2*h
                ucol(2*h*(k - 1) + i) = uxs(k) + i*(uxs(k + 1) - uxs(k))/(2*h)
            end do
        end do
        sep = level(chain%xpoint(:, order(1)))
        do k = 1, m
            if (abs(level(chain%xpoint(:, k)) - sep) > &
                1.e-9_dp*max(1.0_dp, abs(sep))) return
        end do
        opolar = polar(chain%opoint(:, island_opoint(1)))
        do k = 0, 3
            corner(:, k) = chain%center + chain%r_core*[cos(ucol(k*nc)), &
                                                        sin(ucol(k*nc))]
        end do
        ! Closed rows from the core square to the wall. The two separatrix rows
        ! are rsep (inner branch) and so (outer branch).
        nrow = counts%core_rings + counts%inner_rows + 1 + counts%outer_rows + &
               counts%wall_rows
        allocate (rowkind(0:nrow), rowlevel(0:nrow), rowa(0:nrow), rowb(0:nrow), &
                  rowfrac(0:nrow))
        rowkind = 0; rowlevel = 0; rowa = 0; rowb = 0; rowfrac = 0
        rowkind(0) = ROW_SQUARE
        rin = counts%core_rings
        rsep = rin + counts%inner_rows
        so = rsep + 1
        rout = so + counts%outer_rows
        rowkind(rin) = ROW_INNER
        rowlevel(rin) = level(point(chain%r_inner, opolar(2)))
        rowkind(rsep) = ROW_INNER
        rowlevel(rsep) = sep
        rowkind(so) = ROW_OUTER
        rowlevel(so) = sep
        do i = 1, rin - 1
            rowkind(i) = ROW_BLEND; rowa(i) = 0; rowb(i) = rin
            rowfrac(i) = real(i, dp)/rin
        end do
        do i = rin + 1, rsep - 1
            rowkind(i) = ROW_BLEND; rowa(i) = rin; rowb(i) = rsep
            rowfrac(i) = real(i - rin, dp)/counts%inner_rows
        end do
        rowkind(rout) = ROW_OUTER
        rowlevel(rout) = level(point(chain%r_outer, opolar(2)))
        do i = so + 1, rout - 1
            rowkind(i) = ROW_BLEND; rowa(i) = so; rowb(i) = rout
            rowfrac(i) = real(i - so, dp)/counts%outer_rows
        end do
        do i = rout + 1, nrow
            rowkind(i) = ROW_BLEND; rowa(i) = rout; rowb(i) = -1
            rowfrac(i) = real(i - rout, dp)/counts%wall_rows
        end do
        rowkind(nrow) = ROW_CIRCLE
        nmax = (nc + 1)**2 + n*(nrow + 1) + m*(4*h*(nisl + 1) + (h + 1)**2)
        emax = 4*nmax + 16
        cmax = nmax + 16
        hsize = 2*emax + 1
        allocate (vx(2, nmax), samples(2, q + 1, emax), edgev(2, emax), &
                  table(hsize), cells(4, cmax), vid(0:nrow, 0:n))
        table = 0
        vid = 0
        nv = 0; ne = 0; ncell = 0
        ! Core square: transfinite grid whose boundary meets the column rays.
        block
            integer :: grid(0:nc, 0:nc), ia, ja
            real(dp) :: bnd(2, 0:n - 1), al, be, x(2)
            do j = 0, n - 1
                bnd(:, j) = point(row_radius(0, ucol(j)), ucol(j))
            end do
            do ja = 0, nc
                do ia = 0, nc
                    if (ia == 0 .or. ia == nc .or. ja == 0 .or. ja == nc) then
                        j = square_column(ia, ja)
                        if (vid(0, j) == 0) vid(0, j) = add_vertex(bnd(:, j))
                        grid(ia, ja) = vid(0, j)
                        cycle
                    end if
                    al = real(ia, dp)/nc
                    be = real(ja, dp)/nc
                    x = (1 - be)*bnd(:, square_column(ia, 0)) + &
                        be*bnd(:, square_column(ia, nc)) + &
                        (1 - al)*bnd(:, square_column(0, ja)) + &
                        al*bnd(:, square_column(nc, ja)) - &
                        ((1 - al)*(1 - be)*bnd(:, 0) + al*(1 - be)*bnd(:, nc) + &
                         al*be*bnd(:, 2*nc) + (1 - al)*be*bnd(:, 3*nc))
                    grid(ia, ja) = add_vertex(x)
                end do
            end do
            do ja = 0, nc - 1
                do ia = 0, nc - 1
                    call straight_cell([grid(ia, ja), grid(ia + 1, ja), &
                                        grid(ia + 1, ja + 1), grid(ia, ja + 1)])
                end do
            end do
        end block
        vid(0, n) = vid(0, 0)
        ! Closed rows. The outer separatrix row shares the X-point vertices.
        do row = 1, nrow
            do j = 0, n - 1
                if ((row == rsep .or. row == so) .and. modulo(j, 2*h) == 0) then
                    if (row == so) then
                        vid(row, j) = vid(rsep, j)
                    else
                        vid(row, j) = add_vertex(chain%xpoint(:, order(j/(2*h) + 1)))
                    end if
                else
                    vid(row, j) = add_vertex(point(row_radius(row, ucol(j)), ucol(j)))
                end if
            end do
            vid(row, n) = vid(row, 0)
        end do
        ! Row edges first, so islands reuse the separatrix curves.
        do row = 1, nrow
            do j = 0, n - 1
                call row_edge(row, j)
            end do
        end do
        do row = 0, nrow - 1
            if (row == rsep) cycle
            do j = 0, n - 1
                call straight_cell([vid(row, j), vid(row, j + 1), &
                                    vid(row + 1, j + 1), vid(row + 1, j)])
            end do
        end do
        ! Islands: lens loop -> contour rings -> transition -> diamond core.
        allocate (lens(2, 0:4*h), lensid(0:4*h), ring(2, 0:4*h, 0:nisl), &
                  ringid(0:4*h, 0:nisl), diamond(2, 0:h, 0:h), &
                  diamondid(0:h, 0:h))
        do k = 1, m
            j0 = 2*h*(k - 1)
            opolar = polar(chain%opoint(:, island_opoint(k)))
            do i = 0, 4*h - 1
                if (i <= 2*h) then
                    jb = j0 + i
                    lensid(i) = vid(rsep, jb)
                else
                    jb = j0 + 4*h - i
                    lensid(i) = vid(so, jb)
                end if
                lens(:, i) = [norm2(vx(:, lensid(i)) - chain%center), ucol(jb)]
            end do
            lens(:, 4*h) = lens(:, 0)
            lensid(4*h) = lensid(0)
            ring(:, :, 0) = lens
            ringid(:, 0) = lensid
            do kk = 1, nisl
                t = 1 - real(kk, dp)/(nisl + 1)
                do i = 0, 4*h - 1
                    ring(:, i, kk) = opolar + t*(lens(:, i) - opolar)
                    ringid(i, kk) = add_vertex(chart(ring(:, i, kk)))
                end do
                ring(:, 4*h, kk) = ring(:, 0, kk)
                ringid(4*h, kk) = ringid(0, kk)
                do i = 0, 4*h - 1
                    call island_ring_edge(ringid(i, kk), ringid(i + 1, kk), i, t)
                end do
                do i = 0, 4*h - 1
                    call chart_edge(ringid(i, kk - 1), ringid(i, kk), &
                                    ring(:, i, kk - 1), ring(:, i, kk))
                    call chart_cell([ringid(i, kk - 1), ringid(i + 1, kk - 1), &
                                     ringid(i + 1, kk), ringid(i, kk)], &
                                    ring(:, [i, i + 1], kk - 1), ring(:, [i + 1, i], kk))
                end do
            end do
            do i = 0, 3
                c0 = opolar + 0.5_dp*(ring(:, i*h, nisl) - opolar)
                diamond(:, merge(0, h, i == 0 .or. i == 3), &
                        merge(0, h, i <= 1)) = c0
            end do
            block
                real(dp) :: d00(2), d10(2), d11(2), d01(2), al, be
                integer :: ia, ja
                d00 = diamond(:, 0, 0); d10 = diamond(:, h, 0)
                d11 = diamond(:, h, h); d01 = diamond(:, 0, h)
                do ja = 0, h
                    do ia = 0, h
                        al = real(ia, dp)/h
                        be = real(ja, dp)/h
                        diamond(:, ia, ja) = (1 - al)*(1 - be)*d00 + al*(1 - be)*d10 + &
                                             al*be*d11 + (1 - al)*be*d01
                        diamondid(ia, ja) = add_vertex(chart(diamond(:, ia, ja)))
                    end do
                end do
                do ja = 0, h - 1
                    do ia = 0, h - 1
                        call chart_cell([diamondid(ia, ja), diamondid(ia + 1, ja), &
                                         diamondid(ia + 1, ja + 1), &
                                         diamondid(ia, ja + 1)], &
                                        reshape([diamond(:, ia, ja), &
                                                 diamond(:, ia + 1, ja)], [2, 2]), &
                                        reshape([diamond(:, ia + 1, ja + 1), &
                                                 diamond(:, ia, ja + 1)], [2, 2]))
                    end do
                end do
                do i = 0, 4*h - 1
                    call diamond_index(i, ia, ja)
                    call diamond_index(i + 1, ib, jb)
                    call chart_cell([ringid(i, nisl), ringid(i + 1, nisl), &
                                     diamondid(ib, jb), diamondid(ia, ja)], &
                                    ring(:, [i, i + 1], nisl), &
                                    reshape([diamond(:, ib, jb), diamond(:, ia, ja)], &
                                            [2, 2]))
                end do
            end block
        end do
        ! Coons patches over the shared edges, then orientation and topology.
        allocate (geom(2, (q + 1)**2, ncell), corners(4, ncell))
        do i = 1, ncell
            call coons(cells(:, i), geom(:, :, i))
            corners(:, i) = cells([1, 2, 4, 3], i)
        end do
        call orient_tensor_cells(2, q, corners, geom, flipped)
        call build_tensor_mesh(2, p, q, corners, geom, mesh, status)
        if (status /= 0) return
        call mark_boundary_faces(mesh)
    contains
        function point(radius, angle) result(x)
            real(dp), intent(in) :: radius, angle
            real(dp) :: x(2)
            x = chain%center + radius*[cos(angle), sin(angle)]
        end function

        function chart(ru) result(x)
            real(dp), intent(in) :: ru(2)
            real(dp) :: x(2)
            x = point(ru(1), ru(2))
        end function

        function polar(x) result(ru)
            ! Polar chart coordinates with the angle unwrapped onto the columns.
            real(dp), intent(in) :: x(2)
            real(dp) :: ru(2)
            ru(1) = norm2(x - chain%center)
            ru(2) = atan2(x(2) - chain%center(2), x(1) - chain%center(1))
            do while (ru(2) < ucol(0))
                ru(2) = ru(2) + 2*pi
            end do
            do while (ru(2) >= ucol(0) + 2*pi)
                ru(2) = ru(2) - 2*pi
            end do
        end function

        integer function island_opoint(island)
            ! The O-point whose angle lies between X-points island and island + 1.
            integer, intent(in) :: island
            integer :: o
            real(dp) :: ru(2)
            island_opoint = 0
            do o = 1, m
                ru = polar(chain%opoint(:, o))
                if (ru(2) > ucol(2*h*(island - 1)) .and. &
                    ru(2) < ucol(2*h*island)) island_opoint = o
            end do
            if (island_opoint == 0) island_opoint = island
        end function

        integer function square_column(ia, ja)
            integer, intent(in) :: ia, ja
            if (ja == 0) then
                square_column = ia
            else if (ia == nc) then
                square_column = nc + ja
            else if (ja == nc) then
                square_column = 2*nc + (nc - ia)
            else
                square_column = modulo(3*nc + (nc - ja), n)
            end if
        end function

        subroutine diamond_index(i, ia, ja)
            integer, intent(in) :: i
            integer, intent(out) :: ia, ja
            integer :: side, offset
            side = modulo(i, 4*h)/h
            offset = modulo(i, 4*h) - side*h
            select case (side)
            case (0)
                ia = offset; ja = 0
            case (1)
                ia = h; ja = offset
            case (2)
                ia = h - offset; ja = h
            case default
                ia = 0; ja = h - offset
            end select
        end subroutine

        real(dp) function ridge(angle)
            ! Golden-section maximum of the level along the ray.
            real(dp), intent(in) :: angle
            real(dp), parameter :: g = (sqrt(5.0_dp) - 1)/2
            real(dp) :: lo, hi, x1, x2, f1, f2
            integer :: it
            lo = chain%r_core; hi = chain%r_wall
            x1 = hi - g*(hi - lo); x2 = lo + g*(hi - lo)
            f1 = level(point(x1, angle)); f2 = level(point(x2, angle))
            do it = 1, 200
                if (hi - lo < 1.e-14_dp*chain%r_wall) exit
                if (f1 < f2) then
                    lo = x1; x1 = x2; f1 = f2
                    x2 = lo + g*(hi - lo); f2 = level(point(x2, angle))
                else
                    hi = x2; x2 = x1; f2 = f1
                    x1 = hi - g*(hi - lo); f1 = level(point(x1, angle))
                end if
            end do
            ridge = (lo + hi)/2
        end function

        real(dp) function contour(value, angle, inner)
            ! Bisection for the level set on the inner or outer flank of the ridge.
            real(dp), intent(in) :: value, angle
            logical, intent(in) :: inner
            real(dp) :: lo, hi, mid, top
            integer :: it
            top = ridge(angle)
            if (level(point(top, angle)) <= value) then
                contour = top
                return
            end if
            if (inner) then
                lo = chain%r_core; hi = top
            else
                lo = top; hi = chain%r_wall
            end if
            do it = 1, 200
                mid = (lo + hi)/2
                if (mid <= lo .or. mid >= hi) exit
                if ((level(point(mid, angle)) < value) .eqv. inner) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            contour = (lo + hi)/2
        end function

        recursive real(dp) function row_radius(row, angle) result(radius)
            integer, intent(in) :: row
            real(dp), intent(in) :: angle
            real(dp) :: e(2), m2(2, 2), rhs(2), det, other
            integer :: side, kk
            select case (rowkind(row))
            case (ROW_SQUARE)
                side = 3
                do kk = 0, 3
                    if (angle >= ucol(kk*nc) .and. angle <= ucol((kk + 1)*nc)) then
                        side = kk
                        exit
                    end if
                end do
                e = [cos(angle), sin(angle)]
                m2(:, 1) = e
                m2(:, 2) = corner(:, side) - corner(:, modulo(side + 1, 4))
                rhs = corner(:, side) - chain%center
                det = m2(1, 1)*m2(2, 2) - m2(1, 2)*m2(2, 1)
                radius = (rhs(1)*m2(2, 2) - rhs(2)*m2(1, 2))/det
            case (ROW_BLEND)
                if (rowb(row) < 0) then
                    other = chain%r_wall
                else
                    other = row_radius(rowb(row), angle)
                end if
                radius = (1 - rowfrac(row))*row_radius(rowa(row), angle) + &
                         rowfrac(row)*other
            case (ROW_INNER)
                radius = contour(rowlevel(row), angle, .true.)
            case (ROW_OUTER)
                radius = contour(rowlevel(row), angle, .false.)
            case default
                radius = chain%r_wall
            end select
        end function

        integer function add_vertex(x)
            real(dp), intent(in) :: x(2)
            nv = nv + 1
            vx(:, nv) = x
            add_vertex = nv
        end function

        integer function find_edge(va, vb, create)
            integer, intent(in) :: va, vb
            logical, intent(in) :: create
            integer(int64) :: hash
            integer :: slot, lo, hi
            lo = min(va, vb); hi = max(va, vb)
            hash = modulo(int(lo, int64)*1000003_int64 + int(hi, int64), &
                          int(hsize, int64))
            slot = int(hash) + 1
            find_edge = 0
            do
                if (table(slot) == 0) then
                    if (.not. create) return
                    ne = ne + 1
                    edgev(:, ne) = [va, vb]
                    table(slot) = ne
                    find_edge = -ne
                    return
                end if
                if (min(edgev(1, table(slot)), edgev(2, table(slot))) == lo .and. &
                    max(edgev(1, table(slot)), edgev(2, table(slot))) == hi) then
                    find_edge = table(slot)
                    return
                end if
                slot = modulo(slot, hsize) + 1
            end do
        end function

        subroutine row_edge(row, j)
            integer, intent(in) :: row, j
            integer :: id, l
            real(dp) :: angle
            id = find_edge(vid(row, j), vid(row, j + 1), .true.)
            if (id > 0) return
            id = -id
            do l = 1, q + 1
                angle = ucol(j) + (ucol(j + 1) - ucol(j))*(z(l) + 1)/2
                samples(:, l, id) = point(row_radius(row, angle), angle)
            end do
            call pin_ends(id)
        end subroutine

        subroutine island_ring_edge(va, vb, i, t)
            ! Lens edge i (separatrix arc between lens vertices i and i + 1),
            ! scaled by t about the O-point in the polar chart.
            integer, intent(in) :: va, vb, i
            real(dp), intent(in) :: t
            integer :: id, l, jb
            real(dp) :: angle, ru(2)
            id = find_edge(va, vb, .true.)
            if (id > 0) return
            id = -id
            do l = 1, q + 1
                if (i < 2*h) then
                    jb = j0 + i
                    angle = ucol(jb) + (ucol(jb + 1) - ucol(jb))*(z(l) + 1)/2
                    ru = [row_radius(rsep, angle), angle]
                else
                    jb = j0 + 4*h - i
                    angle = ucol(jb) + (ucol(jb - 1) - ucol(jb))*(z(l) + 1)/2
                    ru = [row_radius(so, angle), angle]
                end if
                samples(:, l, id) = chart(opolar + t*(ru - opolar))
            end do
            call pin_ends(id)
        end subroutine

        subroutine chart_edge(va, vb, ra, rb)
            ! Straight segment in the polar chart, curved in the plane.
            integer, intent(in) :: va, vb
            real(dp), intent(in) :: ra(2), rb(2)
            integer :: id, l
            id = find_edge(va, vb, .true.)
            if (id > 0) return
            id = -id
            do l = 1, q + 1
                samples(:, l, id) = chart(ra + (rb - ra)*(z(l) + 1)/2)
            end do
            call pin_ends(id)
        end subroutine

        subroutine straight_edge(va, vb)
            integer, intent(in) :: va, vb
            integer :: id, l
            id = find_edge(va, vb, .true.)
            if (id > 0) return
            id = -id
            do l = 1, q + 1
                samples(:, l, id) = vx(:, va) + (vx(:, vb) - vx(:, va))*(z(l) + 1)/2
            end do
            call pin_ends(id)
        end subroutine

        subroutine pin_ends(id)
            integer, intent(in) :: id
            samples(:, 1, id) = vx(:, edgev(1, id))
            samples(:, q + 1, id) = vx(:, edgev(2, id))
        end subroutine

        subroutine straight_cell(v)
            integer, intent(in) :: v(4)
            integer :: l
            do l = 1, 4
                call straight_edge(v(l), v(modulo(l, 4) + 1))
            end do
            ncell = ncell + 1
            cells(:, ncell) = v
        end subroutine

        subroutine chart_cell(v, bottom, top)
            ! Cell v1 v2 v3 v4 with chart points bottom = (v1, v2), top = (v3, v4).
            integer, intent(in) :: v(4)
            real(dp), intent(in) :: bottom(2, 2), top(2, 2)
            call chart_edge(v(1), v(2), bottom(:, 1), bottom(:, 2))
            call chart_edge(v(2), v(3), bottom(:, 2), top(:, 1))
            call chart_edge(v(3), v(4), top(:, 1), top(:, 2))
            call chart_edge(v(4), v(1), top(:, 2), bottom(:, 1))
            ncell = ncell + 1
            cells(:, ncell) = v
        end subroutine

        subroutine oriented(va, vb, x)
            integer, intent(in) :: va, vb
            real(dp), intent(out) :: x(2, q + 1)
            integer :: id
            id = find_edge(va, vb, .false.)
            if (edgev(1, id) == va) then
                x = samples(:, :, id)
            else
                x = samples(:, q + 1:1:-1, id)
            end if
        end subroutine

        subroutine coons(v, x)
            integer, intent(in) :: v(4)
            real(dp), intent(out) :: x(:, :)
            real(dp) :: bottom(2, q + 1), top(2, q + 1), left(2, q + 1), &
                        right(2, q + 1), sx, sy
            integer :: l, o
            call oriented(v(1), v(2), bottom)
            call oriented(v(4), v(3), top)
            call oriented(v(1), v(4), left)
            call oriented(v(2), v(3), right)
            do o = 1, q + 1
                sy = (z(o) + 1)/2
                do l = 1, q + 1
                    sx = (z(l) + 1)/2
                    x(:, l + (q + 1)*(o - 1)) = (1 - sy)*bottom(:, l) + &
                                                sy*top(:, l) + (1 - sx)*left(:, o) + &
                                                sx*right(:, o) - &
                                                ((1 - sx)*(1 - sy)*vx(:, v(1)) + &
                                                 sx*(1 - sy)*vx(:, v(2)) + &
                                                 sx*sy*vx(:, v(3)) + &
                                                 (1 - sx)*sy*vx(:, v(4)))
                end do
            end do
        end subroutine
    end subroutine
end module
