module fortfem_line_mesh
    ! Tensor meshes built from per-direction line tables (GLL, B-spline, Fourier).
    ! Cells are ordered e = c1 + n1*((c2-1) + n2*((c3-1) + n3*(patch-1))), with
    ! cell(:, e) the 1D cell of every direction. Local dofs are tensor ordered,
    ! a = a1 + nloc1*((a2-1) + nloc2*(a3-1)). Geometry coefficients per cell use
    ! Lagrange nodes of degree q in every local direction and, in single-cell
    ! Fourier directions, values at the line's nq equispaced quadrature points
    ! (trigonometric interpolation, independent of the solution's mode count).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_line_tables, only: tensor_line_tables_t, LINE_FOURIER
    use fortfem_gll_line, only: gll_rule, lagrange_table
    use fortfem_tensor_mesh, only: tensor_mesh_t, geometry_map
    use fortfem_tensor_contractions, only: tensor_contract
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    private
    public :: line_geometry_count, line_geometry_tables, unit_line
    public :: line_tensor_mesh, line_ogrid_mesh, line_mesh_prolongation
    public :: line_matrix_t, ogrid_disk_point
    type :: line_matrix_t
        real(dp), allocatable :: a(:, :)
    end type
    real(dp), parameter :: core = 0.32_dp, pi = acos(-1.0_dp)
contains
    subroutine unit_line(line)
        ! Trivial third direction of a two-dimensional mesh.
        type(tensor_line_tables_t), intent(out) :: line
        line%family = 0
        line%nelem = 1; line%nloc = 1; line%nq = 1; line%ndof = 1
        allocate (line%xq(1, 1), line%wq(1, 1), line%b(1, 1, 1), line%d(1, 1, 1), &
                  line%dof(1, 1))
        line%xq = 0; line%wq = 1; line%b = 1; line%d = 0; line%dof = 1
    end subroutine

    pure integer function line_geometry_count(line, q) result(count)
        type(tensor_line_tables_t), intent(in) :: line
        integer, intent(in) :: q
        if (line%family == 0) then
            count = 1
        else if (line%family == LINE_FOURIER) then
            count = line%nq
        else
            count = q + 1
        end if
    end function

    pure subroutine cell_interval(line, c, mid, h)
        ! Cell length and midpoint from the line's symmetric cell rule.
        type(tensor_line_tables_t), intent(in) :: line
        integer, intent(in) :: c
        real(dp), intent(out) :: mid, h
        h = sum(line%wq(:, c))
        mid = sum(line%wq(:, c)*line%xq(:, c))/h
    end subroutine

    subroutine line_geometry_tables(line, c, q, gb, gd, status, face)
        ! Geometry basis values and parameter derivatives at the line quadrature
        ! points of cell c, or at the cell end `face` (-1 or 1) when present.
        type(tensor_line_tables_t), intent(in) :: line
        integer, intent(in) :: c, q
        real(dp), allocatable, intent(out) :: gb(:, :), gd(:, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: face
        real(dp), allocatable :: nodes(:), w(:), ref(:)
        real(dp) :: mid, h
        status = 1
        if (line%family == 0) then
            allocate (gb(1, 1), gd(1, 1))
            gb = 1; gd = 0
        else if (line%family == LINE_FOURIER) then
            if (present(face)) return
            call spectral_tables(line%nq, line%x1 - line%x0, gb, gd)
        else
            call gll_rule(q, nodes, w, status)
            if (status /= 0) return
            call cell_interval(line, c, mid, h)
            if (present(face)) then
                ref = [face]
            else
                ref = 2*(line%xq(:, c) - mid)/h
            end if
            call lagrange_table(nodes, ref, gb, gd)
            gd = gd*2/h
        end if
        status = 0
    end subroutine

    subroutine geometry_samples(line, c, q, xi, projector, status)
        ! Sample parameters of a cell's geometry coefficients and the matrix that
        ! turns samples into coefficients (identity for Lagrange nodes, discrete
        ! L2 projection with the line quadrature for Fourier directions).
        type(tensor_line_tables_t), intent(in) :: line
        integer, intent(in) :: c, q
        real(dp), allocatable, intent(out) :: xi(:, :), projector(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: nodes(:), w(:), gram(:, :)
        real(dp) :: mid, h
        integer :: i
        status = 1
        if (line%family == 0) then
            allocate (xi(1, 1), projector(1, 1))
            xi = 0
            projector = 1
        else if (line%family == LINE_FOURIER) then
            xi = reshape(line%xq(:, c), [line%nq, 1])
            call spectral_tables(line%nq, line%x1 - line%x0, projector, gram)
        else
            call gll_rule(q, nodes, w, status)
            if (status /= 0) return
            call cell_interval(line, c, mid, h)
            xi = reshape(mid + h*nodes/2, [q + 1, 1])
            ! Cell ends rebuilt from the rule may overshoot by an ulp.
            if (.not. line%periodic) xi = min(max(xi, line%x0), line%x1)
            allocate (projector(q + 1, q + 1))
            projector = 0
            do i = 1, q + 1
                projector(i, i) = 1
            end do
        end if
        status = 0
    end subroutine

    subroutine spectral_tables(n, period, b, d)
        ! Identity values and the periodic spectral differentiation matrix on n
        ! equispaced points x_j = (j-1) period/n (trigonometric interpolation).
        integer, intent(in) :: n
        real(dp), intent(in) :: period
        real(dp), allocatable, intent(out) :: b(:, :), d(:, :)
        real(dp) :: h, angle
        integer :: i, j
        allocate (b(n, n), d(n, n))
        b = 0
        d = 0
        h = 2*pi/n
        do j = 1, n
            b(j, j) = 1
            do i = 1, n
                if (i == j) cycle
                angle = (i - j)*h/2
                if (modulo(n, 2) == 0) then
                    d(i, j) = 0.5_dp*(-1)**(i - j)/tan(angle)
                else
                    d(i, j) = 0.5_dp*(-1)**(i - j)/sin(angle)
                end if
            end do
        end do
        d = d*2*pi/period
    end subroutine

    subroutine evaluate_point(dim, block, x0, length, xi, x, map)
        ! Physical point of parameter xi. block < 0: x = map(xi(1:dim)) (or xi);
        ! block 0-4: O-grid patch, x = map(rx, ry[, phi]) (or (rx, ry)).
        integer, intent(in) :: dim, block
        real(dp), intent(in) :: x0, length, xi(3)
        real(dp), intent(out) :: x(:)
        procedure(geometry_map), optional :: map
        real(dp) :: rx, ry
        if (block < 0) then
            if (present(map)) then
                call map(xi(1:dim), x)
            else
                x = xi(1:dim)
            end if
            return
        end if
        call ogrid_disk_point(block, (xi(1) - x0)/length, 2*(xi(2) - x0)/length - 1, &
                              rx, ry)
        if (present(map)) then
            if (dim == 3) then
                call map([rx, ry, xi(3)], x)
            else
                call map([rx, ry], x)
            end if
        else
            x = [rx, ry]
        end if
    end subroutine

    subroutine cell_geometry(lines, c, q, dim, block, x0, length, geometry, status, map)
        ! Geometry coefficients of one cell from the parameter-space point map.
        type(tensor_line_tables_t), intent(in) :: lines(3)
        integer, intent(in) :: c(3), q, dim, block
        real(dp), intent(in) :: x0, length
        real(dp), intent(out) :: geometry(:, :)
        integer, intent(out) :: status
        procedure(geometry_map), optional :: map
        type(line_matrix_t) :: xi(3), projector(3)
        real(dp), allocatable :: samples(:), values(:), t1(:), t2(:)
        real(dp) :: x(dim)
        integer :: ns(3), ng(3), i, j, k, comp, ws
        do i = 1, 3
            call geometry_samples(lines(i), c(i), q, xi(i)%a, projector(i)%a, &
                                  status)
            if (status /= 0) return
            ns(i) = size(projector(i)%a, 2)
            ng(i) = size(projector(i)%a, 1)
        end do
        ws = maxval([ns, ng])**3
        allocate (samples(product(ns)), values(product(ng)), t1(ws), t2(ws))
        do comp = 1, dim
            do k = 1, ns(3)
                do j = 1, ns(2)
                    do i = 1, ns(1)
                        call evaluate_point(dim, block, x0, length, &
                                        [xi(1)%a(i, 1), xi(2)%a(j, 1), xi(3)%a(k, 1)], &
                                            x, map)
                        samples(i + ns(1)*(j - 1 + ns(2)*(k - 1))) = x(comp)
                    end do
                end do
            end do
            call tensor_contract(ns, ng, projector(1)%a, projector(2)%a, &
                                 projector(3)%a, samples, values, t1, t2)
            geometry(comp, :) = values
        end do
        status = 0
    end subroutine

    subroutine prepare_lines(dim, lines, mesh, status)
        integer, intent(in) :: dim
        type(tensor_line_tables_t), intent(in) :: lines(:)
        type(tensor_mesh_t), intent(inout) :: mesh
        integer, intent(out) :: status
        integer :: i
        status = 1
        if (dim < 2 .or. dim > 3 .or. size(lines) < dim) return
        do i = 1, dim
            if (lines(i)%nelem < 1 .or. lines(i)%nloc < 1) return
            if (.not. allocated(lines(i)%dof)) return
            mesh%line(i) = lines(i)
        end do
        if (dim == 2) call unit_line(mesh%line(3))
        mesh%dim = dim
        mesh%degree = maxval(mesh%line(1:dim)%degree)
        ! Default geometry degree: the highest local (non-Fourier) degree.
        mesh%geometry_degree = 1
        do i = 1, dim
            if (mesh%line(i)%family /= LINE_FOURIER) &
                mesh%geometry_degree = max(mesh%geometry_degree, mesh%line(i)%degree)
        end do
        status = 0
    end subroutine

    subroutine line_tensor_mesh(dim, lines, mesh, status, map, q)
        ! Single structured patch over the tensor product of `lines`. Dirichlet
        ! dofs are the end dofs of non-periodic lines (GLL or open B-spline).
        integer, intent(in) :: dim
        type(tensor_line_tables_t), intent(in) :: lines(:)
        type(tensor_mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        procedure(geometry_map), optional :: map
        integer, intent(in), optional :: q
        type(line_matrix_t) :: anchor(3)
        integer :: nc(3), nd(3), nl(3), e, a, c(3), loc(3), g(3), i, ng
        call prepare_lines(dim, lines, mesh, status)
        if (status /= 0) return
        if (present(q)) mesh%geometry_degree = q
        nc = mesh%line%nelem; nd = mesh%line%ndof; nl = mesh%line%nloc
        mesh%nelem = product(nc)
        mesh%ndof = product(nd)
        ng = 1
        do i = 1, 3
            ng = ng*line_geometry_count(mesh%line(i), mesh%geometry_degree)
        end do
        allocate (mesh%cell(3, mesh%nelem), mesh%patch(mesh%nelem), &
                mesh%dof(product(nl), mesh%nelem), mesh%geometry(dim, ng, mesh%nelem), &
                  mesh%boundary(mesh%ndof), mesh%vertices(0, mesh%nelem))
        mesh%patch = 1
        call line_anchors(mesh, anchor)
        do e = 1, mesh%nelem
            c = [modulo(e - 1, nc(1)), modulo((e - 1)/nc(1), nc(2)), &
                 (e - 1)/(nc(1)*nc(2))] + 1
            mesh%cell(:, e) = c
            do a = 1, product(nl)
                loc = [modulo(a - 1, nl(1)), modulo((a - 1)/nl(1), nl(2)), &
                       (a - 1)/(nl(1)*nl(2))] + 1
                do i = 1, 3
                    g(i) = mesh%line(i)%dof(loc(i), c(i))
                end do
                mesh%dof(a, e) = g(1) + nd(1)*(g(2) - 1 + nd(2)*(g(3) - 1))
            end do
            call cell_geometry(mesh%line, c, mesh%geometry_degree, dim, -1, 0.0_dp, &
                               1.0_dp, mesh%geometry(:, :, e), status, map)
            if (status /= 0) return
            call set_anchors(mesh, e, anchor, -1, 0.0_dp, 1.0_dp, map)
        end do
        mesh%boundary = .false.
        do a = 1, mesh%ndof
            g = [modulo(a - 1, nd(1)), modulo((a - 1)/nd(1), nd(2)), &
                 (a - 1)/(nd(1)*nd(2))] + 1
            do i = 1, dim
                if (mesh%line(i)%periodic) cycle
                if (g(i) == 1 .or. g(i) == nd(i)) mesh%boundary(a) = .true.
            end do
        end do
        status = 0
    end subroutine

    pure subroutine ogrid_disk_point(block, radial, tangent, rx, ry)
        ! MF1's five-block disk O-grid: block 0 is the core square, blocks 1-4
        ! blend its sides to the unit circle. radial in [0,1], tangent in [-1,1].
        integer, intent(in) :: block
        real(dp), intent(in) :: radial, tangent
        real(dp), intent(out) :: rx, ry
        real(dp) :: angle, tmp
        integer :: i
        if (block == 0) then
            rx = core*(2*radial - 1)
            ry = core*tangent
        else
            angle = pi*tangent/4
            rx = (1 - radial)*core + radial*cos(angle)
            ry = (1 - radial)*core*tangent + radial*sin(angle)
            do i = 2, block
                tmp = rx; rx = -ry; ry = tmp
            end do
        end if
    end subroutine

    subroutine line_ogrid_mesh(dim, lines, mesh, status, map, q)
        ! Five conforming (C0 across interfaces) patches through the axis, MF1
        ! O-grid layout. lines(1) is used in both block directions and must have
        ! a symmetric knot/node layout (uniform); lines(3) is the phi direction.
        ! map(xi, x) receives xi = (rx, ry, phi) and is required in 3D.
        integer, intent(in) :: dim
        type(tensor_line_tables_t), intent(in) :: lines(:)
        type(tensor_mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        procedure(geometry_map), optional :: map
        integer, intent(in), optional :: q
        real(dp) :: corner(2, 8), xy(2), length
        integer :: vid(2, 2, 0:4), edge(8, 8), nv, ne, n, nn, n2d, nl(3), nc(3)
        integer :: block, i, j, k, v, e, a, c(3), loc(3), local, id
        integer, allocatable :: patch_dof(:, :, :)
        type(tensor_line_tables_t) :: block_lines(3)
        type(line_matrix_t) :: anchor(3)
        status = 1
        if (dim == 3 .and. .not. present(map)) return
        block_lines(1) = lines(1)
        block_lines(2) = lines(1)
        if (dim == 3) block_lines(3) = lines(size(lines))
        call prepare_lines(dim, block_lines, mesh, status)
        if (status /= 0) return
        status = 1
        if (mesh%line(1)%periodic) return
        if (present(q)) mesh%geometry_degree = q
        n = mesh%line(1)%ndof
        if (n < 2) return
        nl = mesh%line%nloc; nc = mesh%line%nelem
        length = mesh%line(1)%x1 - mesh%line(1)%x0
        ! Patch corners, identified by position (4 core corners, 4 on the circle).
        nv = 0
        do block = 0, 4
            do j = 0, 1
                do i = 0, 1
                    call ogrid_disk_point(block, real(i, dp), real(2*j - 1, dp), &
                                          xy(1), xy(2))
                    do v = 1, nv
                        if (maxval(abs(corner(:, v) - xy)) < 1.e-12_dp) exit
                    end do
                    if (v > nv) then
                        nv = v
                        corner(:, v) = xy
                    end if
                    vid(i + 1, j + 1, block) = v
                end do
            end do
        end do
        ! Global 2D numbering: vertices, then edge interiors, then patch interiors.
        allocate (patch_dof(n, n, 0:4))
        edge = 0
        ne = 0
        do block = 0, 4
            do j = 1, n
                do i = 1, n
                    patch_dof(i, j, block) = dof2d(block, i, j)
                end do
            end do
        end do
        nn = nv + ne*(n - 2)
        n2d = nn + 5*(n - 2)**2
        do block = 0, 4
            do j = 2, n - 1
                do i = 2, n - 1
                    patch_dof(i, j, block) = nn + block*(n - 2)**2 + &
                                             (i - 1) + (n - 2)*(j - 2)
                end do
            end do
        end do
        mesh%nelem = 5*nc(1)*nc(2)*nc(3)
        mesh%ndof = n2d*mesh%line(3)%ndof
        nn = 1
        do i = 1, 3
            nn = nn*line_geometry_count(mesh%line(i), mesh%geometry_degree)
        end do
        allocate (mesh%cell(3, mesh%nelem), mesh%patch(mesh%nelem), &
                mesh%dof(product(nl), mesh%nelem), mesh%geometry(dim, nn, mesh%nelem), &
                  mesh%boundary(mesh%ndof), mesh%vertices(0, mesh%nelem))
        mesh%boundary = .false.
        call line_anchors(mesh, anchor)
        do e = 1, mesh%nelem
            local = modulo(e - 1, product(nc)) + 1
            block = (e - 1)/product(nc)
            c = [modulo(local - 1, nc(1)), modulo((local - 1)/nc(1), nc(2)), &
                 (local - 1)/(nc(1)*nc(2))] + 1
            mesh%cell(:, e) = c
            mesh%patch(e) = block + 1
            do a = 1, product(nl)
                loc = [modulo(a - 1, nl(1)), modulo((a - 1)/nl(1), nl(2)), &
                       (a - 1)/(nl(1)*nl(2))] + 1
                i = mesh%line(1)%dof(loc(1), c(1))
                j = mesh%line(2)%dof(loc(2), c(2))
                k = mesh%line(3)%dof(loc(3), c(3))
                id = patch_dof(i, j, block) + n2d*(k - 1)
                mesh%dof(a, e) = id
                if (block > 0 .and. i == n) mesh%boundary(id) = .true.
            end do
            call cell_geometry(mesh%line, c, mesh%geometry_degree, dim, block, &
                           mesh%line(1)%x0, length, mesh%geometry(:, :, e), status, map)
            if (status /= 0) return
            call set_anchors(mesh, e, anchor, block, mesh%line(1)%x0, length, map)
        end do
        status = 0
    contains
        integer function dof2d(b, i, j) result(id)
            integer, intent(in) :: b, i, j
            integer :: va, vb, pos
            id = 0
            if ((i == 1 .or. i == n) .and. (j == 1 .or. j == n)) then
                id = vid(merge(1, 2, i == 1), merge(1, 2, j == 1), b)
                return
            end if
            if (i > 1 .and. i < n .and. j > 1 .and. j < n) return
            if (i == 1 .or. i == n) then
                va = vid(merge(1, 2, i == 1), 1, b)
                vb = vid(merge(1, 2, i == 1), 2, b)
                pos = j
            else
                va = vid(1, merge(1, 2, j == 1), b)
                vb = vid(2, merge(1, 2, j == 1), b)
                pos = i
            end if
            if (va > vb) then
                pos = n + 1 - pos
                v = va; va = vb; vb = v
            end if
            if (edge(va, vb) == 0) then
                ne = ne + 1
                edge(va, vb) = ne
            end if
            id = nv + (edge(va, vb) - 1)*(n - 2) + pos - 1
        end function
    end subroutine

    subroutine line_anchors(mesh, anchor)
        ! Diagnostic 1D dof anchors: each basis function's centroid (the cell
        ! start for Fourier modes). They are not interpolation nodes.
        type(tensor_mesh_t), intent(inout) :: mesh
        type(line_matrix_t), intent(out) :: anchor(3)
        real(dp), allocatable :: mass(:)
        integer :: d, c, i, j
        do d = 1, 3
            allocate (anchor(d)%a(mesh%line(d)%ndof, 1), mass(mesh%line(d)%ndof))
            anchor(d)%a = 0; mass = 0
            do c = 1, mesh%line(d)%nelem
                do i = 1, mesh%line(d)%nloc
                    j = mesh%line(d)%dof(i, c)
                  mass(j) = mass(j) + sum(mesh%line(d)%wq(:, c)*mesh%line(d)%b(:, i, c))
                    anchor(d)%a(j, 1) = anchor(d)%a(j, 1) + sum(mesh%line(d)%wq(:, c)* &
                                          mesh%line(d)%xq(:, c)*mesh%line(d)%b(:, i, c))
                end do
            end do
            if (mesh%line(d)%family == LINE_FOURIER) then
                anchor(d)%a = mesh%line(d)%x0
            else
             where (abs(mass) > tiny(1.0_dp)) anchor(d)%a(:, 1) = anchor(d)%a(:, 1)/mass
            end if
            deallocate (mass)
        end do
        allocate (mesh%coordinates(mesh%dim, mesh%ndof))
        mesh%coordinates = 0
    end subroutine

    subroutine set_anchors(mesh, e, anchor, block, x0, length, map)
        type(tensor_mesh_t), intent(inout) :: mesh
        integer, intent(in) :: e, block
        type(line_matrix_t), intent(in) :: anchor(3)
        real(dp), intent(in) :: x0, length
        procedure(geometry_map), optional :: map
        real(dp) :: x(mesh%dim), xi(3)
        integer :: a, d, loc(3), nl(3)
        nl = mesh%line%nloc
        do a = 1, size(mesh%dof, 1)
            loc = [modulo(a - 1, nl(1)), modulo((a - 1)/nl(1), nl(2)), &
                   (a - 1)/(nl(1)*nl(2))] + 1
            do d = 1, 3
                xi(d) = anchor(d)%a(mesh%line(d)%dof(loc(d), mesh%cell(d, e)), 1)
            end do
            call evaluate_point(mesh%dim, block, x0, length, xi, x, map)
            mesh%coordinates(:, mesh%dof(a, e)) = x
        end do
    end subroutine

    subroutine line_mesh_prolongation(fine, coarse, transfer, ratio, p, status)
        ! Global nested prolongation P (fine x coarse) from per-direction 1D
        ! matrices transfer(d)%a (fine line dofs x coarse line dofs). Fine cell c
        ! lies in coarse cell (c-1)/ratio+1 of the same patch. Each fine row is
        ! formed once, from its first cell; conforming interface rows agree.
        type(tensor_mesh_t), intent(in) :: fine, coarse
        type(line_matrix_t), intent(in) :: transfer(3)
        integer, intent(in) :: ratio(3)
        type(csr_matrix_t), intent(out) :: p
        integer, intent(out) :: status
        integer, allocatable :: count(:), cols(:, :)
        real(dp), allocatable :: vals(:, :)
        logical, allocatable :: done(:)
        integer :: nf(3), nc(3), lf(3), lc(3), e, ec, a, b, i, k, row, col, width
        integer :: cf(3), cc(3), ia(3), ib(3), d, nnz, pos
        real(dp) :: value
        status = 1
        if (.not. allocated(fine%cell) .or. .not. allocated(coarse%cell)) return
        nf = fine%line%nelem; nc = coarse%line%nelem
        lf = fine%line%nloc; lc = coarse%line%nloc
        do d = 1, 3
            if (nf(d) /= nc(d)*ratio(d)) return
            if (fine%line(d)%family == 0) cycle
            if (size(transfer(d)%a, 1) /= fine%line(d)%ndof) return
            if (size(transfer(d)%a, 2) /= coarse%line(d)%ndof) return
        end do
        width = product(lc)
        allocate (count(fine%ndof), cols(width, fine%ndof), vals(width, fine%ndof), &
                  done(fine%ndof))
        count = 0
        done = .false.
        do e = 1, fine%nelem
            cf = fine%cell(:, e)
            cc = (cf - 1)/ratio + 1
          ec = cc(1) + nc(1)*(cc(2) - 1 + nc(2)*(cc(3) - 1 + nc(3)*(fine%patch(e) - 1)))
            if (coarse%patch(ec) /= fine%patch(e)) return
            do a = 1, product(lf)
                row = fine%dof(a, e)
                if (done(row)) cycle
                done(row) = .true.
                ia = [modulo(a - 1, lf(1)), modulo((a - 1)/lf(1), lf(2)), &
                      (a - 1)/(lf(1)*lf(2))] + 1
                do b = 1, width
                    ib = [modulo(b - 1, lc(1)), modulo((b - 1)/lc(1), lc(2)), &
                          (b - 1)/(lc(1)*lc(2))] + 1
                    value = 1
                    do d = 1, fine%dim
                        value = value*transfer(d)%a(fine%line(d)%dof(ia(d), cf(d)), &
                                                    coarse%line(d)%dof(ib(d), cc(d)))
                    end do
                    if (abs(value) < tiny(1.0_dp)) cycle
                    col = coarse%dof(b, ec)
                    pos = 0
                    do k = 1, count(row)
                        if (cols(k, row) == col) pos = k
                    end do
                    if (pos == 0) then
                        count(row) = count(row) + 1
                        pos = count(row)
                        cols(pos, row) = col
                        vals(pos, row) = 0
                    end if
                    vals(pos, row) = vals(pos, row) + value
                end do
            end do
        end do
        if (.not. all(done)) return
        nnz = sum(count)
        p%n = fine%ndof
        allocate (p%rowptr(fine%ndof + 1), p%colidx(nnz), p%val(nnz))
        p%rowptr(1) = 1
        do i = 1, fine%ndof
            p%rowptr(i + 1) = p%rowptr(i) + count(i)
            p%colidx(p%rowptr(i):p%rowptr(i + 1) - 1) = cols(1:count(i), i)
            p%val(p%rowptr(i):p%rowptr(i + 1) - 1) = vals(1:count(i), i)
        end do
        status = 0
    end subroutine
end module
