module fortfem_verified_bernstein_majorant
    !! Outward-rounded Prager--Synge audit of a sampled RT_k equilibrated flux
    !! (k = 1, 2, 3) on affine triangles, for -div(A grad u) = f with
    !! A = 1/w, w > 0 affine on every cell (vertex values), f = s/w^m with s a
    !! polynomial supplied by its Bernstein coefficients (m = 0, 1, 2), or
    !! s = polynomial + delta with a caller-proved |delta| <= rho per cell, and
    !! u = 0 on the whole boundary. For conforming u_h in H1_0 and any sigma
    !! in H(div),
    !!
    !!   ||u - u_h||_A <= ||w^(1/2)(sigma - grad u_h/w)|| + C_F ||f + div sigma||.
    !!
    !! The audited sigma is not the floating RT_k field itself. On every cell
    !! its degree-(k + 1) Bernstein coefficients are fitted in floating point;
    !! then one exactly represented normal component is selected per global
    !! edge and Bernstein position, and the mathematical field is defined by
    !! replacing the normal components of the edge coefficients by these shared
    !! values (two normals fix a corner coefficient). The restriction of a
    !! Bernstein polynomial to an edge depends only on its edge coefficients,
    !! so the normal trace of this exact BDM_(k+1) field is continuous and the
    !! field lies in H(div). Fitting errors only loosen the bound.
    !!
    !! With q = w sigma - grad u_h and t = s + w^m div sigma (polynomials),
    !! int |q|^2/w <= int |q|^2 l_1 and int t^2/w^(2m) <= int t^2 l_2m, where
    !! l_j is the chord of the convex w^-j between the cell extremes of w, an
    !! affine upper bound. The polynomial integrals are evaluated exactly by
    !! Bernstein moments in outward interval arithmetic, so the result
    !! encloses quadrature and arithmetic error; there is no sampled
    !! quadrature, assumed exactness or projection. Leaves come from the
    !! FortSym-generated module fortfem_verified_bernstein_products.
    !!
    !! The caller guarantees a conforming, nonoverlapping, counter-clockwise
    !! affine mesh, a discrete solution vanishing on the boundary, valid source
    !! coefficient/remainder enclosures and a valid Friedrichs enclosure C_F with
    !! ||v|| <= C_F ||v||_A on H1_0. Mixed or Neumann boundaries are outside
    !! this audit.
    use fortfem_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, operator(+), operator(-), &
                                operator(*), sqrt
    use fortfem_equilibrated_flux_topology, only: topology_t, build_topology, dgesv
    use fortfem_triangle_rt_arbitrary_order, only: triangle_rt_basis_t, &
        initialize_triangle_raviart_thomas, evaluate_triangle_raviart_thomas, &
        triangle_rt_dof_count
    use fortfem_verified_bernstein_products, only: verified_bernstein_geometry, &
        verified_bernstein_edge_normal, verified_bernstein_edge_projection, &
        verified_bernstein_corner, verified_bernstein_affine_product, &
        verified_bernstein_elevate, verified_bernstein_gradient, &
        verified_bernstein_square_integral, verified_bernstein_chord, &
        verified_bernstein_reciprocal, verified_bernstein_lagrange_p2, &
        verified_bernstein_lagrange_p3, verified_bernstein_source_tail_0, &
        verified_bernstein_source_tail_1, verified_bernstein_source_tail_2
    use fortfem_verified_bernstein_products, only: mapped_product, mapped_elevate, &
        mapped_square_integral, verified_mapped_divide, verified_mapped_tail
    implicit none
    private
    public :: verified_bernstein_majorant, verified_bernstein_source
    public :: verified_mapped_bernstein_majorant
    public :: verified_bernstein_source_remainder

    integer, parameter, public :: VERIFIED_BERNSTEIN_OK = 0
    integer, parameter, public :: VERIFIED_BERNSTEIN_BAD_INPUT = 1
    integer, parameter, public :: VERIFIED_BERNSTEIN_BAD_CELL = 2
    integer, parameter, public :: VERIFIED_BERNSTEIN_BAD_SOURCE = 3
    integer, parameter :: max_coefficients = 28

    abstract interface
        subroutine verified_bernstein_source(cell, coefficients, ok)
            !! Outward Bernstein coefficients of the polynomial numerator (or
            !! its approximation when source_remainder is supplied), using
            !! the vertex order of `triangles`. Called
            !! concurrently from OpenMP threads.
            import :: interval_t
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: coefficients(:)
            logical, intent(out) :: ok
        end subroutine verified_bernstein_source
        subroutine verified_bernstein_source_remainder(cell, radius, ok)
            !! Nonnegative outward upper bound radius%hi for the absolute
            !! numerator remainder on the WHOLE cell, not sampled points.
            !! The source polynomial is supplied by verified_bernstein_source.
            !! The caller proves this contract. Called concurrently.
            import :: interval_t
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: radius
            logical, intent(out) :: ok
        end subroutine verified_bernstein_source_remainder
    end interface

    !> Floating fit from RT_k coefficients to Bernstein coefficients.
    type :: fit_t
        integer :: degree = 0, nrt = 0, nb = 0
        real(dp), allocatable :: rt(:, :, :)     !< (2, nrt, nb) at lattice points
        real(dp), allocatable :: inverse(:, :)   !< (nb, nb) Bernstein from values
        integer, allocatable :: index(:, :)      !< (3, nb) multi-indices
    end type fit_t

    type :: mapped_poly_t
        integer :: degree = 0
        type(interval_t) :: c(66)
    end type mapped_poly_t

    type, abstract, public :: bernstein_source_t
    contains
        procedure(bernstein_source_t_binding), deferred :: evaluate
    end type
    abstract interface
        subroutine bernstein_source_t_binding(self, cell, coefficients, ok)
            import :: bernstein_source_t, interval_t
            class(bernstein_source_t), intent(in) :: self
            !! Outward Bernstein coefficients of the polynomial numerator (or
            !! its approximation when source_remainder is supplied), using
            !! the vertex order of `triangles`. Called
            !! concurrently from OpenMP threads.
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: coefficients(:)
            logical, intent(out) :: ok
        end subroutine
    end interface

    type, extends(bernstein_source_t) :: bernstein_source_adapter_t
        procedure(verified_bernstein_source), pointer, nopass :: callback => null()
    contains
        procedure :: evaluate => bernstein_source_adapter_t_evaluate
    end type

    type, abstract, public :: bernstein_remainder_t
    contains
        procedure(bernstein_remainder_t_binding), deferred :: evaluate
    end type
    abstract interface
        subroutine bernstein_remainder_t_binding(self, cell, radius, ok)
            import :: bernstein_remainder_t, interval_t
            class(bernstein_remainder_t), intent(in) :: self
            !! Nonnegative outward upper bound radius%hi for the absolute
            !! numerator remainder on the WHOLE cell, not sampled points.
            !! The source polynomial is supplied by verified_bernstein_source.
            !! The caller proves this contract. Called concurrently.
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: radius
            logical, intent(out) :: ok
        end subroutine
    end interface

    type, extends(bernstein_remainder_t) :: bernstein_remainder_adapter_t
        procedure(verified_bernstein_source_remainder), pointer, nopass :: callback &
            => null()
    contains
        procedure :: evaluate => bernstein_remainder_adapter_t_evaluate
    end type

    interface verified_mapped_bernstein_majorant
        module procedure verified_mapped_bernstein_majorant_object, &
            verified_mapped_bernstein_majorant_procedure
    end interface

    interface verified_bernstein_majorant
        module procedure verified_bernstein_majorant_object, &
            verified_bernstein_majorant_procedure
    end interface

contains

    subroutine verified_mapped_bernstein_majorant_object(vertices, triangles, degree, &
                          flux_coefficients, solution_cells, solution, geometry_nodes, &
                        weight_nodes, source, source_degree, weight_power, friedrichs, &
                        bound, flux_squared, residual_squared, status, source_remainder)
        !! Outward audit on conforming isoparametric P1/P2/P3 triangles.
        !! weight_nodes is cell-major: n_geometry_nodes values for each cell.
        !! It defines w=1/A by the same reference Lagrange interpolation as
        !! geometry. Source supplies a reference Bernstein numerator s in
        !! f=(s+delta)/w**weight_power, |delta|<=source_remainder on each cell.
        !! No quadrature is sampled. Exact polynomial square moments divided
        !! by positive denominator hulls give lower/upper integral budgets.
        !! The caller supplies nonoverlap, a conforming H1_0 candidate and C_F.
        !! Shared mapped edges are checked node for node, including orientation.
        !! A nonpositive Bernstein hull fails closed (a sufficient test only).
        real(dp), intent(in) :: vertices(:, :), flux_coefficients(:, :)
        real(dp), intent(in) :: solution(:), geometry_nodes(:, :, :)
        real(dp), intent(in) :: weight_nodes(:)
        integer, intent(in) :: triangles(:, :), solution_cells(:, :)
        integer, intent(in) :: degree, source_degree, weight_power
        class(bernstein_source_t), intent(in) :: source
        class(bernstein_remainder_t), intent(in), optional :: source_remainder
        type(interval_t), intent(in) :: friedrichs
        type(interval_t), intent(out) :: bound, flux_squared, residual_squared
        integer, intent(out) :: status
        type(topology_t) :: topo
        type(fit_t) :: fit
        real(dp), parameter :: normals(2, 3) = reshape( &
                               [0._dp, -1._dp, 1._dp, 1._dp, -1._dp, 0._dp], [2, 3])
        real(dp), allocatable :: shared(:, :), edge_nodes(:, :, :)
        integer, allocatable :: owner(:), uses(:), turns(:)
        type(interval_t) :: fp, rp, fc, rc
        integer :: nt, nv, ng, pg, pu, t, e, edge, stat, j, k, slot, orient
        integer :: failed, cell_status, thread_failed

        status = VERIFIED_BERNSTEIN_BAD_INPUT
        bound = interval(0._dp, huge(1._dp))
        flux_squared = bound; residual_squared = bound
        if (size(vertices, 1) /= 2 .or. size(triangles, 1) /= 3) return
        nv = size(vertices, 2); nt = size(triangles, 2)
        if (nv < 3 .or. nt < 1) return
        if (size(geometry_nodes, 1) /= 2) return
        ng = size(geometry_nodes, 2)
        pg = lagrange_degree(ng); pu = lagrange_degree(size(solution_cells, 1))
        if (pg == 0 .or. pu == 0) return
        if (degree < 1 .or. degree > 3) return
        if (source_degree < 0 .or. source_degree > 6) return
        if (weight_power < 0 .or. weight_power > 2) return
        if (size(geometry_nodes, 3) /= nt) return
        if (size(weight_nodes) /= ng*nt) return
        if (size(solution_cells, 2) /= nt) return
        if (size(flux_coefficients, 2) /= nt) return
        if (any(triangles < 1) .or. any(triangles > nv)) return
        if (any(solution_cells < 1) .or. any(solution_cells > size(solution))) return
        if (.not. all(ieee_is_finite(vertices))) return
        if (.not. all(ieee_is_finite(geometry_nodes))) return
        if (.not. all(ieee_is_finite(weight_nodes))) return
        if (.not. all(ieee_is_finite(solution))) return
        if (.not. all(ieee_is_finite(flux_coefficients))) return
        if (.not. valid(friedrichs)) return
        if (friedrichs%lo < 0) return
        call build_fit(degree, fit, stat)
        if (stat /= 0) return
        if (size(flux_coefficients, 1) /= fit%nrt) return
        call build_topology(nv, triangles, topo, stat)
        if (stat /= 0) return
        allocate (owner(maxval(topo%cell_edge)), uses(maxval(topo%cell_edge)), &
                  turns(maxval(topo%cell_edge)))
        allocate (edge_nodes(2, pg + 1, size(owner)))
        owner = 0; uses = 0; turns = 0
        do t = 1, nt
            if (any(geometry_nodes(:, :3, t) /= vertices(:, triangles(:, t)))) return
            do e = 1, 3
                edge = topo%cell_edge(e, t)
                orient = topo%cell_orientation(e, t)
                uses(edge) = uses(edge) + 1
                turns(edge) = turns(edge) + orient
                do j = 0, pg
                    k = j
                    if (orient < 0) k = pg - j
                    slot = geometry_edge_slot(pg, e, k)
                    if (owner(edge) == 0) then
                        edge_nodes(:, j + 1, edge) = geometry_nodes(:, slot, t)
                    else
                        if (any(edge_nodes(:, j + 1, edge) &
                                /= geometry_nodes(:, slot, t))) return
                    end if
                end do
                if (owner(edge) == 0) owner(edge) = t
            end do
        end do
        if (any(uses > 2)) return
        if (any(uses == 2 .and. turns /= 0)) return
        deallocate (edge_nodes)
        allocate (shared(fit%degree + 1, size(owner)))
        !$omp parallel do schedule(static)
        do edge = 1, size(owner)
            call select_normal(edge)
        end do
        !$omp end parallel do
        flux_squared = interval(0._dp); residual_squared = interval(0._dp)
        failed = 0
        !$omp parallel private(fp,rp,fc,rc,cell_status,thread_failed)
        fp = interval(0._dp); rp = interval(0._dp); thread_failed = 0
        !$omp do schedule(static)
        do t = 1, nt
            call audit(t, fc, rc, cell_status, source_remainder)
            if (cell_status /= 0) then
                thread_failed = max(thread_failed, cell_status)
            else
                fp = fp + fc; rp = rp + rc
            end if
        end do
        !$omp end do
        !$omp critical (verified_mapped_sum)
        flux_squared = flux_squared + fp
        residual_squared = residual_squared + rp
        failed = max(failed, thread_failed)
        !$omp end critical (verified_mapped_sum)
        !$omp end parallel
        if (failed /= 0) then
            status = failed
            flux_squared = bound; residual_squared = bound
            return
        end if
        flux_squared%lo = max(0._dp, flux_squared%lo)
        residual_squared%lo = max(0._dp, residual_squared%lo)
        bound = sqrt(flux_squared) + friedrichs*sqrt(residual_squared)
        status = VERIFIED_BERNSTEIN_BAD_CELL
        if (.not. valid(bound)) then
            bound = interval(0._dp, huge(1._dp))
            return
        end if
        status = VERIFIED_BERNSTEIN_OK
    contains
        subroutine fitted_reference(cell, c)
            integer, intent(in) :: cell
            real(dp), intent(out) :: c(2, max_coefficients)
            real(dp) :: v(2, max_coefficients)
            integer :: i, q
            v = 0; c = 0
            do q = 1, fit%nb
                do i = 1, fit%nrt
                    v(:, q) = v(:, q) + fit%rt(:, i, q)*flux_coefficients(i, cell)
                end do
            end do
            do i = 1, fit%nb
                do q = 1, fit%nb
                    c(:, i) = c(:, i) + fit%inverse(i, q)*v(:, q)
                end do
            end do
        end subroutine

        integer function edge_position(cell, local, j) result(k)
            integer, intent(in) :: cell, local, j
            integer :: a(3), next, jj
            next = mod(local, 3) + 1
            jj = j
            if (topo%cell_orientation(local, cell) < 0) jj = fit%degree - j
            a = 0; a(local) = fit%degree - jj; a(next) = jj
            k = position(fit, a)
        end function

        subroutine select_normal(edge)
            integer, intent(in) :: edge
            real(dp) :: c(2, max_coefficients)
            integer :: cell, local, j, k
            cell = owner(edge)
            do local = 1, 3
                if (topo%cell_edge(local, cell) == edge) exit
            end do
            call fitted_reference(cell, c)
            do j = 0, fit%degree
                k = edge_position(cell, local, j)
                shared(j + 1, edge) = topo%cell_orientation(local, cell) &
                                      *dot_product(normals(:, local), c(:, k))
            end do
        end subroutine

        subroutine audit(cell, flux_value, residual_value, cell_status, tail)
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: flux_value, residual_value
            integer, intent(out) :: cell_status
            class(bernstein_remainder_t), intent(in), optional :: tail
            type(mapped_poly_t) :: x, y, w, u, j11, j12, j21, j22, det
            type(mapped_poly_t) :: sx, sy, ux, uy, dummy, qx, qy, dx, dy
            type(mapped_poly_t) :: div, s, wm, r
            type(interval_t) :: dh, wh, wmh, integral_x, integral_y, rho, enclosed
            type(interval_t) :: normal1(2), normal2(2), value1, value2
            real(dp) :: c(2, max_coefficients)
            real(dp) :: values(10)
            integer :: k, a(3), zeros, vertex, e1, e2, local, next, slot1, slot2
            logical :: ok

            cell_status = VERIFIED_BERNSTEIN_BAD_CELL
            flux_value = interval(0._dp); residual_value = flux_value
            call mapped_lagrange(pg, geometry_nodes(1, :, cell), x)
            call mapped_lagrange(pg, geometry_nodes(2, :, cell), y)
            call mapped_lagrange(pg, weight_nodes((cell - 1)*ng + 1:cell*ng), w)
            call mapped_derivative(x, j11, j12)
            call mapped_derivative(y, j21, j22)
            det = mapped_difference(mapped_multiply(j11, j22), &
                                    mapped_multiply(j12, j21))
            dh = mapped_hull(det); wh = mapped_hull(w)
            if (.not. valid(dh) .or. .not. valid(wh)) return
            if (dh%lo <= 0 .or. wh%lo <= 0) return
            call fitted_reference(cell, c)
            sx%degree = fit%degree; sy%degree = fit%degree
            sx%c = interval(0._dp); sy%c = interval(0._dp)
            do k = 1, fit%nb
                a = fit%index(:, k); zeros = count(a == 0)
                if (zeros == 2) then
                    vertex = maxloc(a, dim=1)
                    e1 = vertex; e2 = mod(vertex + 1, 3) + 1
                    slot1 = merge(fit%degree + 1, 1, &
                              (vertex == e1 .and. topo%cell_orientation(e1, cell) < 0) &
                                  .or. (vertex /= e1 .and. &
                                        topo%cell_orientation(e1, cell) > 0))
                    slot2 = merge(fit%degree + 1, 1, &
                              (vertex == e2 .and. topo%cell_orientation(e2, cell) < 0) &
                                  .or. (vertex /= e2 .and. &
                                        topo%cell_orientation(e2, cell) > 0))
                    normal1 = interval(normals(:, e1))
                    normal2 = interval(normals(:, e2))
                    value1 = interval(topo%cell_orientation(e1, cell) &
                                      *shared(slot1, topo%cell_edge(e1, cell)))
                    value2 = interval(topo%cell_orientation(e2, cell) &
                                      *shared(slot2, topo%cell_edge(e2, cell)))
                    call verified_bernstein_corner(value1, normal1(1), normal1(2), &
                                       value2, normal2(1), normal2(2), sx%c(k), sy%c(k))
                else if (zeros == 1) then
                    local = mod(minloc(a, dim=1), 3) + 1
                    next = mod(local, 3) + 1
                    slot1 = a(next) + 1
                    if (topo%cell_orientation(local, cell) < 0) &
                        slot1 = a(local) + 1
                    value1 = interval(topo%cell_orientation(local, cell) &
                                      *shared(slot1, topo%cell_edge(local, cell)))
                    normal1 = interval(normals(:, local))
                    call verified_bernstein_edge_projection(interval(c(1, k)), &
                                    interval(c(2, k)), value1, normal1(1), normal1(2), &
                                                            sx%c(k), sy%c(k))
                else
                    sx%c(k) = interval(c(1, k)); sy%c(k) = interval(c(2, k))
                end if
            end do
            do k = 1, count_of(pu)
                values(k) = solution(solution_cells(k, cell))
            end do
            call mapped_lagrange(pu, values, u)
            call mapped_derivative(u, ux, uy)
            ! Q = w J sigma_hat - adj(J)^T grad_hat u.
            qx = mapped_difference(mapped_multiply(w, &
                      mapped_sum(mapped_multiply(j11, sx), mapped_multiply(j12, sy))), &
                                   mapped_difference(mapped_multiply(j22, ux), &
                                                     mapped_multiply(j21, uy)))
            qy = mapped_difference(mapped_multiply(w, &
                      mapped_sum(mapped_multiply(j21, sx), mapped_multiply(j22, sy))), &
                                   mapped_difference(mapped_multiply(j11, uy), &
                                                     mapped_multiply(j12, ux)))
            call mapped_square_integral(qx%degree, qx%c, integral_x)
            call mapped_square_integral(qy%degree, qy%c, integral_y)
            call verified_mapped_divide(integral_x + integral_y, &
                                        wh*dh, flux_value)
            ! T = s det(J) + w^m div_hat sigma_hat.
            call mapped_derivative(sx, dx, dummy)
            call mapped_derivative(sy, dummy, dy)
            div = mapped_sum(dx, dy)
            s%degree = source_degree; s%c = interval(0._dp)
            call source%evaluate(cell, s%c(:count_of(source_degree)), ok)
            cell_status = VERIFIED_BERNSTEIN_BAD_SOURCE
            if (.not. ok) return
            if (.not. all(valid(s%c(:count_of(source_degree))))) return
            wm%degree = 0; wm%c = interval(0._dp); wm%c(1) = interval(1._dp)
            do k = 1, weight_power
                wm = mapped_multiply(wm, w)
            end do
            r = mapped_sum(mapped_multiply(s, det), mapped_multiply(wm, div))
            call mapped_square_integral(r%degree, r%c, integral_x)
            wmh = mapped_hull(wm)
            call verified_mapped_divide(integral_x, &
                                        wmh*wmh*dh, residual_value)
            if (present(tail)) then
                call tail%evaluate(cell, rho, ok)
                if (.not. ok) return
                if (.not. valid(rho)) return
                if (rho%lo < 0) return
                if (.not. valid(residual_value)) return
                residual_value%lo = max(0._dp, residual_value%lo)
                call verified_mapped_tail(residual_value, rho, &
                                          interval(dh%hi), interval(wmh%lo), enclosed)
                if (.not. valid(enclosed)) return
                residual_value = interval(0._dp, enclosed%hi)
            end if
            cell_status = VERIFIED_BERNSTEIN_BAD_CELL
            if (.not. valid(flux_value) .or. .not. valid(residual_value)) return
            if (flux_value%hi < 0 .or. residual_value%hi < 0) return
            flux_value%lo = max(0._dp, flux_value%lo)
            residual_value%lo = max(0._dp, residual_value%lo)
            cell_status = VERIFIED_BERNSTEIN_OK
        end subroutine
    end subroutine verified_mapped_bernstein_majorant_object

    pure integer function lagrange_degree(n) result(p)
        integer, intent(in) :: n
        p = 0
        if (n == 3) p = 1
        if (n == 6) p = 2
        if (n == 10) p = 3
    end function

    pure integer function geometry_edge_slot(p, edge, j) result(k)
        integer, intent(in) :: p, edge, j
        k = 3 + (edge - 1)*(p - 1) + j
        if (j == 0) k = edge
        if (j == p) k = mod(edge, 3) + 1
    end function

    subroutine mapped_lagrange(p, values, b)
        integer, intent(in) :: p
        real(dp), intent(in) :: values(:)
        type(mapped_poly_t), intent(out) :: b
        type(interval_t) :: v(10)
        b%degree = p; b%c = interval(0._dp)
        v(:count_of(p)) = interval(values(:count_of(p)))
        select case (p)
        case (1)
            b%c(:3) = v(:3)
        case (2)
            call verified_bernstein_lagrange_p2(v(1), v(2), v(3), v(4), v(5), &
                                   v(6), b%c(1), b%c(2), b%c(3), b%c(4), b%c(5), b%c(6))
        case (3)
            call verified_bernstein_lagrange_p3(v(1), v(2), v(3), v(4), v(5), &
                                v(6), v(7), v(8), v(9), v(10), b%c(1), b%c(2), b%c(3), &
                                b%c(4), b%c(5), b%c(6), b%c(7), b%c(8), b%c(9), b%c(10))
        end select
    end subroutine

    subroutine mapped_derivative(a, dx, dy)
        type(mapped_poly_t), intent(in) :: a
        type(mapped_poly_t), intent(out) :: dx, dy
        type(interval_t) :: g(4)
        dx%degree = a%degree - 1; dy%degree = a%degree - 1
        dx%c = interval(0._dp); dy%c = interval(0._dp)
        g = interval([1._dp, 0._dp, 0._dp, 1._dp])
        call verified_bernstein_gradient(a%degree, a%c, g, &
                                 dx%c(:count_of(dx%degree)), dy%c(:count_of(dy%degree)))
    end subroutine

    function mapped_multiply(a, b) result(c)
        type(mapped_poly_t), intent(in) :: a, b
        type(mapped_poly_t) :: c
        c%degree = a%degree + b%degree
        call mapped_product(a%degree, b%degree, a%c, b%c, c%c)
    end function

    function mapped_sum(a, b) result(c)
        type(mapped_poly_t), intent(in) :: a, b
        type(mapped_poly_t) :: c, tmp, elevated
        c = a; tmp = b
        do while (c%degree < tmp%degree)
            call mapped_elevate(c%degree, c%c, elevated%c)
            c%c = elevated%c; c%degree = c%degree + 1
        end do
        do while (tmp%degree < c%degree)
            call mapped_elevate(tmp%degree, tmp%c, elevated%c)
            tmp%c = elevated%c; tmp%degree = tmp%degree + 1
        end do
        c%c = c%c + tmp%c
    end function

    function mapped_difference(a, b) result(c)
        type(mapped_poly_t), intent(in) :: a, b
        type(mapped_poly_t) :: c, negative
        negative = b; negative%c = -b%c
        c = mapped_sum(a, negative)
    end function

    function mapped_hull(a) result(h)
        type(mapped_poly_t), intent(in) :: a
        type(interval_t) :: h
        integer :: n
        n = count_of(a%degree)
        h = interval(minval(a%c(:n)%lo), maxval(a%c(:n)%hi))
    end function


    subroutine verified_bernstein_majorant_object(vertices, triangles, degree, &
            flux_coefficients, solution_cells, solution, weight, source, &
            source_degree, weight_power, friedrichs, bound, flux_squared, &
            residual_squared, status, source_remainder)
        !! degree: k of the RT_k flux (equilibrated_flux_estimate_t,
        !! flux_coefficients(dofs, cells)). solution_cells(nloc, nt) indexes the
        !! Lagrange values of u_h in `solution`: nloc = 3 (P1), 6 (P2), or 10
        !! (P3): vertices, then edge nodes along 1-2, 2-3, 3-1, then centroid.
        !! P3 has two edge nodes at 1/3 and 2/3 along each directed edge.
        real(dp), intent(in) :: vertices(:, :), flux_coefficients(:, :)
        real(dp), intent(in) :: solution(:), weight(:)
        integer, intent(in) :: triangles(:, :), solution_cells(:, :)
        integer, intent(in) :: degree, source_degree, weight_power
        class(bernstein_source_t), intent(in) :: source
        class(bernstein_remainder_t), intent(in), optional :: source_remainder
        type(interval_t), intent(in) :: friedrichs
        type(interval_t), intent(out) :: bound, flux_squared, residual_squared
        integer, intent(out) :: status
        type(topology_t) :: topo
        type(fit_t) :: fit
        type(interval_t) :: flux_part, residual_part, cell_flux, cell_residual
        real(dp), allocatable :: shared(:, :)
        integer, allocatable :: owner(:), uses(:), turns(:)
        integer :: nv, nt, t, e, edge, stat, failed, cell_status, thread_failed

        status = VERIFIED_BERNSTEIN_BAD_INPUT
        bound = interval(0.0_dp, huge(1.0_dp))
        flux_squared = bound
        residual_squared = bound
        if (size(vertices, 1) /= 2 .or. size(triangles, 1) /= 3) return
        nv = size(vertices, 2); nt = size(triangles, 2)
        if (nv < 3 .or. nt < 1) return
        if (degree < 1 .or. degree > 3) return
        if (weight_power < 0 .or. weight_power > 2) return
        if (source_degree < 0) return
        if (max(source_degree, degree + weight_power) > 6) return
        if (all(size(solution_cells, 1) /= [3, 6, 10])) return
        if (size(solution_cells, 2) /= nt .or. size(flux_coefficients, 2) /= nt) return
        if (size(weight) /= nv) return
        if (any(triangles < 1) .or. any(triangles > nv)) return
        if (any(solution_cells < 1) .or. any(solution_cells > size(solution))) return
        if (.not. all(ieee_is_finite(vertices))) return
        if (.not. all(ieee_is_finite(weight))) return
        if (any(weight <= 0)) return
        if (.not. all(ieee_is_finite(solution))) return
        if (.not. all(ieee_is_finite(flux_coefficients))) return
        if (.not. valid(friedrichs)) return
        if (friedrichs%lo < 0) return
        call build_fit(degree, fit, stat)
        if (stat /= 0) return
        if (size(flux_coefficients, 1) /= fit%nrt) return
        call build_topology(nv, triangles, topo, stat)
        if (stat /= 0) return
        ! One owner cell per global edge selects the shared normal components.
        allocate (owner(maxval(topo%cell_edge)), uses(maxval(topo%cell_edge)))
        allocate (turns(maxval(topo%cell_edge)))
        owner = 0; uses = 0; turns = 0
        do t = 1, nt
            do e = 1, 3
                edge = topo%cell_edge(e, t)
                uses(edge) = uses(edge) + 1
                turns(edge) = turns(edge) + topo%cell_orientation(e, t)
                if (owner(edge) == 0) owner(edge) = t
            end do
        end do
        if (any(uses > 2)) return
        if (any(uses == 2 .and. turns /= 0)) return
        allocate (shared(fit%degree + 2, size(owner)))
        !$omp parallel do schedule(static)
        do edge = 1, size(owner)
            call owner_values(edge)
        end do
        !$omp end parallel do
        status = VERIFIED_BERNSTEIN_BAD_CELL
        failed = 0
        flux_squared = interval(0.0_dp); residual_squared = interval(0.0_dp)
        ! Interval sums enclose the exact sum in any order, so each thread
        ! accumulates its own partial enclosure.
        !$omp parallel private(flux_part, residual_part, cell_flux, cell_residual, &
        !$omp cell_status, thread_failed)
        flux_part = interval(0.0_dp); residual_part = interval(0.0_dp)
        thread_failed = 0
        !$omp do schedule(static)
        do t = 1, nt
            if (thread_failed /= 0) cycle
            call audit_cell(t, cell_flux, cell_residual, cell_status, source_remainder)
            if (cell_status /= 0) then
                thread_failed = cell_status
                cycle
            end if
            flux_part = flux_part + cell_flux
            residual_part = residual_part + cell_residual
        end do
        !$omp end do
        !$omp critical (verified_bernstein_sum)
        flux_squared = flux_squared + flux_part
        residual_squared = residual_squared + residual_part
        failed = max(failed, thread_failed)
        !$omp end critical (verified_bernstein_sum)
        !$omp end parallel
        if (failed /= 0) then
            status = failed
            flux_squared = interval(0.0_dp, huge(1.0_dp))
            residual_squared = flux_squared
            return
        end if
        flux_squared%lo = max(0.0_dp, flux_squared%lo)
        residual_squared%lo = max(0.0_dp, residual_squared%lo)
        bound = sqrt(flux_squared) + friedrichs*sqrt(residual_squared)
        if (.not. valid(bound)) then
            bound = interval(0.0_dp, huge(1.0_dp))
            return
        end if
        status = VERIFIED_BERNSTEIN_OK

    contains

        subroutine owner_values(edge)
            !! Floating normal components of the owner's fit at the Bernstein
            !! positions j = 0..n counted from the lower global vertex.
            integer, intent(in) :: edge
            real(dp) :: c(2, max_coefficients), normal(2)
            integer :: cell, local, lo, hi, j, a(3), k
            cell = owner(edge)
            do local = 1, 3
                if (topo%cell_edge(local, cell) == edge) exit
            end do
            call fitted_coefficients(cell, c)
            call edge_ends(cell, local, lo, hi)
            lo = triangles(lo, cell); hi = triangles(hi, cell)
            normal = [vertices(2, hi) - vertices(2, lo), &
                      vertices(1, lo) - vertices(1, hi)]
            call edge_ends(cell, local, lo, hi)
            do j = 0, fit%degree
                a = 0; a(hi) = j; a(lo) = fit%degree - j
                k = position(fit, a)
                shared(j + 1, edge) = dot_product(normal, c(:, k))
            end do
        end subroutine owner_values

        subroutine edge_ends(cell, local, lo, hi)
            !! Local vertex slots of the lower and higher global vertex of edge
            !! `local` (from slot local to slot local + 1).
            integer, intent(in) :: cell, local
            integer, intent(out) :: lo, hi
            integer :: next
            next = mod(local, 3) + 1
            if (triangles(local, cell) < triangles(next, cell)) then
                lo = local; hi = next
            else
                lo = next; hi = local
            end if
        end subroutine edge_ends

        subroutine fitted_coefficients(cell, c)
            integer, intent(in) :: cell
            real(dp), intent(out) :: c(:, :)
            real(dp) :: jac(2, 2), det, values(2, max_coefficients), v(2)
            integer :: q, i
            jac(:, 1) = vertices(:, triangles(2, cell)) &
                        - vertices(:, triangles(1, cell))
            jac(:, 2) = vertices(:, triangles(3, cell)) &
                        - vertices(:, triangles(1, cell))
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            c = 0
            if (.not. (det > 0)) return
            do q = 1, fit%nb
                v = 0
                do i = 1, fit%nrt
                    v = v + fit%rt(:, i, q)*flux_coefficients(i, cell)
                end do
                ! Contravariant Piola map of the reference RT_k field.
                values(1, q) = (jac(1, 1)*v(1) + jac(1, 2)*v(2))/det
                values(2, q) = (jac(2, 1)*v(1) + jac(2, 2)*v(2))/det
            end do
            do i = 1, fit%nb
                v = 0
                do q = 1, fit%nb
                    v = v + fit%inverse(i, q)*values(:, q)
                end do
                c(:, i) = v
            end do
        end subroutine fitted_coefficients

        subroutine audit_cell(cell, flux_value, residual_value, cell_status, tail)
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: flux_value, residual_value
            integer, intent(out) :: cell_status
            class(bernstein_remainder_t), intent(in), optional :: tail
            real(dp) :: c(2, max_coefficients)
            type(interval_t) :: xy(2, 3), det, g(4), normal(2, 3), w(3), chord(3)
            type(interval_t) :: sx(max_coefficients), sy(max_coefficients)
            type(interval_t) :: qx(max_coefficients), qy(max_coefficients)
            type(interval_t) :: ux(max_coefficients), uy(max_coefficients)
            type(interval_t) :: u(10), b(10), tmp(max_coefficients)
            type(interval_t) :: div(max_coefficients), dyy(max_coefficients)
            type(interval_t) :: s(max_coefficients), part, rho, enclosed
            real(dp) :: wmin, wmax
            integer :: j, k, n, nq, local, lo, hi, v, e1, e2, j1, j2, a(3)
            integer :: zero, nloc, pu, ns, nt_degree
            logical :: ok

            cell_status = VERIFIED_BERNSTEIN_BAD_CELL
            flux_value = interval(0.0_dp); residual_value = interval(0.0_dp)
            n = fit%degree
            do j = 1, 3
                xy(:, j) = interval(vertices(:, triangles(j, cell)))
                w(j) = interval(weight(triangles(j, cell)))
            end do
            call verified_bernstein_geometry(xy(1, 1), xy(2, 1), xy(1, 2), &
                xy(2, 2), xy(1, 3), xy(2, 3), det, g(1), g(2), g(3), g(4))
            if (.not. valid(det)) return
            if (det%lo <= 0) return
            do local = 1, 3
                call edge_ends(cell, local, lo, hi)
                call verified_bernstein_edge_normal(xy(1, lo), xy(2, lo), &
                    xy(1, hi), xy(2, hi), normal(1, local), normal(2, local))
            end do
            ! Exact H(div) field from the fitted coefficients and shared values.
            call fitted_coefficients(cell, c)
            do k = 1, fit%nb
                a = fit%index(:, k)
                zero = count(a == 0)
                if (zero == 2) then
                    v = maxloc(a, dim=1)
                    e1 = v; e2 = mod(v + 1, 3) + 1
                    j1 = corner_position(cell, e1, v)
                    j2 = corner_position(cell, e2, v)
                    call verified_bernstein_corner( &
                        interval(shared(j1, topo%cell_edge(e1, cell))), &
                        normal(1, e1), normal(2, e1), &
                        interval(shared(j2, topo%cell_edge(e2, cell))), &
                        normal(1, e2), normal(2, e2), sx(k), sy(k))
                else if (zero == 1) then
                    ! Edge local is opposite the vanishing barycentric index.
                    local = mod(minloc(a, dim=1), 3) + 1
                    call edge_ends(cell, local, lo, hi)
                    call verified_bernstein_edge_projection(interval(c(1, k)), &
                        interval(c(2, k)), &
                        interval(shared(a(hi) + 1, topo%cell_edge(local, cell))), &
                        normal(1, local), normal(2, local), sx(k), sy(k))
                else
                    sx(k) = interval(c(1, k)); sy(k) = interval(c(2, k))
                end if
            end do
            ! Discrete solution in Bernstein form and its physical gradient.
            nloc = size(solution_cells, 1)
            do j = 1, nloc
                u(j) = interval(solution(solution_cells(j, cell)))
            end do
            if (nloc == 10) then
                pu = 3
                call verified_bernstein_lagrange_p3(u(1), u(2), u(3), u(4), &
                    u(5), u(6), u(7), u(8), u(9), u(10), b(1), b(2), b(3), &
                    b(4), b(5), b(6), b(7), b(8), b(9), b(10))
            else if (nloc == 6) then
                pu = 2
                call verified_bernstein_lagrange_p2(u(1), u(2), u(3), u(4), &
                    u(5), u(6), b(1), b(2), b(3), b(4), b(5), b(6))
            else
                pu = 1
                b(1:3) = u(1:3)
            end if
            call verified_bernstein_gradient(pu, b, g, ux, uy)
            do j = pu - 1, n
                call verified_bernstein_elevate(j, ux, tmp)
                ux = tmp
                call verified_bernstein_elevate(j, uy, tmp)
                uy = tmp
            end do
            ! q = w sigma - grad u_h, degree n + 1.
            call verified_bernstein_affine_product(n, sx, w, qx)
            call verified_bernstein_affine_product(n, sy, w, qy)
            nq = count_of(n + 1)
            qx(1:nq) = qx(1:nq) - ux(1:nq)
            qy(1:nq) = qy(1:nq) - uy(1:nq)
            wmin = minval(weight(triangles(:, cell)))
            wmax = maxval(weight(triangles(:, cell)))
            call chord_values(1, wmin, wmax, w, chord)
            call verified_bernstein_square_integral(n + 1, qx, chord, part)
            flux_value = part
            call verified_bernstein_square_integral(n + 1, qy, chord, part)
            flux_value = (flux_value + part)*det
            ! t = s + w^m div sigma at the common degree.
            call verified_bernstein_gradient(n, sx, g, div, tmp)
            call verified_bernstein_gradient(n, sy, g, tmp, dyy)
            k = count_of(n - 1)
            div(1:k) = div(1:k) + dyy(1:k)
            nt_degree = n - 1
            do j = 1, weight_power
                call verified_bernstein_affine_product(nt_degree, div, w, tmp)
                nt_degree = nt_degree + 1
                div(1:count_of(nt_degree)) = tmp(1:count_of(nt_degree))
            end do
            ns = source_degree
            call source%evaluate(cell, s(1:count_of(ns)), ok)
            if (.not. ok) then
                cell_status = VERIFIED_BERNSTEIN_BAD_SOURCE
                return
            end if
            if (.not. all(valid(s(1:count_of(ns))))) then
                cell_status = VERIFIED_BERNSTEIN_BAD_SOURCE
                return
            end if
            do while (ns < nt_degree)
                call verified_bernstein_elevate(ns, s, tmp)
                ns = ns + 1
                s(1:count_of(ns)) = tmp(1:count_of(ns))
            end do
            do while (nt_degree < ns)
                call verified_bernstein_elevate(nt_degree, div, tmp)
                nt_degree = nt_degree + 1
                div(1:count_of(nt_degree)) = tmp(1:count_of(nt_degree))
            end do
            k = count_of(ns)
            s(1:k) = s(1:k) + div(1:k)
            if (weight_power > 0) then
                call chord_values(2*weight_power, wmin, wmax, w, chord)
            else
                chord = interval(1.0_dp)
            end if
            call verified_bernstein_square_integral(ns, s, chord, part)
            residual_value = part*det
            if (present(tail)) then
                call tail%evaluate(cell, rho, ok)
                cell_status = VERIFIED_BERNSTEIN_BAD_SOURCE
                if (.not. ok) return
                if (.not. valid(rho)) return
                if (rho%lo < 0) return
                ! Minkowski and |delta/w^m| <= rho/wmin^m. The squared
                ! residual is enclosed by [0, upper], since delta may cancel
                ! the polynomial residual. No claim of a positive lower bound.
                if (.not. valid(residual_value)) return
                residual_value%lo = max(0.0_dp, residual_value%lo)
                select case (weight_power)
                case (0)
                    call verified_bernstein_source_tail_0(residual_value, rho, &
                        det, enclosed)
                case (1)
                    call verified_bernstein_source_tail_1(residual_value, rho, &
                        det, interval(wmin), enclosed)
                case (2)
                    call verified_bernstein_source_tail_2(residual_value, rho, &
                        det, interval(wmin), enclosed)
                end select
                if (.not. valid(enclosed)) return
                if (enclosed%hi < 0) return
                residual_value = interval(0.0_dp, enclosed%hi)
            end if
            if (.not. (valid(flux_value) .and. valid(residual_value))) return
            flux_value%lo = max(0.0_dp, flux_value%lo)
            residual_value%lo = max(0.0_dp, residual_value%lo)
            cell_status = 0
        end subroutine audit_cell

        integer function corner_position(cell, local, v) result(j)
            !! Shared-value slot (1-based) of local vertex v on edge `local`.
            integer, intent(in) :: cell, local, v
            integer :: lo, hi
            call edge_ends(cell, local, lo, hi)
            j = merge(fit%degree + 1, 1, v == hi)
        end function corner_position

    end subroutine verified_bernstein_majorant_object

    subroutine chord_values(power, wmin, wmax, w, chord)
        !! Vertex values of an affine upper bound of w^-power on the cell.
        integer, intent(in) :: power
        real(dp), intent(in) :: wmin, wmax
        type(interval_t), intent(in) :: w(3)
        type(interval_t), intent(out) :: chord(3)
        integer :: j
        do j = 1, 3
            if (wmax > wmin) then
                call verified_bernstein_chord(power, interval(wmin), interval(wmax), &
                                              w(j), chord(j))
            else
                call verified_bernstein_reciprocal(power, w(j), chord(j))
            end if
        end do
    end subroutine chord_values

    subroutine build_fit(degree, fit, status)
        integer, intent(in) :: degree
        type(fit_t), intent(out) :: fit
        integer, intent(out) :: status
        type(triangle_rt_basis_t) :: basis
        real(dp), allocatable :: vandermonde(:, :), divergences(:)
        integer, allocatable :: pivots(:)
        real(dp) :: xi, eta, l(3)
        integer :: q, k, n, info

        call initialize_triangle_raviart_thomas(degree, basis, status)
        if (status /= 0) return
        status = 1
        n = degree + 1
        fit%degree = n
        fit%nrt = triangle_rt_dof_count(basis)
        fit%nb = count_of(n)
        fit%index = indices(n)
        allocate (fit%rt(2, fit%nrt, fit%nb), divergences(fit%nrt))
        allocate (vandermonde(fit%nb, fit%nb), pivots(fit%nb))
        ! Lattice point of index a is the barycentric point a/n.
        do q = 1, fit%nb
            xi = real(fit%index(2, q), dp)/n
            eta = real(fit%index(3, q), dp)/n
            call evaluate_triangle_raviart_thomas(basis, xi, eta, fit%rt(:, :, q), &
                                                  divergences, status)
            if (status /= 0) return
            l = [1.0_dp - xi - eta, xi, eta]
            do k = 1, fit%nb
                vandermonde(q, k) = bernstein_value(n, fit%index(:, k), l)
            end do
        end do
        allocate (fit%inverse(fit%nb, fit%nb))
        fit%inverse = 0
        do k = 1, fit%nb
            fit%inverse(k, k) = 1
        end do
        call dgesv(fit%nb, fit%nb, vandermonde, fit%nb, pivots, fit%inverse, &
                   fit%nb, info)
        status = merge(0, 1, info == 0)
    end subroutine build_fit

    pure integer function count_of(degree)
        integer, intent(in) :: degree
        count_of = (degree + 1)*(degree + 2)/2
    end function count_of

    pure function indices(degree) result(a)
        !! Bernstein multi-indices in the generator's order: descending a1,
        !! then descending a2.
        integer, intent(in) :: degree
        integer :: a(3, (degree + 1)*(degree + 2)/2)
        integer :: i, j, k
        k = 0
        do i = degree, 0, -1
            do j = degree - i, 0, -1
                k = k + 1
                a(:, k) = [i, j, degree - i - j]
            end do
        end do
    end function indices

    pure integer function position(fit, a)
        type(fit_t), intent(in) :: fit
        integer, intent(in) :: a(3)
        integer :: k
        position = 0
        do k = 1, fit%nb
            if (all(fit%index(:, k) == a)) position = k
        end do
    end function position

    pure real(dp) function bernstein_value(n, a, l) result(value)
        integer, intent(in) :: n, a(3)
        real(dp), intent(in) :: l(3)
        value = gamma(real(n + 1, dp))/(gamma(real(a(1) + 1, dp)) &
                *gamma(real(a(2) + 1, dp))*gamma(real(a(3) + 1, dp))) &
                *l(1)**a(1)*l(2)**a(2)*l(3)**a(3)
    end function bernstein_value

    elemental logical function valid(value)
        type(interval_t), intent(in) :: value
        ! Test finiteness before the ordered comparison: Fortran .and. need
        ! not short-circuit, and an ordered comparison with NaN signals invalid.
        valid = .false.
        if (.not. ieee_is_finite(value%lo)) return
        if (.not. ieee_is_finite(value%hi)) return
        valid = value%lo <= value%hi
    end function valid

    subroutine verified_mapped_bernstein_majorant_procedure(vertices, triangles, &
        degree, &
                          flux_coefficients, solution_cells, solution, geometry_nodes, &
                        weight_nodes, source, source_degree, weight_power, friedrichs, &
                        bound, flux_squared, residual_squared, status, source_remainder)
        !! Outward audit on conforming isoparametric P1/P2/P3 triangles.
        !! weight_nodes is cell-major: n_geometry_nodes values for each cell.
        !! It defines w=1/A by the same reference Lagrange interpolation as
        !! geometry. Source supplies a reference Bernstein numerator s in
        !! f=(s+delta)/w**weight_power, |delta|<=source_remainder on each cell.
        !! No quadrature is sampled. Exact polynomial square moments divided
        !! by positive denominator hulls give lower/upper integral budgets.
        !! The caller supplies nonoverlap, a conforming H1_0 candidate and C_F.
        !! Shared mapped edges are checked node for node, including orientation.
        !! A nonpositive Bernstein hull fails closed (a sufficient test only).
        real(dp), intent(in) :: vertices(:, :), flux_coefficients(:, :)
        real(dp), intent(in) :: solution(:), geometry_nodes(:, :, :)
        real(dp), intent(in) :: weight_nodes(:)
        integer, intent(in) :: triangles(:, :), solution_cells(:, :)
        integer, intent(in) :: degree, source_degree, weight_power
        procedure(verified_bernstein_source) :: source
        procedure(verified_bernstein_source_remainder), optional :: source_remainder
        type(interval_t), intent(in) :: friedrichs
        type(interval_t), intent(out) :: bound, flux_squared, residual_squared
        integer, intent(out) :: status
        type(bernstein_source_adapter_t), target :: source_adapter
        type(bernstein_remainder_adapter_t), target :: source_remainder_adapter
        class(bernstein_remainder_t), pointer :: source_remainder_object
        nullify(source_remainder_object)

        source_adapter%callback => source
        if (present(source_remainder)) then
            source_remainder_adapter%callback => source_remainder
            source_remainder_object => source_remainder_adapter
        end if
        call verified_mapped_bernstein_majorant_object(vertices, &
            triangles, &
            degree, &
            flux_coefficients, &
            solution_cells, &
            solution, &
            geometry_nodes, &
            weight_nodes, &
            source_adapter, &
            source_degree, &
            weight_power, &
            friedrichs, &
            bound, &
            flux_squared, &
            residual_squared, &
            status, &
            source_remainder_object)
    end subroutine

    subroutine verified_bernstein_majorant_procedure(vertices, triangles, degree, &
            flux_coefficients, solution_cells, solution, weight, source, &
            source_degree, weight_power, friedrichs, bound, flux_squared, &
            residual_squared, status, source_remainder)
        !! degree: k of the RT_k flux (equilibrated_flux_estimate_t,
        !! flux_coefficients(dofs, cells)). solution_cells(nloc, nt) indexes the
        !! Lagrange values of u_h in `solution`: nloc = 3 (P1), 6 (P2), or 10
        !! (P3): vertices, then edge nodes along 1-2, 2-3, 3-1, then centroid.
        !! P3 has two edge nodes at 1/3 and 2/3 along each directed edge.
        real(dp), intent(in) :: vertices(:, :), flux_coefficients(:, :)
        real(dp), intent(in) :: solution(:), weight(:)
        integer, intent(in) :: triangles(:, :), solution_cells(:, :)
        integer, intent(in) :: degree, source_degree, weight_power
        procedure(verified_bernstein_source) :: source
        procedure(verified_bernstein_source_remainder), optional :: source_remainder
        type(interval_t), intent(in) :: friedrichs
        type(interval_t), intent(out) :: bound, flux_squared, residual_squared
        integer, intent(out) :: status
        type(bernstein_source_adapter_t), target :: source_adapter
        type(bernstein_remainder_adapter_t), target :: source_remainder_adapter
        class(bernstein_remainder_t), pointer :: source_remainder_object
        nullify(source_remainder_object)

        source_adapter%callback => source
        if (present(source_remainder)) then
            source_remainder_adapter%callback => source_remainder
            source_remainder_object => source_remainder_adapter
        end if
        call verified_bernstein_majorant_object(vertices, &
            triangles, &
            degree, &
            flux_coefficients, &
            solution_cells, &
            solution, &
            weight, &
            source_adapter, &
            source_degree, &
            weight_power, &
            friedrichs, &
            bound, &
            flux_squared, &
            residual_squared, &
            status, &
            source_remainder_object)
    end subroutine

    subroutine bernstein_source_adapter_t_evaluate(self, cell, coefficients, ok)
        class(bernstein_source_adapter_t), intent(in) :: self
            !! Outward Bernstein coefficients of the polynomial numerator (or
            !! its approximation when source_remainder is supplied), using
            !! the vertex order of `triangles`. Called
            !! concurrently from OpenMP threads.
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: coefficients(:)
            logical, intent(out) :: ok
        call self%callback(cell, coefficients, ok)
    end subroutine

    subroutine bernstein_remainder_adapter_t_evaluate(self, cell, radius, ok)
        class(bernstein_remainder_adapter_t), intent(in) :: self
            !! Nonnegative outward upper bound radius%hi for the absolute
            !! numerator remainder on the WHOLE cell, not sampled points.
            !! The source polynomial is supplied by verified_bernstein_source.
            !! The caller proves this contract. Called concurrently.
            integer, intent(in) :: cell
            type(interval_t), intent(out) :: radius
            logical, intent(out) :: ok
        call self%callback(cell, radius, ok)
    end subroutine

end module fortfem_verified_bernstein_majorant
