module fortfem_tensor_points
    ! Point location and evaluation on lab tensor meshes (2D/3D, curved
    ! isoparametric cells), plus a local extremum search for a nodal field.
    ! Location inverts each candidate cell's geometry map by Newton's method
    ! inside an inflated bounding box; evaluation returns value and physical
    ! gradient. Intended for sampling, diagnostics and nonlocal functionals
    ! (for example a free magnetic axis); not for quadrature hot loops.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t
    use fortfem_gll_line, only: gll_rule, lagrange_table
    implicit none
    private
    public :: tensor_locator_t, build_tensor_locator, locate_tensor_point
    public :: tensor_point_basis, evaluate_tensor_point, tensor_field_extremum

    type :: tensor_locator_t
        integer :: dim = 0, p = 0, q = 0
        real(dp), allocatable :: field_nodes(:), geometry_nodes(:)
        real(dp), allocatable :: lower(:, :), upper(:, :), size(:)
    end type

contains

    subroutine build_tensor_locator(mesh, locator, status)
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_locator_t), intent(out) :: locator
        integer, intent(out) :: status
        real(dp), allocatable :: w(:)
        real(dp) :: margin
        integer :: e, d
        locator%dim = mesh%dim
        locator%p = mesh%degree
        locator%q = mesh%geometry_degree
        call gll_rule(mesh%degree, locator%field_nodes, w, status)
        if (status /= 0) return
        call gll_rule(mesh%geometry_degree, locator%geometry_nodes, w, status)
        if (status /= 0) return
        allocate (locator%lower(mesh%dim, mesh%nelem), &
                  locator%upper(mesh%dim, mesh%nelem), locator%size(mesh%nelem))
        do e = 1, mesh%nelem
            do d = 1, mesh%dim
                locator%lower(d, e) = minval(mesh%geometry(d, :, e))
                locator%upper(d, e) = maxval(mesh%geometry(d, :, e))
            end do
            locator%size(e) = maxval(locator%upper(:, e) - locator%lower(:, e))
            ! Curved edges may bulge between geometry nodes.
            margin = 0.25_dp*locator%size(e)
            locator%lower(:, e) = locator%lower(:, e) - margin
            locator%upper(:, e) = locator%upper(:, e) + margin
        end do
    end subroutine

    subroutine line_basis(nodes, x, b, d)
        real(dp), intent(in) :: nodes(:), x
        real(dp), intent(out) :: b(:), d(:)
        real(dp), allocatable :: bt(:, :), dt(:, :)
        call lagrange_table(nodes, [x], bt, dt)
        b = bt(1, :)
        d = dt(1, :)
    end subroutine

    subroutine tensor_basis(nodes, dim, xi, values, gradients)
        ! Tensor Lagrange values and reference gradients at xi.
        real(dp), intent(in) :: nodes(:), xi(:)
        integer, intent(in) :: dim
        real(dp), intent(out) :: values(:), gradients(:, :)
        real(dp) :: b(size(nodes), 3), d(size(nodes), 3)
        integer :: n, i, j, k, a, l
        n = size(nodes)
        b = 0; d = 0
        b(1, 3) = 1
        do l = 1, dim
            call line_basis(nodes, xi(l), b(:, l), d(:, l))
        end do
        a = 0
        do k = 1, merge(n, 1, dim == 3)
            do j = 1, n
                do i = 1, n
                    a = a + 1
                    values(a) = b(i, 1)*b(j, 2)*b(k, 3)
                    gradients(1, a) = d(i, 1)*b(j, 2)*b(k, 3)
                    gradients(2, a) = b(i, 1)*d(j, 2)*b(k, 3)
                    if (dim == 3) gradients(3, a) = b(i, 1)*b(j, 2)*d(k, 3)
                end do
            end do
        end do
    end subroutine

    subroutine geometry_at(mesh, locator, e, xi, x, jac)
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_locator_t), intent(in) :: locator
        integer, intent(in) :: e
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:), jac(:, :)
        real(dp), allocatable :: values(:), gradients(:, :)
        integer :: g, d
        allocate (values(size(mesh%geometry, 2)), &
                  gradients(mesh%dim, size(mesh%geometry, 2)))
        call tensor_basis(locator%geometry_nodes, mesh%dim, xi, values, gradients)
        x = 0
        jac = 0
        do g = 1, size(values)
            x = x + values(g)*mesh%geometry(:, g, e)
            do d = 1, mesh%dim
                jac(:, d) = jac(:, d) + gradients(d, g)*mesh%geometry(:, g, e)
            end do
        end do
    end subroutine

    pure subroutine small_solve(a, b, x, status)
        ! Gaussian elimination with partial pivoting for dim <= 3 systems.
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: status
        real(dp) :: m(size(b), size(b) + 1), row(size(b) + 1), factor
        integer :: n, i, k, pivot
        n = size(b)
        m(:, 1:n) = a
        m(:, n + 1) = b
        status = 1
        do k = 1, n
            pivot = k - 1 + maxloc(abs(m(k:n, k)), 1)
            if (m(pivot, k) == 0) return
            if (pivot /= k) then
                row = m(k, :); m(k, :) = m(pivot, :); m(pivot, :) = row
            end if
            do i = k + 1, n
                factor = m(i, k)/m(k, k)
                m(i, k:) = m(i, k:) - factor*m(k, k:)
            end do
        end do
        do i = n, 1, -1
            x(i) = (m(i, n + 1) - dot_product(m(i, i + 1:n), x(i + 1:n)))/m(i, i)
        end do
        status = 0
    end subroutine

    subroutine invert_cell(mesh, locator, e, x, xi, inside)
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_locator_t), intent(in) :: locator
        integer, intent(in) :: e
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: xi(:)
        logical, intent(out) :: inside
        real(dp) :: y(mesh%dim), jac(mesh%dim, mesh%dim), step(mesh%dim)
        integer :: it, status
        inside = .false.
        xi = 0
        do it = 1, 60
            call geometry_at(mesh, locator, e, xi, y, jac)
            call small_solve(jac, x - y, step, status)
            if (status /= 0) return
            xi = max(-2.0_dp, min(2.0_dp, xi + step))
            if (maxval(abs(step)) < 1.e-14_dp) exit
        end do
        if (maxval(abs(step)) > 1.e-10_dp) return
        inside = maxval(abs(xi)) <= 1 + 1.e-10_dp
        xi = max(-1.0_dp, min(1.0_dp, xi))
    end subroutine

    subroutine locate_tensor_point(mesh, locator, x, cell, xi, status)
        ! On entry a positive cell is tried first. status 0: found.
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_locator_t), intent(in) :: locator
        real(dp), intent(in) :: x(:)
        integer, intent(inout) :: cell
        real(dp), intent(out) :: xi(:)
        integer, intent(out) :: status
        logical :: inside
        integer :: e
        status = 1
        if (cell >= 1 .and. cell <= mesh%nelem) then
            call invert_cell(mesh, locator, cell, x, xi, inside)
            if (inside) then
                status = 0
                return
            end if
        end if
        do e = 1, mesh%nelem
            if (any(x(1:mesh%dim) < locator%lower(:, e))) cycle
            if (any(x(1:mesh%dim) > locator%upper(:, e))) cycle
            call invert_cell(mesh, locator, e, x, xi, inside)
            if (inside) then
                cell = e
                status = 0
                return
            end if
        end do
    end subroutine

    subroutine tensor_point_basis(mesh, locator, xi, values)
        ! Field basis values at reference point xi, ordered as mesh%dof(:, e).
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_locator_t), intent(in) :: locator
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: values(:)
        real(dp), allocatable :: gradients(:, :)
        allocate (gradients(mesh%dim, size(values)))
        call tensor_basis(locator%field_nodes, mesh%dim, xi, values, gradients)
    end subroutine

    subroutine evaluate_tensor_point(mesh, locator, u, cell, xi, value, gradient, &
                                     status)
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_locator_t), intent(in) :: locator
        real(dp), intent(in) :: u(:), xi(:)
        integer, intent(in) :: cell
        real(dp), intent(out) :: value, gradient(:)
        integer, intent(out) :: status
        real(dp), allocatable :: values(:), gradients(:, :)
        real(dp) :: x(mesh%dim), jac(mesh%dim, mesh%dim), reference(mesh%dim)
        integer :: a
        allocate (values(size(mesh%dof, 1)), gradients(mesh%dim, size(mesh%dof, 1)))
        call tensor_basis(locator%field_nodes, mesh%dim, xi, values, gradients)
        value = 0
        reference = 0
        do a = 1, size(values)
            value = value + values(a)*u(mesh%dof(a, cell))
            reference = reference + gradients(:, a)*u(mesh%dof(a, cell))
        end do
        call geometry_at(mesh, locator, cell, xi, x, jac)
        call small_solve(transpose(jac), reference, gradient(1:mesh%dim), status)
    end subroutine

    subroutine reference_gradient(mesh, locator, u, cell, xi, value, gradient)
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_locator_t), intent(in) :: locator
        real(dp), intent(in) :: u(:), xi(:)
        integer, intent(in) :: cell
        real(dp), intent(out) :: value, gradient(:)
        real(dp), allocatable :: values(:), gradients(:, :)
        integer :: a
        allocate (values(size(mesh%dof, 1)), gradients(mesh%dim, size(mesh%dof, 1)))
        call tensor_basis(locator%field_nodes, mesh%dim, xi, values, gradients)
        value = 0
        gradient = 0
        do a = 1, size(values)
            value = value + values(a)*u(mesh%dof(a, cell))
            gradient = gradient + gradients(:, a)*u(mesh%dof(a, cell))
        end do
    end subroutine

    subroutine tensor_field_extremum(mesh, locator, u, cell, xi, x, value, status, &
                                     maximum)
        ! Local minimum (default) or maximum of the nodal field u. On entry a
        ! positive cell with its xi is the start; otherwise the extreme node.
        ! Stationarity is coordinate invariant, so Newton runs in the cell's
        ! reference coordinates (exact polynomial gradient, central-difference
        ! Hessian) and relocates only when a step leaves the cell. Steps are
        ! damped so that u never moves away from the extremum.
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_locator_t), intent(in) :: locator
        real(dp), intent(in) :: u(:)
        integer, intent(inout) :: cell
        real(dp), intent(inout) :: xi(:)
        real(dp), intent(out) :: x(:), value
        integer, intent(out) :: status
        logical, intent(in), optional :: maximum
        real(dp) :: sgn, jac(mesh%dim, mesh%dim), gradient(mesh%dim), &
                    hessian(mesh%dim, mesh%dim), gp(mesh%dim), gm(mesh%dim), &
                    step(mesh%dim), txi(mesh%dim), trial(mesh%dim), h, tvalue, vp
        real(dp), allocatable :: nodes(:)
        integer :: it, d, e, a, node, tcell, half, n
        sgn = 1
        if (present(maximum)) then
            if (maximum) sgn = -1
        end if
        status = 1
        if (cell < 1 .or. cell > mesh%nelem) then
            node = minloc(sgn*u, 1)
            n = mesh%degree + 1
            nodes = locator%field_nodes
            outer: do e = 1, mesh%nelem
                do a = 1, size(mesh%dof, 1)
                    if (mesh%dof(a, e) /= node) cycle
                    cell = e
                    xi(1) = nodes(modulo(a - 1, n) + 1)
                    xi(2) = nodes(modulo((a - 1)/n, n) + 1)
                    if (mesh%dim == 3) xi(3) = nodes((a - 1)/n**2 + 1)
                    exit outer
                end do
            end do outer
        end if
        h = 1.e-4_dp
        call reference_gradient(mesh, locator, u, cell, xi(1:mesh%dim), value, &
                                gradient)
        do it = 1, 50
            do d = 1, mesh%dim
                txi = xi(1:mesh%dim)
                txi(d) = txi(d) + h
                call reference_gradient(mesh, locator, u, cell, txi, vp, gp)
                txi(d) = txi(d) - 2*h
                call reference_gradient(mesh, locator, u, cell, txi, vp, gm)
                hessian(:, d) = (gp - gm)/(2*h)
            end do
            hessian = 0.5_dp*(hessian + transpose(hessian))
            call small_solve(hessian, -gradient, step, status)
            if (status /= 0) return
            if (maxval(abs(step)) > 0.5_dp) step = 0.5_dp*step/maxval(abs(step))
            do half = 1, 40
                txi = xi(1:mesh%dim) + step
                tcell = cell
                status = 0
                if (maxval(abs(txi)) > 1) then
                    ! Leave the cell through the extrapolated geometry map.
                    call geometry_at(mesh, locator, cell, txi, trial, jac)
                    call locate_tensor_point(mesh, locator, trial, tcell, txi, &
                                             status)
                end if
                if (status == 0) then
                    call reference_gradient(mesh, locator, u, tcell, txi, tvalue, gp)
                    if (sgn*tvalue <= sgn*value + 1.e-15_dp*abs(value)) exit
                end if
                step = 0.5_dp*step
            end do
            if (half > 40) exit
            if (tcell /= cell) then
                ! Reference coordinates changed; restart Newton in the new cell.
                cell = tcell
                xi(1:mesh%dim) = txi
                value = tvalue
                gradient = gp
                cycle
            end if
            xi(1:mesh%dim) = txi
            value = tvalue
            gradient = gp
            if (maxval(abs(step)) < 1.e-11_dp) exit
        end do
        call geometry_at(mesh, locator, cell, xi(1:mesh%dim), x, jac)
        status = 0
    end subroutine
end module fortfem_tensor_points
