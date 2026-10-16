program test_verified_bernstein_majorant
    !! Verified Bernstein audit of sampled RT_k equilibrated fluxes (k = 1, 2, 3).
    !! Independent oracles: Duffy-quadrature energy errors and the sampled
    !! estimator's own exactly integrated flux defect; shared-edge repair,
    !! exact zero data and rejections.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortnum_interval, only: interval_t, interval, operator(+), operator(-), &
                                operator(*)
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_t, EQUILIBRATION_OK
    use fortfem_verified_bernstein_majorant, only: verified_bernstein_majorant, &
        VERIFIED_BERNSTEIN_OK, VERIFIED_BERNSTEIN_BAD_SOURCE
    use fortfem_verified_bernstein_products, only: verified_bernstein_affine_product
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: vertices(:, :), xi(:), eta(:), w(:), weight(:)
    real(dp), allocatable :: nodal(:)
    integer, allocatable :: triangles(:, :), cells(:, :)
    type(equilibrated_flux_estimate_t) :: est
    type(interval_t) :: bound, flux, residual, cf
    real(dp) :: error, sampled
    integer :: status, k, mode, tail_mode, power

    ! C_F = 1/(sqrt(2) pi sqrt(A_min)) with A_min = 1/2; 0.3184 > 1/pi.
    cf = interval(0.3184_dp)
    call square_mesh(8, vertices, triangles)
    call triangle_duffy_quadrature(14, xi, eta, w, status)
    call check_condition(status == 0, "quadrature")
    do k = 1, 2
        ! A = 1/(1+x): certified bound above the exact interpolation error and
        ! close to the sampled estimate (chord and global C_F only).
        mode = 1
        call interpolant_case(k, est, error, status)
        call check_condition(status == EQUILIBRATION_OK, "sampled equilibration")
        call verified_bernstein_majorant(vertices, triangles, k, &
            est%flux_coefficients, cells, nodal, weight, numerator, 4, 2, cf, &
            bound, flux, residual, status)
        call check_condition(status == VERIFIED_BERNSTEIN_OK, "weighted audit")
        call check_condition(bound%hi >= error, "certified bound above error")
        call check_condition(bound%hi <= 1.25_dp*est%bound, "certified bound narrow")
        call check_condition(abs(sqrt(flux%hi) - est%flux_defect) &
                             <= 0.02_dp*est%flux_defect, "chord flux defect")
        print '(a,i1,3es12.4)', 'k', k, error, est%bound, bound%hi
        ! A = 1, f = 2 in P_k: exact polynomial integrals reproduce the sampled
        ! flux defect, and the residual is the sampled divergence defect only.
        mode = 2
        call interpolant_case(k, est, error, status)
        call check_condition(status == EQUILIBRATION_OK, "unit equilibration")
        call verified_bernstein_majorant(vertices, triangles, k, &
            est%flux_coefficients, cells, nodal, weight, numerator, 0, 0, cf, &
            bound, flux, residual, status)
        call check_condition(status == VERIFIED_BERNSTEIN_OK, "unit audit")
        sampled = est%flux_defect**2
        call check_condition(flux%lo <= sampled*(1 + 1.0e-9_dp) .and. &
                             flux%hi >= sampled*(1 - 1.0e-9_dp), &
                             "exact flux defect integral")
        call check_condition(abs(sqrt(residual%hi) - est%divergence_defect/0.3184_dp) &
                             <= 1.0e-8_dp*est%bound, "residual is divergence defect")
        ! Deliberately inconsistent neighbour coefficients are canonicalized.
        sampled = bound%hi
        est%flux_coefficients(:, 2) = est%flux_coefficients(:, 2)*(1 + 1.0e-9_dp)
        call verified_bernstein_majorant(vertices, triangles, k, &
            est%flux_coefficients, cells, nodal, weight, numerator, 0, 0, cf, &
            bound, flux, residual, status)
        call check_condition(status == VERIFIED_BERNSTEIN_OK, "shared edges")
        call check_condition(abs(bound%hi - sampled) <= 1.0e-6_dp*sampled, &
                             "shared edge repair keeps the bound")
        ! Exact zero data: only outward rounding remains.
        est%flux_coefficients = 0; nodal = 0; mode = 3
        call verified_bernstein_majorant(vertices, triangles, k, &
            est%flux_coefficients, cells, nodal, weight, numerator, 0, 0, cf, &
            bound, flux, residual, status)
        call check_condition(status == VERIFIED_BERNSTEIN_OK .and. &
                             bound%hi <= 1.0e-100_dp, "zero containment")
        tail_mode = 0
        call verified_bernstein_majorant(vertices, triangles, k, &
            est%flux_coefficients, cells, nodal, weight, numerator, 0, 0, cf, &
            bound, flux, residual, status, source_remainder=tail)
        call check_condition(status == VERIFIED_BERNSTEIN_OK .and. &
                             bound%hi <= 1.0e-100_dp, "zero source tail")
        tail_mode = 1
        weight = 2
        do power = 0, 2
            call verified_bernstein_majorant(vertices, triangles, k, &
                est%flux_coefficients, cells, nodal, weight, numerator, 0, power, cf, &
                bound, flux, residual, status, source_remainder=tail)
            call check_condition(status == VERIFIED_BERNSTEIN_OK, "source tail audit")
            sampled = 4.0_dp/4.0_dp**power
            call check_condition(residual%hi >= sampled .and. &
                residual%hi < sampled + 1.0e-10_dp, "weighted constant tail L2 norm")
        end do
        do tail_mode = 2, 4
            call verified_bernstein_majorant(vertices, triangles, k, &
                est%flux_coefficients, cells, nodal, weight, numerator, 0, 0, cf, &
                bound, flux, residual, status, source_remainder=tail)
            call check_condition(status == VERIFIED_BERNSTEIN_BAD_SOURCE, &
                                 "invalid source tail rejected")
        end do
        mode = 4
        call verified_bernstein_majorant(vertices, triangles, k, &
            est%flux_coefficients, cells, nodal, weight, numerator, 0, 0, cf, &
            bound, flux, residual, status)
        call check_condition(status == VERIFIED_BERNSTEIN_BAD_SOURCE, &
                             "rejected source")
        weight(1) = 0
        mode = 3
        call verified_bernstein_majorant(vertices, triangles, k, &
            est%flux_coefficients, cells, nodal, weight, numerator, 0, 0, cf, &
            bound, flux, residual, status)
        call check_condition(status /= VERIFIED_BERNSTEIN_OK, "nonpositive weight")
    end do
    call verified_bernstein_majorant(vertices, triangles, 4, &
        est%flux_coefficients, cells, nodal, weight, numerator, 0, 0, cf, &
        bound, flux, residual, status)
    call check_condition(status /= VERIFIED_BERNSTEIN_OK, "unsupported degree")
    call cubic_recovery()
    call check_summary("verified_bernstein_majorant")

contains

    subroutine cubic_recovery()
        !! Nonzero exact P3 state on one triangle: u = xy(1-x-y),
        !! f = 2(x+y), A = 1. Its quadratic gradient is inside RT_3.
        real(dp) :: points(2, 3), values(10), a(size(w), 1), f(size(w), 1)
        real(dp) :: dflux(2, size(w), 1), lower(1), xx, yy
        integer :: tri(3, 1), nodes(10, 1), q, stat, j
        type(equilibrated_flux_estimate_t) :: recovered
        type(interval_t) :: radius, flux2, residual2

        points = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        tri(:, 1) = [1, 2, 3]
        nodes(:, 1) = [(j, j = 1, 10)]
        values = 0; values(10) = 1.0_dp/27
        a = 1; lower = 1
        do q = 1, size(w)
            xx = xi(q); yy = eta(q)
            f(q, 1) = 2*(xx + yy)
            dflux(:, q, 1) = [yy*(1 - 2*xx - yy), xx*(1 - xx - 2*yy)]
        end do
        call equilibrated_flux_estimate(points, tri, 3, xi, eta, w, a, f, &
            dflux, lower, 1.0_dp, recovered, stat)
        call check_condition(stat == EQUILIBRATION_OK, "cubic exact reconstruction")
        call verified_bernstein_majorant(points, tri, 3, &
            recovered%flux_coefficients, nodes, values, [1.0_dp, 1.0_dp, 1.0_dp], &
            cubic_numerator, 1, 0, interval(1.0_dp), radius, flux2, residual2, stat)
        call check_condition(stat == VERIFIED_BERNSTEIN_OK, "cubic exact audit")
        call check_condition(radius%hi < 1.0e-11_dp, "nonzero cubic exact recovery")
        print '(a,es13.5)', 'cubic exact radius ', radius%hi
    end subroutine cubic_recovery

    subroutine cubic_numerator(cell, coefficients, ok)
        integer, intent(in) :: cell
        type(interval_t), intent(out) :: coefficients(:)
        logical, intent(out) :: ok
        coefficients = interval([0.0_dp, 2.0_dp, 2.0_dp])
        ok = cell == 1
    end subroutine cubic_numerator

    subroutine tail(cell, radius, ok)
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        integer, intent(in) :: cell
        type(interval_t), intent(out) :: radius
        logical, intent(out) :: ok
        ok = tail_mode /= 3 .and. cell > 0
        radius = interval(0.0_dp)
        if (tail_mode == 1) radius = interval(2.0_dp)
        if (tail_mode == 2) radius = interval(-1.0_dp)
        if (tail_mode == 4) radius%hi = ieee_value(0.0_dp, ieee_quiet_nan)
    end subroutine tail

    subroutine square_mesh(n, points, cells)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: points(:, :)
        integer, allocatable, intent(out) :: cells(:, :)
        integer :: i, j, v, c

        allocate (points(2, (n + 1)**2), cells(3, 2*n*n))
        do j = 0, n
            do i = 0, n
                points(:, j*(n + 1) + i + 1) = [real(i, dp), real(j, dp)]/n
            end do
        end do
        c = 0
        do j = 0, n - 1
            do i = 0, n - 1
                v = j*(n + 1) + i + 1
                ! Alternate diagonals so both edge orientations occur.
                if (mod(i + j, 2) == 0) then
                    cells(:, c + 1) = [v, v + 1, v + n + 2]
                    cells(:, c + 2) = [v, v + n + 2, v + n + 1]
                else
                    cells(:, c + 1) = [v, v + 1, v + n + 1]
                    cells(:, c + 2) = [v + 1, v + n + 2, v + n + 1]
                end if
                c = c + 2
            end do
        end do
    end subroutine square_mesh

    pure subroutine exact(x, u, grad)
        !! u = x(1-x) y(1-y) + x y (1-x)(1-y)(x-y) vanishes on the boundary.
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: u, grad(2)
        real(dp) :: b, bx, by
        b = x(1)*(1 - x(1))*x(2)*(1 - x(2))
        bx = (1 - 2*x(1))*x(2)*(1 - x(2))
        by = x(1)*(1 - x(1))*(1 - 2*x(2))
        u = b*(1 + x(1) - x(2))
        grad = [bx*(1 + x(1) - x(2)) + b, by*(1 + x(1) - x(2)) - b]
    end subroutine exact

    subroutine interpolant_case(degree, res, energy_error, stat)
        !! Lagrange interpolant of degree `degree` of u; mode 1: A = 1/(1+x)
        !! and f = -div(A grad u); mode 2: A = 1 and f = 2.
        integer, intent(in) :: degree
        type(equilibrated_flux_estimate_t), intent(out) :: res
        real(dp), intent(out) :: energy_error
        integer, intent(out) :: stat
        real(dp), allocatable :: a(:, :), f(:, :), dflux(:, :, :), lower(:)
        real(dp) :: x(2), p(2, 3), jac(2, 2), det, u, grad(2), lap, coef, g(2)
        real(dp) :: values(6), shape_grad(2, 6)
        integer :: t, q, nt, j, nloc

        nt = size(triangles, 2)
        nloc = merge(3, 6, degree == 1)
        allocate (a(size(w), nt), f(size(w), nt), dflux(2, size(w), nt), lower(nt))
        if (allocated(cells)) deallocate (cells, nodal)
        call number_nodes(nloc, cells, nodal)
        weight = 1 + vertices(1, :)
        if (mode == 2) weight = 1
        energy_error = 0
        do t = 1, nt
            p = vertices(:, triangles(:, t))
            jac(:, 1) = p(:, 2) - p(:, 1)
            jac(:, 2) = p(:, 3) - p(:, 1)
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            values(1:nloc) = nodal(cells(:, t))
            lower(t) = 1/maxval(1 + p(1, :))
            if (mode == 2) lower(t) = 1
            do q = 1, size(w)
                x = p(:, 1) + xi(q)*jac(:, 1) + eta(q)*jac(:, 2)
                call exact(x, u, grad)
                call shape_gradients(degree, xi(q), eta(q), jac, det, shape_grad)
                g = matmul(shape_grad(:, 1:nloc), values(1:nloc))
                coef = 1/(1 + x(1))
                lap = -2*x(2)*(1 - x(2))*(1 + x(1) - x(2)) &
                      - 2*x(1)*(1 - x(1))*(1 + x(1) - x(2)) &
                      + 2*(1 - 2*x(1))*x(2)*(1 - x(2)) &
                      - 2*x(1)*(1 - x(1))*(1 - 2*x(2))
                if (mode == 2) then
                    coef = 1; f(q, t) = 2
                else
                    f(q, t) = -(coef*lap - coef*coef*grad(1))
                end if
                a(q, t) = coef
                dflux(:, q, t) = coef*g
                energy_error = energy_error + w(q)*det*coef*sum((grad - g)**2)
            end do
        end do
        energy_error = sqrt(energy_error)
        call equilibrated_flux_estimate(vertices, triangles, degree, xi, eta, w, a, &
            f, dflux, lower, 0.3184_dp, res, stat)
    end subroutine interpolant_case

    subroutine number_nodes(nloc, cells, nodal)
        !! Global P1/P2 Lagrange nodes (edge midpoints from the lower global
        !! vertex) and the interpolated values of u, so u_h is continuous.
        integer, intent(in) :: nloc
        integer, allocatable, intent(out) :: cells(:, :)
        real(dp), allocatable, intent(out) :: nodal(:)
        integer, allocatable :: ends(:, :)
        real(dp) :: g(2)
        integer :: t, j, a, b, e, n, nv
        nv = size(vertices, 2)
        allocate (cells(nloc, size(triangles, 2)), ends(2, 3*size(triangles, 2)))
        allocate (nodal(nv + 3*size(triangles, 2)))
        do j = 1, nv
            call exact(vertices(:, j), nodal(j), g)
        end do
        n = 0
        do t = 1, size(triangles, 2)
            cells(1:3, t) = triangles(:, t)
            if (nloc == 3) cycle
            do j = 1, 3
                a = minval(triangles([j, mod(j, 3) + 1], t))
                b = maxval(triangles([j, mod(j, 3) + 1], t))
                do e = 1, n
                    if (ends(1, e) == a .and. ends(2, e) == b) exit
                end do
                if (e > n) then
                    n = n + 1
                    ends(:, n) = [a, b]
                    call exact(0.5_dp*(vertices(:, a) + vertices(:, b)), &
                               nodal(nv + n), g)
                end if
                cells(3 + j, t) = nv + e
            end do
        end do
        ! Keep only defined values: the audit rejects any non-finite entry.
        nodal = nodal(1:nv + n)
    end subroutine number_nodes

    pure subroutine shape_gradients(degree, s, r, jac, det, grads)
        integer, intent(in) :: degree
        real(dp), intent(in) :: s, r, jac(2, 2), det
        real(dp), intent(out) :: grads(2, 6)
        real(dp) :: ref(2, 6), l(3), dl(2, 3), inv_t(2, 2)
        integer :: j, next
        l = [1 - s - r, s, r]
        dl = reshape([-1.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
        if (degree == 1) then
            ref(:, 1:3) = dl
        else
            do j = 1, 3
                next = mod(j, 3) + 1
                ref(:, j) = (4*l(j) - 1)*dl(:, j)
                ref(:, 3 + j) = 4*(l(j)*dl(:, next) + l(next)*dl(:, j))
            end do
        end if
        inv_t = reshape([jac(2, 2), -jac(1, 2), -jac(2, 1), jac(1, 1)], [2, 2])/det
        grads = 0
        do j = 1, merge(3, 6, degree == 1)
            grads(:, j) = matmul(inv_t, ref(:, j))
        end do
    end subroutine shape_gradients

    subroutine numerator(cell, coefficients, ok)
        !! Bernstein coefficients of s = (1+x)^2 f: mode 1 degree 4,
        !! s = -(1+x) lap u + u_x; mode 2 the constant 2; mode 3 zero.
        integer, intent(in) :: cell
        type(interval_t), intent(out) :: coefficients(:)
        logical, intent(out) :: ok
        type(interval_t) :: x(3), y(3), one(3), mx(3), my(3), s(28)
        integer :: j
        ok = mode /= 4
        coefficients = interval(0.0_dp)
        if (mode == 2) coefficients = interval(2.0_dp)
        if (mode /= 1) return
        do j = 1, 3
            x(j) = interval(vertices(1, triangles(j, cell)))
            y(j) = interval(vertices(2, triangles(j, cell)))
        end do
        one = interval(1.0_dp); mx = one - x; my = one - y
        ! s = (1+x)(2 y(1-y)(1+x-y) + 2 x(1-x)(1+x-y)
        !     - 2(1-2x) y(1-y) + 2 x(1-x)(1-2y))
        !     + (1-2x) y(1-y)(1+x-y) + x(1-x) y(1-y).
        s = interval(0.0_dp)
        s = s + 2.0_dp*product4([y, my, one + x - y, one + x])
        s = s + 2.0_dp*product4([x, mx, one + x - y, one + x])
        s = s - 2.0_dp*product4([one - 2.0_dp*x, y, my, one + x])
        s = s + 2.0_dp*product4([x, mx, one - 2.0_dp*y, one + x])
        s = s + product4([one - 2.0_dp*x, y, my, one + x - y])
        s = s + product4([x, mx, y, my])
        ok = size(coefficients) == 15
        if (ok) coefficients = s(1:15)
    end subroutine numerator

    function product4(factors) result(r)
        !! Degree-4 Bernstein coefficients (first 15 entries) of the product of
        !! four affine factors given by their vertex values.
        type(interval_t), intent(in) :: factors(12)
        type(interval_t) :: r(28), c(28)
        integer :: f
        c = interval(1.0_dp)
        do f = 0, 3
            call verified_bernstein_affine_product(f, c, factors(3*f + 1:3*f + 3), r)
            c = r
        end do
    end function product4

end program test_verified_bernstein_majorant
