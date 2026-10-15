program test_equilibrated_flux_mixed_2d
    !! Equilibrated-flux bound for mixed Dirichlet/Neumann data on the unit
    !! square (affine cells, P1 and P2 geometry): (1) the normal trace of an exact
    !! RT_(k-1) flux prescribed as Neumann data on part or all of the boundary is
    !! recovered by the RT_k patch reconstruction, so the bound vanishes to
    !! rounding (a wrong edge sign, moment order or normal-trace scale gives O(1));
    !! (2) for the P_m Galerkin solution of -div(A grad u) = f, A = 1/(1 + x),
    !! u = sin(pi x/3) sin(pi y/3), Dirichlet data on x = 0 and y = 0 and flux
    !! data A grad u . n on x = 1 and y = 1, the bound exceeds the exact energy
    !! error, the efficiency index stays below a measured cap and the bound
    !! converges at the energy rate.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_rt_arbitrary_order, only: triangle_rt_basis_t, &
        evaluate_triangle_raviart_thomas
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_t, EQUILIBRATION_OK
    use equilibrated_test_mesh, only: test_mesh_t, build_test_mesh, cell_point, &
        galerkin_solve, random_rt_flux, edge_flux
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: third_pi = pi/3.0_dp
    type(equilibrated_flux_estimate_t) :: est
    real(dp) :: error, bound(3), efficiency
    integer :: status, m, k, level, n, mode, ns(3)

    do mode = 1, 3
        do m = 1, 2
            call recovery(m, 2, mode, est, status)
            call check_condition(status == EQUILIBRATION_OK, "Neumann recovery solve")
            call check_condition(est%bound < 1.0e-9_dp, &
                "RT_k recovers an RT_(k-1) flux from its Neumann trace")
            print '(a, i2, a, i2, 3(a, es10.2))', "recovery mode", mode, " P", m, &
                " bound", est%bound, " flux", est%flux_defect, " neumann", &
                est%neumann_oscillation
        end do
    end do

    do m = 1, 2
        ns = [4, 8, 16]
        do k = m, m + 1
            do level = 1, 3
                n = ns(level)
                call manufactured(m, k, n, est, error, status)
                call check_condition(status == EQUILIBRATION_OK, "mixed solve")
                efficiency = est%bound/error
                bound(level) = est%bound
                call check_condition(est%bound >= error, &
                    "mixed bound is above the energy error")
                call check_condition(efficiency < 1.2_dp, &
                    "mixed efficiency index is below its cap")
                print '(a, i2, a, i2, a, i3, 6(a, es10.2))', "mixed P", m, " k", k, &
                    " n", n, " error", error, " bound", est%bound, " eff", &
                    efficiency, " flux", est%flux_defect, " div", &
                    est%divergence_defect, " neumann", est%neumann_oscillation
            end do
            call check_condition(log(bound(2)/bound(3))/log(2.0_dp) > m - 0.4_dp, &
                "mixed bound converges at the energy rate")
        end do
    end do
    call check_summary("equilibrated_flux_mixed_2d")

contains

    subroutine rule(degree, xi, eta, w)
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: xi(:), eta(:), w(:)
        integer :: stat

        call triangle_duffy_quadrature(degree, xi, eta, w, stat)
        if (stat /= 0) error stop "quadrature failed"
    end subroutine rule

    logical function neumann_side(mesh, t, e, mode) result(flag)
        !! Local edge e of cell t lies on a Neumann side of the unit square:
        !! mode 1: x = 1 or y = 1, mode 2: the whole boundary, mode 3: x = 1.
        type(test_mesh_t), intent(in) :: mesh
        integer, intent(in) :: t, e, mode
        real(dp) :: p(2), q(2)

        p = mesh%vertices(:, mesh%triangles(e, t))
        q = mesh%vertices(:, mesh%triangles(mod(e, 3) + 1, t))
        flag = .false.
        if (all([p(1), q(1)] == 1.0_dp)) flag = .true.
        if (mode == 3) return
        if (all([p(2), q(2)] == 1.0_dp)) flag = .true.
        if (mode == 1) return
        if (all([p(1), q(1)] == 0.0_dp)) flag = .true.
        if (all([p(2), q(2)] == 0.0_dp)) flag = .true.
    end function neumann_side

    subroutine recovery(m, k0, mode, res, stat)
        !! Random global RT_k0 flux on an affine mesh; A = 1, f = -div sigma and
        !! g = sigma . n on the Neumann edges (mode as in neumann_side).
        integer, intent(in) :: m, k0, mode
        type(equilibrated_flux_estimate_t), intent(out) :: res
        integer, intent(out) :: stat
        type(test_mesh_t) :: mesh
        type(triangle_rt_basis_t) :: rt
        real(dp), allocatable :: xi(:), eta(:), w(:), dofs(:, :), edge_pts(:)
        real(dp), allocatable :: edge_w(:), coef(:, :), f(:, :), flux(:, :, :)
        real(dp), allocatable :: lower(:), gdata(:, :, :), phi(:), gref(:, :)
        real(dp), allocatable :: rt_val(:, :), rt_div(:)
        logical, allocatable :: flags(:, :)
        real(dp) :: x(2), jac(2, 2), det, speed
        integer :: t, e, q, nt, nq1

        call build_test_mesh(4, m, 0.0_dp, mesh, stat)
        if (stat /= 0) error stop "mesh"
        call random_rt_flux(mesh, k0, rt, dofs)
        nt = mesh%nt
        nq1 = k0 + 4
        allocate (edge_pts(nq1), edge_w(nq1))
        call gauss_legendre_ab(nq1, 0.0_dp, 1.0_dp, edge_pts, edge_w)
        call rule(2*k0 + 6, xi, eta, w)
        allocate (coef(size(w), nt), f(size(w), nt), flux(2, size(w), nt))
        allocate (lower(nt), flags(3, nt), gdata(nq1, 3, nt))
        allocate (phi(mesh%nb), gref(2, mesh%nb))
        allocate (rt_val(2, size(dofs, 1)), rt_div(size(dofs, 1)))
        coef = 1.0_dp
        lower = 1.0_dp
        do t = 1, nt
            do q = 1, size(w)
                call cell_point(mesh, t, xi(q), eta(q), x, jac, det, phi, gref)
                call evaluate_triangle_raviart_thomas(rt, xi(q), eta(q), rt_val, &
                    rt_div, stat)
                flux(:, q, t) = matmul(jac, matmul(rt_val, dofs(:, t)))/det
                f(q, t) = -dot_product(rt_div, dofs(:, t))/det
            end do
            do e = 1, 3
                flags(e, t) = neumann_side(mesh, t, e, mode)
                speed = norm2(mesh%vertices(:, mesh%triangles(mod(e, 3) + 1, t)) &
                    - mesh%vertices(:, mesh%triangles(e, t)))
                do q = 1, nq1
                    gdata(q, e, t) = edge_flux(mesh, rt, dofs, t, e, edge_pts(q))/speed
                end do
            end do
        end do
        call equilibrated_flux_estimate(mesh%vertices, mesh%triangles, k0 + 1, xi, &
            eta, w, coef, f, flux, lower, 1.0_dp/pi, res, stat, mesh%nodes, flags, &
            gdata, edge_pts, edge_w)
    end subroutine recovery

    pure function exact_gradient(x) result(g)
        real(dp), intent(in) :: x(2)
        real(dp) :: g(2)

        g = third_pi*[cos(third_pi*x(1))*sin(third_pi*x(2)), &
            sin(third_pi*x(1))*cos(third_pi*x(2))]
    end function exact_gradient

    pure real(dp) function exact_source(x) result(f)
        real(dp), intent(in) :: x(2)
        real(dp) :: c, lap, grad(2)

        c = 1.0_dp/(1.0_dp + x(1))
        grad = exact_gradient(x)
        lap = -2*third_pi**2*sin(third_pi*x(1))*sin(third_pi*x(2))
        f = -(-c*c*grad(1) + c*lap)
    end function exact_source

    subroutine manufactured(m, degree, n, res, energy_error, stat)
        integer, intent(in) :: m, degree, n
        type(equilibrated_flux_estimate_t), intent(out) :: res
        real(dp), intent(out) :: energy_error
        integer, intent(out) :: stat
        type(test_mesh_t) :: mesh
        real(dp), allocatable :: xi(:), eta(:), w(:), coef(:, :), f(:, :)
        real(dp), allocatable :: flux(:, :, :), lower(:), u(:), extra(:)
        real(dp), allocatable :: phi(:), gref(:, :), grad(:, :), edge_pts(:)
        real(dp), allocatable :: edge_w(:), gdata(:, :, :), xe(:), ee(:), we(:)
        logical, allocatable :: fixed(:), flags(:, :)
        real(dp) :: x(2), jac(2, 2), det, inv_t(2, 2), gu(2), c, speed, normal(2)
        real(dp) :: p(2), tangent(2), s, cf
        integer :: t, q, i, e, info, nq1

        call build_test_mesh(n, m, 0.0_dp, mesh, stat)
        if (stat /= 0) error stop "mesh"
        allocate (phi(mesh%nb), gref(2, mesh%nb), grad(2, mesh%nb))
        allocate (fixed(mesh%ndof), extra(mesh%ndof), u(mesh%ndof))
        allocate (flags(3, mesh%nt), lower(mesh%nt))
        nq1 = degree + 4
        allocate (edge_pts(nq1), edge_w(nq1), gdata(nq1, 3, mesh%nt))
        call gauss_legendre_ab(nq1, 0.0_dp, 1.0_dp, edge_pts, edge_w)
        do i = 1, mesh%ndof
            fixed(i) = any(mesh%dof_position(:, i) == 0.0_dp)
        end do

        ! Neumann data g = A grad u . n on x = 1 and y = 1 (mode 1) and its
        ! load int g phi_i, assembled with the basis of the cells.
        extra = 0.0_dp
        do t = 1, mesh%nt
            do e = 1, 3
                flags(e, t) = neumann_side(mesh, t, e, 1)
                p = mesh%vertices(:, mesh%triangles(e, t))
                tangent = mesh%vertices(:, mesh%triangles(mod(e, 3) + 1, t)) - p
                speed = norm2(tangent)
                normal = [tangent(2), -tangent(1)]/speed
                do q = 1, nq1
                    s = edge_pts(q)
                    x = p + s*tangent
                    c = 1.0_dp/(1.0_dp + x(1))
                    gdata(q, e, t) = c*dot_product(exact_gradient(x), normal)
                    if (.not. flags(e, t)) cycle
                    call edge_reference_point(e, s, gref(:, 1))
                    call cell_point(mesh, t, gref(1, 1), gref(2, 1), x, jac, det, &
                        phi, grad)
                    do i = 1, mesh%nb
                        extra(mesh%gdof(i, t)) = extra(mesh%gdof(i, t)) &
                            + edge_w(q)*speed*gdata(q, e, t)*phi(i)
                    end do
                end do
            end do
        end do

        call rule(2*m + 4, xi, eta, w)
        allocate (coef(size(w), mesh%nt), f(size(w), mesh%nt))
        do t = 1, mesh%nt
            do q = 1, size(w)
                call cell_point(mesh, t, xi(q), eta(q), x, jac, det, phi, gref)
                coef(q, t) = 1.0_dp/(1.0_dp + x(1))
                f(q, t) = exact_source(x)
            end do
        end do
        call galerkin_solve(mesh, xi, eta, w, coef, f, fixed, extra, u, info)
        if (info /= 0) error stop "galerkin"

        call rule(2*degree + 2, xi, eta, w)
        deallocate (coef, f)
        allocate (coef(size(w), mesh%nt), f(size(w), mesh%nt))
        allocate (flux(2, size(w), mesh%nt))
        do t = 1, mesh%nt
            do q = 1, size(w)
                call cell_point(mesh, t, xi(q), eta(q), x, jac, det, phi, gref)
                inv_t(1, :) = [jac(2, 2), -jac(2, 1)]/det
                inv_t(2, :) = [-jac(1, 2), jac(1, 1)]/det
                grad = matmul(inv_t, gref)
                gu = matmul(grad, u(mesh%gdof(:, t)))
                c = 1.0_dp/(1.0_dp + x(1))
                coef(q, t) = c
                f(q, t) = exact_source(x)
                flux(:, q, t) = c*gu
            end do
            ! A >= 1/(1 + max x) on the cell.
            lower(t) = 1.0_dp/(1.0_dp + maxval(mesh%vertices(1, mesh%triangles(:, t))))
        end do
        ! Friedrichs constant for v = 0 on x = 0 only (1D, first Dirichlet-
        ! Neumann eigenvalue (pi/2)^2) and A >= 1/2: ||v|| <= (2/pi) sqrt(2) |||v|||.
        cf = 2.0_dp*sqrt(2.0_dp)/pi
        call equilibrated_flux_estimate(mesh%vertices, mesh%triangles, degree, xi, &
            eta, w, coef, f, flux, lower, cf, res, stat, mesh%nodes, flags, gdata, &
            edge_pts, edge_w)

        call rule(24, xe, ee, we)
        energy_error = 0.0_dp
        do t = 1, mesh%nt
            do q = 1, size(we)
                call cell_point(mesh, t, xe(q), ee(q), x, jac, det, phi, gref)
                inv_t(1, :) = [jac(2, 2), -jac(2, 1)]/det
                inv_t(2, :) = [-jac(1, 2), jac(1, 1)]/det
                grad = matmul(inv_t, gref)
                gu = matmul(grad, u(mesh%gdof(:, t)))
                energy_error = energy_error + we(q)*det/(1.0_dp + x(1)) &
                    *sum((exact_gradient(x) - gu)**2)
            end do
        end do
        energy_error = sqrt(energy_error)
    end subroutine manufactured

    pure subroutine edge_reference_point(e, s, p)
        integer, intent(in) :: e
        real(dp), intent(in) :: s
        real(dp), intent(out) :: p(2)

        select case (e)
        case (1)
            p = [s, 0.0_dp]
        case (2)
            p = [1.0_dp - s, s]
        case default
            p = [0.0_dp, 1.0_dp - s]
        end select
    end subroutine edge_reference_point

end program test_equilibrated_flux_mixed_2d
