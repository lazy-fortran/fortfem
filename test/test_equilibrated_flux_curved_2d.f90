program test_equilibrated_flux_curved_2d
    !! Equilibrated-flux bound on isoparametric P2/P3 meshes of the unit square
    !! (interior cells and edges curved, boundary straight): (1) a random global
    !! Piola-mapped RT_(k-1) flux is recovered exactly by the RT_k patch
    !! reconstruction (the test first checks the normal-flux continuity of that
    !! flux across the curved shared edges with its own geometry); (2) for the
    !! P_m Galerkin solution of -div(A grad u) = f, A = 1/(1 + x), with the exact
    !! solution u = sin(pi x) sin(pi y), the bound exceeds the exact energy error
    !! on every mesh, the efficiency index stays below a measured cap and the
    !! bound converges at the energy rate m.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_rt_arbitrary_order, only: triangle_rt_basis_t, &
        evaluate_triangle_raviart_thomas
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_t, EQUILIBRATION_OK
    use equilibrated_test_mesh, only: test_mesh_t, build_test_mesh, cell_point, &
        galerkin_solve, random_rt_flux, edge_jump
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: amplitude = 0.06_dp
    type(equilibrated_flux_estimate_t) :: est
    real(dp) :: jump, error, bound(3), efficiency
    integer :: status, m, k, level, n, ns(3), kmin, kmax

    do m = 2, 3
        call recovery(m, m - 1, 3, est, jump, status)
        call check_condition(status == EQUILIBRATION_OK, "recovery solve succeeds")
        call check_condition(jump < 1.0e-12_dp, &
            "test flux has continuous normal trace across curved edges")
        call check_condition(est%bound < 1.0e-9_dp, &
            "RT_k recovers a Piola-mapped RT_(k-1) flux on a curved mesh")
        print '(a, i2, a, es10.2, a, es10.2)', "recovery P", m, " normal jump", jump, &
            " bound", est%bound
    end do

    do m = 2, 3
        ns = [4, 8, 16]
        if (m == 3) ns = [3, 6, 12]
        kmin = m
        kmax = m + 1
        do k = kmin, kmax
            do level = 1, 3
                n = ns(level)
                call manufactured(m, k, n, est, error, status)
                call check_condition(status == EQUILIBRATION_OK, &
                    "curved solve succeeds")
                efficiency = est%bound/error
                bound(level) = est%bound
                call check_condition(est%bound >= error, &
                    "curved bound is above the energy error")
                call check_condition(efficiency < 1.25_dp, &
                    "curved efficiency index is below its cap")
                print '(a, i2, a, i2, a, i3, 5(a, es10.2))', "curved P", m, " k", k, &
                    " n", n, " error", error, " bound", est%bound, " eff", &
                    efficiency, " flux", est%flux_defect, " osc", est%oscillation
            end do
            call check_condition(log(bound(2)/bound(3))/log(2.0_dp) > m - 0.4_dp, &
                "curved bound converges at the energy rate")
        end do
    end do
    call check_summary("equilibrated_flux_curved_2d")

contains

    subroutine rule(degree, xi, eta, w)
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: xi(:), eta(:), w(:)
        integer :: stat

        call triangle_duffy_quadrature(degree, xi, eta, w, stat)
        if (stat /= 0) error stop "quadrature failed"
    end subroutine rule

    subroutine recovery(m, k0, n, res, normal_jump, stat)
        !! Random global RT_k0 flux (see random_rt_flux), mapped by the pointwise
        !! Piola transformation of the isoparametric cells. A = 1, f = -div sigma.
        integer, intent(in) :: m, k0, n
        type(equilibrated_flux_estimate_t), intent(out) :: res
        real(dp), intent(out) :: normal_jump
        integer, intent(out) :: stat
        type(test_mesh_t) :: mesh
        type(triangle_rt_basis_t) :: rt
        real(dp), allocatable :: xi(:), eta(:), w(:), dofs(:, :)
        real(dp), allocatable :: coef(:, :), f(:, :), flux(:, :, :), lower(:)
        real(dp), allocatable :: phi(:), gref(:, :), rt_val(:, :), rt_div(:)
        integer :: nrt0, nt, t, e, a, b, j, i, q
        real(dp) :: x(2), jac(2, 2), det, sigma(2), divergence, curvature

        call build_test_mesh(n, m, amplitude, mesh, stat)
        if (stat /= 0) error stop "mesh"
        call random_rt_flux(mesh, k0, rt, dofs)
        nrt0 = size(dofs, 1)
        nt = mesh%nt

        ! Normal-flux continuity across every interior edge, own geometry.
        normal_jump = 0.0_dp
        do t = 1, nt
            do e = 1, 3
                a = mesh%triangles(e, t)
                b = mesh%triangles(mod(e, 3) + 1, t)
                if (a < b) then
                    do i = 1, nt
                        do j = 1, 3
                            if (mesh%triangles(j, i) == b .and. &
                                mesh%triangles(mod(j, 3) + 1, i) == a) then
                                normal_jump = max(normal_jump, &
                                    edge_jump(mesh, rt, dofs, i, j, t, e))
                            end if
                        end do
                    end do
                end if
            end do
        end do

        ! The mesh must really be curved: an edge midpoint node off its chord.
        curvature = 0.0_dp
        do t = 1, nt
            curvature = max(curvature, norm2(mesh%nodes(:, 4, t) &
                - 0.5_dp*(mesh%nodes(:, 1, t) + mesh%nodes(:, 2, t))))
        end do
        if (curvature < 1.0e-3_dp) error stop "test mesh is not curved"

        call rule(14, xi, eta, w)
        allocate (coef(size(w), nt), f(size(w), nt), flux(2, size(w), nt), lower(nt))
        allocate (phi(mesh%nb), gref(2, mesh%nb), rt_val(2, nrt0), rt_div(nrt0))
        coef = 1.0_dp
        lower = 1.0_dp
        do t = 1, nt
            do q = 1, size(w)
                call cell_point(mesh, t, xi(q), eta(q), x, jac, det, phi, gref)
                call evaluate_triangle_raviart_thomas(rt, xi(q), eta(q), rt_val, &
                    rt_div, stat)
                sigma = matmul(jac, matmul(rt_val, dofs(:, t)))/det
                divergence = dot_product(rt_div, dofs(:, t))/det
                flux(:, q, t) = sigma
                f(q, t) = -divergence
            end do
        end do
        call equilibrated_flux_estimate(mesh%vertices, mesh%triangles, k0 + 1, xi, &
            eta, w, coef, f, flux, lower, 1.0_dp/pi, res, stat, mesh%nodes)

    end subroutine recovery


    subroutine manufactured(m, degree, n, res, energy_error, stat)
        !! Galerkin P_m solution on the curved mesh; the exact energy error uses
        !! a higher-order rule than the estimator's.
        integer, intent(in) :: m, degree, n
        type(equilibrated_flux_estimate_t), intent(out) :: res
        real(dp), intent(out) :: energy_error
        integer, intent(out) :: stat
        type(test_mesh_t) :: mesh
        type(equilibrated_flux_estimate_t) :: wider
        real(dp), allocatable :: xi(:), eta(:), w(:), coef(:, :), f(:, :)
        real(dp), allocatable :: flux(:, :, :), lower(:), u(:), extra(:)
        real(dp), allocatable :: phi(:), gref(:, :), grad(:, :), xe(:), ee(:), we(:)
        logical, allocatable :: fixed(:)
        real(dp) :: x(2), jac(2, 2), det, inv_t(2, 2), gu(2), ge(2), c, lap
        integer :: t, q, i, info

        call build_test_mesh(n, m, amplitude, mesh, stat)
        if (stat /= 0) error stop "mesh"
        allocate (phi(mesh%nb), gref(2, mesh%nb), grad(2, mesh%nb))
        allocate (fixed(mesh%ndof), extra(mesh%ndof), u(mesh%ndof))
        extra = 0.0_dp
        do i = 1, mesh%ndof
            fixed(i) = any(abs(mesh%dof_position(:, i) - 0.5_dp) > 0.5_dp - 1.0e-12_dp)
        end do

        call rule(2*m + 4, xi, eta, w)
        allocate (coef(size(w), mesh%nt), f(size(w), mesh%nt), lower(mesh%nt))
        do t = 1, mesh%nt
            do q = 1, size(w)
                call cell_point(mesh, t, xi(q), eta(q), x, jac, det, phi, gref)
                c = 1.0_dp/(1.0_dp + x(1))
                coef(q, t) = c
                lap = -2*pi*pi*sin(pi*x(1))*sin(pi*x(2))
                f(q, t) = -(-c*c*pi*cos(pi*x(1))*sin(pi*x(2)) + c*lap)
            end do
        end do
        call galerkin_solve(mesh, xi, eta, w, coef, f, fixed, extra, u, info)
        if (info /= 0) error stop "galerkin"
        lower = 0.5_dp

        ! Estimator data: flux A grad u_h at the estimator's points.
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
                lap = -2*pi*pi*sin(pi*x(1))*sin(pi*x(2))
                f(q, t) = -(-c*c*pi*cos(pi*x(1))*sin(pi*x(2)) + c*lap)
                flux(:, q, t) = c*gu
            end do
        end do
        call equilibrated_flux_estimate(mesh%vertices, mesh%triangles, degree, xi, &
            eta, w, coef, f, flux, lower, 1.0_dp/pi, res, stat, mesh%nodes)

        if (m == 3 .and. degree == 3 .and. n == 12) then
            ! A valid but looser global C_F must not inflate the local source
            ! oscillation. This fails when curved cells always use C_F.
            call equilibrated_flux_estimate(mesh%vertices, mesh%triangles, degree, &
                xi, eta, w, coef, f, flux, lower, 2.0_dp/pi, wider, stat, mesh%nodes)
            call check_condition(stat == EQUILIBRATION_OK, "looser C_F solve")
            call check_condition(res%oscillation > 0.0_dp, "nonzero oscillation")
            call check_condition(abs(wider%oscillation - res%oscillation) &
                < 1.0e-12_dp*res%oscillation, "mapped oscillation uses local constant")
        end if

        call rule(24, xe, ee, we)
        energy_error = 0.0_dp
        do t = 1, mesh%nt
            do q = 1, size(we)
                call cell_point(mesh, t, xe(q), ee(q), x, jac, det, phi, gref)
                inv_t(1, :) = [jac(2, 2), -jac(2, 1)]/det
                inv_t(2, :) = [-jac(1, 2), jac(1, 1)]/det
                grad = matmul(inv_t, gref)
                gu = matmul(grad, u(mesh%gdof(:, t)))
                ge = pi*[cos(pi*x(1))*sin(pi*x(2)), sin(pi*x(1))*cos(pi*x(2))]
                energy_error = energy_error + we(q)*det/(1.0_dp + x(1)) &
                    *sum((ge - gu)**2)
            end do
        end do
        energy_error = sqrt(energy_error)
    end subroutine manufactured

end program test_equilibrated_flux_curved_2d
