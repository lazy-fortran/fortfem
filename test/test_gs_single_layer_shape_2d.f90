program test_gs_single_layer_shape_2d
    !! Single-layer operator of the axisymmetric Grad-Shafranov kernel on
    !! curved P_k panels and its shape derivatives.  Oracles:
    !!  - for P2 panels, ops%v of gs_bem_assemble (separately written
    !!    kernel and assembly);
    !!  - the same curve written as P3 panels (sampled from the P2 panels,
    !!    so the geometry is identical) gives the same V;
    !!  - d(V tau)/dX, the gradient and the Hessian of tau'V tau/2, and the
    !!    load derivatives against central differences at the frozen plan.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: &
        gs_bem_curve_t, gs_bem_operators_t, gs_bem_curve_from_nodes, &
        gs_bem_assemble, gs_bem_panel_point
    use fortfem_gs_single_layer_shape_2d, only: gs_sl_curve_t, gs_sl_plan_t, &
        gs_sl_curve_create, gs_sl_plan_create, gs_sl_assemble, gs_sl_load, &
        radial_n, radial_q, near_order, far_order, sub_scale
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: np = 10, flux_degree = 1
    logical :: all_passed

    all_passed = .true.
    call test_against_p2_operator()
    call test_shape_derivatives(3)
    call test_shape_derivatives(2)
    call test_load_derivatives()
    call test_quadrature_converged()
    call check_summary("GS single-layer shape derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine test_quadrature_converged()
        !! The default rules against rules refined in every parameter: V,
        !! d(V tau)/dX and the shape gradient agree to quadrature accuracy.
        type(gs_sl_curve_t) :: curve
        type(gs_sl_plan_t) :: plan
        real(dp), allocatable :: v0(:, :), v1(:, :), d0(:, :), d1(:, :), g0(:), g1(:)
        real(dp), allocatable :: tau(:)
        integer :: st, nt, nx, i, saved(4)
        real(dp) :: saved_sub

        call gs_sl_curve_create(3, sampled_nodes(3), curve, st)
        nt = (flux_degree + 1)*np
        nx = 2*curve%n_nodes
        allocate (v0(nt, nt), v1(nt, nt), d0(nt, nx), d1(nt, nx), g0(nx), g1(nx))
        tau = [(1.0_dp + 0.3_dp*sin(0.7_dp*i), i=1, nt)]
        call gs_sl_plan_create(curve, flux_degree, plan, st)
        call gs_sl_assemble(curve, plan, v0, tau, d0, g0)
        saved = [radial_n, radial_q, near_order, far_order]
        saved_sub = sub_scale
        radial_n = 40; radial_q = 7; near_order = 12; far_order = 10; sub_scale = 3.5_dp
        call gs_sl_plan_create(curve, flux_degree, plan, st)
        call gs_sl_assemble(curve, plan, v1, tau, d1, g1)
        radial_n = saved(1); radial_q = saved(2); near_order = saved(3)
        far_order = saved(4); sub_scale = saved_sub
        call record(maxval(abs(v0 - v1)) <= 1.0e-10_dp*maxval(abs(v1)) .and. &
            maxval(abs(d0 - d1)) <= 1.0e-9_dp*maxval(abs(d1)) .and. &
            maxval(abs(g0 - g1)) <= 1.0e-9_dp*maxval(abs(g1)), &
            "default quadrature agrees with refined rules", &
            maxval(abs(v0 - v1))/maxval(abs(v1)))
    end subroutine test_quadrature_converged

    subroutine record(condition, description, err)
        logical, intent(in) :: condition
        character(*), intent(in) :: description
        real(dp), intent(in), optional :: err

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
        if (present(err)) write (*, "(a,es10.3)") "   relative error ", err
    end subroutine record

    pure function d_shape(theta) result(x)
        !! R = 3 + cos(theta + 0.4 sin theta), Z = 0.5 + 1.4 sin theta.
        real(dp), intent(in) :: theta
        real(dp) :: x(2)

        x = [3.0_dp + cos(theta + 0.4_dp*sin(theta)), 0.5_dp + 1.4_dp*sin(theta)]
    end function d_shape

    function sampled_nodes(k) result(nodes)
        !! Nodes of np P_k panels on the D shape, equispaced in theta.
        integer, intent(in) :: k
        real(dp) :: nodes(2, k*np)
        integer :: n

        do n = 1, k*np
            nodes(:, n) = d_shape(2.0_dp*pi*real(n - 1, dp)/real(k*np, dp))
        end do
    end function sampled_nodes

    subroutine test_against_p2_operator()
        type(gs_bem_curve_t) :: c2
        type(gs_bem_operators_t) :: ops
        type(gs_sl_curve_t) :: s2, s3
        type(gs_sl_plan_t) :: plan2, plan3
        real(dp) :: nodes2(2, 2*np), nodes3(2, 3*np), x(2), tau(2), sp
        real(dp), allocatable :: v2(:, :), v3(:, :)
        integer :: status, p, j

        nodes2 = sampled_nodes(2)
        call gs_bem_curve_from_nodes(nodes2, c2, status)
        call gs_bem_assemble(c2, flux_degree, ops, status, hypersingular=.false.)
        call record(status == 0, "reference P2 operator assembles")
        call gs_sl_curve_create(2, nodes2, s2, status)
        call gs_sl_plan_create(s2, flux_degree, plan2, status)
        allocate (v2(ops%ny, ops%ny), v3(ops%ny, ops%ny))
        call gs_sl_assemble(s2, plan2, v2)
        call record(maxval(abs(v2 - ops%v)) <= 1.0e-13_dp*maxval(abs(ops%v)), &
            "P2 single layer equals gs_bem_assemble", &
            maxval(abs(v2 - ops%v))/maxval(abs(ops%v)))

        ! P3 nodes on the P2 panels: the same curve and the same V.
        do p = 1, np
            do j = 0, 2
                call gs_bem_panel_point(c2, p, real(j, dp)/3.0_dp, x, tau, sp)
                nodes3(:, 3*(p - 1) + j + 1) = x
            end do
        end do
        call gs_sl_curve_create(3, nodes3, s3, status)
        call record(status == 0, "P3 curve is created")
        call gs_sl_plan_create(s3, flux_degree, plan3, status)
        call gs_sl_assemble(s3, plan3, v3)
        call record(maxval(abs(v3 - ops%v)) <= 1.0e-12_dp*maxval(abs(ops%v)), &
            "P3 panels of a P2 curve give the same single layer", &
            maxval(abs(v3 - ops%v))/maxval(abs(ops%v)))
    end subroutine test_against_p2_operator

    subroutine test_shape_derivatives(k)
        integer, intent(in) :: k
        type(gs_sl_curve_t) :: curve, moved
        type(gs_sl_plan_t) :: plan
        real(dp), allocatable :: nodes(:, :), v(:, :), vp(:, :), vm(:, :)
        real(dp), allocatable :: dvtau(:, :), grad(:), hess(:, :), tau(:)
        real(dp), allocatable :: dvp(:, :), gp(:), hp(:, :), gm(:), dir(:)
        real(dp), allocatable :: fd_vtau(:), fd_grad(:)
        real(dp) :: h, fd_energy, err
        integer :: ny, nx, status, i
        character(len=8) :: label

        write (label, "(a,i0)") "P", k
        nodes = sampled_nodes(k)
        call gs_sl_curve_create(k, nodes, curve, status)
        call gs_sl_plan_create(curve, flux_degree, plan, status)
        ny = (flux_degree + 1)*np
        nx = 2*curve%n_nodes
        allocate (v(ny, ny), vp(ny, ny), vm(ny, ny), dvtau(ny, nx), grad(nx), &
            hess(nx, nx), tau(ny), dvp(ny, nx), gp(nx), hp(nx, nx), gm(nx), &
            dir(nx), fd_vtau(ny), fd_grad(nx))
        do i = 1, ny
            tau(i) = 1.0_dp + 0.3_dp*sin(1.7_dp*i)
        end do
        do i = 1, nx
            dir(i) = 0.05_dp*cos(2.3_dp*i + 0.4_dp)
        end do
        call gs_sl_assemble(curve, plan, v, tau, dvtau, grad, hess)
        ! The value and gradient leaves give what the Hessian leaf gives.
        call gs_sl_assemble(curve, plan, vp)
        call gs_sl_assemble(curve, plan, vm, tau, dvp, gp)
        call record(maxval(abs(vp - v)) <= 1.0e-13_dp*maxval(abs(v)) .and. &
            maxval(abs(vm - v)) <= 1.0e-13_dp*maxval(abs(v)) .and. &
            maxval(abs(dvp - dvtau)) <= 1.0e-12_dp*maxval(abs(dvtau)) .and. &
            maxval(abs(gp - grad)) <= 1.0e-12_dp*maxval(abs(grad)), &
            trim(label)//" value and gradient leaves agree with the Hessian leaf")
        call record(maxval(abs(hess - transpose(hess))) <= &
            1.0e-12_dp*maxval(abs(hess)), trim(label)//" shape Hessian is symmetric")

        h = 1.0e-5_dp
        moved = curve
        moved%nodes = curve%nodes + h*reshape(dir, shape(curve%nodes))
        call gs_sl_assemble(moved, plan, vp, tau, dvp, gp, hp)
        moved%nodes = curve%nodes - h*reshape(dir, shape(curve%nodes))
        call gs_sl_assemble(moved, plan, vm, tau, dvp, gm, hp)
        fd_vtau = (matmul(vp, tau) - matmul(vm, tau))/(2*h)
        err = maxval(abs(fd_vtau - matmul(dvtau, dir)))/maxval(abs(fd_vtau))
        call record(err <= 1.0e-7_dp, trim(label)//" d(V tau)/dX matches differences", err)
        fd_energy = (0.5_dp*dot_product(tau, matmul(vp, tau)) &
            - 0.5_dp*dot_product(tau, matmul(vm, tau)))/(2*h)
        err = abs(fd_energy - dot_product(grad, dir))/abs(fd_energy)
        call record(err <= 1.0e-7_dp, trim(label)//" energy gradient matches differences", err)
        call record(maxval(abs(grad - 0.5_dp*matmul(tau, dvtau))) <= &
            1.0e-12_dp*maxval(abs(grad)), trim(label)//" gradient is tau'dV tau/2")
        fd_grad = (gp - gm)/(2*h)
        err = maxval(abs(fd_grad - matmul(hess, dir)))/maxval(abs(fd_grad))
        call record(err <= 1.0e-6_dp, trim(label)//" shape Hessian matches differences", &
            err)
    end subroutine test_shape_derivatives

    subroutine field(x, value, gradient, hessian)
        !! psi = R**2 Z + sin(Z) + R
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: value, gradient(2), hessian(3)

        value = x(1)**2*x(2) + sin(x(2)) + x(1)
        gradient = [2*x(1)*x(2) + 1.0_dp, x(1)**2 + cos(x(2))]
        hessian = [2*x(2), 2*x(1), -sin(x(2))]
    end subroutine field

    subroutine test_load_derivatives()
        type(gs_sl_curve_t) :: curve, moved
        real(dp), allocatable :: g(:), gp(:), gm(:), dg(:, :), grad(:), hess(:, :)
        real(dp), allocatable :: w(:), dir(:), d2(:, :), g2p(:), g2m(:), h2(:, :)
        real(dp) :: h, err
        integer :: ny, nx, status, i

        call gs_sl_curve_create(3, sampled_nodes(3), curve, status)
        ny = (flux_degree + 1)*np
        nx = 2*curve%n_nodes
        allocate (g(ny), gp(ny), gm(ny), dg(ny, nx), grad(nx), hess(nx, nx), &
            w(ny), dir(nx), d2(ny, nx), g2p(nx), g2m(nx), h2(nx, nx))
        do i = 1, ny
            w(i) = 0.5_dp + cos(0.9_dp*i)
        end do
        do i = 1, nx
            dir(i) = 0.05_dp*sin(1.3_dp*i + 0.2_dp)
        end do
        call gs_sl_load(curve, flux_degree, field, g, w, dg, grad, hess)
        h = 1.0e-5_dp
        moved = curve
        moved%nodes = curve%nodes + h*reshape(dir, shape(curve%nodes))
        call gs_sl_load(moved, flux_degree, field, gp, w, d2, g2p, h2)
        moved%nodes = curve%nodes - h*reshape(dir, shape(curve%nodes))
        call gs_sl_load(moved, flux_degree, field, gm, w, d2, g2m, h2)
        err = maxval(abs((gp - gm)/(2*h) - matmul(dg, dir)))/maxval(abs(gp - gm)/(2*h))
        call record(err <= 1.0e-8_dp, "load dg/dX matches differences", err)
        call record(maxval(abs(grad - matmul(w, dg))) <= 1.0e-12_dp*maxval(abs(grad)), &
            "load gradient is w'dg")
        err = maxval(abs((g2p - g2m)/(2*h) - matmul(hess, dir))) &
            /maxval(abs(g2p - g2m)/(2*h))
        call record(err <= 1.0e-7_dp, "load Hessian matches differences", err)
    end subroutine test_load_derivatives

end program test_gs_single_layer_shape_2d
