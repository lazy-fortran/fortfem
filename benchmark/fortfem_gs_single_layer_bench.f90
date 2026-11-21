program fortfem_gs_single_layer_bench
    !! Wall time of gs_sl_assemble at orders 0, 1, 2 on a DIII-D-sized
    !! D-shaped curve (np P3 panels, flux degree 2), with checksums of V,
    !! the gradient and the Hessian for before/after comparisons.
    use fortfem_kinds, only: dp
    use fortfem_gs_single_layer_shape_2d, only: gs_sl_curve_t, gs_sl_plan_t, &
        gs_sl_curve_create, gs_sl_plan_create, gs_sl_assemble, radial_n, radial_q, near_order, far_order, far_ratio, sub_scale
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: degree = 3, flux_degree = 2
    type(gs_sl_curve_t) :: curve
    type(gs_sl_plan_t) :: plan
    real(dp), allocatable :: nodes(:, :), v(:, :), tau(:), dvt(:, :), g(:), h(:, :)
    real(dp) :: t0, t1, theta
    character(len=16) :: arg
    integer :: np, n, nt, nx, st, reps, r, order
    np = 48
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg); read (arg, *) np
    end if
    reps = 3
    if (command_argument_count() > 2) then
        call get_command_argument(2, arg); read (arg, *) radial_n
        call get_command_argument(3, arg); read (arg, *) radial_q
    end if
    if (command_argument_count() > 6) then
        call get_command_argument(4, arg); read (arg, *) near_order
        call get_command_argument(5, arg); read (arg, *) far_order
        call get_command_argument(6, arg); read (arg, *) far_ratio
        call get_command_argument(7, arg); read (arg, *) sub_scale
    end if
    allocate (nodes(2, degree*np))
    do n = 1, degree*np
        theta = 2*pi*real(n - 1, dp)/real(degree*np, dp)
        nodes(:, n) = [1.66_dp + 0.6_dp*cos(theta + 0.45_dp*sin(theta)), &
            1.1_dp*sin(theta)]
    end do
    call gs_sl_curve_create(degree, nodes, curve, st)
    if (st /= 0) error stop 'curve'
    call gs_sl_plan_create(curve, flux_degree, plan, st)
    if (st /= 0) error stop 'plan'
    block
        integer(8) :: pts(0:3), a, b
        pts = 0
        do b = 1, np
            do a = 1, np
                select case (plan%kind(a, b))
                case (1); pts(1) = pts(1) + size(plan%self_xi)
                case (2); pts(2) = pts(2) + size(plan%nxt_xi)
                case (3); pts(3) = pts(3) + size(plan%prv_xi)
                case default
                    pts(0) = pts(0) + int(plan%sub(a, b)*plan%order(a, b), 8)**2
                end select
            end do
        end do
        print "(a,4i10,a,i10)", "points separated, self, next, prev:", pts, " total", sum(pts)
        print '(a,i6,a,f8.2)', 'separated pairs with sub > 1:', count_sub(), &
            ', mean sub', real(sum(plan%sub), dp)/max(1, count(plan%kind == 0))
    end block
    nt = (flux_degree + 1)*np
    nx = 2*degree*np
    allocate (v(nt, nt), tau(nt), dvt(nt, nx), g(nx), h(nx, nx))
    tau = [(1 + 0.3_dp*sin(0.7_dp*n), n = 1, nt)]
    do order = 0, 2
        call cpu_time(t0)
        do r = 1, reps
            select case (order)
            case (0)
                call gs_sl_assemble(curve, plan, v)
            case (1)
                call gs_sl_assemble(curve, plan, v, tau, dvt, g)
            case (2)
                call gs_sl_assemble(curve, plan, v, tau, dvt, g, h)
            end select
        end do
        call cpu_time(t1)
        print '(a,i1,a,f10.4,a)', 'order ', order, ': ', (t1 - t0)/reps, ' s CPU per assembly'
    end do
    print '(a,3es24.16)', 'checksum V, g, H:', sum(v*v), sum(g*g), sum(h*h)
contains
    integer function count_sub()
        count_sub = count(plan%kind == 0 .and. plan%sub > 1)
    end function
end program fortfem_gs_single_layer_bench
