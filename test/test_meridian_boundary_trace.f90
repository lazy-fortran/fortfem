program test_meridian_boundary_trace
    !! Boundary trace of a physical field in V1 (meridian_boundary_trace).
    !!  1. Straight edges: a polynomial A_pol of degree k - 1 and R A_phi of
    !!     degree k are reproduced exactly on the boundary.
    !!  2. Curved edges (geometry order 2): the tangential trace error of a
    !!     smooth field falls at least like h^k under refinement.
    !!  3. Interior dofs are zero.
    !!  4. Second-kind edges reproduce polynomial A_pol of degree k.
    !! The reference tangent dx/dtau is a fourth-order central difference of the mapped
    !! edge point, independent of the routine's Jacobian.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init, &
        meridian_point_basis, meridian_boundary_trace
    use meridian_test_meshes, only: make_warped_rectangle_mesh
    implicit none

    logical :: all_passed
    integer :: k, kind2 = 1, pdeg = 1
    real(dp) :: e1, e2, order

    all_passed = .true.
    do k = 1, 3
        pdeg = k
        call run(k, 1, 4, polynomial, e1)
        print '(a,i0,a,es12.4)', "order ", k, " polynomial trace error ", e1
        call record(e1 < 1.0e-11_dp, "polynomial trace exact, order " &
            //char(48 + k))
    end do
    do k = 1, 2
        ! Second-kind N2_k contains every degree-k field; reuse the
        ! order-(k + 1) polynomial, whose R A_phi has degree k + 1.
        kind2 = 2
        pdeg = k + 1
        call run(k, 1, 4, polynomial, e1)
        kind2 = 1
        print '(a,i0,a,es12.4)', "second kind ", k, " polynomial error ", e1
        call record(e1 < 1.0e-11_dp, "second-kind polynomial trace exact, " &
            //"order "//char(48 + k))
    end do
    do k = 2, 3
        call run(k, 2, 4, smooth, e1)
        call run(k, 2, 8, smooth, e2)
        order = log(e1/e2)/log(2.0_dp)
        print '(a,i0,a,2es12.4,a,f7.3)', "order ", k, " trace errors ", e1, &
            e2, " observed ", order
        call record(order > k - 0.3_dp, "curved trace converges, order " &
            //char(48 + k))
    end do
    call check_summary("Meridian boundary trace")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    subroutine polynomial(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        real(dp) :: r, z
        r = x(1); z = x(2)
        select case (pdeg)
        case (1)
            v = [cmplx(0.7_dp, -0.2_dp, dp), cmplx(1.3_dp + 0.4_dp*r &
                - 0.3_dp*z, 0.1_dp*z, dp)/r, cmplx(-0.5_dp, 0.9_dp, dp)]
        case (2)
            v = [cmplx(0.7_dp + 0.3_dp*r - 0.2_dp*z, 0.5_dp*z, dp), &
                cmplx(1.3_dp + 0.4_dp*r*z - 0.3_dp*z*z, r*r, dp)/r, &
                cmplx(-0.5_dp + 0.6_dp*z, 0.9_dp*r, dp)]
        case default
            v = [cmplx(0.7_dp + 0.3_dp*r*z - 0.2_dp*z*z, 0.5_dp*r*r, dp), &
                cmplx(1.3_dp + 0.4_dp*r*z*z - 0.3_dp*z**3, r**3, dp)/r, &
                cmplx(-0.5_dp + 0.6_dp*r*r, 0.9_dp*r*z, dp)]
        end select
    end subroutine polynomial

    subroutine smooth(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        v = [cmplx(sin(2*x(2)) + x(1)**2, cos(x(1)*x(2)), dp), &
            cmplx(exp(0.5_dp*x(2))*cos(x(1)), 0.0_dp, dp), &
            cmplx(cos(3*x(1) - x(2)), sin(x(1) + x(2)), dp)]
    end subroutine smooth

    subroutine run(k, g, n, fun, err)
        integer, intent(in) :: k, g, n
        procedure(polynomial) :: fun
        real(dp), intent(out) :: err
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        complex(dp), allocatable :: values(:)
        real(dp), allocatable :: phi(:, :)
        real(dp), parameter :: stencil(-2:2) = [1.0_dp/12, -2.0_dp/3, 0.0_dp, &
            2.0_dp/3, -1.0_dp/12]
        real(dp) :: tau, xi, eta, x(2), xp(2), t(2), h, scale
        complex(dp) :: v(3), ah(3), c
        integer :: status, tri, e, q, l, j0, nln, ke, s
        logical :: interior_zero

        call make_warped_rectangle_mesh(n, g, mm)
        call meridian_feec_build(mm, k, fs, status, edge_kind=kind2)
        call record(status == 0, "build")
        call meridian_point_basis_init(fs, pb, status)
        call meridian_boundary_trace(fs, fun, values, status)
        call record(status == 0, "trace computed")
        interior_zero = all(abs(pack(values(:fs%nn), &
            .not. fs%sp(1)%bnd(:fs%nn))) == 0.0_dp)
        interior_zero = interior_zero .and. all(abs(pack(values(fs%nn + 1:), &
            .not. fs%sp(0)%bnd)) == 0.0_dp)
        call record(interior_zero, "interior dofs zero")
        nln = fs%nln; ke = fs%ke
        allocate (phi(3, fs%sp(1)%nloc))
        h = 1.0e-3_dp
        err = 0.0_dp; scale = 0.0_dp
        do tri = 1, fs%mm%topology%n_triangles
            do e = 1, 3
                j0 = (e - 1)*ke
                if (.not. fs%sp(1)%bnd(fs%sp(1)%dof(j0 + 1, tri))) cycle
                do q = 1, 7
                    tau = (q - 0.5_dp)/7
                    t = 0
                    do s = -2, 2
                        if (s == 0) cycle
                        call edge_point(e, tau + s*h, xi, eta)
                        call meridian_point_basis(fs, pb, 1, tri, xi, eta, xp, &
                            phi, status)
                        t = t + stencil(s)*xp/h
                    end do
                    call edge_point(e, tau, xi, eta)
                    call meridian_point_basis(fs, pb, 1, tri, xi, eta, x, phi, &
                        status)
                    ah = 0
                    do l = 1, size(phi, 2)
                        c = values(fs%sp(1)%dof(l, tri))
                        if (l <= nln) c = c*fs%sp(1)%sgn(l, tri)
                        ah = ah + c*phi(:, l)
                    end do
                    call fun(x, v)
                    err = max(err, abs(dot_product([real(ah(1)), real(ah(3))], &
                        t) - dot_product([real(v(1)), real(v(3))], t)), &
                        abs(dot_product([aimag(ah(1)), aimag(ah(3))], t) &
                        - dot_product([aimag(v(1)), aimag(v(3))], t)), &
                        x(1)*abs(ah(2) - v(2)))
                    scale = max(scale, abs(v(1))*norm2(t), abs(v(3))*norm2(t))
                end do
            end do
        end do
        err = err/scale
    end subroutine run

    pure subroutine edge_point(e, tau, xi, eta)
        integer, intent(in) :: e
        real(dp), intent(in) :: tau
        real(dp), intent(out) :: xi, eta
        select case (e)
        case (1)
            xi = tau; eta = 0
        case (2)
            xi = 1 - tau; eta = tau
        case default
            xi = 0; eta = 1 - tau
        end select
    end subroutine edge_point
end program test_meridian_boundary_trace
