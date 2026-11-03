program test_meridian_second_gradient
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_lagrange_nodes, meridian_qpoint
    use fortsparse, only: csc_matvec
    implicit none
    type(mesh_2d_t) :: mesh
    type(meridian_mesh_t) :: mm
    type(meridian_feec_t) :: fs
    real(dp), allocatable :: nodes(:, :), phi(:, :)
    complex(dp), allocatable :: scalar(:), gradient(:)
    complex(dp) :: field(3), exact(3)
    real(dp) :: x(2), weight, error, scale
    integer :: k, degree, status, t, q, j
    logical :: passed

    passed = .true.
    call mesh%create_rectangular(4, 4, 1.6_dp, 2.6_dp, -0.3_dp, 0.7_dp)
    call meridian_mesh_straight(mesh, 1, mm, status)
    if (status /= 0) error stop 'mesh'
    do k = 2, 3
        call meridian_feec_build(mm, k, fs, status, quad_degree=4*k + 4, edge_kind=2)
        if (status /= 0) error stop 'complex'
        call meridian_lagrange_nodes(fs, nodes)
        allocate (phi(3, fs%sp(1)%nloc))
        do degree = 1, 2
            scalar = polynomial(nodes, degree)
            gradient = cmplx(csc_matvec(fs%grad%a, real(scalar, dp)), &
                csc_matvec(fs%grad%a, aimag(scalar)), dp)
            error = 0; scale = 0
            do t = 1, mm%topology%n_triangles
                do q = 1, fs%nq
                    call meridian_qpoint(fs, 1, t, q, x, weight, phi, status)
                    if (status /= 0) error stop 'basis'
                    field = 0
                    do j = 1, fs%sp(1)%nloc
                        field = field + phi(:, j)*fs%sp(1)%sgn(j, t) &
                            *gradient(fs%sp(1)%dof(j, t))
                    end do
                    exact = polynomial_gradient(x, degree)
                    error = max(error, maxval(abs(field - exact)))
                    scale = max(scale, maxval(abs(exact)))
                end do
            end do
            print '(a,2i3,es20.12)', 'edge/scalar degree, gradient relative error: ', &
                k, degree, error/scale
            call check_condition(error/scale < 1.0e-13_dp, &
                'Second-kind gradient reproduces analytic polynomial')
            passed = passed .and. error/scale < 1.0e-13_dp
        end do
        deallocate (phi)
    end do
    call check_summary('Second-kind analytic gradient')
    if (.not. passed) error stop 1
contains
    function polynomial(x, degree) result(f)
        real(dp), intent(in) :: x(:, :)
        integer, intent(in) :: degree
        complex(dp) :: f(size(x, 2))
        real(dp) :: d, re(size(x, 2)), im(size(x, 2))
        d = degree - 1
        re = 0.3_dp + 0.7_dp*x(1, :) - 0.4_dp*x(2, :) &
            + d*(0.5_dp*x(1, :)**2 - 0.2_dp*x(1, :)*x(2, :) + 0.6_dp*x(2, :)**2)
        im = -0.2_dp + 0.1_dp*x(1, :) + 0.9_dp*x(2, :) &
            + d*(-0.3_dp*x(1, :)**2 + 0.8_dp*x(1, :)*x(2, :) - 0.1_dp*x(2, :)**2)
        f = cmplx(re, im, dp)
    end function
    function polynomial_gradient(x, degree) result(g)
        real(dp), intent(in) :: x(2)
        integer, intent(in) :: degree
        complex(dp) :: g(3)
        real(dp) :: d
        d = degree - 1
        g(1) = cmplx(0.7_dp + d*(x(1) - 0.2_dp*x(2)), &
            0.1_dp + d*(-0.6_dp*x(1) + 0.8_dp*x(2)), dp)
        g(2) = 0
        g(3) = cmplx(-0.4_dp + d*(-0.2_dp*x(1) + 1.2_dp*x(2)), &
            0.9_dp + d*(0.8_dp*x(1) - 0.2_dp*x(2)), dp)
    end function
end program
