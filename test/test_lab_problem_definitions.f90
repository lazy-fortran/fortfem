program test_lab_problem_definitions
    use check, only: check_condition, check_summary
    use fortfem_generated_lab_scalar, only: lab_scalar
    use fortfem_generated_lab_vector, only: lab_vector
    use fortfem_generated_lab_geometry, only: lab_geometry
    use fortfem_generated_lab_transport, only: lab_transport
    use fortfem_generated_lab_guenter, only: lab_guenter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, parameter :: qp = selected_real_kind(30)
    real(dp) :: q(3), p(3), v(8), vp(8), vm(8), g(12), t(9)
    real(dp) :: step, lap, point(3), fluxp(3), fluxm(3), div, ratio
    real(dp) :: vector(9), vecp(9), vecm(9), dcurl(3, 3), curlcurl(3)
    integer :: i, j, k
    step = 2.0e-4_dp
    do k = 1, 4
        q = [0.1_dp*k, -0.07_dp*k, 0.19_dp*k]
        call lab_geometry(q(1), q(2), q(3), g)
        p = [g(1)*cos(q(3)), g(1)*sin(q(3)), g(2)]
        call lab_vector(p(1), p(2), p(3), vector)
        do i = 1, 3
            point = p
            point(i) = point(i) + step
            call lab_vector(point(1), point(2), point(3), vecp)
            point(i) = point(i) - 2*step
            call lab_vector(point(1), point(2), point(3), vecm)
            dcurl(:, i) = (vecp(4:6) - vecm(4:6))/(2*step)
        end do
        curlcurl = [dcurl(3, 2) - dcurl(2, 3), dcurl(1, 3) - dcurl(3, 1), &
                    dcurl(2, 1) - dcurl(1, 2)]
     call check_condition(maxval(abs(curlcurl + vector(:3) - vector(7:9))) < 3.e-8_dp, &
                             'L2 exact curl-curl plus mass agrees with forcing')
        call lab_scalar(p(1), p(2), p(3), v)
        lap = 0
        do i = 1, 3
            point = p
            point(i) = p(i) + step
            call lab_scalar(point(1), point(2), point(3), vp)
            point(i) = p(i) - step
            call lab_scalar(point(1), point(2), point(3), vm)
            ! First differences of independently sampled exact gradients.
            lap = lap + (vp(i + 1) - vm(i + 1))/(2*step)
        end do
        call check_condition(abs(lap + v(5)) < 3.0e-8_dp, &
                             'L1/L5 exact gradient divergence agrees with forcing')
        call lab_transport(g(1), q(3), g(2), t)
        call check_quad_source(p, t(5))
        call check_condition(abs(t(1) - (1 - sum(q(:2)**2))) < 2.e-14_dp, &
                          'L3 inverse geometry reproduces flux surfaces including axis')
        call check_condition(abs(dot_product(t(2:4), t(6:8))) < 2.e-13_dp, &
                             'L3 field direction is tangent to exact temperature')
        do j = 0, 3
            ratio = 10.0_dp**(2*j)
            div = 0
            do i = 1, 3
                point = p
                point(i) = p(i) + step
                call transport_flux(point, ratio, fluxp)
                point(i) = p(i) - step
                call transport_flux(point, ratio, fluxm)
                div = div + (fluxp(i) - fluxm(i))/(2*step)
            end do
            call check_condition(abs(div + t(5)) < 2.e-5_dp, &
                                 'L3 Cartesian flux divergence agrees with forcing')
        end do
    end do
    call lab_geometry(0.0_dp, 0.0_dp, 0.4_dp, g)
    call lab_scalar(g(1)*cos(0.4_dp), g(1)*sin(0.4_dp), g(2), v)
    call check_condition(sum(v(2:4)**2) > 0.01_dp, &
                         'L5 nonzero axis gradient')
    call lab_guenter(0.13_dp, -0.21_dp, v)
    lap = 0
    do i = 1, 2
        point = [0.13_dp, -0.21_dp, 0.0_dp]
        point(i) = point(i) + step
        call lab_guenter(point(1), point(2), vp)
        point(i) = point(i) - 2*step
        call lab_guenter(point(1), point(2), vm)
        lap = lap + (vp(i + 1) - vm(i + 1))/(2*step)
    end do
    call check_condition(abs(lap + v(4)) < 1.e-6_dp, &
                         'L3a Guenter cosine closed-line forcing')
    call check_condition(abs(dot_product(v(2:3), v(5:6))) < 1.e-14_dp, &
                         'L3a parallel exact derivative vanishes')
    call check_summary('laboratory problem definitions')
contains
    subroutine check_quad_source(point, forcing)
        real(dp), intent(in) :: point(3), forcing
        real(qp) :: x(3), shift(3), lap, h, center, values(4)
        integer :: i, j
        integer, parameter :: offsets(4) = [-2, -1, 1, 2]
        x = real(point, qp)
        h = 1.e-5_qp
        center = reference_temperature(x)
        lap = 0
        do i = 1, 3
            do j = 1, 4
                shift = x
                shift(i) = shift(i) + offsets(j)*h
                values(j) = reference_temperature(shift)
            end do
            lap = lap + (-values(1) + 16*values(2) - 30*center + &
                         16*values(3) - values(4))/(12*h*h)
        end do
        call check_condition(abs(real(-lap, dp) - forcing) < &
                             256*epsilon(1._dp)*max(1._dp, abs(forcing)), &
                   'L3 independent quad-precision FD forcing agrees at double roundoff')
    end subroutine
    pure function reference_temperature(x) result(value)
        real(qp), intent(in) :: x(3)
        real(qp) :: value, r, phi, u, v
        r = sqrt(x(1)*x(1) + x(2)*x(2))
        phi = atan2(x(2), x(1))
        u = (r - 3._qp - cos(5*phi)/10)/(1._qp + cos(5*phi)/10)
        v = (x(3) - sin(5*phi)/10)*10/13
        value = 1 - u*u - v*v
    end function
    subroutine transport_flux(x, ratio, flux)
        real(dp), intent(in) :: x(3), ratio
        real(dp), intent(out) :: flux(3)
        real(dp) :: data(9), r, phi, b(3), grad(3)
        r = sqrt(sum(x(:2)**2))
        phi = atan2(x(2), x(1))
        call lab_transport(r, phi, x(3), data)
        grad = [data(2)*cos(phi) - data(3)*sin(phi), &
                data(2)*sin(phi) + data(3)*cos(phi), data(4)]
        b = [data(6)*cos(phi) - data(7)*sin(phi), &
             data(6)*sin(phi) + data(7)*cos(phi), data(8)]
        flux = grad + (ratio - 1)*b*dot_product(b, grad)
    end subroutine
end program
