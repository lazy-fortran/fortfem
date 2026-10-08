program test_tetra_modal_radial_products
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    use fortfem_generated_tetra_modal_radial_products, only: &
        evaluate_tetra_modal_radial_products
    use fortfem_generated_tetra_modal_radial_products_jvp, only: &
        evaluate_tetra_modal_radial_products_jvp
    use fortfem_generated_tetra_modal_component_curls_jvp, only: &
        evaluate_tetra_modal_component_curls_jvp
    implicit none
    real(dp), parameter :: points(3, 4) = reshape([ &
        0._dp, 0._dp, 0._dp, 1._dp, 0._dp, 0._dp, &
        .17_dp, .21_dp, .13_dp, -.7_dp, .2_dp, .9_dp], [3, 4])
    real(dp), parameter :: directions(3, 3) = reshape([ &
        1._dp, 0._dp, 0._dp, 0._dp, 0._dp, -1._dp, &
        .17_dp, -.11_dp, .09_dp], [3, 3])
    real(dp), parameter :: step = .002_dp
    real(dp) :: p(3), d(3), phi, gradient(3), hessian(3, 3), gradient_dot(3)
    real(dp) :: values(3), values_dot(3), divergence, divergence_dot, phi_dot
    real(dp) :: curl_dot(3, 3), expected_values_dot(3), expected_curl_dot(3, 3)
    real(dp) :: expected_divergence_dot
    integer :: point_id, direction_id

    do point_id = 1, size(points, 2)
        p = points(:, point_id)
        call polynomial_jet(p, phi, gradient, hessian)
        call evaluate_tetra_modal_radial_products( &
            p(1), p(2), p(3), phi, gradient(1), gradient(2), gradient(3), &
            values, divergence)
        call check_condition(maxval(abs(values - radial_field(p))) < 2e-13_dp, &
            'Radial value reproduces an independent quartic scalar field')
        call check_condition(abs(divergence - divergence_fd(p)) < 2e-10_dp, &
            'Radial divergence matches independent spatial finite differences')
        do direction_id = 1, size(directions, 2)
            d = directions(:, direction_id)
            phi_dot = dot_product(gradient, d)
            gradient_dot = matmul(hessian, d)
            call evaluate_tetra_modal_radial_products_jvp( &
                p(1), p(2), p(3), phi, gradient(1), gradient(2), gradient(3), &
                d(1), d(2), d(3), phi_dot, &
                gradient_dot(1), gradient_dot(2), gradient_dot(3), &
                values_dot, divergence_dot)
            expected_values_dot = (-radial_field(p + 2*step*d) + &
                8*radial_field(p + step*d) - 8*radial_field(p - step*d) + &
                radial_field(p - 2*step*d))/(12*step)
            expected_divergence_dot = (-divergence_fd(p + 2*step*d) + &
                8*divergence_fd(p + step*d) - 8*divergence_fd(p - step*d) + &
                divergence_fd(p - 2*step*d))/(12*step)
            call check_condition(maxval(abs(values_dot - expected_values_dot)) &
                < 3e-8_dp, 'Radial value JVP matches independent field differences')
            call check_condition(abs(divergence_dot - expected_divergence_dot) &
                < 5e-8_dp, 'Radial divergence JVP matches spatial differences')
            call evaluate_tetra_modal_component_curls_jvp( &
                gradient_dot(1), gradient_dot(2), gradient_dot(3), curl_dot)
            expected_curl_dot = (-component_curls_fd(p + 2*step*d) + &
                8*component_curls_fd(p + step*d) - &
                8*component_curls_fd(p - step*d) + &
                component_curls_fd(p - 2*step*d))/(12*step)
            call check_condition(maxval(abs(curl_dot - expected_curl_dot)) &
                < 5e-8_dp, 'Component curl JVP matches spatial differences')
        end do
    end do
    call check_summary('Tetrahedral generated modal differential products')
contains
    pure real(dp) function polynomial(p) result(value)
        real(dp), intent(in) :: p(3)
        value = 1 + .7_dp*p(1) - .5_dp*p(2) + .2_dp*p(3) + &
            p(1)*p(2) - 2*p(1)*p(3) + 3*p(2)*p(3) + &
            p(1)**2*p(2)**2 + .4_dp*p(3)**3
    end function polynomial

    pure subroutine polynomial_jet(p, value, gradient, hessian)
        real(dp), intent(in) :: p(3)
        real(dp), intent(out) :: value, gradient(3), hessian(3, 3)
        value = polynomial(p)
        gradient = [.7_dp + p(2) - 2*p(3) + 2*p(1)*p(2)**2, &
            -.5_dp + p(1) + 3*p(3) + 2*p(1)**2*p(2), &
            .2_dp - 2*p(1) + 3*p(2) + 1.2_dp*p(3)**2]
        hessian(1, :) = [2*p(2)**2, 1 + 4*p(1)*p(2), -2._dp]
        hessian(2, :) = [1 + 4*p(1)*p(2), 2*p(1)**2, 3._dp]
        hessian(3, :) = [-2._dp, 3._dp, 2.4_dp*p(3)]
    end subroutine polynomial_jet

    pure function radial_field(p) result(value)
        real(dp), intent(in) :: p(3)
        real(dp) :: value(3)
        value = p*polynomial(p)
    end function radial_field

    pure function component_field(p, component) result(value)
        real(dp), intent(in) :: p(3)
        integer, intent(in) :: component
        real(dp) :: value(3)
        value = 0
        value(component) = polynomial(p)
    end function component_field

    pure function spatial_derivative(p, axis, component) result(value)
        real(dp), intent(in) :: p(3)
        integer, intent(in) :: axis, component
        real(dp), parameter :: h = .0005_dp
        real(dp) :: value(3), delta(3)
        delta = 0
        delta(axis) = h
        if (component == 0) then
            value = (-radial_field(p + 2*delta) + 8*radial_field(p + delta) - &
                8*radial_field(p - delta) + radial_field(p - 2*delta))/(12*h)
        else
            value = (-component_field(p + 2*delta, component) + &
                8*component_field(p + delta, component) - &
                8*component_field(p - delta, component) + &
                component_field(p - 2*delta, component))/(12*h)
        end if
    end function spatial_derivative

    pure real(dp) function divergence_fd(p) result(value)
        real(dp), intent(in) :: p(3)
        real(dp) :: derivative(3)
        integer :: axis
        value = 0
        do axis = 1, 3
            derivative = spatial_derivative(p, axis, 0)
            value = value + derivative(axis)
        end do
    end function divergence_fd

    pure function component_curls_fd(p) result(value)
        real(dp), intent(in) :: p(3)
        real(dp) :: value(3, 3), dx(3), dy(3), dz(3)
        integer :: component
        do component = 1, 3
            dx = spatial_derivative(p, 1, component)
            dy = spatial_derivative(p, 2, component)
            dz = spatial_derivative(p, 3, component)
            value(:, component) = [dy(3) - dz(2), dz(1) - dx(3), dx(2) - dy(1)]
        end do
    end function component_curls_fd
end program test_tetra_modal_radial_products
