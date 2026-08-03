module fortfem_cut_polygon_fifth_moments_2d
    !! Exact degree-five raw moments of a fixed-topology polygon.
    !!
    !! This neutral cut-cell primitive is intended for a polygon produced by a
    !! level-set/XIGA clip.  Green-theorem edge polynomials provide exact value,
    !! tangent, and transpose products while the polygon topology is fixed.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    real(dp), parameter :: topology_tolerance = 64.0_dp*epsilon(1.0_dp)

    public :: evaluate_cut_polygon_fifth_moments_2d
    public :: evaluate_cut_polygon_fifth_moments_2d_jvp
    public :: evaluate_cut_polygon_fifth_moments_2d_vjp

contains

    subroutine evaluate_cut_polygon_fifth_moments_2d(polygon, moment, status)
        real(dp), intent(in) :: polygon(:, :)
        real(dp), intent(out) :: moment(2, 2, 2, 2, 2)
        type(fortsparse_status_t), intent(out) :: status
        integer :: edge, next_edge, first, second, third, fourth, fifth
        integer :: x_degree, y_degree
        real(dp) :: value

        moment = 0.0_dp
        call validate_polygon(polygon, status)
        if (status%code /= FORTSPARSE_OK) return
        do edge = 1, size(polygon, 2)
            next_edge = 1 + mod(edge, size(polygon, 2))
            do first = 1, 2
                do second = 1, 2
                    do third = 1, 2
                        do fourth = 1, 2
                            do fifth = 1, 2
                                x_degree = count([first, second, third, fourth, fifth] == 1)
                                y_degree = 5 - x_degree
                                value = edge_monomial_moment( &
                                    polygon(:, edge), polygon(:, next_edge), x_degree, y_degree)
                                moment(first, second, third, fourth, fifth) = &
                                    moment(first, second, third, fourth, fifth) + value
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_cut_polygon_fifth_moments_2d

    subroutine evaluate_cut_polygon_fifth_moments_2d_jvp( &
            polygon, polygon_dot, moment_dot, status)
        real(dp), intent(in) :: polygon(:, :), polygon_dot(:, :)
        real(dp), intent(out) :: moment_dot(2, 2, 2, 2, 2)
        type(fortsparse_status_t), intent(out) :: status
        integer :: edge, next_edge, first, second, third, fourth, fifth
        integer :: x_degree, y_degree
        real(dp) :: value, value_dot

        moment_dot = 0.0_dp
        call validate_polygon(polygon, status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(polygon_dot) /= shape(polygon)) .or. &
            any(.not. ieee_is_finite(polygon_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "fifth cut moment JVP received incompatible vertices")
            return
        end if
        do edge = 1, size(polygon, 2)
            next_edge = 1 + mod(edge, size(polygon, 2))
            do first = 1, 2
                do second = 1, 2
                    do third = 1, 2
                        do fourth = 1, 2
                            do fifth = 1, 2
                                x_degree = count([first, second, third, fourth, fifth] == 1)
                                y_degree = 5 - x_degree
                                call edge_monomial_moment_jvp( &
                                    polygon(:, edge), polygon(:, next_edge), polygon_dot(:, edge), &
                                    polygon_dot(:, next_edge), x_degree, y_degree, value, value_dot)
                                moment_dot(first, second, third, fourth, fifth) = &
                                    moment_dot(first, second, third, fourth, fifth) + value_dot
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_cut_polygon_fifth_moments_2d_jvp

    subroutine evaluate_cut_polygon_fifth_moments_2d_vjp( &
            polygon, moment_bar, polygon_bar, status)
        real(dp), intent(in) :: polygon(:, :), moment_bar(2, 2, 2, 2, 2)
        real(dp), intent(out) :: polygon_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: edge, next_edge, first, second, third, fourth, fifth
        integer :: x_degree, y_degree
        real(dp) :: value, value_bar

        polygon_bar = 0.0_dp
        call validate_polygon(polygon, status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(polygon_bar) /= shape(polygon)) .or. &
            any(.not. ieee_is_finite(moment_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "fifth cut moment VJP received incompatible cotangents")
            return
        end if
        do edge = 1, size(polygon, 2)
            next_edge = 1 + mod(edge, size(polygon, 2))
            do first = 1, 2
                do second = 1, 2
                    do third = 1, 2
                        do fourth = 1, 2
                            do fifth = 1, 2
                                x_degree = count([first, second, third, fourth, fifth] == 1)
                                y_degree = 5 - x_degree
                                value = moment_bar(first, second, third, fourth, fifth)
                                call edge_monomial_moment_vjp( &
                                    polygon(:, edge), polygon(:, next_edge), x_degree, y_degree, value, &
                                    polygon_bar(:, edge), polygon_bar(:, next_edge))
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_cut_polygon_fifth_moments_2d_vjp

    subroutine validate_polygon(polygon, status)
        real(dp), intent(in) :: polygon(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: signed_area
        integer :: edge, next_edge

        call status_set(status, FORTSPARSE_INVALID_MATRIX, "fifth cut moment polygon is invalid")
        if (size(polygon, 1) /= 2 .or. size(polygon, 2) < 3 .or. size(polygon, 2) > 6 .or. &
            any(.not. ieee_is_finite(polygon))) return
        signed_area = 0.0_dp
        do edge = 1, size(polygon, 2)
            next_edge = 1 + mod(edge, size(polygon, 2))
            signed_area = signed_area + polygon(1, edge)*polygon(2, next_edge) - &
                polygon(2, edge)*polygon(1, next_edge)
        end do
        if (abs(signed_area) <= topology_tolerance) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_polygon

    pure function edge_monomial_moment(point, next_point, x_degree, y_degree) result(value)
        real(dp), intent(in) :: point(2), next_point(2)
        integer, intent(in) :: x_degree, y_degree
        real(dp) :: value, ax, ay, dx, dy
        integer :: x_power, y_power

        ax = point(1); ay = point(2); dx = next_point(1) - ax; dy = next_point(2) - ay
        value = 0.0_dp
        do x_power = 0, x_degree + 1
            do y_power = 0, y_degree
                value = value + dy*binomial(x_degree + 1, x_power)*binomial(y_degree, y_power)* &
                    ax**(x_degree + 1 - x_power)*dx**x_power*ay**(y_degree - y_power)*dy**y_power / &
                    real((x_degree + 1)*(x_power + y_power + 1), dp)
            end do
        end do
    end function edge_monomial_moment

    pure subroutine edge_monomial_moment_jvp(point, next_point, point_dot, next_point_dot, &
            x_degree, y_degree, value, value_dot)
        real(dp), intent(in) :: point(2), next_point(2), point_dot(2), next_point_dot(2)
        integer, intent(in) :: x_degree, y_degree
        real(dp), intent(out) :: value, value_dot
        real(dp) :: ax, ay, dx, dy, ax_dot, ay_dot, dx_dot, dy_dot
        real(dp) :: fa, fb, fc, fd, fa_dot, fb_dot, fc_dot, fd_dot, product, product_dot
        real(dp) :: coefficient
        integer :: x_power, y_power

        ax = point(1); ay = point(2); dx = next_point(1) - ax; dy = next_point(2) - ay
        ax_dot = point_dot(1); ay_dot = point_dot(2)
        dx_dot = next_point_dot(1) - ax_dot; dy_dot = next_point_dot(2) - ay_dot
        value = 0.0_dp; value_dot = 0.0_dp
        do x_power = 0, x_degree + 1
            do y_power = 0, y_degree
                call power_with_derivative(ax, ax_dot, x_degree + 1 - x_power, fa, fa_dot)
                call power_with_derivative(dx, dx_dot, x_power, fb, fb_dot)
                call power_with_derivative(ay, ay_dot, y_degree - y_power, fc, fc_dot)
                call power_with_derivative(dy, dy_dot, y_power, fd, fd_dot)
                product = fa*fb*fc*fd
                product_dot = fa_dot*fb*fc*fd + fa*fb_dot*fc*fd + fa*fb*fc_dot*fd + fa*fb*fc*fd_dot
                coefficient = binomial(x_degree + 1, x_power)*binomial(y_degree, y_power) / &
                    real((x_degree + 1)*(x_power + y_power + 1), dp)
                value = value + coefficient*dy*product
                value_dot = value_dot + coefficient*(dy_dot*product + dy*product_dot)
            end do
        end do
    end subroutine edge_monomial_moment_jvp

    pure subroutine edge_monomial_moment_vjp(point, next_point, x_degree, y_degree, value_bar, &
            point_bar, next_point_bar)
        real(dp), intent(in) :: point(2), next_point(2), value_bar
        integer, intent(in) :: x_degree, y_degree
        real(dp), intent(inout) :: point_bar(2), next_point_bar(2)
        real(dp) :: ax, ay, dx, dy, ax_bar, ay_bar, dx_bar, dy_bar
        real(dp) :: fa, fb, fc, fd, fa_bar, fb_bar, fc_bar, fd_bar, coefficient, seed
        integer :: x_power, y_power, exponent

        ax = point(1); ay = point(2); dx = next_point(1) - ax; dy = next_point(2) - ay
        ax_bar = 0.0_dp; ay_bar = 0.0_dp; dx_bar = 0.0_dp; dy_bar = 0.0_dp
        do x_power = 0, x_degree + 1
            do y_power = 0, y_degree
                fa = ax**(x_degree + 1 - x_power); fb = dx**x_power
                fc = ay**(y_degree - y_power); fd = dy**y_power
                coefficient = binomial(x_degree + 1, x_power)*binomial(y_degree, y_power) / &
                    real((x_degree + 1)*(x_power + y_power + 1), dp)
                seed = value_bar*coefficient
                fa_bar = seed*dy*fb*fc*fd
                fb_bar = seed*dy*fa*fc*fd
                fc_bar = seed*dy*fa*fb*fd
                fd_bar = seed*dy*fa*fb*fc
                dy_bar = dy_bar + seed*fa*fb*fc*fd
                if (y_power > 0) dy_bar = dy_bar + fd_bar*real(y_power, dp)*dy**(y_power - 1)
                exponent = x_degree + 1 - x_power
                if (exponent > 0) ax_bar = ax_bar + fa_bar*real(exponent, dp)*ax**(exponent - 1)
                if (x_power > 0) dx_bar = dx_bar + fb_bar*real(x_power, dp)*dx**(x_power - 1)
                exponent = y_degree - y_power
                if (exponent > 0) ay_bar = ay_bar + fc_bar*real(exponent, dp)*ay**(exponent - 1)
            end do
        end do
        point_bar(1) = point_bar(1) + ax_bar - dx_bar
        next_point_bar(1) = next_point_bar(1) + dx_bar
        point_bar(2) = point_bar(2) + ay_bar - dy_bar
        next_point_bar(2) = next_point_bar(2) + dy_bar
    end subroutine edge_monomial_moment_vjp

    pure subroutine power_with_derivative(base, base_dot, exponent, value, value_dot)
        real(dp), intent(in) :: base, base_dot
        integer, intent(in) :: exponent
        real(dp), intent(out) :: value, value_dot

        value = base**exponent
        value_dot = 0.0_dp
        if (exponent > 0) value_dot = real(exponent, dp)*base**(exponent - 1)*base_dot
    end subroutine power_with_derivative

    pure function binomial(number, choose) result(value)
        integer, intent(in) :: number, choose
        real(dp) :: value
        integer :: factor

        value = 0.0_dp
        if (choose < 0 .or. choose > number) return
        value = 1.0_dp
        do factor = 1, choose
            value = value*real(number - choose + factor, dp)/real(factor, dp)
        end do
    end function binomial

end module fortfem_cut_polygon_fifth_moments_2d
