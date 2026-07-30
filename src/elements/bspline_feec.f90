module fortfem_bspline_feec
    !! Tensor-product B-spline de Rham complexes on the parametric domain.
    !!
    !! The degree-reduced component spaces and commuting coefficient
    !! derivatives follow Buffa, Rivas, Sangalli, and Vázquez,
    !! SIAM J. Numer. Anal. 49 (2011), doi:10.1137/100786708.
    use fortfem_kinds, only: dp
    use fortfem_generated_nurbs_geometry_jvp, only: &
        generated_nurbs_geometry_quotient_jvp
    use fortfem_generated_nurbs_geometry_vjp, only: &
        generated_nurbs_geometry_quotient_vjp
    use fortfem_generated_nurbs_volume_jvp, only: &
        generated_nurbs_volume_quotient_jvp
    use fortfem_generated_nurbs_volume_vjp, only: &
        generated_nurbs_volume_quotient_vjp
    implicit none
    private

    public :: build_bspline_derivative_matrix
    public :: build_bspline_feec_2d_operators
    public :: build_bspline_feec_3d_operators
    public :: evaluate_bspline_basis
    public :: evaluate_nurbs_surface_geometry
    public :: evaluate_nurbs_surface_geometry_jvp
    public :: evaluate_nurbs_surface_geometry_vjp
    public :: evaluate_nurbs_volume_geometry
    public :: evaluate_nurbs_volume_geometry_jvp
    public :: evaluate_nurbs_volume_geometry_vjp
    public :: map_isogeometric_h1_gradient
    public :: map_isogeometric_hcurl
    public :: map_isogeometric_hdiv
    public :: map_isogeometric_l2

contains

    subroutine evaluate_bspline_basis( &
            knots, degree, coordinate, values, derivatives, status)
        real(dp), intent(in) :: knots(:), coordinate
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: values(:), derivatives(:)
        integer, intent(out) :: status

        real(dp), allocatable :: lower(:), work(:)
        real(dp) :: denominator, evaluation_point, left, right
        integer :: basis, basis_count, level

        status = 1
        if (.not. valid_knot_vector(knots, degree)) return
        basis_count = size(knots) - degree - 1
        left = knots(degree + 1)
        right = knots(basis_count + 1)
        if (coordinate < left .or. coordinate > right) return
        evaluation_point = coordinate
        if (coordinate == right) evaluation_point = nearest(right, -1.0_dp)
        allocate(work(size(knots) - 1))
        work = 0.0_dp
        do basis = 1, size(work)
            if (knots(basis) <= evaluation_point .and. &
                evaluation_point < knots(basis + 1)) work(basis) = 1.0_dp
        end do
        do level = 1, degree
            if (level == degree) then
                allocate(lower(basis_count + 1))
                lower = work(:basis_count + 1)
            end if
            do basis = 1, size(knots) - level - 1
                left = 0.0_dp
                denominator = knots(basis + level) - knots(basis)
                if (denominator > 0.0_dp) left = &
                    (evaluation_point - knots(basis))*work(basis)/denominator
                right = 0.0_dp
                denominator = knots(basis + level + 1) - knots(basis + 1)
                if (denominator > 0.0_dp) right = &
                    (knots(basis + level + 1) - evaluation_point)* &
                    work(basis + 1)/denominator
                work(basis) = left + right
            end do
        end do
        allocate(values(basis_count), derivatives(basis_count))
        values = work(:basis_count)
        if (degree == 0) then
            derivatives = 0.0_dp
        else
            do basis = 1, basis_count
                derivatives(basis) = 0.0_dp
                denominator = knots(basis + degree) - knots(basis)
                if (denominator > 0.0_dp) derivatives(basis) = &
                    real(degree, dp)*lower(basis)/denominator
                denominator = &
                    knots(basis + degree + 1) - knots(basis + 1)
                if (denominator > 0.0_dp) derivatives(basis) = &
                    derivatives(basis) - &
                    real(degree, dp)*lower(basis + 1)/denominator
            end do
        end if
        status = 0
    end subroutine evaluate_bspline_basis

    subroutine build_bspline_derivative_matrix( &
            knots, degree, matrix, status)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp) :: scale
        integer :: basis, basis_count

        status = 1
        if (.not. valid_knot_vector(knots, degree) .or. degree < 1) return
        basis_count = size(knots) - degree - 1
        allocate(matrix(basis_count - 1, basis_count))
        matrix = 0.0_dp
        do basis = 1, basis_count - 1
            if (knots(basis + degree + 1) <= knots(basis + 1)) return
            scale = real(degree, dp)/ &
                (knots(basis + degree + 1) - knots(basis + 1))
            matrix(basis, basis) = -scale
            matrix(basis, basis + 1) = scale
        end do
        status = 0
    end subroutine build_bspline_derivative_matrix

    subroutine evaluate_nurbs_surface_geometry( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            coordinate_x, coordinate_y, point, jacobian, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        real(dp), intent(in) :: coordinate_x, coordinate_y
        real(dp), intent(out) :: point(:), jacobian(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: derivative_x(:), derivative_y(:)
        real(dp), allocatable :: value_x(:), value_y(:)
        real(dp) :: denominator, derivative_denominator(2), factor
        real(dp) :: numerator(size(point)), derivative_numerator(size(point), 2)
        integer :: ix, iy

        status = 1
        point = 0.0_dp
        jacobian = 0.0_dp
        call evaluate_bspline_basis( &
            knots_x, degree_x, coordinate_x, value_x, derivative_x, status)
        if (status /= 0) return
        call evaluate_bspline_basis( &
            knots_y, degree_y, coordinate_y, value_y, derivative_y, status)
        if (status /= 0) return
        if (size(control_points, 1) /= size(point)) return
        if (size(control_points, 2) /= size(value_x) .or. &
            size(control_points, 3) /= size(value_y)) return
        if (any(shape(weights) /= [size(value_x), size(value_y)])) return
        if (size(jacobian, 1) /= size(point) .or. &
            size(jacobian, 2) /= 2 .or. any(weights <= 0.0_dp)) return
        denominator = 0.0_dp
        derivative_denominator = 0.0_dp
        numerator = 0.0_dp
        derivative_numerator = 0.0_dp
        do iy = 1, size(value_y)
            do ix = 1, size(value_x)
                factor = weights(ix, iy)*value_x(ix)*value_y(iy)
                denominator = denominator + factor
                numerator = numerator + factor*control_points(:, ix, iy)
                factor = weights(ix, iy)*derivative_x(ix)*value_y(iy)
                derivative_denominator(1) = derivative_denominator(1) + factor
                derivative_numerator(:, 1) = derivative_numerator(:, 1) + &
                    factor*control_points(:, ix, iy)
                factor = weights(ix, iy)*value_x(ix)*derivative_y(iy)
                derivative_denominator(2) = derivative_denominator(2) + factor
                derivative_numerator(:, 2) = derivative_numerator(:, 2) + &
                    factor*control_points(:, ix, iy)
            end do
        end do
        if (denominator <= tiny(1.0_dp)) return
        point = numerator/denominator
        do ix = 1, 2
            jacobian(:, ix) = ( &
                derivative_numerator(:, ix) - &
                point*derivative_denominator(ix))/denominator
        end do
        status = 0
    end subroutine evaluate_nurbs_surface_geometry

    subroutine evaluate_nurbs_surface_geometry_jvp( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            control_points_dot, weights_dot, coordinate_x, coordinate_y, &
            point_dot, jacobian_dot, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        real(dp), intent(in) :: control_points_dot(:, :, :), weights_dot(:, :)
        real(dp), intent(in) :: coordinate_x, coordinate_y
        real(dp), intent(out) :: point_dot(:), jacobian_dot(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: derivative_x(:), derivative_y(:)
        real(dp), allocatable :: value_x(:), value_y(:)
        real(dp) :: denominator, denominator_dot
        real(dp) :: derivative_denominator(2), derivative_denominator_dot(2)
        real(dp) :: derivative_numerator(size(point_dot), 2)
        real(dp) :: derivative_numerator_dot(size(point_dot), 2)
        real(dp) :: factor, factor_dot, numerator(size(point_dot))
        real(dp) :: numerator_dot(size(point_dot))
        integer :: component, ix, iy

        status = 1
        point_dot = 0.0_dp
        jacobian_dot = 0.0_dp
        call evaluate_bspline_basis( &
            knots_x, degree_x, coordinate_x, value_x, derivative_x, status)
        if (status /= 0) return
        call evaluate_bspline_basis( &
            knots_y, degree_y, coordinate_y, value_y, derivative_y, status)
        if (status /= 0) return
        if (size(control_points, 1) /= size(point_dot)) return
        if (any(shape(control_points_dot) /= shape(control_points))) return
        if (any(shape(weights_dot) /= shape(weights))) return
        if (size(control_points, 2) /= size(value_x) .or. &
            size(control_points, 3) /= size(value_y)) return
        if (any(shape(weights) /= [size(value_x), size(value_y)])) return
        if (size(jacobian_dot, 1) /= size(point_dot) .or. &
            size(jacobian_dot, 2) /= 2 .or. any(weights <= 0.0_dp)) return

        denominator = 0.0_dp
        denominator_dot = 0.0_dp
        derivative_denominator = 0.0_dp
        derivative_denominator_dot = 0.0_dp
        numerator = 0.0_dp
        numerator_dot = 0.0_dp
        derivative_numerator = 0.0_dp
        derivative_numerator_dot = 0.0_dp
        do iy = 1, size(value_y)
            do ix = 1, size(value_x)
                factor = weights(ix, iy)*value_x(ix)*value_y(iy)
                factor_dot = weights_dot(ix, iy)*value_x(ix)*value_y(iy)
                denominator = denominator + factor
                denominator_dot = denominator_dot + factor_dot
                numerator = numerator + factor*control_points(:, ix, iy)
                numerator_dot = numerator_dot + &
                    factor_dot*control_points(:, ix, iy) + &
                    factor*control_points_dot(:, ix, iy)

                factor = weights(ix, iy)*derivative_x(ix)*value_y(iy)
                factor_dot = weights_dot(ix, iy)*derivative_x(ix)*value_y(iy)
                derivative_denominator(1) = derivative_denominator(1) + factor
                derivative_denominator_dot(1) = &
                    derivative_denominator_dot(1) + factor_dot
                derivative_numerator(:, 1) = derivative_numerator(:, 1) + &
                    factor*control_points(:, ix, iy)
                derivative_numerator_dot(:, 1) = &
                    derivative_numerator_dot(:, 1) + &
                    factor_dot*control_points(:, ix, iy) + &
                    factor*control_points_dot(:, ix, iy)

                factor = weights(ix, iy)*value_x(ix)*derivative_y(iy)
                factor_dot = weights_dot(ix, iy)*value_x(ix)*derivative_y(iy)
                derivative_denominator(2) = derivative_denominator(2) + factor
                derivative_denominator_dot(2) = &
                    derivative_denominator_dot(2) + factor_dot
                derivative_numerator(:, 2) = derivative_numerator(:, 2) + &
                    factor*control_points(:, ix, iy)
                derivative_numerator_dot(:, 2) = &
                    derivative_numerator_dot(:, 2) + &
                    factor_dot*control_points(:, ix, iy) + &
                    factor*control_points_dot(:, ix, iy)
            end do
        end do
        if (denominator <= tiny(1.0_dp)) return
        do component = 1, size(point_dot)
            call generated_nurbs_geometry_quotient_jvp( &
                numerator(component), derivative_numerator(component, 1), &
                derivative_numerator(component, 2), denominator, &
                derivative_denominator(1), derivative_denominator(2), &
                numerator_dot(component), &
                derivative_numerator_dot(component, 1), &
                derivative_numerator_dot(component, 2), denominator_dot, &
                derivative_denominator_dot(1), &
                derivative_denominator_dot(2), point_dot(component), &
                jacobian_dot(component, 1), jacobian_dot(component, 2))
        end do
        status = 0
    end subroutine evaluate_nurbs_surface_geometry_jvp

    subroutine evaluate_nurbs_surface_geometry_vjp( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            coordinate_x, coordinate_y, point_bar, jacobian_bar, &
            control_points_bar, weights_bar, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        real(dp), intent(in) :: coordinate_x, coordinate_y
        real(dp), intent(in) :: point_bar(:), jacobian_bar(:, :)
        real(dp), intent(out) :: control_points_bar(:, :, :), weights_bar(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: derivative_x(:), derivative_y(:)
        real(dp), allocatable :: value_x(:), value_y(:)
        real(dp) :: denominator, denominator_bar, denominator_bar_component
        real(dp) :: derivative_denominator(2), derivative_denominator_bar(2)
        real(dp) :: derivative_denominator_bar_component_x
        real(dp) :: derivative_denominator_bar_component_y
        real(dp) :: derivative_numerator(size(point_bar), 2)
        real(dp) :: derivative_numerator_bar(size(point_bar), 2)
        real(dp) :: factor, numerator(size(point_bar))
        real(dp) :: numerator_bar(size(point_bar))
        real(dp) :: weighted_basis_bar(size(point_bar))
        integer :: component, ix, iy

        status = 1
        control_points_bar = 0.0_dp
        weights_bar = 0.0_dp
        call evaluate_bspline_basis( &
            knots_x, degree_x, coordinate_x, value_x, derivative_x, status)
        if (status /= 0) return
        call evaluate_bspline_basis( &
            knots_y, degree_y, coordinate_y, value_y, derivative_y, status)
        if (status /= 0) return
        if (size(control_points, 1) /= size(point_bar)) return
        if (any(shape(control_points_bar) /= shape(control_points))) return
        if (any(shape(weights_bar) /= shape(weights))) return
        if (size(control_points, 2) /= size(value_x) .or. &
            size(control_points, 3) /= size(value_y)) return
        if (any(shape(weights) /= [size(value_x), size(value_y)])) return
        if (size(jacobian_bar, 1) /= size(point_bar) .or. &
            size(jacobian_bar, 2) /= 2 .or. any(weights <= 0.0_dp)) return

        denominator = 0.0_dp
        derivative_denominator = 0.0_dp
        numerator = 0.0_dp
        derivative_numerator = 0.0_dp
        do iy = 1, size(value_y)
            do ix = 1, size(value_x)
                factor = weights(ix, iy)*value_x(ix)*value_y(iy)
                denominator = denominator + factor
                numerator = numerator + factor*control_points(:, ix, iy)
                factor = weights(ix, iy)*derivative_x(ix)*value_y(iy)
                derivative_denominator(1) = derivative_denominator(1) + factor
                derivative_numerator(:, 1) = derivative_numerator(:, 1) + &
                    factor*control_points(:, ix, iy)
                factor = weights(ix, iy)*value_x(ix)*derivative_y(iy)
                derivative_denominator(2) = derivative_denominator(2) + factor
                derivative_numerator(:, 2) = derivative_numerator(:, 2) + &
                    factor*control_points(:, ix, iy)
            end do
        end do
        if (denominator <= tiny(1.0_dp)) return

        denominator_bar = 0.0_dp
        derivative_denominator_bar = 0.0_dp
        do component = 1, size(point_bar)
            call generated_nurbs_geometry_quotient_vjp( &
                numerator(component), derivative_numerator(component, 1), &
                derivative_numerator(component, 2), denominator, &
                derivative_denominator(1), derivative_denominator(2), &
                point_bar(component), jacobian_bar(component, 1), &
                jacobian_bar(component, 2), numerator_bar(component), &
                derivative_numerator_bar(component, 1), &
                derivative_numerator_bar(component, 2), &
                denominator_bar_component, &
                derivative_denominator_bar_component_x, &
                derivative_denominator_bar_component_y)
            denominator_bar = denominator_bar + denominator_bar_component
            derivative_denominator_bar(1) = &
                derivative_denominator_bar(1) + &
                derivative_denominator_bar_component_x
            derivative_denominator_bar(2) = &
                derivative_denominator_bar(2) + &
                derivative_denominator_bar_component_y
        end do

        do iy = 1, size(value_y)
            do ix = 1, size(value_x)
                weighted_basis_bar = &
                    value_x(ix)*value_y(iy)*numerator_bar + &
                    derivative_x(ix)*value_y(iy)* &
                    derivative_numerator_bar(:, 1) + &
                    value_x(ix)*derivative_y(iy)* &
                    derivative_numerator_bar(:, 2)
                control_points_bar(:, ix, iy) = &
                    weights(ix, iy)*weighted_basis_bar
                weights_bar(ix, iy) = &
                    value_x(ix)*value_y(iy)*denominator_bar + &
                    derivative_x(ix)*value_y(iy)* &
                    derivative_denominator_bar(1) + &
                    value_x(ix)*derivative_y(iy)* &
                    derivative_denominator_bar(2) + &
                    dot_product(control_points(:, ix, iy), weighted_basis_bar)
            end do
        end do
        status = 0
    end subroutine evaluate_nurbs_surface_geometry_vjp

    subroutine evaluate_nurbs_volume_geometry( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            control_points, weights, coordinate_x, coordinate_y, coordinate_z, &
            point, jacobian, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:), knots_z(:)
        integer, intent(in) :: degree_x, degree_y, degree_z
        real(dp), intent(in) :: control_points(:, :, :, :), weights(:, :, :)
        real(dp), intent(in) :: coordinate_x, coordinate_y, coordinate_z
        real(dp), intent(out) :: point(:), jacobian(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: dx(:), dy(:), dz(:), vx(:), vy(:), vz(:)
        real(dp) :: denominator, derivative_denominator(3), factor
        real(dp) :: numerator(size(point)), derivative_numerator(size(point), 3)
        integer :: ix, iy, iz

        status = 1
        point = 0.0_dp
        jacobian = 0.0_dp
        call evaluate_bspline_basis( &
            knots_x, degree_x, coordinate_x, vx, dx, status)
        if (status /= 0) return
        call evaluate_bspline_basis( &
            knots_y, degree_y, coordinate_y, vy, dy, status)
        if (status /= 0) return
        call evaluate_bspline_basis( &
            knots_z, degree_z, coordinate_z, vz, dz, status)
        if (status /= 0) return
        if (size(control_points, 1) /= size(point)) return
        if (any(shape(control_points(1, :, :, :)) /= &
            [size(vx), size(vy), size(vz)])) return
        if (any(shape(weights) /= [size(vx), size(vy), size(vz)])) return
        if (size(jacobian, 1) /= size(point) .or. &
            size(jacobian, 2) /= 3 .or. any(weights <= 0.0_dp)) return
        denominator = 0.0_dp
        derivative_denominator = 0.0_dp
        numerator = 0.0_dp
        derivative_numerator = 0.0_dp
        do iz = 1, size(vz)
            do iy = 1, size(vy)
                do ix = 1, size(vx)
                    factor = weights(ix, iy, iz)*vx(ix)*vy(iy)*vz(iz)
                    denominator = denominator + factor
                    numerator = numerator + &
                        factor*control_points(:, ix, iy, iz)
                    factor = weights(ix, iy, iz)*dx(ix)*vy(iy)*vz(iz)
                    derivative_denominator(1) = &
                        derivative_denominator(1) + factor
                    derivative_numerator(:, 1) = &
                        derivative_numerator(:, 1) + &
                        factor*control_points(:, ix, iy, iz)
                    factor = weights(ix, iy, iz)*vx(ix)*dy(iy)*vz(iz)
                    derivative_denominator(2) = &
                        derivative_denominator(2) + factor
                    derivative_numerator(:, 2) = &
                        derivative_numerator(:, 2) + &
                        factor*control_points(:, ix, iy, iz)
                    factor = weights(ix, iy, iz)*vx(ix)*vy(iy)*dz(iz)
                    derivative_denominator(3) = &
                        derivative_denominator(3) + factor
                    derivative_numerator(:, 3) = &
                        derivative_numerator(:, 3) + &
                        factor*control_points(:, ix, iy, iz)
                end do
            end do
        end do
        if (denominator <= tiny(1.0_dp)) return
        point = numerator/denominator
        do ix = 1, 3
            jacobian(:, ix) = ( &
                derivative_numerator(:, ix) - &
                point*derivative_denominator(ix))/denominator
        end do
        status = 0
    end subroutine evaluate_nurbs_volume_geometry

    subroutine evaluate_nurbs_volume_geometry_jvp( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            control_points, weights, control_points_dot, weights_dot, &
            coordinate_x, coordinate_y, coordinate_z, point_dot, &
            jacobian_dot, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:), knots_z(:)
        integer, intent(in) :: degree_x, degree_y, degree_z
        real(dp), intent(in) :: control_points(:, :, :, :), weights(:, :, :)
        real(dp), intent(in) :: control_points_dot(:, :, :, :)
        real(dp), intent(in) :: weights_dot(:, :, :)
        real(dp), intent(in) :: coordinate_x, coordinate_y, coordinate_z
        real(dp), intent(out) :: point_dot(:), jacobian_dot(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: dx(:), dy(:), dz(:), vx(:), vy(:), vz(:)
        real(dp) :: denominator, denominator_dot
        real(dp) :: derivative_denominator(3), derivative_denominator_dot(3)
        real(dp) :: derivative_numerator(size(point_dot), 3)
        real(dp) :: derivative_numerator_dot(size(point_dot), 3)
        real(dp) :: factor, factor_dot, numerator(size(point_dot))
        real(dp) :: numerator_dot(size(point_dot))
        integer :: component, ix, iy, iz

        status = 1
        point_dot = 0.0_dp
        jacobian_dot = 0.0_dp
        call evaluate_bspline_basis( &
            knots_x, degree_x, coordinate_x, vx, dx, status)
        if (status /= 0) return
        call evaluate_bspline_basis( &
            knots_y, degree_y, coordinate_y, vy, dy, status)
        if (status /= 0) return
        call evaluate_bspline_basis( &
            knots_z, degree_z, coordinate_z, vz, dz, status)
        if (status /= 0) return
        if (size(control_points, 1) /= size(point_dot)) return
        if (any(shape(control_points_dot) /= shape(control_points))) return
        if (any(shape(weights_dot) /= shape(weights))) return
        if (any(shape(control_points(1, :, :, :)) /= &
            [size(vx), size(vy), size(vz)])) return
        if (any(shape(weights) /= [size(vx), size(vy), size(vz)])) return
        if (size(jacobian_dot, 1) /= size(point_dot) .or. &
            size(jacobian_dot, 2) /= 3 .or. any(weights <= 0.0_dp)) return

        denominator = 0.0_dp
        denominator_dot = 0.0_dp
        derivative_denominator = 0.0_dp
        derivative_denominator_dot = 0.0_dp
        numerator = 0.0_dp
        numerator_dot = 0.0_dp
        derivative_numerator = 0.0_dp
        derivative_numerator_dot = 0.0_dp
        do iz = 1, size(vz)
            do iy = 1, size(vy)
                do ix = 1, size(vx)
                    factor = weights(ix, iy, iz)*vx(ix)*vy(iy)*vz(iz)
                    factor_dot = weights_dot(ix, iy, iz)*vx(ix)*vy(iy)*vz(iz)
                    denominator = denominator + factor
                    denominator_dot = denominator_dot + factor_dot
                    numerator = numerator + &
                        factor*control_points(:, ix, iy, iz)
                    numerator_dot = numerator_dot + &
                        factor_dot*control_points(:, ix, iy, iz) + &
                        factor*control_points_dot(:, ix, iy, iz)

                    factor = weights(ix, iy, iz)*dx(ix)*vy(iy)*vz(iz)
                    factor_dot = weights_dot(ix, iy, iz)*dx(ix)*vy(iy)*vz(iz)
                    derivative_denominator(1) = &
                        derivative_denominator(1) + factor
                    derivative_denominator_dot(1) = &
                        derivative_denominator_dot(1) + factor_dot
                    derivative_numerator(:, 1) = &
                        derivative_numerator(:, 1) + &
                        factor*control_points(:, ix, iy, iz)
                    derivative_numerator_dot(:, 1) = &
                        derivative_numerator_dot(:, 1) + &
                        factor_dot*control_points(:, ix, iy, iz) + &
                        factor*control_points_dot(:, ix, iy, iz)

                    factor = weights(ix, iy, iz)*vx(ix)*dy(iy)*vz(iz)
                    factor_dot = weights_dot(ix, iy, iz)*vx(ix)*dy(iy)*vz(iz)
                    derivative_denominator(2) = &
                        derivative_denominator(2) + factor
                    derivative_denominator_dot(2) = &
                        derivative_denominator_dot(2) + factor_dot
                    derivative_numerator(:, 2) = &
                        derivative_numerator(:, 2) + &
                        factor*control_points(:, ix, iy, iz)
                    derivative_numerator_dot(:, 2) = &
                        derivative_numerator_dot(:, 2) + &
                        factor_dot*control_points(:, ix, iy, iz) + &
                        factor*control_points_dot(:, ix, iy, iz)

                    factor = weights(ix, iy, iz)*vx(ix)*vy(iy)*dz(iz)
                    factor_dot = weights_dot(ix, iy, iz)*vx(ix)*vy(iy)*dz(iz)
                    derivative_denominator(3) = &
                        derivative_denominator(3) + factor
                    derivative_denominator_dot(3) = &
                        derivative_denominator_dot(3) + factor_dot
                    derivative_numerator(:, 3) = &
                        derivative_numerator(:, 3) + &
                        factor*control_points(:, ix, iy, iz)
                    derivative_numerator_dot(:, 3) = &
                        derivative_numerator_dot(:, 3) + &
                        factor_dot*control_points(:, ix, iy, iz) + &
                        factor*control_points_dot(:, ix, iy, iz)
                end do
            end do
        end do
        if (denominator <= tiny(1.0_dp)) return
        do component = 1, size(point_dot)
            call generated_nurbs_volume_quotient_jvp( &
                numerator(component), derivative_numerator(component, 1), &
                derivative_numerator(component, 2), &
                derivative_numerator(component, 3), denominator, &
                derivative_denominator(1), derivative_denominator(2), &
                derivative_denominator(3), numerator_dot(component), &
                derivative_numerator_dot(component, 1), &
                derivative_numerator_dot(component, 2), &
                derivative_numerator_dot(component, 3), denominator_dot, &
                derivative_denominator_dot(1), &
                derivative_denominator_dot(2), &
                derivative_denominator_dot(3), point_dot(component), &
                jacobian_dot(component, 1), jacobian_dot(component, 2), &
                jacobian_dot(component, 3))
        end do
        status = 0
    end subroutine evaluate_nurbs_volume_geometry_jvp

    subroutine evaluate_nurbs_volume_geometry_vjp( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            control_points, weights, coordinate_x, coordinate_y, coordinate_z, &
            point_bar, jacobian_bar, control_points_bar, weights_bar, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:), knots_z(:)
        integer, intent(in) :: degree_x, degree_y, degree_z
        real(dp), intent(in) :: control_points(:, :, :, :), weights(:, :, :)
        real(dp), intent(in) :: coordinate_x, coordinate_y, coordinate_z
        real(dp), intent(in) :: point_bar(:), jacobian_bar(:, :)
        real(dp), intent(out) :: control_points_bar(:, :, :, :)
        real(dp), intent(out) :: weights_bar(:, :, :)
        integer, intent(out) :: status

        real(dp), allocatable :: dx(:), dy(:), dz(:), vx(:), vy(:), vz(:)
        real(dp) :: denominator, denominator_bar, denominator_bar_component
        real(dp) :: derivative_denominator(3), derivative_denominator_bar(3)
        real(dp) :: derivative_denominator_bar_component(3)
        real(dp) :: derivative_numerator(size(point_bar), 3)
        real(dp) :: derivative_numerator_bar(size(point_bar), 3)
        real(dp) :: factor, numerator(size(point_bar))
        real(dp) :: numerator_bar(size(point_bar))
        real(dp) :: weighted_basis_bar(size(point_bar))
        integer :: component, ix, iy, iz

        status = 1
        control_points_bar = 0.0_dp
        weights_bar = 0.0_dp
        call evaluate_bspline_basis( &
            knots_x, degree_x, coordinate_x, vx, dx, status)
        if (status /= 0) return
        call evaluate_bspline_basis( &
            knots_y, degree_y, coordinate_y, vy, dy, status)
        if (status /= 0) return
        call evaluate_bspline_basis( &
            knots_z, degree_z, coordinate_z, vz, dz, status)
        if (status /= 0) return
        if (size(control_points, 1) /= size(point_bar)) return
        if (any(shape(control_points_bar) /= shape(control_points))) return
        if (any(shape(weights_bar) /= shape(weights))) return
        if (any(shape(control_points(1, :, :, :)) /= &
            [size(vx), size(vy), size(vz)])) return
        if (any(shape(weights) /= [size(vx), size(vy), size(vz)])) return
        if (size(jacobian_bar, 1) /= size(point_bar) .or. &
            size(jacobian_bar, 2) /= 3 .or. any(weights <= 0.0_dp)) return

        denominator = 0.0_dp
        derivative_denominator = 0.0_dp
        numerator = 0.0_dp
        derivative_numerator = 0.0_dp
        do iz = 1, size(vz)
            do iy = 1, size(vy)
                do ix = 1, size(vx)
                    factor = weights(ix, iy, iz)*vx(ix)*vy(iy)*vz(iz)
                    denominator = denominator + factor
                    numerator = numerator + &
                        factor*control_points(:, ix, iy, iz)
                    factor = weights(ix, iy, iz)*dx(ix)*vy(iy)*vz(iz)
                    derivative_denominator(1) = &
                        derivative_denominator(1) + factor
                    derivative_numerator(:, 1) = &
                        derivative_numerator(:, 1) + &
                        factor*control_points(:, ix, iy, iz)
                    factor = weights(ix, iy, iz)*vx(ix)*dy(iy)*vz(iz)
                    derivative_denominator(2) = &
                        derivative_denominator(2) + factor
                    derivative_numerator(:, 2) = &
                        derivative_numerator(:, 2) + &
                        factor*control_points(:, ix, iy, iz)
                    factor = weights(ix, iy, iz)*vx(ix)*vy(iy)*dz(iz)
                    derivative_denominator(3) = &
                        derivative_denominator(3) + factor
                    derivative_numerator(:, 3) = &
                        derivative_numerator(:, 3) + &
                        factor*control_points(:, ix, iy, iz)
                end do
            end do
        end do
        if (denominator <= tiny(1.0_dp)) return

        denominator_bar = 0.0_dp
        derivative_denominator_bar = 0.0_dp
        do component = 1, size(point_bar)
            call generated_nurbs_volume_quotient_vjp( &
                numerator(component), derivative_numerator(component, 1), &
                derivative_numerator(component, 2), &
                derivative_numerator(component, 3), denominator, &
                derivative_denominator(1), derivative_denominator(2), &
                derivative_denominator(3), point_bar(component), &
                jacobian_bar(component, 1), jacobian_bar(component, 2), &
                jacobian_bar(component, 3), numerator_bar(component), &
                derivative_numerator_bar(component, 1), &
                derivative_numerator_bar(component, 2), &
                derivative_numerator_bar(component, 3), &
                denominator_bar_component, &
                derivative_denominator_bar_component(1), &
                derivative_denominator_bar_component(2), &
                derivative_denominator_bar_component(3))
            denominator_bar = denominator_bar + denominator_bar_component
            derivative_denominator_bar = derivative_denominator_bar + &
                derivative_denominator_bar_component
        end do

        do iz = 1, size(vz)
            do iy = 1, size(vy)
                do ix = 1, size(vx)
                    weighted_basis_bar = &
                        vx(ix)*vy(iy)*vz(iz)*numerator_bar + &
                        dx(ix)*vy(iy)*vz(iz)* &
                        derivative_numerator_bar(:, 1) + &
                        vx(ix)*dy(iy)*vz(iz)* &
                        derivative_numerator_bar(:, 2) + &
                        vx(ix)*vy(iy)*dz(iz)* &
                        derivative_numerator_bar(:, 3)
                    control_points_bar(:, ix, iy, iz) = &
                        weights(ix, iy, iz)*weighted_basis_bar
                    weights_bar(ix, iy, iz) = &
                        vx(ix)*vy(iy)*vz(iz)*denominator_bar + &
                        dx(ix)*vy(iy)*vz(iz)* &
                        derivative_denominator_bar(1) + &
                        vx(ix)*dy(iy)*vz(iz)* &
                        derivative_denominator_bar(2) + &
                        vx(ix)*vy(iy)*dz(iz)* &
                        derivative_denominator_bar(3) + &
                        dot_product( &
                        control_points(:, ix, iy, iz), weighted_basis_bar)
                end do
            end do
        end do
        status = 0
    end subroutine evaluate_nurbs_volume_geometry_vjp

    subroutine map_isogeometric_h1_gradient( &
            jacobian, reference_gradients, physical_gradients, determinant, &
            status)
        real(dp), intent(in) :: jacobian(:, :), reference_gradients(:, :)
        real(dp), allocatable, intent(out) :: physical_gradients(:, :)
        real(dp), intent(out) :: determinant
        integer, intent(out) :: status

        call map_covariant( &
            jacobian, reference_gradients, physical_gradients, determinant, &
            status)
    end subroutine map_isogeometric_h1_gradient

    subroutine map_isogeometric_hcurl( &
            jacobian, reference_values, physical_values, determinant, status)
        real(dp), intent(in) :: jacobian(:, :), reference_values(:, :)
        real(dp), allocatable, intent(out) :: physical_values(:, :)
        real(dp), intent(out) :: determinant
        integer, intent(out) :: status

        call map_covariant( &
            jacobian, reference_values, physical_values, determinant, status)
    end subroutine map_isogeometric_hcurl

    subroutine map_covariant( &
            jacobian, reference_values, physical_values, determinant, status)
        real(dp), intent(in) :: jacobian(:, :), reference_values(:, :)
        real(dp), allocatable, intent(out) :: physical_values(:, :)
        real(dp), intent(out) :: determinant
        integer, intent(out) :: status

        real(dp) :: inverse(size(jacobian, 1), size(jacobian, 2))

        status = 1
        determinant = 0.0_dp
        if (size(jacobian, 1) /= size(jacobian, 2)) return
        if (size(reference_values, 1) /= size(jacobian, 1)) return
        call small_inverse(jacobian, inverse, determinant, status)
        if (status /= 0) return
        allocate(physical_values, mold=reference_values)
        physical_values = matmul(transpose(inverse), reference_values)
        status = 0
    end subroutine map_covariant

    subroutine map_isogeometric_hdiv( &
            jacobian, reference_values, physical_values, determinant, status)
        real(dp), intent(in) :: jacobian(:, :), reference_values(:, :)
        real(dp), allocatable, intent(out) :: physical_values(:, :)
        real(dp), intent(out) :: determinant
        integer, intent(out) :: status

        real(dp) :: inverse(size(jacobian, 1), size(jacobian, 2))

        status = 1
        determinant = 0.0_dp
        if (size(jacobian, 1) /= size(jacobian, 2)) return
        if (size(reference_values, 1) /= size(jacobian, 1)) return
        call small_inverse(jacobian, inverse, determinant, status)
        if (status /= 0) return
        allocate(physical_values, mold=reference_values)
        physical_values = matmul(jacobian, reference_values)/determinant
        status = 0
    end subroutine map_isogeometric_hdiv

    subroutine map_isogeometric_l2( &
            jacobian, reference_values, physical_values, determinant, status)
        real(dp), intent(in) :: jacobian(:, :), reference_values(:)
        real(dp), allocatable, intent(out) :: physical_values(:)
        real(dp), intent(out) :: determinant
        integer, intent(out) :: status

        real(dp) :: inverse(size(jacobian, 1), size(jacobian, 2))

        status = 1
        determinant = 0.0_dp
        if (size(jacobian, 1) /= size(jacobian, 2)) return
        call small_inverse(jacobian, inverse, determinant, status)
        if (status /= 0) return
        allocate(physical_values, mold=reference_values)
        physical_values = reference_values/determinant
        status = 0
    end subroutine map_isogeometric_l2

    pure subroutine small_inverse(matrix, inverse, determinant, status)
        real(dp), intent(in) :: matrix(:, :)
        real(dp), intent(out) :: inverse(:, :), determinant
        integer, intent(out) :: status

        real(dp) :: scale

        inverse = 0.0_dp
        determinant = 0.0_dp
        status = 1
        if (any(shape(matrix) /= shape(inverse))) return
        select case (size(matrix, 1))
        case (2)
            determinant = matrix(1, 1)*matrix(2, 2) - &
                matrix(1, 2)*matrix(2, 1)
            scale = max(1.0_dp, maxval(abs(matrix))**2)
            if (abs(determinant) <= 64.0_dp*epsilon(1.0_dp)*scale) return
            inverse = reshape([ &
                matrix(2, 2), -matrix(2, 1), &
                -matrix(1, 2), matrix(1, 1)], [2, 2])/determinant
        case (3)
            determinant = &
                matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - &
                matrix(2, 3)*matrix(3, 2)) - &
                matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - &
                matrix(2, 3)*matrix(3, 1)) + &
                matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - &
                matrix(2, 2)*matrix(3, 1))
            scale = max(1.0_dp, maxval(abs(matrix))**3)
            if (abs(determinant) <= 128.0_dp*epsilon(1.0_dp)*scale) return
            inverse(1, 1) = matrix(2, 2)*matrix(3, 3) - &
                matrix(2, 3)*matrix(3, 2)
            inverse(1, 2) = matrix(1, 3)*matrix(3, 2) - &
                matrix(1, 2)*matrix(3, 3)
            inverse(1, 3) = matrix(1, 2)*matrix(2, 3) - &
                matrix(1, 3)*matrix(2, 2)
            inverse(2, 1) = matrix(2, 3)*matrix(3, 1) - &
                matrix(2, 1)*matrix(3, 3)
            inverse(2, 2) = matrix(1, 1)*matrix(3, 3) - &
                matrix(1, 3)*matrix(3, 1)
            inverse(2, 3) = matrix(1, 3)*matrix(2, 1) - &
                matrix(1, 1)*matrix(2, 3)
            inverse(3, 1) = matrix(2, 1)*matrix(3, 2) - &
                matrix(2, 2)*matrix(3, 1)
            inverse(3, 2) = matrix(1, 2)*matrix(3, 1) - &
                matrix(1, 1)*matrix(3, 2)
            inverse(3, 3) = matrix(1, 1)*matrix(2, 2) - &
                matrix(1, 2)*matrix(2, 1)
            inverse = inverse/determinant
        case default
            return
        end select
        status = 0
    end subroutine small_inverse

    subroutine build_bspline_feec_2d_operators( &
            knots_x, knots_y, degree_x, degree_y, gradient, curl, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y
        real(dp), allocatable, intent(out) :: gradient(:, :), curl(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: derivative_x(:, :), derivative_y(:, :)
        integer :: column, ix, iy, nx, ny, row, x_component_count

        status = 1
        call build_bspline_derivative_matrix( &
            knots_x, degree_x, derivative_x, status)
        if (status /= 0) return
        call build_bspline_derivative_matrix( &
            knots_y, degree_y, derivative_y, status)
        if (status /= 0) return
        nx = size(derivative_x, 2)
        ny = size(derivative_y, 2)
        x_component_count = (nx - 1)*ny
        allocate( &
            gradient(x_component_count + nx*(ny - 1), nx*ny), &
            curl((nx - 1)*(ny - 1), x_component_count + nx*(ny - 1)))
        gradient = 0.0_dp
        curl = 0.0_dp
        do iy = 1, ny
            do ix = 1, nx - 1
                row = ix + (iy - 1)*(nx - 1)
                do column = 1, nx
                    gradient(row, column + (iy - 1)*nx) = &
                        derivative_x(ix, column)
                end do
            end do
        end do
        do iy = 1, ny - 1
            do ix = 1, nx
                row = x_component_count + ix + (iy - 1)*nx
                do column = 1, ny
                    gradient(row, ix + (column - 1)*nx) = &
                        derivative_y(iy, column)
                end do
            end do
        end do
        do iy = 1, ny - 1
            do ix = 1, nx - 1
                row = ix + (iy - 1)*(nx - 1)
                do column = 1, ny
                    curl(row, ix + (column - 1)*(nx - 1)) = &
                        -derivative_y(iy, column)
                end do
                do column = 1, nx
                    curl(row, x_component_count + &
                        column + (iy - 1)*nx) = derivative_x(ix, column)
                end do
            end do
        end do
        status = 0
    end subroutine build_bspline_feec_2d_operators

    subroutine build_bspline_feec_3d_operators( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            gradient, curl, divergence, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:), knots_z(:)
        integer, intent(in) :: degree_x, degree_y, degree_z
        real(dp), allocatable, intent(out) :: gradient(:, :), curl(:, :)
        real(dp), allocatable, intent(out) :: divergence(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: dx(:, :), dy(:, :), dz(:, :)
        integer :: bx_count, by_count, column, ex_count, ey_count
        integer :: ix, iy, iz, nx, ny, nz, row

        status = 1
        call build_bspline_derivative_matrix(knots_x, degree_x, dx, status)
        if (status /= 0) return
        call build_bspline_derivative_matrix(knots_y, degree_y, dy, status)
        if (status /= 0) return
        call build_bspline_derivative_matrix(knots_z, degree_z, dz, status)
        if (status /= 0) return
        nx = size(dx, 2)
        ny = size(dy, 2)
        nz = size(dz, 2)
        ex_count = (nx - 1)*ny*nz
        ey_count = nx*(ny - 1)*nz
        bx_count = nx*(ny - 1)*(nz - 1)
        by_count = (nx - 1)*ny*(nz - 1)
        allocate( &
            gradient(ex_count + ey_count + nx*ny*(nz - 1), nx*ny*nz), &
            curl(bx_count + by_count + (nx - 1)*(ny - 1)*nz, &
            ex_count + ey_count + nx*ny*(nz - 1)), &
            divergence((nx - 1)*(ny - 1)*(nz - 1), &
            bx_count + by_count + (nx - 1)*(ny - 1)*nz))
        gradient = 0.0_dp
        curl = 0.0_dp
        divergence = 0.0_dp

        do iz = 1, nz
            do iy = 1, ny
                do ix = 1, nx - 1
                    row = index_3d(ix, iy, iz, nx - 1, ny)
                    do column = 1, nx
                        gradient(row, index_3d(column, iy, iz, nx, ny)) = &
                            dx(ix, column)
                    end do
                end do
            end do
        end do
        do iz = 1, nz
            do iy = 1, ny - 1
                do ix = 1, nx
                    row = ex_count + index_3d(ix, iy, iz, nx, ny - 1)
                    do column = 1, ny
                        gradient(row, index_3d(ix, column, iz, nx, ny)) = &
                            dy(iy, column)
                    end do
                end do
            end do
        end do
        do iz = 1, nz - 1
            do iy = 1, ny
                do ix = 1, nx
                    row = ex_count + ey_count + index_3d(ix, iy, iz, nx, ny)
                    do column = 1, nz
                        gradient(row, index_3d(ix, iy, column, nx, ny)) = &
                            dz(iz, column)
                    end do
                end do
            end do
        end do

        do iz = 1, nz - 1
            do iy = 1, ny - 1
                do ix = 1, nx
                    row = index_3d(ix, iy, iz, nx, ny - 1)
                    do column = 1, ny
                        curl(row, ex_count + ey_count + &
                            index_3d(ix, column, iz, nx, ny)) = dy(iy, column)
                    end do
                    do column = 1, nz
                        curl(row, ex_count + &
                            index_3d(ix, iy, column, nx, ny - 1)) = &
                            -dz(iz, column)
                    end do
                end do
            end do
        end do
        do iz = 1, nz - 1
            do iy = 1, ny
                do ix = 1, nx - 1
                    row = bx_count + index_3d(ix, iy, iz, nx - 1, ny)
                    do column = 1, nz
                        curl(row, index_3d(ix, iy, column, nx - 1, ny)) = &
                            dz(iz, column)
                    end do
                    do column = 1, nx
                        curl(row, ex_count + ey_count + &
                            index_3d(column, iy, iz, nx, ny)) = -dx(ix, column)
                    end do
                end do
            end do
        end do
        do iz = 1, nz
            do iy = 1, ny - 1
                do ix = 1, nx - 1
                    row = bx_count + by_count + &
                        index_3d(ix, iy, iz, nx - 1, ny - 1)
                    do column = 1, nx
                        curl(row, ex_count + &
                            index_3d(column, iy, iz, nx, ny - 1)) = &
                            dx(ix, column)
                    end do
                    do column = 1, ny
                        curl(row, index_3d(ix, column, iz, nx - 1, ny)) = &
                            -dy(iy, column)
                    end do
                end do
            end do
        end do

        do iz = 1, nz - 1
            do iy = 1, ny - 1
                do ix = 1, nx - 1
                    row = index_3d(ix, iy, iz, nx - 1, ny - 1)
                    do column = 1, nx
                        divergence(row, &
                            index_3d(column, iy, iz, nx, ny - 1)) = &
                            dx(ix, column)
                    end do
                    do column = 1, ny
                        divergence(row, bx_count + &
                            index_3d(ix, column, iz, nx - 1, ny)) = &
                            dy(iy, column)
                    end do
                    do column = 1, nz
                        divergence(row, bx_count + by_count + &
                            index_3d(ix, iy, column, nx - 1, ny - 1)) = &
                            dz(iz, column)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine build_bspline_feec_3d_operators

    pure integer function index_3d( &
            ix, iy, iz, count_x, count_y) result(index)
        integer, intent(in) :: ix, iy, iz, count_x, count_y

        index = ix + (iy - 1)*count_x + (iz - 1)*count_x*count_y
    end function index_3d

    pure logical function valid_knot_vector(knots, degree) result(valid)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree

        valid = degree >= 0 .and. size(knots) >= 2*degree + 2
        if (.not. valid) return
        valid = all(knots(2:) >= knots(:size(knots) - 1))
        if (.not. valid) return
        valid = knots(degree + 1) < knots(size(knots) - degree)
    end function valid_knot_vector

end module fortfem_bspline_feec
