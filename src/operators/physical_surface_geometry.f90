module fortfem_physical_surface_geometry
    !! Fixed-topology sampler for a two-parameter surface in three dimensions.
    !!
    !! A geometry provider evaluates its map and the two physical tangent
    !! columns at each reference quadrature coordinate.  This small layer
    !! copies the physical map values, normalizes the oriented cross product,
    !! and exposes one positive surface measure for mortar, FEM, BEM, and DtN
    !! consumers.  The map evaluation itself remains provider-owned: a NURBS,
    !! Fourier, panel, or cut-cell provider can feed the same contract.
    !!
    !! The orientation sign is a fixed discrete choice (+1 or -1).  Changing
    !! it, changing quadrature ownership, or creating a degenerate tangent
    !! pair is a topology event and is rejected rather than differentiated.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: sample_physical_surface_geometry
    public :: sample_physical_surface_geometry_jvp
    public :: sample_physical_surface_geometry_vjp

contains

    subroutine sample_physical_surface_geometry( &
            reference_coordinates, map_points, map_tangents, orientation_sign, &
            physical_points, normals, surface_jacobian, status)
        !! Sample physical points, oriented unit normals, and surface measure.
        real(dp), intent(in) :: reference_coordinates(:, :)
        real(dp), intent(in) :: map_points(:, :)
        real(dp), intent(in) :: map_tangents(:, :, :)
        integer, intent(in) :: orientation_sign
        real(dp), intent(out) :: physical_points(:, :), normals(:, :)
        real(dp), intent(out) :: surface_jacobian(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, quadrature_count
        real(dp) :: cross_product(3), jacobian

        physical_points = 0.0_dp
        normals = 0.0_dp
        surface_jacobian = 0.0_dp
        call validate_surface_geometry_inputs( &
            reference_coordinates, map_points, map_tangents, orientation_sign, &
            physical_points, normals, surface_jacobian, status)
        if (status%code /= FORTSPARSE_OK) return

        quadrature_count = size(reference_coordinates, 2)
        physical_points = map_points
        do quadrature = 1, quadrature_count
            cross_product = cross3( &
                map_tangents(:, 1, quadrature), map_tangents(:, 2, quadrature))
            jacobian = sqrt(dot_product(cross_product, cross_product))
            surface_jacobian(quadrature) = jacobian
            normals(:, quadrature) = real(orientation_sign, dp)* &
                cross_product/jacobian
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine sample_physical_surface_geometry

    subroutine sample_physical_surface_geometry_jvp( &
            reference_coordinates, map_points, map_tangents, orientation_sign, &
            map_points_dot, map_tangents_dot, physical_points_dot, normals_dot, &
            surface_jacobian_dot, status)
        !! Apply the fixed-topology tangent product of the surface sampler.
        real(dp), intent(in) :: reference_coordinates(:, :)
        real(dp), intent(in) :: map_points(:, :)
        real(dp), intent(in) :: map_tangents(:, :, :)
        integer, intent(in) :: orientation_sign
        real(dp), intent(in) :: map_points_dot(:, :)
        real(dp), intent(in) :: map_tangents_dot(:, :, :)
        real(dp), intent(out) :: physical_points_dot(:, :), normals_dot(:, :)
        real(dp), intent(out) :: surface_jacobian_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: cross_product(3), cross_product_dot(3)
        real(dp) :: jacobian, jacobian_dot
        integer :: quadrature, quadrature_count

        physical_points_dot = 0.0_dp
        normals_dot = 0.0_dp
        surface_jacobian_dot = 0.0_dp
        call validate_surface_geometry_inputs( &
            reference_coordinates, map_points, map_tangents, orientation_sign, &
            physical_points_dot, normals_dot, surface_jacobian_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(reference_coordinates, 2)
        if (any(shape(map_points_dot) /= shape(map_points)) .or. &
            any(shape(map_tangents_dot) /= shape(map_tangents))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "physical surface geometry JVP has incompatible increments")
            return
        end if
        if (any(.not. ieee_is_finite(map_points_dot)) .or. &
            any(.not. ieee_is_finite(map_tangents_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "physical surface geometry JVP received non-finite increments")
            return
        end if

        physical_points_dot = map_points_dot
        do quadrature = 1, quadrature_count
            cross_product = cross3( &
                map_tangents(:, 1, quadrature), map_tangents(:, 2, quadrature))
            cross_product_dot = cross3( &
                map_tangents_dot(:, 1, quadrature), &
                map_tangents(:, 2, quadrature)) + cross3( &
                map_tangents(:, 1, quadrature), &
                map_tangents_dot(:, 2, quadrature))
            jacobian = sqrt(dot_product(cross_product, cross_product))
            jacobian_dot = dot_product(cross_product, cross_product_dot)/jacobian
            surface_jacobian_dot(quadrature) = jacobian_dot
            normals_dot(:, quadrature) = real(orientation_sign, dp)* &
                (cross_product_dot/jacobian - &
                cross_product*jacobian_dot/(jacobian*jacobian))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine sample_physical_surface_geometry_jvp

    subroutine sample_physical_surface_geometry_vjp( &
            reference_coordinates, map_points, map_tangents, orientation_sign, &
            physical_points_bar, normals_bar, surface_jacobian_bar, &
            reference_coordinates_bar, map_points_bar, map_tangents_bar, status)
        !! Apply the real reverse product of the physical surface sampler.
        real(dp), intent(in) :: reference_coordinates(:, :)
        real(dp), intent(in) :: map_points(:, :)
        real(dp), intent(in) :: map_tangents(:, :, :)
        integer, intent(in) :: orientation_sign
        real(dp), intent(in) :: physical_points_bar(:, :), normals_bar(:, :)
        real(dp), intent(in) :: surface_jacobian_bar(:)
        real(dp), intent(out) :: reference_coordinates_bar(:, :)
        real(dp), intent(out) :: map_points_bar(:, :), map_tangents_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: cross_product(3), cross_product_bar(3)
        real(dp) :: jacobian
        real(dp) :: normal_bar_dot_cross
        integer :: quadrature, quadrature_count

        reference_coordinates_bar = 0.0_dp
        map_points_bar = 0.0_dp
        map_tangents_bar = 0.0_dp
        call validate_surface_geometry_inputs( &
            reference_coordinates, map_points, map_tangents, orientation_sign, &
            physical_points_bar, normals_bar, surface_jacobian_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(reference_coordinates, 2)
        if (any(shape(physical_points_bar) /= shape(map_points)) .or. &
            any(shape(normals_bar) /= [3, quadrature_count]) .or. &
            size(surface_jacobian_bar) /= quadrature_count .or. &
            any(shape(reference_coordinates_bar) /= shape(reference_coordinates))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "physical surface geometry VJP has incompatible cotangents")
            return
        end if
        if (any(.not. ieee_is_finite(physical_points_bar)) .or. &
            any(.not. ieee_is_finite(normals_bar)) .or. &
            any(.not. ieee_is_finite(surface_jacobian_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "physical surface geometry VJP received non-finite cotangents")
            return
        end if

        map_points_bar = physical_points_bar
        do quadrature = 1, quadrature_count
            cross_product = cross3( &
                map_tangents(:, 1, quadrature), map_tangents(:, 2, quadrature))
            jacobian = sqrt(dot_product(cross_product, cross_product))
            normal_bar_dot_cross = dot_product( &
                normals_bar(:, quadrature), cross_product)
            cross_product_bar = real(orientation_sign, dp)* &
                normals_bar(:, quadrature)/jacobian - &
                real(orientation_sign, dp)*cross_product* &
                normal_bar_dot_cross/(jacobian*jacobian*jacobian) + &
                surface_jacobian_bar(quadrature)*cross_product/jacobian
            map_tangents_bar(:, 1, quadrature) = cross3( &
                map_tangents(:, 2, quadrature), cross_product_bar)
            map_tangents_bar(:, 2, quadrature) = cross3( &
                cross_product_bar, map_tangents(:, 1, quadrature))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine sample_physical_surface_geometry_vjp

    subroutine validate_surface_geometry_inputs( &
            reference_coordinates, map_points, map_tangents, orientation_sign, &
            physical_points, normals, surface_jacobian, status)
        real(dp), intent(in) :: reference_coordinates(:, :), map_points(:, :)
        real(dp), intent(in) :: map_tangents(:, :, :)
        integer, intent(in) :: orientation_sign
        real(dp), intent(in) :: physical_points(:, :), normals(:, :)
        real(dp), intent(in) :: surface_jacobian(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: quadrature, quadrature_count
        real(dp) :: cross_product(3), jacobian

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "physical surface geometry received incompatible arrays")
        quadrature_count = size(reference_coordinates, 2)
        if (size(reference_coordinates, 1) /= 2 .or. quadrature_count < 1) return
        if (size(map_points, 1) /= 3 .or. size(map_points, 2) /= quadrature_count) return
        if (size(map_tangents, 1) /= 3 .or. size(map_tangents, 2) /= 2 .or. &
            size(map_tangents, 3) /= quadrature_count) return
        if (size(physical_points, 1) /= 3 .or. &
            size(physical_points, 2) /= quadrature_count) return
        if (size(normals, 1) /= 3 .or. size(normals, 2) /= quadrature_count) return
        if (size(surface_jacobian) /= quadrature_count) return
        if (orientation_sign /= 1 .and. orientation_sign /= -1) return
        if (any(.not. ieee_is_finite(reference_coordinates)) .or. &
            any(.not. ieee_is_finite(map_points)) .or. &
            any(.not. ieee_is_finite(map_tangents))) return
        if (any(.not. ieee_is_finite(physical_points)) .or. &
            any(.not. ieee_is_finite(normals)) .or. &
            any(.not. ieee_is_finite(surface_jacobian))) return
        do quadrature = 1, quadrature_count
            cross_product = cross3( &
                map_tangents(:, 1, quadrature), map_tangents(:, 2, quadrature))
            jacobian = sqrt(dot_product(cross_product, cross_product))
            if (.not. ieee_is_finite(jacobian) .or. &
                jacobian <= 64.0_dp*epsilon(1.0_dp)) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_surface_geometry_inputs

    pure function cross3(first, second) result(cross_product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: cross_product(3)

        cross_product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross3

end module fortfem_physical_surface_geometry
