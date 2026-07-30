module fortfem_sphere_curved_panel
    !! Exact radial map from an affine surface triangle to a sphere.
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: evaluate_sphere_curved_panel

contains

    pure subroutine evaluate_sphere_curved_panel( &
            vertices, radius, xi, eta, point, tangent_xi, tangent_eta, &
            surface_jacobian, status)
        real(dp), intent(in) :: vertices(3, 3), radius, xi, eta
        real(dp), intent(out) :: point(3), tangent_xi(3), tangent_eta(3)
        real(dp), intent(out) :: surface_jacobian
        integer, intent(out) :: status

        real(dp) :: affine_point(3), derivative_eta(3), derivative_xi(3)
        real(dp) :: norm_affine

        point = 0.0_dp
        tangent_xi = 0.0_dp
        tangent_eta = 0.0_dp
        surface_jacobian = 0.0_dp
        status = 1
        if (radius <= 0.0_dp) return
        if (xi < 0.0_dp .or. eta < 0.0_dp .or. xi + eta > 1.0_dp) return
        derivative_xi = vertices(:, 2) - vertices(:, 1)
        derivative_eta = vertices(:, 3) - vertices(:, 1)
        affine_point = vertices(:, 1) + &
            xi*derivative_xi + eta*derivative_eta
        norm_affine = norm2(affine_point)
        if (norm_affine <= tiny(1.0_dp)) return
        point = radius*affine_point/norm_affine
        tangent_xi = radius*( &
            derivative_xi/norm_affine - affine_point* &
            dot_product(affine_point, derivative_xi)/norm_affine**3)
        tangent_eta = radius*( &
            derivative_eta/norm_affine - affine_point* &
            dot_product(affine_point, derivative_eta)/norm_affine**3)
        surface_jacobian = norm2(cross_product(tangent_xi, tangent_eta))
        if (surface_jacobian <= tiny(1.0_dp)) return
        status = 0
    end subroutine evaluate_sphere_curved_panel

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_sphere_curved_panel
