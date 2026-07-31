module fortfem_sphere_curved_panel
    !! Exact radial map from an affine surface triangle to a sphere.
    use fortfem_generated_sphere_curved_panel, only: &
        generated_sphere_curved_panel
    use fortfem_generated_sphere_curved_panel_jvp, only: &
        generated_sphere_curved_panel_jvp
    use fortfem_generated_sphere_curved_panel_vjp, only: &
        generated_sphere_curved_panel_vjp
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: inv3
    implicit none
    private

    public :: evaluate_sphere_curved_panel
    public :: evaluate_sphere_curved_panel_jvp
    public :: evaluate_sphere_curved_panel_vjp
    public :: invert_sphere_curved_panel

contains

    pure subroutine invert_sphere_curved_panel( &
            vertices, radius, point, xi, eta, status)
        real(dp), intent(in) :: vertices(3, 3), radius, point(3)
        real(dp), intent(out) :: xi, eta
        integer, intent(out) :: status

        real(dp) :: inverse_matrix(3, 3), matrix(3, 3), solution(3)
        real(dp) :: tolerance
        integer :: inverse_status

        xi = 0.0_dp
        eta = 0.0_dp
        status = 1
        if (radius <= 0.0_dp) return
        tolerance = 512.0_dp*epsilon(1.0_dp)*max(1.0_dp, radius)
        if (abs(norm2(point) - radius) > tolerance) return
        matrix(:, 1) = vertices(:, 2) - vertices(:, 1)
        matrix(:, 2) = vertices(:, 3) - vertices(:, 1)
        matrix(:, 3) = -point/radius
        call inv3(matrix, inverse_matrix, inverse_status)
        if (inverse_status /= 0) return
        solution = matmul(inverse_matrix, -vertices(:, 1))
        if (solution(3) <= 0.0_dp) return
        if (solution(1) < -tolerance .or. solution(2) < -tolerance .or. &
            solution(1) + solution(2) > 1.0_dp + tolerance) return
        xi = solution(1)
        eta = solution(2)
        status = 0
    end subroutine invert_sphere_curved_panel

    pure subroutine evaluate_sphere_curved_panel( &
            vertices, radius, xi, eta, point, tangent_xi, tangent_eta, &
            surface_jacobian, status)
        real(dp), intent(in) :: vertices(3, 3), radius, xi, eta
        real(dp), intent(out) :: point(3), tangent_xi(3), tangent_eta(3)
        real(dp), intent(out) :: surface_jacobian
        integer, intent(out) :: status

        real(dp) :: affine_point(3)

        point = 0.0_dp
        tangent_xi = 0.0_dp
        tangent_eta = 0.0_dp
        surface_jacobian = 0.0_dp
        status = 1
        if (radius <= 0.0_dp) return
        if (xi < 0.0_dp) return
        if (eta < 0.0_dp) return
        if (xi + eta > 1.0_dp) return
        affine_point = vertices(:, 1) + &
            xi*(vertices(:, 2) - vertices(:, 1)) + &
            eta*(vertices(:, 3) - vertices(:, 1))
        if (norm2(affine_point) <= tiny(1.0_dp)) return
        call generated_sphere_curved_panel( &
            vertices(1, 1), vertices(2, 1), vertices(3, 1), &
            vertices(1, 2), vertices(2, 2), vertices(3, 2), &
            vertices(1, 3), vertices(2, 3), vertices(3, 3), radius, xi, eta, &
            point(1), point(2), point(3), tangent_xi(1), tangent_xi(2), &
            tangent_xi(3), tangent_eta(1), tangent_eta(2), tangent_eta(3), &
            surface_jacobian)
        if (surface_jacobian <= tiny(1.0_dp)) return
        status = 0
    end subroutine evaluate_sphere_curved_panel

    pure subroutine evaluate_sphere_curved_panel_jvp( &
            vertices, radius, xi, eta, vertices_dot, radius_dot, xi_dot, &
            eta_dot, point_dot, tangent_xi_dot, tangent_eta_dot, &
            surface_jacobian_dot, status)
        real(dp), intent(in) :: vertices(3, 3), radius, xi, eta
        real(dp), intent(in) :: vertices_dot(3, 3), radius_dot, xi_dot, eta_dot
        real(dp), intent(out) :: point_dot(3), tangent_xi_dot(3)
        real(dp), intent(out) :: tangent_eta_dot(3), surface_jacobian_dot
        integer, intent(out) :: status

        real(dp) :: point(3), tangent_xi(3), tangent_eta(3), jacobian

        point_dot = 0.0_dp
        tangent_xi_dot = 0.0_dp
        tangent_eta_dot = 0.0_dp
        surface_jacobian_dot = 0.0_dp
        call evaluate_sphere_curved_panel( &
            vertices, radius, xi, eta, point, tangent_xi, tangent_eta, &
            jacobian, status)
        if (status /= 0) return
        call generated_sphere_curved_panel_jvp( &
            vertices(1, 1), vertices(2, 1), vertices(3, 1), &
            vertices(1, 2), vertices(2, 2), vertices(3, 2), &
            vertices(1, 3), vertices(2, 3), vertices(3, 3), radius, xi, eta, &
            vertices_dot(1, 1), vertices_dot(2, 1), vertices_dot(3, 1), &
            vertices_dot(1, 2), vertices_dot(2, 2), vertices_dot(3, 2), &
            vertices_dot(1, 3), vertices_dot(2, 3), vertices_dot(3, 3), &
            radius_dot, xi_dot, eta_dot, point_dot(1), point_dot(2), &
            point_dot(3), tangent_xi_dot(1), tangent_xi_dot(2), &
            tangent_xi_dot(3), tangent_eta_dot(1), tangent_eta_dot(2), &
            tangent_eta_dot(3), surface_jacobian_dot)
    end subroutine evaluate_sphere_curved_panel_jvp

    pure subroutine evaluate_sphere_curved_panel_vjp( &
            vertices, radius, xi, eta, point_bar, tangent_xi_bar, &
            tangent_eta_bar, surface_jacobian_bar, vertices_bar, radius_bar, &
            xi_bar, eta_bar, status)
        real(dp), intent(in) :: vertices(3, 3), radius, xi, eta
        real(dp), intent(in) :: point_bar(3), tangent_xi_bar(3)
        real(dp), intent(in) :: tangent_eta_bar(3), surface_jacobian_bar
        real(dp), intent(out) :: vertices_bar(3, 3), radius_bar, xi_bar, eta_bar
        integer, intent(out) :: status

        real(dp) :: point(3), tangent_xi(3), tangent_eta(3), jacobian

        vertices_bar = 0.0_dp
        radius_bar = 0.0_dp
        xi_bar = 0.0_dp
        eta_bar = 0.0_dp
        call evaluate_sphere_curved_panel( &
            vertices, radius, xi, eta, point, tangent_xi, tangent_eta, &
            jacobian, status)
        if (status /= 0) return
        call generated_sphere_curved_panel_vjp( &
            vertices(1, 1), vertices(2, 1), vertices(3, 1), &
            vertices(1, 2), vertices(2, 2), vertices(3, 2), &
            vertices(1, 3), vertices(2, 3), vertices(3, 3), radius, xi, eta, &
            point_bar(1), point_bar(2), point_bar(3), tangent_xi_bar(1), &
            tangent_xi_bar(2), tangent_xi_bar(3), tangent_eta_bar(1), &
            tangent_eta_bar(2), tangent_eta_bar(3), surface_jacobian_bar, &
            vertices_bar(1, 1), vertices_bar(2, 1), vertices_bar(3, 1), &
            vertices_bar(1, 2), vertices_bar(2, 2), vertices_bar(3, 2), &
            vertices_bar(1, 3), vertices_bar(2, 3), vertices_bar(3, 3), &
            radius_bar, xi_bar, eta_bar)
    end subroutine evaluate_sphere_curved_panel_vjp

end module fortfem_sphere_curved_panel
