module fortfem_helmholtz_panel_pair_3d
    use fortfem_generated_helmholtz_single_layer_integrand, only: &
        generated_helmholtz_single_layer_integrand
    use fortfem_generated_helmholtz_single_layer_integrand_jvp, only: &
        generated_helmholtz_single_layer_integrand_jvp
    use fortfem_generated_helmholtz_single_layer_integrand_vjp, only: &
        generated_helmholtz_single_layer_integrand_vjp
    use fortfem_kinds, only: dp
    use fortfem_surface_triangle_geometry_3d, only: &
        evaluate_surface_triangle_geometry_3d, &
        evaluate_surface_triangle_geometry_3d_jvp, &
        evaluate_surface_triangle_geometry_3d_vjp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    private

    public :: integrate_helmholtz_single_layer_regular_panel_pair_p0_3d
    public :: integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_jvp
    public :: integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_vjp

contains

    subroutine integrate_helmholtz_single_layer_regular_panel_pair_p0_3d( &
            first_vertices, second_vertices, wave_number, quadrature_degree, &
            value, status)
        real(dp), intent(in) :: first_vertices(3, 3), second_vertices(3, 3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: quadrature_degree
        complex(dp), intent(out) :: value
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, first_point(3), kernel_scale
        real(dp) :: normal(3), second_jacobian, second_point(3)
        real(dp) :: term_real, term_imag
        integer :: first_node, quadrature_status, second_node

        value = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (wave_number < 0.0_dp) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do first_node = 1, size(weights)
            call evaluate_surface_triangle_geometry_3d( &
                first_vertices, xi(first_node), eta(first_node), first_point, &
                first_jacobian, normal, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_surface_triangle_geometry_3d( &
                    second_vertices, xi(second_node), eta(second_node), &
                    second_point, second_jacobian, normal, status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    64.0_dp*epsilon(1.0_dp)) return
                kernel_scale = weights(first_node)*weights(second_node)/ &
                    (4.0_dp*acos(-1.0_dp))
                call generated_helmholtz_single_layer_integrand( &
                    first_point(1), first_point(2), first_point(3), &
                    second_point(1), second_point(2), second_point(3), &
                    first_jacobian, second_jacobian, wave_number, kernel_scale, &
                    term_real, term_imag)
                value = value + cmplx(term_real, term_imag, dp)
            end do
        end do
        status = 0
    end subroutine integrate_helmholtz_single_layer_regular_panel_pair_p0_3d

    subroutine integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_jvp( &
            first_vertices, second_vertices, wave_number, quadrature_degree, &
            first_vertices_dot, second_vertices_dot, wave_number_dot, &
            value_dot, status)
        real(dp), intent(in) :: first_vertices(3, 3), second_vertices(3, 3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: quadrature_degree
        real(dp), intent(in) :: first_vertices_dot(3, 3)
        real(dp), intent(in) :: second_vertices_dot(3, 3), wave_number_dot
        complex(dp), intent(out) :: value_dot
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, first_jacobian_dot, first_point(3)
        real(dp) :: first_point_dot(3), kernel_scale, normal(3), normal_dot(3)
        real(dp) :: second_jacobian, second_jacobian_dot, second_point(3)
        real(dp) :: second_point_dot(3), term_real_dot, term_imag_dot
        integer :: first_node, quadrature_status, second_node

        value_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (wave_number < 0.0_dp) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        do first_node = 1, size(weights)
            call evaluate_surface_triangle_geometry_3d( &
                first_vertices, xi(first_node), eta(first_node), first_point, &
                first_jacobian, normal, status)
            if (status /= 0) return
            call evaluate_surface_triangle_geometry_3d_jvp( &
                first_vertices, xi(first_node), eta(first_node), &
                first_vertices_dot, 0.0_dp, 0.0_dp, first_point_dot, &
                first_jacobian_dot, normal_dot, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_surface_triangle_geometry_3d( &
                    second_vertices, xi(second_node), eta(second_node), &
                    second_point, second_jacobian, normal, status)
                if (status /= 0) return
                call evaluate_surface_triangle_geometry_3d_jvp( &
                    second_vertices, xi(second_node), eta(second_node), &
                    second_vertices_dot, 0.0_dp, 0.0_dp, second_point_dot, &
                    second_jacobian_dot, normal_dot, status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    64.0_dp*epsilon(1.0_dp)) return
                kernel_scale = weights(first_node)*weights(second_node)/ &
                    (4.0_dp*acos(-1.0_dp))
                call generated_helmholtz_single_layer_integrand_jvp( &
                    first_point(1), first_point(2), first_point(3), &
                    second_point(1), second_point(2), second_point(3), &
                    first_jacobian, second_jacobian, wave_number, kernel_scale, &
                    first_point_dot(1), first_point_dot(2), &
                    first_point_dot(3), second_point_dot(1), &
                    second_point_dot(2), second_point_dot(3), &
                    first_jacobian_dot, second_jacobian_dot, wave_number_dot, &
                    term_real_dot, term_imag_dot)
                value_dot = value_dot + &
                    cmplx(term_real_dot, term_imag_dot, dp)
            end do
        end do
        status = 0
    end subroutine integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_jvp

    subroutine integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_vjp( &
            first_vertices, second_vertices, wave_number, quadrature_degree, &
            value_bar, first_vertices_bar, second_vertices_bar, &
            wave_number_bar, status)
        real(dp), intent(in) :: first_vertices(3, 3), second_vertices(3, 3)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: quadrature_degree
        complex(dp), intent(in) :: value_bar
        real(dp), intent(out) :: first_vertices_bar(3, 3)
        real(dp), intent(out) :: second_vertices_bar(3, 3), wave_number_bar
        integer, intent(out) :: status

        real(dp), allocatable :: eta(:), first_jacobian_bar(:)
        real(dp), allocatable :: first_point_bar(:, :)
        real(dp), allocatable :: second_jacobian_bar(:)
        real(dp), allocatable :: second_point_bar(:, :), weights(:), xi(:)
        real(dp) :: first_jacobian, first_point(3), kernel_scale
        real(dp) :: normal(3), second_jacobian, second_point(3)
        real(dp) :: local_first_jacobian_bar, local_first_point_bar(3)
        real(dp) :: local_second_jacobian_bar, local_second_point_bar(3)
        real(dp) :: local_vertices_bar(3, 3), local_wave_number_bar
        real(dp) :: xi_bar, eta_bar
        integer :: first_node, quadrature_status, second_node

        first_vertices_bar = 0.0_dp
        second_vertices_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        status = 1
        if (wave_number < 0.0_dp) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        allocate(first_point_bar(3, size(weights)), source=0.0_dp)
        allocate(second_point_bar(3, size(weights)), source=0.0_dp)
        allocate(first_jacobian_bar(size(weights)), source=0.0_dp)
        allocate(second_jacobian_bar(size(weights)), source=0.0_dp)
        do first_node = 1, size(weights)
            call evaluate_surface_triangle_geometry_3d( &
                first_vertices, xi(first_node), eta(first_node), first_point, &
                first_jacobian, normal, status)
            if (status /= 0) return
            do second_node = 1, size(weights)
                call evaluate_surface_triangle_geometry_3d( &
                    second_vertices, xi(second_node), eta(second_node), &
                    second_point, second_jacobian, normal, status)
                if (status /= 0) return
                if (norm2(first_point - second_point) <= &
                    64.0_dp*epsilon(1.0_dp)) return
                kernel_scale = weights(first_node)*weights(second_node)/ &
                    (4.0_dp*acos(-1.0_dp))
                call generated_helmholtz_single_layer_integrand_vjp( &
                    first_point(1), first_point(2), first_point(3), &
                    second_point(1), second_point(2), second_point(3), &
                    first_jacobian, second_jacobian, wave_number, kernel_scale, &
                    real(value_bar, dp), aimag(value_bar), &
                    local_first_point_bar(1), local_first_point_bar(2), &
                    local_first_point_bar(3), local_second_point_bar(1), &
                    local_second_point_bar(2), local_second_point_bar(3), &
                    local_first_jacobian_bar, local_second_jacobian_bar, &
                    local_wave_number_bar)
                first_point_bar(:, first_node) = &
                    first_point_bar(:, first_node) + local_first_point_bar
                second_point_bar(:, second_node) = &
                    second_point_bar(:, second_node) + local_second_point_bar
                first_jacobian_bar(first_node) = &
                    first_jacobian_bar(first_node) + local_first_jacobian_bar
                second_jacobian_bar(second_node) = &
                    second_jacobian_bar(second_node) + local_second_jacobian_bar
                wave_number_bar = wave_number_bar + local_wave_number_bar
            end do
        end do
        do first_node = 1, size(weights)
            call evaluate_surface_triangle_geometry_3d_vjp( &
                first_vertices, xi(first_node), eta(first_node), &
                first_point_bar(:, first_node), &
                first_jacobian_bar(first_node), [0.0_dp, 0.0_dp, 0.0_dp], &
                local_vertices_bar, xi_bar, eta_bar, status)
            if (status /= 0) return
            first_vertices_bar = first_vertices_bar + local_vertices_bar
            call evaluate_surface_triangle_geometry_3d_vjp( &
                second_vertices, xi(first_node), eta(first_node), &
                second_point_bar(:, first_node), &
                second_jacobian_bar(first_node), [0.0_dp, 0.0_dp, 0.0_dp], &
                local_vertices_bar, xi_bar, eta_bar, status)
            if (status /= 0) return
            second_vertices_bar = second_vertices_bar + local_vertices_bar
        end do
        status = 0
    end subroutine integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_vjp

end module fortfem_helmholtz_panel_pair_3d
