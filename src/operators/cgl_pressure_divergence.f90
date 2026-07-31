module fortfem_cgl_pressure_divergence
    !! Generated product-rule contraction for the CGL pressure force.
    !!
    !! The input `direction_gradient(i,j)` is the physical derivative
    !! d(direction_i)/d(x_j).  The generated expression evaluates
    !! div(P), where P is the gyrotropic CGL pressure tensor.  This keeps the
    !! constitutive tensor and its force contraction as separate contracts.
    use fortfem_generated_cgl_pressure_divergence, only: &
        generated_cgl_pressure_divergence
    use fortfem_generated_cgl_pressure_divergence_jvp, only: &
        generated_cgl_pressure_divergence_jvp
    use fortfem_generated_cgl_pressure_divergence_vjp, only: &
        generated_cgl_pressure_divergence_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    real(dp), parameter :: unit_tolerance = 1.0e-12_dp

    public :: evaluate_cgl_pressure_divergence
    public :: evaluate_cgl_pressure_divergence_jvp
    public :: evaluate_cgl_pressure_divergence_vjp

contains

    pure subroutine evaluate_cgl_pressure_divergence( &
            p_parallel, p_perpendicular, unit_direction, parallel_gradient, &
            perpendicular_gradient, direction_gradient, force_divergence, &
            status)
        real(dp), intent(in) :: p_parallel, p_perpendicular
        real(dp), intent(in) :: unit_direction(:)
        real(dp), intent(in) :: parallel_gradient(:), perpendicular_gradient(:)
        real(dp), intent(in) :: direction_gradient(:, :)
        real(dp), intent(out) :: force_divergence(:)
        type(fortsparse_status_t), intent(out) :: status

        force_divergence = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "CGL pressure divergence received incompatible arrays")
        if (.not. valid_shapes( &
            unit_direction, parallel_gradient, perpendicular_gradient, &
            direction_gradient, force_divergence)) return
        if (.not. is_unit_direction(unit_direction)) return
        call generated_cgl_pressure_divergence( &
            p_parallel, p_perpendicular, unit_direction(1), unit_direction(2), &
            unit_direction(3), parallel_gradient(1), parallel_gradient(2), &
            parallel_gradient(3), perpendicular_gradient(1), &
            perpendicular_gradient(2), perpendicular_gradient(3), &
            direction_gradient(1, 1), direction_gradient(1, 2), &
            direction_gradient(1, 3), direction_gradient(2, 1), &
            direction_gradient(2, 2), direction_gradient(2, 3), &
            direction_gradient(3, 1), direction_gradient(3, 2), &
            direction_gradient(3, 3), force_divergence(1), force_divergence(2), &
            force_divergence(3))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_cgl_pressure_divergence

    pure subroutine evaluate_cgl_pressure_divergence_jvp( &
            p_parallel, p_perpendicular, unit_direction, parallel_gradient, &
            perpendicular_gradient, direction_gradient, p_parallel_dot, &
            p_perpendicular_dot, direction_dot, parallel_gradient_dot, &
            perpendicular_gradient_dot, direction_gradient_dot, force_dot, &
            status)
        real(dp), intent(in) :: p_parallel, p_perpendicular
        real(dp), intent(in) :: unit_direction(:)
        real(dp), intent(in) :: parallel_gradient(:), perpendicular_gradient(:)
        real(dp), intent(in) :: direction_gradient(:, :)
        real(dp), intent(in) :: p_parallel_dot, p_perpendicular_dot
        real(dp), intent(in) :: direction_dot(:), parallel_gradient_dot(:)
        real(dp), intent(in) :: perpendicular_gradient_dot(:)
        real(dp), intent(in) :: direction_gradient_dot(:, :)
        real(dp), intent(out) :: force_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        force_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "CGL pressure divergence JVP received incompatible arrays")
        if (.not. valid_shapes( &
            unit_direction, parallel_gradient, perpendicular_gradient, &
            direction_gradient, force_dot)) return
        if (size(direction_dot) /= 3) return
        if (size(parallel_gradient_dot) /= 3) return
        if (size(perpendicular_gradient_dot) /= 3) return
        if (size(direction_gradient_dot, 1) /= 3) return
        if (size(direction_gradient_dot, 2) /= 3) return
        if (.not. is_unit_direction(unit_direction)) return
        call generated_cgl_pressure_divergence_jvp( &
            p_parallel, p_perpendicular, unit_direction(1), unit_direction(2), &
            unit_direction(3), parallel_gradient(1), parallel_gradient(2), &
            parallel_gradient(3), perpendicular_gradient(1), &
            perpendicular_gradient(2), perpendicular_gradient(3), &
            direction_gradient(1, 1), direction_gradient(1, 2), &
            direction_gradient(1, 3), direction_gradient(2, 1), &
            direction_gradient(2, 2), direction_gradient(2, 3), &
            direction_gradient(3, 1), direction_gradient(3, 2), &
            direction_gradient(3, 3), p_parallel_dot, p_perpendicular_dot, &
            direction_dot(1), direction_dot(2), direction_dot(3), &
            parallel_gradient_dot(1), parallel_gradient_dot(2), &
            parallel_gradient_dot(3), perpendicular_gradient_dot(1), &
            perpendicular_gradient_dot(2), perpendicular_gradient_dot(3), &
            direction_gradient_dot(1, 1), direction_gradient_dot(1, 2), &
            direction_gradient_dot(1, 3), direction_gradient_dot(2, 1), &
            direction_gradient_dot(2, 2), direction_gradient_dot(2, 3), &
            direction_gradient_dot(3, 1), direction_gradient_dot(3, 2), &
            direction_gradient_dot(3, 3), force_dot(1), force_dot(2), &
            force_dot(3))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_cgl_pressure_divergence_jvp

    pure subroutine evaluate_cgl_pressure_divergence_vjp( &
            p_parallel, p_perpendicular, unit_direction, parallel_gradient, &
            perpendicular_gradient, direction_gradient, force_bar, &
            p_parallel_bar, p_perpendicular_bar, direction_bar, &
            parallel_gradient_bar, perpendicular_gradient_bar, &
            direction_gradient_bar, status)
        real(dp), intent(in) :: p_parallel, p_perpendicular
        real(dp), intent(in) :: unit_direction(:)
        real(dp), intent(in) :: parallel_gradient(:), perpendicular_gradient(:)
        real(dp), intent(in) :: direction_gradient(:, :), force_bar(:)
        real(dp), intent(out) :: p_parallel_bar, p_perpendicular_bar
        real(dp), intent(out) :: direction_bar(:), parallel_gradient_bar(:)
        real(dp), intent(out) :: perpendicular_gradient_bar(:)
        real(dp), intent(out) :: direction_gradient_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        p_parallel_bar = 0.0_dp
        p_perpendicular_bar = 0.0_dp
        direction_bar = 0.0_dp
        parallel_gradient_bar = 0.0_dp
        perpendicular_gradient_bar = 0.0_dp
        direction_gradient_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "CGL pressure divergence VJP received incompatible arrays")
        if (.not. valid_shapes( &
            unit_direction, parallel_gradient, perpendicular_gradient, &
            direction_gradient, force_bar)) return
        if (size(direction_bar) /= 3) return
        if (size(parallel_gradient_bar) /= 3) return
        if (size(perpendicular_gradient_bar) /= 3) return
        if (size(direction_gradient_bar, 1) /= 3) return
        if (size(direction_gradient_bar, 2) /= 3) return
        if (.not. is_unit_direction(unit_direction)) return
        call generated_cgl_pressure_divergence_vjp( &
            p_parallel, p_perpendicular, unit_direction(1), unit_direction(2), &
            unit_direction(3), parallel_gradient(1), parallel_gradient(2), &
            parallel_gradient(3), perpendicular_gradient(1), &
            perpendicular_gradient(2), perpendicular_gradient(3), &
            direction_gradient(1, 1), direction_gradient(1, 2), &
            direction_gradient(1, 3), direction_gradient(2, 1), &
            direction_gradient(2, 2), direction_gradient(2, 3), &
            direction_gradient(3, 1), direction_gradient(3, 2), &
            direction_gradient(3, 3), force_bar(1), force_bar(2), force_bar(3), &
            p_parallel_bar, p_perpendicular_bar, direction_bar(1), &
            direction_bar(2), direction_bar(3), parallel_gradient_bar(1), &
            parallel_gradient_bar(2), parallel_gradient_bar(3), &
            perpendicular_gradient_bar(1), perpendicular_gradient_bar(2), &
            perpendicular_gradient_bar(3), direction_gradient_bar(1, 1), &
            direction_gradient_bar(1, 2), direction_gradient_bar(1, 3), &
            direction_gradient_bar(2, 1), direction_gradient_bar(2, 2), &
            direction_gradient_bar(2, 3), direction_gradient_bar(3, 1), &
            direction_gradient_bar(3, 2), direction_gradient_bar(3, 3))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_cgl_pressure_divergence_vjp

    pure logical function valid_shapes( &
            unit_direction, parallel_gradient, perpendicular_gradient, &
            direction_gradient, force)
        real(dp), intent(in) :: unit_direction(:), parallel_gradient(:)
        real(dp), intent(in) :: perpendicular_gradient(:)
        real(dp), intent(in) :: direction_gradient(:, :), force(:)

        valid_shapes = size(unit_direction) == 3 .and. &
            size(parallel_gradient) == 3 .and. &
            size(perpendicular_gradient) == 3 .and. &
            size(direction_gradient, 1) == 3 .and. &
            size(direction_gradient, 2) == 3 .and. size(force) == 3
    end function valid_shapes

    pure logical function is_unit_direction(direction)
        real(dp), intent(in) :: direction(:)

        is_unit_direction = size(direction) == 3 .and. &
            abs(dot_product(direction, direction) - 1.0_dp) <= unit_tolerance
    end function is_unit_direction

end module fortfem_cgl_pressure_divergence
