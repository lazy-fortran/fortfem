module fortfem_force_balance_metric_residual
    !! Metric-aware collocation composition for force-balance samples.
    !!
    !! A caller supplies positive reference quadrature weights and the positive
    !! volume Jacobian returned by a geometry map.  The volume contribution is
    !! assembled with ``weight = quadrature_weight * volume_jacobian`` while
    !! boundary and explicitly represented sheet terms retain their own
    !! measures.  This is deliberately only a geometric composition: force
    !! laws, profiles, coordinate conventions, and plasma formats remain
    !! external to FortFEM.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_force_balance_residual, only: &
        assemble_force_balance_residual, &
        assemble_force_balance_residual_jvp, &
        assemble_force_balance_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_metric_force_balance_residual
    public :: assemble_metric_force_balance_residual_jvp
    public :: assemble_metric_force_balance_residual_vjp

contains

    subroutine assemble_metric_force_balance_residual( &
            volume_test, volume_force, quadrature_weights, volume_jacobian, &
            boundary_test, boundary_force, boundary_weights, sheet_test, &
            sheet_force, sheet_weights, residual, status)
        real(dp), intent(in) :: volume_test(:, :, :), volume_force(:, :)
        real(dp), intent(in) :: quadrature_weights(:), volume_jacobian(:)
        real(dp), intent(in) :: boundary_test(:, :, :), boundary_force(:, :)
        real(dp), intent(in) :: boundary_weights(:)
        real(dp), intent(in) :: sheet_test(:, :, :), sheet_force(:, :)
        real(dp), intent(in) :: sheet_weights(:)
        real(dp), intent(out) :: residual(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: volume_weights(:)

        residual = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "metric force residual received incompatible volume measures")
        if (.not. valid_measures(quadrature_weights, volume_jacobian)) return
        allocate(volume_weights(size(quadrature_weights)))
        volume_weights = quadrature_weights*volume_jacobian
        call assemble_force_balance_residual( &
            volume_test, volume_force, volume_weights, boundary_test, boundary_force, &
            boundary_weights, sheet_test, sheet_force, sheet_weights, residual, status)
    end subroutine assemble_metric_force_balance_residual

    subroutine assemble_metric_force_balance_residual_jvp( &
            volume_test, volume_force, quadrature_weights, volume_jacobian, &
            boundary_test, boundary_force, boundary_weights, sheet_test, sheet_force, &
            sheet_weights, volume_test_dot, volume_force_dot, quadrature_weights_dot, &
            volume_jacobian_dot, boundary_test_dot, boundary_force_dot, &
            boundary_weights_dot, &
            sheet_test_dot, sheet_force_dot, sheet_weights_dot, residual_dot, status)
        real(dp), intent(in) :: volume_test(:, :, :), volume_force(:, :)
        real(dp), intent(in) :: quadrature_weights(:), volume_jacobian(:)
        real(dp), intent(in) :: boundary_test(:, :, :), boundary_force(:, :)
        real(dp), intent(in) :: boundary_weights(:)
        real(dp), intent(in) :: sheet_test(:, :, :), sheet_force(:, :)
        real(dp), intent(in) :: sheet_weights(:)
        real(dp), intent(in) :: volume_test_dot(:, :, :), volume_force_dot(:, :)
        real(dp), intent(in) :: quadrature_weights_dot(:), volume_jacobian_dot(:)
        real(dp), intent(in) :: boundary_test_dot(:, :, :), boundary_force_dot(:, :)
        real(dp), intent(in) :: boundary_weights_dot(:)
        real(dp), intent(in) :: sheet_test_dot(:, :, :), sheet_force_dot(:, :)
        real(dp), intent(in) :: sheet_weights_dot(:)
        real(dp), intent(out) :: residual_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: volume_weights(:), volume_weights_dot(:)

        residual_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "metric force residual JVP received incompatible volume measures")
        if (.not. valid_measures(quadrature_weights, volume_jacobian) .or. &
            .not. valid_measure_directions( &
                quadrature_weights_dot, volume_jacobian_dot, size(quadrature_weights))) return
        allocate(volume_weights(size(quadrature_weights)), &
            volume_weights_dot(size(quadrature_weights)))
        volume_weights = quadrature_weights*volume_jacobian
        volume_weights_dot = quadrature_weights_dot*volume_jacobian + &
            quadrature_weights*volume_jacobian_dot
        call assemble_force_balance_residual_jvp( &
            volume_test, volume_force, volume_weights, boundary_test, boundary_force, &
            boundary_weights, sheet_test, sheet_force, sheet_weights, volume_test_dot, &
            volume_force_dot, volume_weights_dot, boundary_test_dot, &
            boundary_force_dot, &
            boundary_weights_dot, sheet_test_dot, sheet_force_dot, sheet_weights_dot, &
            residual_dot, status)
    end subroutine assemble_metric_force_balance_residual_jvp

    subroutine assemble_metric_force_balance_residual_vjp( &
            volume_test, volume_force, quadrature_weights, volume_jacobian, &
            boundary_test, boundary_force, boundary_weights, sheet_test, sheet_force, &
            sheet_weights, residual_bar, volume_test_bar, volume_force_bar, &
            quadrature_weights_bar, volume_jacobian_bar, boundary_test_bar, &
            boundary_force_bar, boundary_weights_bar, sheet_test_bar, sheet_force_bar, &
            sheet_weights_bar, status)
        real(dp), intent(in) :: volume_test(:, :, :), volume_force(:, :)
        real(dp), intent(in) :: quadrature_weights(:), volume_jacobian(:)
        real(dp), intent(in) :: boundary_test(:, :, :), boundary_force(:, :)
        real(dp), intent(in) :: boundary_weights(:)
        real(dp), intent(in) :: sheet_test(:, :, :), sheet_force(:, :)
        real(dp), intent(in) :: sheet_weights(:), residual_bar(:, :)
        real(dp), intent(out) :: volume_test_bar(:, :, :), volume_force_bar(:, :)
        real(dp), intent(out) :: quadrature_weights_bar(:), volume_jacobian_bar(:)
        real(dp), intent(out) :: boundary_test_bar(:, :, :), boundary_force_bar(:, :)
        real(dp), intent(out) :: boundary_weights_bar(:)
        real(dp), intent(out) :: sheet_test_bar(:, :, :), sheet_force_bar(:, :)
        real(dp), intent(out) :: sheet_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: volume_weights(:), volume_weights_bar(:)

        volume_test_bar = 0.0_dp
        volume_force_bar = 0.0_dp
        quadrature_weights_bar = 0.0_dp
        volume_jacobian_bar = 0.0_dp
        boundary_test_bar = 0.0_dp
        boundary_force_bar = 0.0_dp
        boundary_weights_bar = 0.0_dp
        sheet_test_bar = 0.0_dp
        sheet_force_bar = 0.0_dp
        sheet_weights_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "metric force residual VJP received incompatible volume measures")
        if (.not. valid_measures(quadrature_weights, volume_jacobian)) return
        if (size(quadrature_weights_bar) /= size(quadrature_weights) .or. &
            size(volume_jacobian_bar) /= size(volume_jacobian)) return
        allocate(volume_weights(size(quadrature_weights)), &
            volume_weights_bar(size(quadrature_weights)))
        volume_weights = quadrature_weights*volume_jacobian
        call assemble_force_balance_residual_vjp( &
            volume_test, volume_force, volume_weights, boundary_test, boundary_force, &
            boundary_weights, sheet_test, sheet_force, sheet_weights, residual_bar, &
            volume_test_bar, volume_force_bar, volume_weights_bar, boundary_test_bar, &
            boundary_force_bar, boundary_weights_bar, sheet_test_bar, sheet_force_bar, &
            sheet_weights_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_weights_bar = volume_weights_bar*volume_jacobian
        volume_jacobian_bar = volume_weights_bar*quadrature_weights
    end subroutine assemble_metric_force_balance_residual_vjp

    logical function valid_measures(quadrature_weights, volume_jacobian) result(valid)
        real(dp), intent(in) :: quadrature_weights(:), volume_jacobian(:)

        valid = size(quadrature_weights) > 0 .and. &
            size(volume_jacobian) == size(quadrature_weights) .and. &
            all(ieee_is_finite(quadrature_weights)) .and. &
            all(ieee_is_finite(volume_jacobian)) .and. &
            all(quadrature_weights > 0.0_dp) .and. all(volume_jacobian > 0.0_dp)
    end function valid_measures

    logical function valid_measure_directions( &
            quadrature_weights_dot, volume_jacobian_dot, sample_count) result(valid)
        real(dp), intent(in) :: quadrature_weights_dot(:), volume_jacobian_dot(:)
        integer, intent(in) :: sample_count

        valid = size(quadrature_weights_dot) == sample_count .and. &
            size(volume_jacobian_dot) == sample_count .and. &
            all(ieee_is_finite(quadrature_weights_dot)) .and. &
            all(ieee_is_finite(volume_jacobian_dot))
    end function valid_measure_directions

end module fortfem_force_balance_metric_residual
