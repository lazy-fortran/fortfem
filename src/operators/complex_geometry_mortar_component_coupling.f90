module fortfem_complex_geometry_mortar_component_coupling
    !! Neutral complex cross-mass for vector- and tensor-valued traces.
    !!
    !! Trace and component-metric values may be complex, while reference
    !! quadrature weights and surface Jacobians remain real geometric data.
    !! The bilinear contraction intentionally does not conjugate either trace;
    !! callers choose the appropriate test representation before assembly.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_complex_geometry_mortar_component_coupling
    public :: assemble_complex_geometry_mortar_component_coupling_jvp
    public :: assemble_complex_geometry_mortar_component_coupling_vjp

    interface finite_complex
        module procedure finite_complex_matrix
        module procedure finite_complex_rank_three
    end interface finite_complex

contains

    subroutine assemble_complex_geometry_mortar_component_coupling( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix, physical_weights, status)
        complex(dp), intent(in) :: test_trace(:, :, :), trial_trace(:, :, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        complex(dp), intent(in) :: component_metric(:, :, :)
        complex(dp), intent(out) :: matrix(:, :)
        real(dp), intent(out) :: physical_weights(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: q, i, j, a, b

        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        physical_weights = 0.0_dp
        call validate_primal( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix, physical_weights, status)
        if (status%code /= FORTSPARSE_OK) return

        physical_weights = reference_weights*surface_jacobian
        do j = 1, size(trial_trace, 3)
            do i = 1, size(test_trace, 3)
                do q = 1, size(reference_weights)
                    do b = 1, size(trial_trace, 1)
                        do a = 1, size(test_trace, 1)
                            matrix(i, j) = matrix(i, j) + physical_weights(q)* &
                                test_trace(a, q, i)*component_metric(a, b, q)* &
                                trial_trace(b, q, j)
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_complex_geometry_mortar_component_coupling

    subroutine assemble_complex_geometry_mortar_component_coupling_jvp( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, test_trace_dot, trial_trace_dot, &
            reference_weights_dot, surface_jacobian_dot, component_metric_dot, &
            matrix_dot, physical_weights_dot, status)
        !! Apply the complete fixed-topology product-rule JVP.
        complex(dp), intent(in) :: test_trace(:, :, :), trial_trace(:, :, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        complex(dp), intent(in) :: component_metric(:, :, :)
        complex(dp), intent(in) :: test_trace_dot(:, :, :)
        complex(dp), intent(in) :: trial_trace_dot(:, :, :)
        real(dp), intent(in) :: reference_weights_dot(:)
        real(dp), intent(in) :: surface_jacobian_dot(:)
        complex(dp), intent(in) :: component_metric_dot(:, :, :)
        complex(dp), intent(out) :: matrix_dot(:, :)
        real(dp), intent(out) :: physical_weights_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: physical_weights(size(reference_weights))
        integer :: q, i, j, a, b

        matrix_dot = cmplx(0.0_dp, 0.0_dp, dp)
        physical_weights_dot = 0.0_dp
        call validate_primal( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix_dot, physical_weights_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex component mortar JVP received incompatible increments")
        if (.not. valid_direction( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, test_trace_dot, trial_trace_dot, &
            reference_weights_dot, surface_jacobian_dot, &
            component_metric_dot)) return

        physical_weights = reference_weights*surface_jacobian
        physical_weights_dot = reference_weights_dot*surface_jacobian + &
            reference_weights*surface_jacobian_dot
        do j = 1, size(trial_trace, 3)
            do i = 1, size(test_trace, 3)
                do q = 1, size(reference_weights)
                    do b = 1, size(trial_trace, 1)
                        do a = 1, size(test_trace, 1)
                            matrix_dot(i, j) = matrix_dot(i, j) + &
                                physical_weights_dot(q)*test_trace(a, q, i)* &
                                component_metric(a, b, q)*trial_trace(b, q, j) + &
                                physical_weights(q)*( &
                                test_trace_dot(a, q, i)* &
                                component_metric(a, b, q)*trial_trace(b, q, j) + &
                                test_trace(a, q, i)* &
                                component_metric_dot(a, b, q)* &
                                trial_trace(b, q, j) + test_trace(a, q, i)* &
                                component_metric(a, b, q)* &
                                trial_trace_dot(b, q, j))
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_complex_geometry_mortar_component_coupling_jvp

    subroutine assemble_complex_geometry_mortar_component_coupling_vjp( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix_bar, physical_weights_bar, test_trace_bar, &
            trial_trace_bar, reference_weights_bar, surface_jacobian_bar, &
            component_metric_bar, status)
        !! Apply the VJP under complex real-part and real Euclidean pairings.
        complex(dp), intent(in) :: test_trace(:, :, :), trial_trace(:, :, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        complex(dp), intent(in) :: component_metric(:, :, :)
        complex(dp), intent(in) :: matrix_bar(:, :)
        real(dp), intent(in) :: physical_weights_bar(:)
        complex(dp), intent(out) :: test_trace_bar(:, :, :)
        complex(dp), intent(out) :: trial_trace_bar(:, :, :)
        real(dp), intent(out) :: reference_weights_bar(:)
        real(dp), intent(out) :: surface_jacobian_bar(:)
        complex(dp), intent(out) :: component_metric_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: physical_weights(size(reference_weights))
        real(dp) :: total_weight_bar(size(reference_weights))
        complex(dp) :: kernel
        integer :: q, i, j, a, b

        test_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        trial_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        reference_weights_bar = 0.0_dp
        surface_jacobian_bar = 0.0_dp
        component_metric_bar = cmplx(0.0_dp, 0.0_dp, dp)
        call validate_primal( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix_bar, physical_weights_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex component mortar VJP received incompatible cotangents")
        if (.not. valid_adjoint_outputs( &
            test_trace, trial_trace, reference_weights, component_metric, &
            test_trace_bar, trial_trace_bar, reference_weights_bar, &
            surface_jacobian_bar, component_metric_bar)) return
        if (.not. finite_complex(matrix_bar)) return
        if (.not. all(ieee_is_finite(physical_weights_bar))) return

        physical_weights = reference_weights*surface_jacobian
        total_weight_bar = physical_weights_bar
        do j = 1, size(trial_trace, 3)
            do i = 1, size(test_trace, 3)
                do q = 1, size(reference_weights)
                    do b = 1, size(trial_trace, 1)
                        do a = 1, size(test_trace, 1)
                            kernel = test_trace(a, q, i)* &
                                component_metric(a, b, q)*trial_trace(b, q, j)
                            test_trace_bar(a, q, i) = &
                                test_trace_bar(a, q, i) + physical_weights(q)* &
                                matrix_bar(i, j)*conjg( &
                                component_metric(a, b, q)*trial_trace(b, q, j))
                            trial_trace_bar(b, q, j) = &
                                trial_trace_bar(b, q, j) + physical_weights(q)* &
                                matrix_bar(i, j)*conjg( &
                                test_trace(a, q, i)*component_metric(a, b, q))
                            component_metric_bar(a, b, q) = &
                                component_metric_bar(a, b, q) + &
                                physical_weights(q)*matrix_bar(i, j)*conjg( &
                                test_trace(a, q, i)*trial_trace(b, q, j))
                            total_weight_bar(q) = total_weight_bar(q) + real( &
                                conjg(matrix_bar(i, j))*kernel, dp)
                        end do
                    end do
                end do
            end do
        end do
        reference_weights_bar = total_weight_bar*surface_jacobian
        surface_jacobian_bar = total_weight_bar*reference_weights
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_complex_geometry_mortar_component_coupling_vjp

    subroutine validate_primal( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix, physical_weights, status)
        complex(dp), intent(in) :: test_trace(:, :, :), trial_trace(:, :, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        complex(dp), intent(in) :: component_metric(:, :, :), matrix(:, :)
        real(dp), intent(in) :: physical_weights(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: quadrature_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex component mortar received incompatible arrays")
        if (size(test_trace, 1) < 1 .or. size(test_trace, 2) < 1 .or. &
            size(test_trace, 3) < 1) return
        if (size(trial_trace, 1) /= size(test_trace, 1)) return
        quadrature_count = size(test_trace, 2)
        if (size(trial_trace, 2) /= quadrature_count .or. &
            size(trial_trace, 3) < 1) return
        if (size(reference_weights) /= quadrature_count .or. &
            size(surface_jacobian) /= quadrature_count) return
        if (size(component_metric, 1) /= size(test_trace, 1) .or. &
            size(component_metric, 2) /= size(trial_trace, 1) .or. &
            size(component_metric, 3) /= quadrature_count) return
        if (size(matrix, 1) /= size(test_trace, 3) .or. &
            size(matrix, 2) /= size(trial_trace, 3)) return
        if (size(physical_weights) /= quadrature_count) return
        if (.not. finite_complex(test_trace) .or. &
            .not. finite_complex(trial_trace) .or. &
            .not. finite_complex(component_metric)) return
        if (.not. all(ieee_is_finite(reference_weights)) .or. &
            .not. all(ieee_is_finite(surface_jacobian))) return
        if (any(reference_weights <= 0.0_dp) .or. &
            any(surface_jacobian <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_primal

    logical function valid_direction( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, test_trace_dot, trial_trace_dot, &
            reference_weights_dot, surface_jacobian_dot, &
            component_metric_dot) result(valid)
        complex(dp), intent(in) :: test_trace(:, :, :), trial_trace(:, :, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        complex(dp), intent(in) :: component_metric(:, :, :)
        complex(dp), intent(in) :: test_trace_dot(:, :, :)
        complex(dp), intent(in) :: trial_trace_dot(:, :, :)
        real(dp), intent(in) :: reference_weights_dot(:)
        real(dp), intent(in) :: surface_jacobian_dot(:)
        complex(dp), intent(in) :: component_metric_dot(:, :, :)

        valid = .false.
        if (any(shape(test_trace_dot) /= shape(test_trace))) return
        if (any(shape(trial_trace_dot) /= shape(trial_trace))) return
        if (size(reference_weights_dot) /= size(reference_weights)) return
        if (size(surface_jacobian_dot) /= size(surface_jacobian)) return
        if (any(shape(component_metric_dot) /= shape(component_metric))) return
        if (.not. finite_complex(test_trace_dot) .or. &
            .not. finite_complex(trial_trace_dot) .or. &
            .not. finite_complex(component_metric_dot)) return
        if (.not. all(ieee_is_finite(reference_weights_dot)) .or. &
            .not. all(ieee_is_finite(surface_jacobian_dot))) return
        valid = .true.
    end function valid_direction

    logical function valid_adjoint_outputs( &
            test_trace, trial_trace, reference_weights, component_metric, &
            test_trace_bar, trial_trace_bar, reference_weights_bar, &
            surface_jacobian_bar, component_metric_bar) result(valid)
        complex(dp), intent(in) :: test_trace(:, :, :), trial_trace(:, :, :)
        real(dp), intent(in) :: reference_weights(:)
        complex(dp), intent(in) :: component_metric(:, :, :)
        complex(dp), intent(in) :: test_trace_bar(:, :, :)
        complex(dp), intent(in) :: trial_trace_bar(:, :, :)
        real(dp), intent(in) :: reference_weights_bar(:)
        real(dp), intent(in) :: surface_jacobian_bar(:)
        complex(dp), intent(in) :: component_metric_bar(:, :, :)

        valid = all(shape(test_trace_bar) == shape(test_trace)) .and. &
            all(shape(trial_trace_bar) == shape(trial_trace)) .and. &
            size(reference_weights_bar) == size(reference_weights) .and. &
            size(surface_jacobian_bar) == size(reference_weights) .and. &
            all(shape(component_metric_bar) == shape(component_metric))
    end function valid_adjoint_outputs

    logical function finite_complex_matrix(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_matrix

    logical function finite_complex_rank_three(values) result(valid)
        complex(dp), intent(in) :: values(:, :, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_rank_three

end module fortfem_complex_geometry_mortar_component_coupling
