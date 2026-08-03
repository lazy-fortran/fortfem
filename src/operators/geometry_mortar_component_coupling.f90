module fortfem_geometry_mortar_component_coupling
    !! Geometry-aware cross-mass for vector- and tensor-valued traces.
    !!
    !! The first index of each trace is a physical component (Cartesian,
    !! covariant, contravariant, or a flattened tensor component), the second
    !! index is a fixed reference quadrature row, and the third index is a
    !! local trace degree of freedom.  A caller supplies the component metric
    !! at each quadrature row.  This keeps Piola, NURBS, Fourier, and material
    !! maps interchangeable while providing one differentiable mortar
    !! contraction for vector and tensor IGA patches.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_geometry_mortar_component_coupling
    public :: assemble_geometry_mortar_component_coupling_jvp
    public :: assemble_geometry_mortar_component_coupling_vjp

contains

    subroutine assemble_geometry_mortar_component_coupling( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix, physical_weights, status)
        !! Assemble a physical-measure component cross-mass.
        real(dp), intent(in) :: test_trace(:, :, :), trial_trace(:, :, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        real(dp), intent(in) :: component_metric(:, :, :)
        real(dp), intent(out) :: matrix(:, :), physical_weights(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof
        integer :: test_component, trial_component

        matrix = 0.0_dp
        physical_weights = 0.0_dp
        call validate_component_inputs( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix, physical_weights, status)
        if (status%code /= FORTSPARSE_OK) return

        quadrature_count = size(test_trace, 2)
        test_count = size(test_trace, 3)
        trial_count = size(trial_trace, 3)
        physical_weights = reference_weights*surface_jacobian
        do quadrature = 1, quadrature_count
            do test_dof = 1, test_count
                do trial_dof = 1, trial_count
                    do test_component = 1, size(test_trace, 1)
                        do trial_component = 1, size(trial_trace, 1)
                            matrix(test_dof, trial_dof) = &
                                matrix(test_dof, trial_dof) + &
                                physical_weights(quadrature)* &
                                test_trace(test_component, quadrature, test_dof)* &
                                component_metric(test_component, trial_component, &
                                quadrature)* &
                                trial_trace(trial_component, quadrature, trial_dof)
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_geometry_mortar_component_coupling

    subroutine assemble_geometry_mortar_component_coupling_jvp( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, test_trace_dot, trial_trace_dot, &
            reference_weights_dot, surface_jacobian_dot, component_metric_dot, &
            matrix_dot, physical_weights_dot, status)
        !! Apply the fixed-topology product-rule JVP.
        real(dp), intent(in) :: test_trace(:, :, :), trial_trace(:, :, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        real(dp), intent(in) :: component_metric(:, :, :)
        real(dp), intent(in) :: test_trace_dot(:, :, :), trial_trace_dot(:, :, :)
        real(dp), intent(in) :: reference_weights_dot(:), surface_jacobian_dot(:)
        real(dp), intent(in) :: component_metric_dot(:, :, :)
        real(dp), intent(out) :: matrix_dot(:, :), physical_weights_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: physical_weights(:)
        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof
        integer :: test_component, trial_component

        matrix_dot = 0.0_dp
        physical_weights_dot = 0.0_dp
        call validate_component_inputs( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix_dot, physical_weights_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_trace, 2)
        test_count = size(test_trace, 3)
        trial_count = size(trial_trace, 3)
        if (size(test_trace_dot, 1) /= size(test_trace, 1) .or. &
            size(test_trace_dot, 2) /= quadrature_count .or. &
            size(test_trace_dot, 3) /= test_count .or. &
            size(trial_trace_dot, 1) /= size(trial_trace, 1) .or. &
            size(trial_trace_dot, 2) /= quadrature_count .or. &
            size(trial_trace_dot, 3) /= trial_count .or. &
            size(reference_weights_dot) /= quadrature_count .or. &
            size(surface_jacobian_dot) /= quadrature_count .or. &
            size(component_metric_dot, 1) /= size(component_metric, 1) .or. &
            size(component_metric_dot, 2) /= size(component_metric, 2) .or. &
            size(component_metric_dot, 3) /= quadrature_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "component mortar JVP has incompatible increment shapes")
            return
        end if
        if (any(.not. ieee_is_finite(test_trace_dot)) .or. &
            any(.not. ieee_is_finite(trial_trace_dot)) .or. &
            any(.not. ieee_is_finite(reference_weights_dot)) .or. &
            any(.not. ieee_is_finite(surface_jacobian_dot)) .or. &
            any(.not. ieee_is_finite(component_metric_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "component mortar JVP received non-finite increments")
            return
        end if

        allocate(physical_weights(quadrature_count))
        physical_weights = reference_weights*surface_jacobian
        physical_weights_dot = reference_weights_dot*surface_jacobian + &
            reference_weights*surface_jacobian_dot
        do quadrature = 1, quadrature_count
            do test_dof = 1, test_count
                do trial_dof = 1, trial_count
                    do test_component = 1, size(test_trace, 1)
                        do trial_component = 1, size(trial_trace, 1)
                            matrix_dot(test_dof, trial_dof) = &
                                matrix_dot(test_dof, trial_dof) + &
                                physical_weights_dot(quadrature)* &
                                test_trace(test_component, quadrature, test_dof)* &
                                component_metric(test_component, trial_component, &
                                quadrature)* &
                                trial_trace(trial_component, quadrature, trial_dof) + &
                                physical_weights(quadrature)* &
                                test_trace_dot(test_component, quadrature, test_dof)* &
                                component_metric(test_component, trial_component, &
                                quadrature)* &
                                trial_trace(trial_component, quadrature, trial_dof) + &
                                physical_weights(quadrature)* &
                                test_trace(test_component, quadrature, test_dof)* &
                                component_metric_dot(test_component, trial_component, &
                                quadrature)* &
                                trial_trace(trial_component, quadrature, trial_dof) + &
                                physical_weights(quadrature)* &
                                test_trace(test_component, quadrature, test_dof)* &
                                component_metric(test_component, trial_component, &
                                quadrature)* &
                                trial_trace_dot(trial_component, quadrature, trial_dof)
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_geometry_mortar_component_coupling_jvp

    subroutine assemble_geometry_mortar_component_coupling_vjp( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix_bar, physical_weights_bar, test_trace_bar, &
            trial_trace_bar, reference_weights_bar, surface_jacobian_bar, &
            component_metric_bar, status)
        !! Apply the real reverse product of the component cross-mass.
        real(dp), intent(in) :: test_trace(:, :, :), trial_trace(:, :, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        real(dp), intent(in) :: component_metric(:, :, :)
        real(dp), intent(in) :: matrix_bar(:, :), physical_weights_bar(:)
        real(dp), intent(out) :: test_trace_bar(:, :, :), trial_trace_bar(:, :, :)
        real(dp), intent(out) :: reference_weights_bar(:), surface_jacobian_bar(:)
        real(dp), intent(out) :: component_metric_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: physical_weights(:), total_weight_bar(:)
        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof
        integer :: test_component, trial_component

        test_trace_bar = 0.0_dp
        trial_trace_bar = 0.0_dp
        reference_weights_bar = 0.0_dp
        surface_jacobian_bar = 0.0_dp
        component_metric_bar = 0.0_dp
        call validate_component_inputs( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix_bar, physical_weights_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_trace, 2)
        test_count = size(test_trace, 3)
        trial_count = size(trial_trace, 3)
        if (size(test_trace_bar, 1) /= size(test_trace, 1) .or. &
            size(test_trace_bar, 2) /= quadrature_count .or. &
            size(test_trace_bar, 3) /= test_count .or. &
            size(trial_trace_bar, 1) /= size(trial_trace, 1) .or. &
            size(trial_trace_bar, 2) /= quadrature_count .or. &
            size(trial_trace_bar, 3) /= trial_count .or. &
            size(reference_weights_bar) /= quadrature_count .or. &
            size(surface_jacobian_bar) /= quadrature_count .or. &
            size(component_metric_bar, 1) /= size(component_metric, 1) .or. &
            size(component_metric_bar, 2) /= size(component_metric, 2) .or. &
            size(component_metric_bar, 3) /= quadrature_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "component mortar VJP has incompatible cotangent shapes")
            return
        end if

        allocate(physical_weights(quadrature_count), total_weight_bar(quadrature_count))
        physical_weights = reference_weights*surface_jacobian
        total_weight_bar = physical_weights_bar
        do quadrature = 1, quadrature_count
            do test_dof = 1, test_count
                do trial_dof = 1, trial_count
                    do test_component = 1, size(test_trace, 1)
                        do trial_component = 1, size(trial_trace, 1)
                            test_trace_bar(test_component, quadrature, test_dof) = &
                                test_trace_bar(test_component, quadrature, test_dof) + &
                                physical_weights(quadrature)* &
                                matrix_bar(test_dof, trial_dof)* &
                                component_metric(test_component, trial_component, &
                                quadrature)* &
                                trial_trace(trial_component, quadrature, trial_dof)
                            trial_trace_bar(trial_component, quadrature, trial_dof) = &
                                trial_trace_bar( &
                                trial_component, quadrature, trial_dof) + &
                                physical_weights(quadrature)* &
                                matrix_bar(test_dof, trial_dof)* &
                                test_trace(test_component, quadrature, test_dof)* &
                                component_metric(test_component, trial_component, &
                                quadrature)
                            component_metric_bar(test_component, trial_component, &
                                quadrature) = &
                                component_metric_bar(test_component, trial_component, &
                                quadrature) + physical_weights(quadrature)* &
                                matrix_bar(test_dof, trial_dof)* &
                                test_trace(test_component, quadrature, test_dof)* &
                                trial_trace(trial_component, quadrature, trial_dof)
                            total_weight_bar(quadrature) = &
                                total_weight_bar(quadrature) + &
                                matrix_bar(test_dof, trial_dof)* &
                                test_trace(test_component, quadrature, test_dof)* &
                                component_metric(test_component, trial_component, &
                                quadrature)* &
                                trial_trace(trial_component, quadrature, trial_dof)
                        end do
                    end do
                end do
            end do
        end do
        reference_weights_bar = total_weight_bar*surface_jacobian
        surface_jacobian_bar = total_weight_bar*reference_weights
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_geometry_mortar_component_coupling_vjp

    subroutine validate_component_inputs( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            component_metric, matrix, physical_weights, status)
        real(dp), intent(in) :: test_trace(:, :, :), trial_trace(:, :, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        real(dp), intent(in) :: component_metric(:, :, :)
        real(dp), intent(in) :: matrix(:, :), physical_weights(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_count, trial_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "component mortar received incompatible arrays")
        if (size(test_trace, 1) < 1 .or. size(test_trace, 2) < 1 .or. &
            size(test_trace, 3) < 1) return
        if (size(trial_trace, 1) /= size(test_trace, 1)) return
        quadrature_count = size(test_trace, 2)
        test_count = size(test_trace, 3)
        trial_count = size(trial_trace, 3)
        if (size(trial_trace, 2) /= quadrature_count) return
        if (size(reference_weights) /= quadrature_count) return
        if (size(surface_jacobian) /= quadrature_count) return
        if (size(component_metric, 1) /= size(test_trace, 1) .or. &
            size(component_metric, 2) /= size(trial_trace, 1) .or. &
            size(component_metric, 3) /= quadrature_count) return
        if (size(matrix, 1) /= test_count .or. size(matrix, 2) /= trial_count) return
        if (size(physical_weights) /= quadrature_count) return
        if (any(.not. ieee_is_finite(test_trace)) .or. &
            any(.not. ieee_is_finite(trial_trace)) .or. &
            any(.not. ieee_is_finite(reference_weights)) .or. &
            any(.not. ieee_is_finite(surface_jacobian)) .or. &
            any(.not. ieee_is_finite(component_metric)) .or. &
            any(.not. ieee_is_finite(matrix)) .or. &
            any(.not. ieee_is_finite(physical_weights))) return
        if (any(reference_weights <= 0.0_dp)) return
        if (any(surface_jacobian <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_component_inputs

end module fortfem_geometry_mortar_component_coupling
