module fortfem_geometry_mortar_trace_coupling
    !! Geometry-aware cross-mass for independently discretized traces.
    !!
    !! A caller supplies trace values on a fixed reference quadrature and the
    !! positive surface metric at each row.  The product
    !!
    !!     physical_weight(q) = reference_weight(q)*surface_jacobian(q)
    !!
    !! is the only geometry operation owned here.  This keeps NURBS, Fourier,
    !! cut-surface, and panel geometry providers interchangeable while giving
    !! all mortar clients one physical-measure and derivative contract.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_geometry_mortar_trace_coupling
    public :: assemble_geometry_mortar_trace_coupling_jvp
    public :: assemble_geometry_mortar_trace_coupling_vjp

contains

    subroutine assemble_geometry_mortar_trace_coupling( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            matrix, physical_weights, status)
        !! Assemble a physical-measure test/trial trace cross-mass.
        real(dp), intent(in) :: test_trace(:, :), trial_trace(:, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        real(dp), intent(out) :: matrix(:, :), physical_weights(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof

        matrix = 0.0_dp
        physical_weights = 0.0_dp
        call validate_geometry_mortar_inputs( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            matrix, physical_weights, status)
        if (status%code /= FORTSPARSE_OK) return

        quadrature_count = size(test_trace, 1)
        test_count = size(test_trace, 2)
        trial_count = size(trial_trace, 2)
        physical_weights = reference_weights*surface_jacobian
        do quadrature = 1, quadrature_count
            do test_dof = 1, test_count
                do trial_dof = 1, trial_count
                    matrix(test_dof, trial_dof) = &
                        matrix(test_dof, trial_dof) + &
                        physical_weights(quadrature)* &
                        test_trace(quadrature, test_dof)* &
                        trial_trace(quadrature, trial_dof)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_geometry_mortar_trace_coupling

    subroutine assemble_geometry_mortar_trace_coupling_jvp( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            test_trace_dot, trial_trace_dot, reference_weights_dot, &
            surface_jacobian_dot, matrix_dot, physical_weights_dot, status)
        !! Apply the fixed-topology product-rule JVP.
        real(dp), intent(in) :: test_trace(:, :), trial_trace(:, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        real(dp), intent(in) :: test_trace_dot(:, :), trial_trace_dot(:, :)
        real(dp), intent(in) :: reference_weights_dot(:)
        real(dp), intent(in) :: surface_jacobian_dot(:)
        real(dp), intent(out) :: matrix_dot(:, :), physical_weights_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: physical_weights(:)
        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof

        matrix_dot = 0.0_dp
        physical_weights_dot = 0.0_dp
        call validate_geometry_mortar_inputs( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            matrix_dot, physical_weights_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_trace, 1)
        test_count = size(test_trace, 2)
        trial_count = size(trial_trace, 2)
        if (size(test_trace_dot, 1) /= quadrature_count .or. &
            size(test_trace_dot, 2) /= test_count .or. &
            size(trial_trace_dot, 1) /= quadrature_count .or. &
            size(trial_trace_dot, 2) /= trial_count .or. &
            size(reference_weights_dot) /= quadrature_count .or. &
            size(surface_jacobian_dot) /= quadrature_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "geometry mortar JVP has incompatible increment shapes")
            return
        end if
        if (any(.not. ieee_is_finite(test_trace_dot)) .or. &
            any(.not. ieee_is_finite(trial_trace_dot)) .or. &
            any(.not. ieee_is_finite(reference_weights_dot)) .or. &
            any(.not. ieee_is_finite(surface_jacobian_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "geometry mortar JVP received non-finite increments")
            return
        end if

        allocate(physical_weights(quadrature_count))
        physical_weights = reference_weights*surface_jacobian
        physical_weights_dot = reference_weights_dot*surface_jacobian + &
            reference_weights*surface_jacobian_dot
        do quadrature = 1, quadrature_count
            do test_dof = 1, test_count
                do trial_dof = 1, trial_count
                    matrix_dot(test_dof, trial_dof) = &
                        matrix_dot(test_dof, trial_dof) + &
                        physical_weights_dot(quadrature)* &
                        test_trace(quadrature, test_dof)* &
                        trial_trace(quadrature, trial_dof) + &
                        physical_weights(quadrature)* &
                        test_trace_dot(quadrature, test_dof)* &
                        trial_trace(quadrature, trial_dof) + &
                        physical_weights(quadrature)* &
                        test_trace(quadrature, test_dof)* &
                        trial_trace_dot(quadrature, trial_dof)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_geometry_mortar_trace_coupling_jvp

    subroutine assemble_geometry_mortar_trace_coupling_vjp( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            matrix_bar, physical_weights_bar, test_trace_bar, trial_trace_bar, &
            reference_weights_bar, surface_jacobian_bar, status)
        !! Apply the real reverse product of the physical cross-mass.
        real(dp), intent(in) :: test_trace(:, :), trial_trace(:, :)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        real(dp), intent(in) :: matrix_bar(:, :), physical_weights_bar(:)
        real(dp), intent(out) :: test_trace_bar(:, :), trial_trace_bar(:, :)
        real(dp), intent(out) :: reference_weights_bar(:)
        real(dp), intent(out) :: surface_jacobian_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: physical_weights(:), physical_weights_total_bar(:)
        integer :: quadrature_count, test_count, trial_count
        integer :: quadrature, test_dof, trial_dof

        test_trace_bar = 0.0_dp
        trial_trace_bar = 0.0_dp
        reference_weights_bar = 0.0_dp
        surface_jacobian_bar = 0.0_dp
        call validate_geometry_mortar_inputs( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            matrix_bar, physical_weights_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_trace, 1)
        test_count = size(test_trace, 2)
        trial_count = size(trial_trace, 2)
        if (size(test_trace_bar, 1) /= quadrature_count .or. &
            size(test_trace_bar, 2) /= test_count .or. &
            size(trial_trace_bar, 1) /= quadrature_count .or. &
            size(trial_trace_bar, 2) /= trial_count .or. &
            size(reference_weights_bar) /= quadrature_count .or. &
            size(surface_jacobian_bar) /= quadrature_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "geometry mortar VJP has incompatible cotangent shapes")
            return
        end if

        allocate(physical_weights(quadrature_count))
        allocate(physical_weights_total_bar(quadrature_count))
        physical_weights = reference_weights*surface_jacobian
        physical_weights_total_bar = physical_weights_bar
        do quadrature = 1, quadrature_count
            do test_dof = 1, test_count
                do trial_dof = 1, trial_count
                    test_trace_bar(quadrature, test_dof) = &
                        test_trace_bar(quadrature, test_dof) + &
                        physical_weights(quadrature)* &
                        matrix_bar(test_dof, trial_dof)* &
                        trial_trace(quadrature, trial_dof)
                    trial_trace_bar(quadrature, trial_dof) = &
                        trial_trace_bar(quadrature, trial_dof) + &
                        physical_weights(quadrature)* &
                        matrix_bar(test_dof, trial_dof)* &
                        test_trace(quadrature, test_dof)
                    physical_weights_total_bar(quadrature) = &
                        physical_weights_total_bar(quadrature) + &
                        matrix_bar(test_dof, trial_dof)* &
                        test_trace(quadrature, test_dof)* &
                        trial_trace(quadrature, trial_dof)
                end do
            end do
        end do
        reference_weights_bar = physical_weights_total_bar*surface_jacobian
        surface_jacobian_bar = physical_weights_total_bar*reference_weights
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_geometry_mortar_trace_coupling_vjp

    subroutine validate_geometry_mortar_inputs( &
            test_trace, trial_trace, reference_weights, surface_jacobian, &
            matrix, physical_weights, status)
        real(dp), intent(in) :: test_trace(:, :), trial_trace(:,:)
        real(dp), intent(in) :: reference_weights(:), surface_jacobian(:)
        real(dp), intent(in) :: matrix(:, :), physical_weights(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_count, trial_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "geometry mortar received incompatible arrays")
        quadrature_count = size(test_trace, 1)
        test_count = size(test_trace, 2)
        trial_count = size(trial_trace, 2)
        if (quadrature_count < 1) return
        if (test_count < 1) return
        if (trial_count < 1) return
        if (size(trial_trace, 1) /= quadrature_count) return
        if (size(reference_weights) /= quadrature_count) return
        if (size(surface_jacobian) /= quadrature_count) return
        if (size(matrix, 1) /= test_count) return
        if (size(matrix, 2) /= trial_count) return
        if (size(physical_weights) /= quadrature_count) return
        if (any(.not. ieee_is_finite(test_trace))) return
        if (any(.not. ieee_is_finite(trial_trace))) return
        if (any(.not. ieee_is_finite(reference_weights))) return
        if (any(.not. ieee_is_finite(surface_jacobian))) return
        if (any(.not. ieee_is_finite(matrix))) return
        if (any(.not. ieee_is_finite(physical_weights))) return
        if (any(reference_weights <= 0.0_dp)) return
        if (any(surface_jacobian <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_geometry_mortar_inputs

end module fortfem_geometry_mortar_trace_coupling
