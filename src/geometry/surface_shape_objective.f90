module fortfem_surface_shape_objective
    !! Weighted fixed-topology mismatch of sampled physical surface coordinates.
    !!
    !! The objective is
    !!
    !!   1/2 sum_q w_q ||x_candidate(q) - x_target(q)||^2.
    !!
    !! Coordinates, point ordering, and quadrature weights are caller-owned.
    !! This is deliberately a geometry-only contract: topology, parameter
    !! maps, profiles, and free-boundary laws remain external to FortFEM.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: evaluate_surface_shape_objective
    public :: evaluate_surface_shape_objective_jvp
    public :: evaluate_surface_shape_objective_vjp

contains

    subroutine evaluate_surface_shape_objective( &
            candidate_coordinates, target_coordinates, weights, objective, &
            status)
        real(dp), intent(in) :: candidate_coordinates(:, :)
        real(dp), intent(in) :: target_coordinates(:, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: objective
        type(fortsparse_status_t), intent(out) :: status
        integer :: component, point
        real(dp) :: difference

        objective = 0.0_dp
        if (.not. validate_surface_samples( &
            candidate_coordinates, target_coordinates, weights, status)) return
        do point = 1, size(weights)
            do component = 1, size(candidate_coordinates, 1)
                difference = candidate_coordinates(component, point) - &
                    target_coordinates(component, point)
                objective = objective + 0.5_dp*weights(point)*difference**2
            end do
        end do
        if (.not. ieee_is_finite(objective)) then
            objective = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape objective is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_surface_shape_objective

    subroutine evaluate_surface_shape_objective_jvp( &
            candidate_coordinates, target_coordinates, weights, &
            candidate_dot, target_dot, weights_dot, objective_dot, status)
        real(dp), intent(in) :: candidate_coordinates(:, :)
        real(dp), intent(in) :: target_coordinates(:, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: candidate_dot(:, :)
        real(dp), intent(in) :: target_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: objective_dot
        type(fortsparse_status_t), intent(out) :: status
        integer :: component, point
        real(dp) :: difference, difference_dot

        objective_dot = 0.0_dp
        if (.not. validate_surface_samples( &
            candidate_coordinates, target_coordinates, weights, status)) return
        if (.not. same_shape(candidate_dot, candidate_coordinates)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape JVP has an incompatible candidate tangent")
            return
        end if
        if (.not. same_shape(target_dot, target_coordinates)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape JVP has an incompatible target tangent")
            return
        end if
        if (size(weights_dot) /= size(weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape JVP has an incompatible weight tangent")
            return
        end if
        if (.not. all(ieee_is_finite(candidate_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape JVP candidate tangent is non-finite")
            return
        end if
        if (.not. all(ieee_is_finite(target_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape JVP target tangent is non-finite")
            return
        end if
        if (.not. all(ieee_is_finite(weights_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape JVP weight tangent is non-finite")
            return
        end if
        do point = 1, size(weights)
            do component = 1, size(candidate_coordinates, 1)
                difference = candidate_coordinates(component, point) - &
                    target_coordinates(component, point)
                difference_dot = candidate_dot(component, point) - &
                    target_dot(component, point)
                objective_dot = objective_dot + weights(point)*difference* &
                    difference_dot + 0.5_dp*weights_dot(point)*difference**2
            end do
        end do
        if (.not. ieee_is_finite(objective_dot)) then
            objective_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape objective JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_surface_shape_objective_jvp

    subroutine evaluate_surface_shape_objective_vjp( &
            candidate_coordinates, target_coordinates, weights, objective_bar, &
            candidate_bar, target_bar, weights_bar, status)
        real(dp), intent(in) :: candidate_coordinates(:, :)
        real(dp), intent(in) :: target_coordinates(:, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: objective_bar
        real(dp), intent(out) :: candidate_bar(:, :)
        real(dp), intent(out) :: target_bar(:, :)
        real(dp), intent(out) :: weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: component, point
        real(dp) :: difference, squared_norm

        candidate_bar = 0.0_dp
        target_bar = 0.0_dp
        weights_bar = 0.0_dp
        if (.not. validate_surface_samples( &
            candidate_coordinates, target_coordinates, weights, status)) return
        if (.not. same_shape(candidate_bar, candidate_coordinates)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape VJP has an incompatible candidate cotangent")
            return
        end if
        if (.not. same_shape(target_bar, target_coordinates)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape VJP has an incompatible target cotangent")
            return
        end if
        if (size(weights_bar) /= size(weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape VJP has an incompatible weight cotangent")
            return
        end if
        if (.not. ieee_is_finite(objective_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape VJP objective cotangent is non-finite")
            return
        end if
        do point = 1, size(weights)
            squared_norm = 0.0_dp
            do component = 1, size(candidate_coordinates, 1)
                difference = candidate_coordinates(component, point) - &
                    target_coordinates(component, point)
                candidate_bar(component, point) = objective_bar*weights(point)* &
                    difference
                target_bar(component, point) = -candidate_bar(component, point)
                squared_norm = squared_norm + difference**2
            end do
            weights_bar(point) = 0.5_dp*objective_bar*squared_norm
        end do
        if (.not. all(ieee_is_finite(candidate_bar))) then
            candidate_bar = 0.0_dp
            target_bar = 0.0_dp
            weights_bar = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape objective candidate VJP is non-finite")
            return
        end if
        if (.not. all(ieee_is_finite(target_bar))) then
            candidate_bar = 0.0_dp
            target_bar = 0.0_dp
            weights_bar = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape objective target VJP is non-finite")
            return
        end if
        if (.not. all(ieee_is_finite(weights_bar))) then
            candidate_bar = 0.0_dp
            target_bar = 0.0_dp
            weights_bar = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape objective weight VJP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_surface_shape_objective_vjp

    logical function validate_surface_samples( &
            candidate_coordinates, target_coordinates, weights, status) result(valid)
        real(dp), intent(in) :: candidate_coordinates(:, :)
        real(dp), intent(in) :: target_coordinates(:, :)
        real(dp), intent(in) :: weights(:)
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        if (size(candidate_coordinates, 1) < 1 .or. &
            size(candidate_coordinates, 2) < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape objective requires non-empty coordinates")
            return
        end if
        if (.not. same_shape(target_coordinates, candidate_coordinates)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape objective requires matching coordinate topology")
            return
        end if
        if (size(weights) /= size(candidate_coordinates, 2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape objective has incompatible quadrature weights")
            return
        end if
        if (.not. all(ieee_is_finite(candidate_coordinates))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape candidate coordinates are non-finite")
            return
        end if
        if (.not. all(ieee_is_finite(target_coordinates))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape target coordinates are non-finite")
            return
        end if
        if (.not. all(ieee_is_finite(weights))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape quadrature weights are non-finite")
            return
        end if
        if (.not. all(weights > 0.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface shape quadrature weights must be positive")
            return
        end if
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_surface_samples

    logical function same_shape(left, right) result(same)
        real(dp), intent(in) :: left(:, :), right(:, :)

        same = .false.
        if (size(left, 1) /= size(right, 1)) return
        if (size(left, 2) /= size(right, 2)) return
        same = .true.
    end function same_shape

end module fortfem_surface_shape_objective
