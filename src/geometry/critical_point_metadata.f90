module fortfem_critical_point_metadata
    !! Neutral metadata for caller-supplied two-dimensional critical points.
    !!
    !! The candidates are supplied by the caller; this module does not locate
    !! nulls or assign an equilibrium interpretation.  It classifies the
    !! supplied scalar-field gradient/Hessian, records an event margin, and
    !! carries an optional limiter/separatrix label.  Classification labels
    !! are fixed discrete metadata for the JVP/VJP actions.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
                          FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    integer, parameter, public :: CRITICAL_POINT_REGULAR = 0
    integer, parameter, public :: CRITICAL_POINT_O_POINT = 1
    integer, parameter, public :: CRITICAL_POINT_X_POINT = 2
    integer, parameter, public :: CRITICAL_POINT_DEGENERATE = 3
    integer, parameter, public :: CRITICAL_CONTOUR_NONE = 0
    integer, parameter, public :: CRITICAL_CONTOUR_LIMITER = 1
    integer, parameter, public :: CRITICAL_CONTOUR_SEPARATRIX = 2

    type, public :: critical_point_metadata_t
        integer :: candidate_count = 0
        real(dp), allocatable :: points(:, :)
        real(dp), allocatable :: scalar_values(:)
        real(dp), allocatable :: gradient_norm(:)
        real(dp), allocatable :: hessian_determinant(:)
        real(dp), allocatable :: event_margin(:)
        integer, allocatable :: classification(:)
        integer, allocatable :: contour_kind(:)
    end type critical_point_metadata_t

    public :: evaluate_critical_point_metadata
    public :: evaluate_critical_point_metadata_jvp
    public :: evaluate_critical_point_metadata_vjp
    public :: validate_critical_point_metadata

contains

    subroutine evaluate_critical_point_metadata( &
        points, scalar_values, gradients, hessians, gradient_tolerance, &
        hessian_tolerance, contour_kind, metadata, status)
        !! Classify caller-owned scalar-field critical-point candidates.
        real(dp), intent(in) :: points(:, :), scalar_values(:), gradients(:, :)
        real(dp), intent(in) :: hessians(:, :, :)
        real(dp), intent(in) :: gradient_tolerance, hessian_tolerance
        integer, intent(in) :: contour_kind(:)
        type(critical_point_metadata_t), intent(inout) :: metadata
        type(fortsparse_status_t), intent(out) :: status

        integer :: candidate
        call clear_metadata(metadata)
        call validate_raw_inputs( &
            points, scalar_values, gradients, hessians, gradient_tolerance, &
            hessian_tolerance, contour_kind, status)
        if (status%code /= FORTSPARSE_OK) return

        metadata%candidate_count = size(points, 2)
        allocate (metadata%points, source=points)
        allocate (metadata%scalar_values, source=scalar_values)
        allocate (metadata%gradient_norm(metadata%candidate_count))
        allocate (metadata%hessian_determinant(metadata%candidate_count))
        allocate (metadata%event_margin(metadata%candidate_count))
        allocate (metadata%classification(metadata%candidate_count))
        allocate (metadata%contour_kind, source=contour_kind)
        do candidate = 1, metadata%candidate_count
            call classify_candidate( &
                gradients(:, candidate), hessians(:, :, candidate), &
                gradient_tolerance, hessian_tolerance, &
                metadata%gradient_norm(candidate), &
                metadata%hessian_determinant(candidate), &
                metadata%event_margin(candidate), metadata%classification(candidate))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_critical_point_metadata

    subroutine evaluate_critical_point_metadata_jvp( &
        points, scalar_values, gradients, hessians, gradient_tolerance, &
        hessian_tolerance, contour_kind, gradients_dot, hessians_dot, &
        event_margin_dot, status)
        !! Apply the fixed-classification tangent of the criticality margins.
        real(dp), intent(in) :: points(:, :), scalar_values(:), gradients(:, :)
        real(dp), intent(in) :: hessians(:, :, :)
        real(dp), intent(in) :: gradient_tolerance, hessian_tolerance
        integer, intent(in) :: contour_kind(:)
        real(dp), intent(in) :: gradients_dot(:, :), hessians_dot(:, :, :)
        real(dp), intent(out) :: event_margin_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        type(critical_point_metadata_t) :: metadata
        integer :: candidate
        real(dp) :: gradient_norm_dot, determinant_dot
        event_margin_dot = 0.0_dp
        call evaluate_critical_point_metadata( &
            points, scalar_values, gradients, hessians, gradient_tolerance, &
            hessian_tolerance, contour_kind, metadata, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(gradients_dot, 1) /= 2 .or. size(gradients_dot, 2) /= &
            metadata%candidate_count .or. size(hessians_dot, 1) /= 2 .or. &
            size(hessians_dot, 2) /= 2 .or. size(hessians_dot, 3) /= &
            metadata%candidate_count .or. size(event_margin_dot) /= &
            metadata%candidate_count .or. .not. all_finite(gradients_dot) .or. &
            .not. all_finite(hessians_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "critical-point metadata JVP has incompatible increments")
            return
        end if
        do candidate = 1, metadata%candidate_count
            select case (metadata%classification(candidate))
            case (CRITICAL_POINT_REGULAR)
                if (metadata%gradient_norm(candidate) > 0.0_dp) then
                    gradient_norm_dot = dot_product(gradients(:, candidate), &
                          gradients_dot(:, candidate))/metadata%gradient_norm(candidate)
                else
                    gradient_norm_dot = 0.0_dp
                end if
                event_margin_dot(candidate) = gradient_norm_dot
            case default
           determinant_dot = hessians(2, 2, candidate)*hessians_dot(1, 1, candidate) + &
                             hessians(1, 1, candidate)*hessians_dot(2, 2, candidate) - &
                             hessians(2, 1, candidate)*hessians_dot(1, 2, candidate) - &
                                 hessians(1, 2, candidate)*hessians_dot(2, 1, candidate)
                if (metadata%hessian_determinant(candidate) < 0.0_dp) then
                    event_margin_dot(candidate) = -determinant_dot
                else
                    event_margin_dot(candidate) = determinant_dot
                end if
            end select
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_critical_point_metadata_jvp

    subroutine evaluate_critical_point_metadata_vjp( &
        points, scalar_values, gradients, hessians, gradient_tolerance, &
        hessian_tolerance, contour_kind, event_margin_bar, points_bar, &
        scalar_values_bar, gradients_bar, hessians_bar, status)
        !! Apply the real transpose of the margin action.
        real(dp), intent(in) :: points(:, :), scalar_values(:), gradients(:, :)
        real(dp), intent(in) :: hessians(:, :, :)
        real(dp), intent(in) :: gradient_tolerance, hessian_tolerance
        integer, intent(in) :: contour_kind(:)
        real(dp), intent(in) :: event_margin_bar(:)
    real(dp), intent(out) :: points_bar(:, :), scalar_values_bar(:), gradients_bar(:, :)
        real(dp), intent(out) :: hessians_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        type(critical_point_metadata_t) :: metadata
        integer :: candidate
        real(dp) :: scale
        points_bar = 0.0_dp
        scalar_values_bar = 0.0_dp
        gradients_bar = 0.0_dp
        hessians_bar = 0.0_dp
        call evaluate_critical_point_metadata( &
            points, scalar_values, gradients, hessians, gradient_tolerance, &
            hessian_tolerance, contour_kind, metadata, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(event_margin_bar) /= metadata%candidate_count .or. &
            any(shape(points_bar) /= shape(points)) .or. size(scalar_values_bar) /= &
           size(scalar_values) .or. any(shape(gradients_bar) /= shape(gradients)) .or. &
            any(shape(hessians_bar) /= shape(hessians)) .or. &
            .not. all_finite(event_margin_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "critical-point metadata VJP has incompatible cotangents")
            return
        end if
        do candidate = 1, metadata%candidate_count
            scale = event_margin_bar(candidate)
            select case (metadata%classification(candidate))
            case (CRITICAL_POINT_REGULAR)
                if (metadata%gradient_norm(candidate) > 0.0_dp) then
                    gradients_bar(:, candidate) = scale*gradients(:, candidate)/ &
                                                  metadata%gradient_norm(candidate)
                end if
            case default
                if (metadata%hessian_determinant(candidate) < 0.0_dp) scale = -scale
                hessians_bar(1, 1, candidate) = scale*hessians(2, 2, candidate)
                hessians_bar(2, 2, candidate) = scale*hessians(1, 1, candidate)
                hessians_bar(1, 2, candidate) = -scale*hessians(2, 1, candidate)
                hessians_bar(2, 1, candidate) = -scale*hessians(1, 2, candidate)
            end select
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_critical_point_metadata_vjp

    subroutine validate_critical_point_metadata(metadata, status)
        !! Validate a previously evaluated descriptor.
        type(critical_point_metadata_t), intent(in) :: metadata
        type(fortsparse_status_t), intent(out) :: status
        if (metadata%candidate_count < 1 .or. .not. allocated(metadata%points) .or. &
            .not. allocated(metadata%scalar_values) .or. &
            .not. allocated(metadata%gradient_norm) .or. &
            .not. allocated(metadata%hessian_determinant) .or. &
            .not. allocated(metadata%event_margin) .or. &
            .not. allocated(metadata%classification) .or. &
            .not. allocated(metadata%contour_kind)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "critical-point metadata is uninitialized")
            return
        end if
        if (size(metadata%points, 1) /= 2 .or. size(metadata%points, 2) /= &
            metadata%candidate_count .or. size(metadata%scalar_values) /= &
            metadata%candidate_count .or. size(metadata%gradient_norm) /= &
            metadata%candidate_count .or. size(metadata%hessian_determinant) /= &
            metadata%candidate_count .or. size(metadata%event_margin) /= &
            metadata%candidate_count .or. size(metadata%classification) /= &
            metadata%candidate_count .or. size(metadata%contour_kind) /= &
            metadata%candidate_count .or. any(metadata%classification < &
                            CRITICAL_POINT_REGULAR) .or. any(metadata%classification > &
                           CRITICAL_POINT_DEGENERATE) .or. any(metadata%contour_kind < &
                               CRITICAL_CONTOUR_NONE) .or. any(metadata%contour_kind > &
                                                      CRITICAL_CONTOUR_SEPARATRIX)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "critical-point metadata has inconsistent dimensions or labels")
            return
        end if
        if (.not. all_finite(metadata%points) .or. &
            .not. all_finite(metadata%scalar_values) .or. &
            .not. all_finite(metadata%gradient_norm) .or. &
            .not. all_finite(metadata%hessian_determinant) .or. &
            .not. all_finite(metadata%event_margin)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "critical-point metadata contains non-finite values")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_critical_point_metadata

    subroutine classify_candidate(gradient, hessian, gradient_tolerance, &
                  hessian_tolerance, gradient_norm, hessian_determinant, event_margin, &
                                  classification)
        real(dp), intent(in) :: gradient(2), hessian(2, 2)
        real(dp), intent(in) :: gradient_tolerance, hessian_tolerance
        real(dp), intent(out) :: gradient_norm, hessian_determinant, event_margin
        integer, intent(out) :: classification
        gradient_norm = sqrt(dot_product(gradient, gradient))
        hessian_determinant = hessian(1, 1)*hessian(2, 2) - &
                              hessian(1, 2)*hessian(2, 1)
        if (gradient_norm > gradient_tolerance) then
            classification = CRITICAL_POINT_REGULAR
            event_margin = gradient_norm - gradient_tolerance
        else if (hessian_determinant > hessian_tolerance) then
            classification = CRITICAL_POINT_O_POINT
            event_margin = hessian_determinant - hessian_tolerance
        else if (hessian_determinant < -hessian_tolerance) then
            classification = CRITICAL_POINT_X_POINT
            event_margin = -hessian_determinant - hessian_tolerance
        else
            classification = CRITICAL_POINT_DEGENERATE
            event_margin = abs(hessian_determinant) - hessian_tolerance
        end if
    end subroutine classify_candidate

    subroutine validate_raw_inputs(points, scalar_values, gradients, hessians, &
                            gradient_tolerance, hessian_tolerance, contour_kind, status)
        real(dp), intent(in) :: points(:, :), scalar_values(:), gradients(:, :)
        real(dp), intent(in) :: hessians(:, :, :)
        real(dp), intent(in) :: gradient_tolerance, hessian_tolerance
        integer, intent(in) :: contour_kind(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: candidate_count
        candidate_count = size(points, 2)
        if (size(points, 1) /= 2 .or. candidate_count < 1 .or. &
            size(scalar_values) /= candidate_count .or. size(gradients, 1) /= 2 .or. &
            size(gradients, 2) /= candidate_count .or. size(hessians, 1) /= 2 .or. &
            size(hessians, 2) /= 2 .or. size(hessians, 3) /= candidate_count .or. &
            size(contour_kind) /= candidate_count .or. &
        .not. ieee_is_finite(gradient_tolerance) .or. gradient_tolerance < 0.0_dp .or. &
          .not. ieee_is_finite(hessian_tolerance) .or. hessian_tolerance < 0.0_dp .or. &
            any(contour_kind < CRITICAL_CONTOUR_NONE) .or. &
            any(contour_kind > CRITICAL_CONTOUR_SEPARATRIX)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "critical-point candidates have incompatible metadata")
            return
        end if
        if (.not. all_finite(points) .or. .not. all_finite(scalar_values) .or. &
            .not. all_finite(gradients) .or. .not. all_finite(hessians)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "critical-point candidates contain non-finite values")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_raw_inputs

    subroutine clear_metadata(metadata)
        type(critical_point_metadata_t), intent(inout) :: metadata
        metadata%candidate_count = 0
        if (allocated(metadata%points)) deallocate (metadata%points)
        if (allocated(metadata%scalar_values)) deallocate (metadata%scalar_values)
        if (allocated(metadata%gradient_norm)) deallocate (metadata%gradient_norm)
  if (allocated(metadata%hessian_determinant)) deallocate (metadata%hessian_determinant)
        if (allocated(metadata%event_margin)) deallocate (metadata%event_margin)
        if (allocated(metadata%classification)) deallocate (metadata%classification)
        if (allocated(metadata%contour_kind)) deallocate (metadata%contour_kind)
    end subroutine clear_metadata

    logical function all_finite(values) result(valid)
        real(dp), intent(in) :: values(..)
        select rank (values)
        rank (1)
            valid = all(ieee_is_finite(values))
        rank (2)
            valid = all(ieee_is_finite(values))
        rank (3)
            valid = all(ieee_is_finite(values))
        rank default
            valid = .false.
        end select
    end function all_finite

end module fortfem_critical_point_metadata
