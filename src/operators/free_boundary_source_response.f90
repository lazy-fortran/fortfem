module fortfem_free_boundary_source_response
    !! Neutral source-basis to free-boundary trace composition.
    !!
    !! A caller-owned source basis is evaluated on a sampled boundary and fed
    !! directly into the weighted free-boundary port residual:
    !!
    !!   t(q,c) = sum_s B(q,c,s) a(s),
    !!   r(q,c) = w(q) (t(q,c) - g(q,c) - k(q,c)).
    !!
    !! This is the small composition shared by source-to-trace maps,
    !! NESTOR/BIEST-like modal adapters, virtual casing, FEM--BEM/DtN, and
    !! wall response clients.  Geometry, kernels, units, and physical
    !! conventions remain caller-owned.  The fixed-topology real JVP/VJP
    !! actions include basis, source, target, weight, and optional sheet data.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_free_boundary_source_response
    public :: assemble_free_boundary_source_response_jvp
    public :: assemble_free_boundary_source_response_vjp

contains

    subroutine assemble_free_boundary_source_response( &
            source_basis, source_coefficients, external_target, weights, &
            physical_trace, residual, status, sheet_current)
        real(dp), intent(in) :: source_basis(:, :, :), source_coefficients(:)
        real(dp), intent(in) :: external_target(:, :), weights(:)
        real(dp), allocatable, intent(out) :: physical_trace(:, :), residual(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: sheet_current(:, :)

        integer :: q, component, source

        call allocate_empty(physical_trace, residual)
        call validate_base(source_basis, source_coefficients, external_target, weights, &
            sheet_current, status)
        if (status%code /= FORTSPARSE_OK) return
        deallocate(physical_trace, residual)
        allocate(physical_trace(size(external_target, 1), size(external_target, 2)), &
            residual(size(external_target, 1), size(external_target, 2)))
        physical_trace = 0.0_dp
        do q = 1, size(external_target, 1)
            do component = 1, size(external_target, 2)
                do source = 1, size(source_coefficients)
                    physical_trace(q, component) = physical_trace(q, component) + &
                        source_basis(q, component, source)*source_coefficients(source)
                end do
                residual(q, component) = weights(q)*(physical_trace(q, component) - &
                    external_target(q, component) - optional_value(sheet_current, q, component))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_free_boundary_source_response

    subroutine assemble_free_boundary_source_response_jvp( &
            source_basis, source_coefficients, external_target, weights, &
            source_basis_dot, source_coefficients_dot, external_target_dot, weights_dot, &
            physical_trace_dot, residual_dot, status, sheet_current, sheet_current_dot)
        real(dp), intent(in) :: source_basis(:, :, :), source_coefficients(:)
        real(dp), intent(in) :: external_target(:, :), weights(:)
        real(dp), intent(in) :: source_basis_dot(:, :, :), source_coefficients_dot(:)
        real(dp), intent(in) :: external_target_dot(:, :), weights_dot(:)
        real(dp), allocatable, intent(out) :: physical_trace_dot(:, :), residual_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: sheet_current(:, :), sheet_current_dot(:, :)

        integer :: q, component, source
        real(dp) :: raw, raw_dot

        call allocate_empty(physical_trace_dot, residual_dot)
        call validate_base(source_basis, source_coefficients, external_target, weights, &
            sheet_current, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_direction( &
            source_basis, source_coefficients, external_target, weights, source_basis_dot, &
            source_coefficients_dot, external_target_dot, weights_dot, sheet_current, &
            sheet_current_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "free-boundary source response JVP has incompatible increments")
            return
        end if
        deallocate(physical_trace_dot, residual_dot)
        allocate(physical_trace_dot(size(external_target, 1), size(external_target, 2)), &
            residual_dot(size(external_target, 1), size(external_target, 2)))
        physical_trace_dot = 0.0_dp
        residual_dot = 0.0_dp
        do q = 1, size(external_target, 1)
            do component = 1, size(external_target, 2)
                do source = 1, size(source_coefficients)
                    physical_trace_dot(q, component) = physical_trace_dot(q, component) + &
                        source_basis_dot(q, component, source)*source_coefficients(source) + &
                        source_basis(q, component, source)*source_coefficients_dot(source)
                end do
                raw = -external_target(q, component) - &
                    optional_value(sheet_current, q, component)
                raw_dot = -external_target_dot(q, component) - &
                    optional_value(sheet_current_dot, q, component)
                do source = 1, size(source_coefficients)
                    raw = raw + source_basis(q, component, source)*source_coefficients(source)
                end do
                raw_dot = raw_dot + physical_trace_dot(q, component)
                residual_dot(q, component) = weights_dot(q)*raw + weights(q)*raw_dot
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_free_boundary_source_response_jvp

    subroutine assemble_free_boundary_source_response_vjp( &
            source_basis, source_coefficients, external_target, weights, &
            physical_trace_bar, residual_bar, source_basis_bar, source_coefficients_bar, &
            external_target_bar, weights_bar, status, sheet_current, sheet_current_bar)
        real(dp), intent(in) :: source_basis(:, :, :), source_coefficients(:)
        real(dp), intent(in) :: external_target(:, :), weights(:)
        real(dp), intent(in) :: physical_trace_bar(:, :), residual_bar(:, :)
        real(dp), allocatable, intent(out) :: source_basis_bar(:, :, :)
        real(dp), allocatable, intent(out) :: source_coefficients_bar(:)
        real(dp), allocatable, intent(out) :: external_target_bar(:, :), weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: sheet_current(:, :)
        real(dp), allocatable, intent(out), optional :: sheet_current_bar(:, :)

        integer :: q, component, source
        real(dp) :: trace_bar_total, raw, raw_bar

        call validate_base(source_basis, source_coefficients, external_target, weights, &
            sheet_current, status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(physical_trace_bar) /= shape(external_target)) .or. &
            any(shape(residual_bar) /= shape(external_target)) .or. &
            .not. all(ieee_is_finite(physical_trace_bar)) .or. &
            .not. all(ieee_is_finite(residual_bar)) .or. &
            present(sheet_current) .neqv. present(sheet_current_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "free-boundary source response VJP has incompatible cotangents")
            return
        end if
        if (allocated(source_basis_bar)) deallocate(source_basis_bar)
        if (allocated(source_coefficients_bar)) deallocate(source_coefficients_bar)
        if (allocated(external_target_bar)) deallocate(external_target_bar)
        if (allocated(weights_bar)) deallocate(weights_bar)
        allocate(source_basis_bar(size(source_basis, 1), size(source_basis, 2), &
            size(source_basis, 3)), source_coefficients_bar(size(source_coefficients)), &
            external_target_bar(size(external_target, 1), size(external_target, 2)), &
            weights_bar(size(weights)))
        source_basis_bar = 0.0_dp
        source_coefficients_bar = 0.0_dp
        external_target_bar = 0.0_dp
        weights_bar = 0.0_dp
        if (present(sheet_current_bar)) then
            if (allocated(sheet_current_bar)) deallocate(sheet_current_bar)
            allocate(sheet_current_bar(size(external_target, 1), size(external_target, 2)))
            sheet_current_bar = 0.0_dp
        end if
        do q = 1, size(external_target, 1)
            do component = 1, size(external_target, 2)
                raw = -external_target(q, component) - &
                    optional_value(sheet_current, q, component)
                do source = 1, size(source_coefficients)
                    raw = raw + source_basis(q, component, source)*source_coefficients(source)
                end do
                raw_bar = weights(q)*residual_bar(q, component)
                trace_bar_total = physical_trace_bar(q, component) + raw_bar
                weights_bar(q) = weights_bar(q) + residual_bar(q, component)*raw
                external_target_bar(q, component) = -raw_bar
                if (present(sheet_current_bar)) sheet_current_bar(q, component) = -raw_bar
                do source = 1, size(source_coefficients)
                    source_basis_bar(q, component, source) = &
                        trace_bar_total*source_coefficients(source)
                    source_coefficients_bar(source) = source_coefficients_bar(source) + &
                        source_basis(q, component, source)*trace_bar_total
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_free_boundary_source_response_vjp

    subroutine validate_base(source_basis, source_coefficients, external_target, weights, &
            sheet_current, status)
        real(dp), intent(in) :: source_basis(:, :, :), source_coefficients(:)
        real(dp), intent(in) :: external_target(:, :), weights(:)
        real(dp), intent(in), optional :: sheet_current(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "free-boundary source response received incompatible arrays")
        if (size(source_basis, 1) < 1 .or. size(source_basis, 2) < 1 .or. &
            size(source_basis, 3) < 1 .or. size(source_coefficients) /= size(source_basis, 3) .or. &
            any(shape(external_target) /= shape(source_basis(:, :, 1))) .or. &
            size(weights) /= size(external_target, 1)) return
        if (.not. all(ieee_is_finite(source_basis)) .or. &
            .not. all(ieee_is_finite(source_coefficients)) .or. &
            .not. all(ieee_is_finite(external_target)) .or. &
            .not. all(ieee_is_finite(weights)) .or. any(weights <= 0.0_dp)) return
        if (present(sheet_current)) then
            if (any(shape(sheet_current) /= shape(external_target)) .or. &
                .not. all(ieee_is_finite(sheet_current))) return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_base

    logical function valid_direction( &
            source_basis, source_coefficients, external_target, weights, source_basis_dot, &
            source_coefficients_dot, external_target_dot, weights_dot, sheet_current, &
            sheet_current_dot) result(valid)
        real(dp), intent(in) :: source_basis(:, :, :), source_coefficients(:)
        real(dp), intent(in) :: external_target(:, :), weights(:)
        real(dp), intent(in) :: source_basis_dot(:, :, :), source_coefficients_dot(:)
        real(dp), intent(in) :: external_target_dot(:, :), weights_dot(:)
        real(dp), intent(in), optional :: sheet_current(:, :), sheet_current_dot(:, :)

        valid = .not. any(shape(source_basis_dot) /= shape(source_basis)) .and. &
            size(source_coefficients_dot) == size(source_coefficients) .and. &
            all(shape(external_target_dot) == shape(external_target)) .and. &
            size(weights_dot) == size(weights) .and. &
            finite(source_basis_dot) .and. finite(source_coefficients_dot) .and. &
            finite(external_target_dot) .and. finite(weights_dot) .and. &
            present(sheet_current) .eqv. present(sheet_current_dot)
        if (valid .and. present(sheet_current_dot)) valid = &
            all(shape(sheet_current_dot) == shape(sheet_current)) .and. finite(sheet_current_dot)
    end function valid_direction

    real(dp) function optional_value(values, q, component) result(value)
        real(dp), intent(in), optional :: values(:, :)
        integer, intent(in) :: q, component

        value = 0.0_dp
        if (present(values)) value = values(q, component)
    end function optional_value

    subroutine allocate_empty(physical_trace, residual)
        real(dp), allocatable, intent(out) :: physical_trace(:, :), residual(:, :)

        allocate(physical_trace(0, 0), residual(0, 0))
    end subroutine allocate_empty

    pure logical function finite(values) result(valid)
        real(dp), intent(in) :: values(..)

        valid = .true.
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
    end function finite

end module fortfem_free_boundary_source_response
