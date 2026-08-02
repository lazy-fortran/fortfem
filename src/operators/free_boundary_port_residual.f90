module fortfem_free_boundary_port_residual
    !! Neutral fixed/free-boundary external-field port residual.
    !!
    !! A caller evaluates a physical trace on a fixed-topology boundary and
    !! supplies the trace expected from an exterior/vacuum block.  In the
    !! trace coordinates selected by that caller, the weighted residual is
    !!
    !!   r(q,c) = w(q) ( t(q,c) - g(q,c) - k(q,c) ),
    !!
    !! where `k` is an optional supplied sheet-current (or jump) target.  The
    !! module does not construct a field, project a current, or choose an
    !! Ampere, BEM, DtN, PML, or equilibrium convention.  Those choices stay
    !! with the external adapter.  This is therefore usable for both fixed
    !! and free-boundary clients, provided their trace coordinates have the
    !! same work pairing.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_free_boundary_port_residual
    public :: assemble_free_boundary_port_residual_jvp
    public :: assemble_free_boundary_port_residual_vjp

contains

    subroutine assemble_free_boundary_port_residual( &
            physical_trace, external_target, weights, residual, status, &
            sheet_current)
        !! Evaluate a positive-weighted physical trace mismatch.
        real(dp), intent(in) :: physical_trace(:, :), external_target(:, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: residual(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: sheet_current(:, :)

        integer :: quadrature, component
        real(dp) :: raw

        residual = 0.0_dp
        call validate_base(physical_trace, external_target, weights, residual, &
            sheet_current, status)
        if (status%code /= FORTSPARSE_OK) return

        do quadrature = 1, size(weights)
            do component = 1, size(physical_trace, 2)
                raw = physical_trace(quadrature, component) - &
                    external_target(quadrature, component)
                if (present(sheet_current)) raw = raw - &
                    sheet_current(quadrature, component)
                residual(quadrature, component) = weights(quadrature)*raw
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_free_boundary_port_residual

    subroutine assemble_free_boundary_port_residual_jvp( &
            physical_trace, external_target, weights, physical_trace_dot, &
            external_target_dot, weights_dot, residual_dot, status, &
            sheet_current, sheet_current_dot)
        !! Apply the full fixed-topology product-rule derivative.
        real(dp), intent(in) :: physical_trace(:, :), external_target(:, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: physical_trace_dot(:, :), external_target_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: residual_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: sheet_current(:, :), sheet_current_dot(:, :)

        integer :: quadrature, component
        real(dp) :: raw, raw_dot

        residual_dot = 0.0_dp
        call validate_base(physical_trace, external_target, weights, residual_dot, &
            sheet_current, status)
        if (status%code /= FORTSPARSE_OK) return
        call validate_direction(physical_trace, external_target, weights, &
            physical_trace_dot, external_target_dot, weights_dot, residual_dot, &
            sheet_current, sheet_current_dot, status)
        if (status%code /= FORTSPARSE_OK) return

        do quadrature = 1, size(weights)
            do component = 1, size(physical_trace, 2)
                raw = physical_trace(quadrature, component) - &
                    external_target(quadrature, component)
                raw_dot = physical_trace_dot(quadrature, component) - &
                    external_target_dot(quadrature, component)
                if (present(sheet_current)) then
                    raw = raw - sheet_current(quadrature, component)
                    raw_dot = raw_dot - sheet_current_dot(quadrature, component)
                end if
                residual_dot(quadrature, component) = &
                    weights_dot(quadrature)*raw + weights(quadrature)*raw_dot
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_free_boundary_port_residual_jvp

    subroutine assemble_free_boundary_port_residual_vjp( &
            physical_trace, external_target, weights, residual_bar, &
            physical_trace_bar, external_target_bar, weights_bar, status, &
            sheet_current, sheet_current_bar)
        !! Apply the real transpose of the fixed-topology port residual.
        real(dp), intent(in) :: physical_trace(:, :), external_target(:, :)
        real(dp), intent(in) :: weights(:), residual_bar(:, :)
        real(dp), intent(out) :: physical_trace_bar(:, :), external_target_bar(:, :)
        real(dp), intent(out) :: weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: sheet_current(:, :)
        real(dp), intent(out), optional :: sheet_current_bar(:, :)

        integer :: quadrature, component
        real(dp) :: raw, weighted_bar

        physical_trace_bar = 0.0_dp
        external_target_bar = 0.0_dp
        weights_bar = 0.0_dp
        if (present(sheet_current_bar)) sheet_current_bar = 0.0_dp
        call validate_base(physical_trace, external_target, weights, residual_bar, &
            sheet_current, status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(physical_trace_bar) /= shape(physical_trace)) .or. &
            any(shape(external_target_bar) /= shape(external_target)) .or. &
            size(weights_bar) /= size(weights) .or. &
            present(sheet_current) .neqv. present(sheet_current_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "free-boundary port VJP has incompatible cotangent arrays")
            return
        end if
        if (present(sheet_current_bar)) then
            if (any(shape(sheet_current_bar) /= shape(sheet_current))) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "free-boundary port VJP has incompatible sheet cotangent")
                return
            end if
        end if
        if (.not. all(ieee_is_finite(residual_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "free-boundary port VJP received non-finite cotangents")
            return
        end if

        do quadrature = 1, size(weights)
            do component = 1, size(physical_trace, 2)
                raw = physical_trace(quadrature, component) - &
                    external_target(quadrature, component)
                if (present(sheet_current)) raw = raw - &
                    sheet_current(quadrature, component)
                weighted_bar = weights(quadrature)* &
                    residual_bar(quadrature, component)
                physical_trace_bar(quadrature, component) = weighted_bar
                external_target_bar(quadrature, component) = -weighted_bar
                weights_bar(quadrature) = weights_bar(quadrature) + &
                    residual_bar(quadrature, component)*raw
                if (present(sheet_current_bar)) sheet_current_bar( &
                    quadrature, component) = -weighted_bar
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_free_boundary_port_residual_vjp

    subroutine validate_base(physical_trace, external_target, weights, residual, &
            sheet_current, status)
        real(dp), intent(in) :: physical_trace(:, :), external_target(:, :)
        real(dp), intent(in) :: weights(:), residual(:, :)
        real(dp), intent(in), optional :: sheet_current(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "free-boundary port received incompatible arrays")
        if (size(physical_trace, 1) < 1 .or. size(physical_trace, 2) < 1) return
        if (any(shape(external_target) /= shape(physical_trace)) .or. &
            any(shape(residual) /= shape(physical_trace)) .or. &
            size(weights) /= size(physical_trace, 1)) return
        if (present(sheet_current)) then
            if (any(shape(sheet_current) /= shape(physical_trace))) return
            if (.not. all(ieee_is_finite(sheet_current))) return
        end if
        if (any(weights <= 0.0_dp) .or. &
            .not. all(ieee_is_finite(weights)) .or. &
            .not. all(ieee_is_finite(physical_trace)) .or. &
            .not. all(ieee_is_finite(external_target))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_base

    subroutine validate_direction(physical_trace, external_target, weights, &
            physical_trace_dot, external_target_dot, weights_dot, residual_dot, &
            sheet_current, sheet_current_dot, status)
        real(dp), intent(in) :: physical_trace(:, :), external_target(:, :)
        real(dp), intent(in) :: weights(:), physical_trace_dot(:, :)
        real(dp), intent(in) :: external_target_dot(:, :), weights_dot(:)
        real(dp), intent(in) :: residual_dot(:, :)
        real(dp), intent(in), optional :: sheet_current(:, :), sheet_current_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "free-boundary port JVP received incompatible increments")
        if (any(shape(physical_trace_dot) /= shape(physical_trace)) .or. &
            any(shape(external_target_dot) /= shape(external_target)) .or. &
            size(weights_dot) /= size(weights) .or. &
            any(shape(residual_dot) /= shape(physical_trace))) return
        if (present(sheet_current) .neqv. present(sheet_current_dot)) return
        if (present(sheet_current_dot)) then
            if (any(shape(sheet_current_dot) /= shape(sheet_current))) return
            if (.not. all(ieee_is_finite(sheet_current_dot))) return
        end if
        if (.not. all(ieee_is_finite(physical_trace_dot)) .or. &
            .not. all(ieee_is_finite(external_target_dot)) .or. &
            .not. all(ieee_is_finite(weights_dot))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_direction

end module fortfem_free_boundary_port_residual
