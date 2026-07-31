module fortfem_nonlinear_surface_flux
    !! Physics-agnostic material-surface residual and derivative contract.
    !!
    !! Applications provide the local trace flux law.  FortFEM owns only the
    !! oriented quadrature, finite-element trace pairing, per-tag surface
    !! ledger, and value/JVP/VJP bookkeeping.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none

    private

    abstract interface
        subroutine surface_flux_value_law(state, normal, tag, flux, status)
            import dp
            real(dp), intent(in) :: state(:), normal(3)
            integer, intent(in) :: tag
            real(dp), intent(out) :: flux(:)
            integer, intent(out) :: status
        end subroutine surface_flux_value_law

        subroutine surface_flux_jvp_law( &
                state, normal, tag, state_dot, normal_dot, flux_dot, status)
            import dp
            real(dp), intent(in) :: state(:), normal(3), state_dot(:), normal_dot(3)
            integer, intent(in) :: tag
            real(dp), intent(out) :: flux_dot(:)
            integer, intent(out) :: status
        end subroutine surface_flux_jvp_law

        subroutine surface_flux_vjp_law( &
                state, normal, tag, flux_bar, state_bar, normal_bar, status)
            import dp
            real(dp), intent(in) :: state(:), normal(3), flux_bar(:)
            integer, intent(in) :: tag
            real(dp), intent(out) :: state_bar(:), normal_bar(3)
            integer, intent(out) :: status
        end subroutine surface_flux_vjp_law
    end interface

    public :: assemble_nonlinear_surface_flux
    public :: assemble_nonlinear_surface_flux_jvp
    public :: assemble_nonlinear_surface_flux_vjp

contains

    subroutine assemble_nonlinear_surface_flux( &
            trace_basis, surface_weights, surface_normals, surface_tags, &
            trace_state, residual, surface_ledger, flux_law, status)
        !! Assemble an application-supplied nonlinear surface flux law.
        real(dp), intent(in) :: trace_basis(:, :, :), surface_weights(:)
        real(dp), intent(in) :: surface_normals(:, :), trace_state(:, :)
        integer, intent(in) :: surface_tags(:)
        real(dp), intent(out) :: residual(:, :), surface_ledger(:, :)
        procedure(surface_flux_value_law) :: flux_law
        integer, intent(out) :: status

        integer :: quadrature, dof, component, law_status
        integer :: quadrature_count, dof_count, component_count
        real(dp), allocatable :: flux(:)

        residual = 0.0_dp
        surface_ledger = 0.0_dp
        status = 1
        if (.not. valid_value_inputs( &
            trace_basis, surface_weights, surface_normals, surface_tags, &
            trace_state, residual, surface_ledger)) return
        quadrature_count = size(trace_basis, 1)
        dof_count = size(trace_basis, 2)
        component_count = size(trace_basis, 3)
        allocate(flux(component_count))
        do quadrature = 1, quadrature_count
            call flux_law( &
                trace_state(quadrature, :), surface_normals(:, quadrature), &
                surface_tags(quadrature), flux, law_status)
            if (law_status /= 0 .or. any(.not. ieee_is_finite(flux))) return
            do component = 1, component_count
                surface_ledger(surface_tags(quadrature), component) = &
                    surface_ledger(surface_tags(quadrature), component) + &
                    surface_weights(quadrature)*flux(component)
                do dof = 1, dof_count
                    residual(dof, component) = residual(dof, component) + &
                        trace_basis(quadrature, dof, component)* &
                        surface_weights(quadrature)*flux(component)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_nonlinear_surface_flux

    subroutine assemble_nonlinear_surface_flux_jvp( &
            trace_basis, surface_weights, surface_normals, surface_tags, &
            trace_state, trace_basis_dot, surface_weights_dot, surface_normals_dot, &
            trace_state_dot, residual_dot, surface_ledger_dot, flux_law, &
            flux_jvp_law, status)
        !! Apply the fixed-tag JVP of a nonlinear surface flux assembly.
        real(dp), intent(in) :: trace_basis(:, :, :), surface_weights(:)
        real(dp), intent(in) :: surface_normals(:, :), trace_state(:, :)
        integer, intent(in) :: surface_tags(:)
        real(dp), intent(in) :: trace_basis_dot(:, :, :), surface_weights_dot(:)
        real(dp), intent(in) :: surface_normals_dot(:, :), trace_state_dot(:, :)
        real(dp), intent(out) :: residual_dot(:, :), surface_ledger_dot(:, :)
        procedure(surface_flux_value_law) :: flux_law
        procedure(surface_flux_jvp_law) :: flux_jvp_law
        integer, intent(out) :: status

        integer :: quadrature, dof, component, law_status
        integer :: quadrature_count, dof_count, component_count
        real(dp), allocatable :: flux(:), flux_dot(:), weighted_flux_dot(:)

        residual_dot = 0.0_dp
        surface_ledger_dot = 0.0_dp
        status = 1
        if (.not. valid_value_inputs( &
            trace_basis, surface_weights, surface_normals, surface_tags, &
            trace_state, residual_dot, surface_ledger_dot)) return
        if (.not. valid_jvp_inputs( &
            trace_basis, surface_weights, surface_normals, trace_state, &
            trace_basis_dot, surface_weights_dot, surface_normals_dot, &
            trace_state_dot)) return
        quadrature_count = size(trace_basis, 1)
        dof_count = size(trace_basis, 2)
        component_count = size(trace_basis, 3)
        allocate(flux(component_count), flux_dot(component_count), &
            weighted_flux_dot(component_count))
        do quadrature = 1, quadrature_count
            call flux_law( &
                trace_state(quadrature, :), surface_normals(:, quadrature), &
                surface_tags(quadrature), flux, law_status)
            if (law_status /= 0 .or. any(.not. ieee_is_finite(flux))) return
            call flux_jvp_law( &
                trace_state(quadrature, :), surface_normals(:, quadrature), &
                surface_tags(quadrature), trace_state_dot(quadrature, :), &
                surface_normals_dot(:, quadrature), flux_dot, law_status)
            if (law_status /= 0 .or. any(.not. ieee_is_finite(flux_dot))) return
            weighted_flux_dot = surface_weights_dot(quadrature)*flux + &
                surface_weights(quadrature)*flux_dot
            do component = 1, component_count
                surface_ledger_dot(surface_tags(quadrature), component) = &
                    surface_ledger_dot(surface_tags(quadrature), component) + &
                    weighted_flux_dot(component)
                do dof = 1, dof_count
                    residual_dot(dof, component) = residual_dot(dof, component) + &
                        trace_basis_dot(quadrature, dof, component)* &
                        surface_weights(quadrature)*flux(component) + &
                        trace_basis(quadrature, dof, component)* &
                        weighted_flux_dot(component)
                end do
            end do
        end do
        status = 0
    end subroutine assemble_nonlinear_surface_flux_jvp

    subroutine assemble_nonlinear_surface_flux_vjp( &
            trace_basis, surface_weights, surface_normals, surface_tags, &
            trace_state, residual_bar, surface_ledger_bar, trace_basis_bar, &
            surface_weights_bar, surface_normals_bar, trace_state_bar, flux_law, &
            flux_vjp_law, status)
        !! Apply the real fixed-tag VJP of a nonlinear surface flux assembly.
        real(dp), intent(in) :: trace_basis(:, :, :), surface_weights(:)
        real(dp), intent(in) :: surface_normals(:, :), trace_state(:, :)
        integer, intent(in) :: surface_tags(:)
        real(dp), intent(in) :: residual_bar(:, :), surface_ledger_bar(:, :)
        real(dp), intent(out) :: trace_basis_bar(:, :, :), surface_weights_bar(:)
        real(dp), intent(out) :: surface_normals_bar(:, :), trace_state_bar(:, :)
        procedure(surface_flux_value_law) :: flux_law
        procedure(surface_flux_vjp_law) :: flux_vjp_law
        integer, intent(out) :: status

        integer :: quadrature, dof, component, law_status
        integer :: quadrature_count, dof_count, component_count
        real(dp), allocatable :: flux(:), flux_bar(:), state_bar(:)
        real(dp) :: normal_bar(3), flux_bar_component

        trace_basis_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        surface_normals_bar = 0.0_dp
        trace_state_bar = 0.0_dp
        status = 1
        if (.not. valid_value_inputs( &
            trace_basis, surface_weights, surface_normals, surface_tags, &
            trace_state, residual_bar, surface_ledger_bar)) return
        if (any(.not. ieee_is_finite(residual_bar)) .or. &
            any(.not. ieee_is_finite(surface_ledger_bar)) .or. &
            any(shape(trace_basis_bar) /= shape(trace_basis)) .or. &
            any(shape(surface_normals_bar) /= shape(surface_normals)) .or. &
            any(shape(trace_state_bar) /= shape(trace_state)) .or. &
            size(surface_weights_bar) /= size(surface_weights)) return
        quadrature_count = size(trace_basis, 1)
        dof_count = size(trace_basis, 2)
        component_count = size(trace_basis, 3)
        allocate(flux(component_count), flux_bar(component_count), &
            state_bar(component_count))
        do quadrature = 1, quadrature_count
            call flux_law( &
                trace_state(quadrature, :), surface_normals(:, quadrature), &
                surface_tags(quadrature), flux, law_status)
            if (law_status /= 0 .or. any(.not. ieee_is_finite(flux))) return
            flux_bar = 0.0_dp
            do component = 1, component_count
                flux_bar_component = surface_ledger_bar( &
                    surface_tags(quadrature), component)
                do dof = 1, dof_count
                    flux_bar_component = flux_bar_component + &
                        trace_basis(quadrature, dof, component)*residual_bar(dof, component)
                    trace_basis_bar(quadrature, dof, component) = &
                        trace_basis_bar(quadrature, dof, component) + &
                        residual_bar(dof, component)*surface_weights(quadrature)*flux(component)
                end do
                flux_bar(component) = surface_weights(quadrature)*flux_bar_component
                surface_weights_bar(quadrature) = surface_weights_bar(quadrature) + &
                    flux(component)*flux_bar_component
            end do
            call flux_vjp_law( &
                trace_state(quadrature, :), surface_normals(:, quadrature), &
                surface_tags(quadrature), flux_bar, state_bar, normal_bar, law_status)
            if (law_status /= 0 .or. any(.not. ieee_is_finite(state_bar)) .or. &
                any(.not. ieee_is_finite(normal_bar))) return
            trace_state_bar(quadrature, :) = trace_state_bar(quadrature, :) + state_bar
            surface_normals_bar(:, quadrature) = &
                surface_normals_bar(:, quadrature) + normal_bar
        end do
        status = 0
    end subroutine assemble_nonlinear_surface_flux_vjp

    pure logical function valid_value_inputs( &
            trace_basis, surface_weights, surface_normals, surface_tags, &
            trace_state, residual, surface_ledger) result(valid)
        real(dp), intent(in) :: trace_basis(:, :, :), surface_weights(:)
        real(dp), intent(in) :: surface_normals(:, :), trace_state(:, :)
        integer, intent(in) :: surface_tags(:)
        real(dp), intent(in) :: residual(:, :), surface_ledger(:, :)

        integer :: quadrature_count, dof_count, component_count

        valid = .false.
        quadrature_count = size(trace_basis, 1)
        dof_count = size(trace_basis, 2)
        component_count = size(trace_basis, 3)
        if (quadrature_count < 1 .or. dof_count < 1 .or. component_count < 1) return
        if (size(surface_weights) /= quadrature_count .or. &
            size(surface_normals, 1) /= 3 .or. size(surface_normals, 2) /= quadrature_count .or. &
            size(surface_tags) /= quadrature_count .or. &
            size(trace_state, 1) /= quadrature_count .or. &
            size(trace_state, 2) /= component_count .or. &
            size(residual, 1) /= dof_count .or. size(residual, 2) /= component_count .or. &
            size(surface_ledger, 1) < 1 .or. size(surface_ledger, 2) /= component_count) return
        if (any(surface_tags < 1) .or. any(surface_tags > size(surface_ledger, 1))) return
        if (any(surface_weights <= 0.0_dp)) return
        if (any(.not. ieee_is_finite(trace_basis)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(surface_normals)) .or. &
            any(.not. ieee_is_finite(trace_state))) return
        valid = .true.
    end function valid_value_inputs

    pure logical function valid_jvp_inputs( &
            trace_basis, surface_weights, surface_normals, trace_state, &
            trace_basis_dot, surface_weights_dot, surface_normals_dot, &
            trace_state_dot) result(valid)
        real(dp), intent(in) :: trace_basis(:, :, :), surface_weights(:)
        real(dp), intent(in) :: surface_normals(:, :), trace_state(:, :)
        real(dp), intent(in) :: trace_basis_dot(:, :, :), surface_weights_dot(:)
        real(dp), intent(in) :: surface_normals_dot(:, :), trace_state_dot(:, :)

        valid = all(shape(trace_basis_dot) == shape(trace_basis)) .and. &
            all(shape(surface_normals_dot) == shape(surface_normals)) .and. &
            all(shape(trace_state_dot) == shape(trace_state)) .and. &
            size(surface_weights_dot) == size(surface_weights)
        if (.not. valid) return
        if (any(.not. ieee_is_finite(trace_basis_dot)) .or. &
            any(.not. ieee_is_finite(surface_weights_dot)) .or. &
            any(.not. ieee_is_finite(surface_normals_dot)) .or. &
            any(.not. ieee_is_finite(trace_state_dot))) valid = .false.
    end function valid_jvp_inputs

end module fortfem_nonlinear_surface_flux
