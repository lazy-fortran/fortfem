module fortfem_fci_power_flux_ledger
    !! Conservative parallel/perpendicular FCI power and flux ledger.
    !!
    !! The parallel support flux is `-K_parallel Q u` on staggered FCI
    !! supports.  Its signed power is the conservative pairing
    !! `-(Q u)^T W_s K_parallel (Q u)`.  The perpendicular action is supplied
    !! by the caller as a canonical-cell action; its signed power is
    !! `u^T W_c A_perpendicular u`.  Geometry, interpolation, material traces,
    !! and boundary conditions remain caller-owned.  Thus this contract can be
    !! used with FEM, IGA, Fourier, or plane-local DG blocks without choosing
    !! a plasma closure.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_generated_fci_parallel_flux_power, only: &
        generated_fci_parallel_flux_power
    use fortfem_generated_fci_parallel_flux_power_jvp, only: &
        generated_fci_parallel_flux_power_jvp
    use fortfem_generated_fci_parallel_flux_power_vjp, only: &
        generated_fci_parallel_flux_power_vjp
    use fortfem_generated_fci_perpendicular_power, only: &
        generated_fci_perpendicular_power
    use fortfem_generated_fci_perpendicular_power_jvp, only: &
        generated_fci_perpendicular_power_jvp
    use fortfem_generated_fci_perpendicular_power_vjp, only: &
        generated_fci_perpendicular_power_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_fci_power_flux_ledger
    public :: evaluate_fci_power_flux_ledger_jvp
    public :: evaluate_fci_power_flux_ledger_vjp

contains

    subroutine evaluate_fci_power_flux_ledger( &
            gradient, parallel_coefficient, staggered_volumes, field, &
            perpendicular_action, canonical_volumes, parallel_flux, &
            parallel_power, perpendicular_power, total_power, status)
        real(dp), intent(in) :: gradient(:), parallel_coefficient(:)
        real(dp), intent(in) :: staggered_volumes(:), field(:)
        real(dp), intent(in) :: perpendicular_action(:), canonical_volumes(:)
        real(dp), intent(out) :: parallel_flux(:)
        real(dp), intent(out) :: parallel_power, perpendicular_power, total_power
        type(fortsparse_status_t), intent(out) :: status

        integer :: support, cell
        real(dp) :: flux, power

        parallel_flux = 0.0_dp
        parallel_power = 0.0_dp
        perpendicular_power = 0.0_dp
        total_power = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI power/flux ledger received incompatible arrays")
        if (.not. valid_value_inputs( &
            gradient, parallel_coefficient, staggered_volumes, field, &
            perpendicular_action, canonical_volumes, parallel_flux)) return

        do support = 1, size(gradient)
            call generated_fci_parallel_flux_power( &
                gradient(support), parallel_coefficient(support), &
                staggered_volumes(support), flux, power)
            parallel_flux(support) = flux
            parallel_power = parallel_power + power
        end do
        do cell = 1, size(field)
            call generated_fci_perpendicular_power( &
                field(cell), perpendicular_action(cell), canonical_volumes(cell), &
                power)
            perpendicular_power = perpendicular_power + power
        end do
        total_power = parallel_power + perpendicular_power
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fci_power_flux_ledger

    subroutine evaluate_fci_power_flux_ledger_jvp( &
            gradient, parallel_coefficient, staggered_volumes, field, &
            perpendicular_action, canonical_volumes, gradient_dot, &
            parallel_coefficient_dot, staggered_volumes_dot, field_dot, &
            perpendicular_action_dot, canonical_volumes_dot, parallel_flux_dot, &
            parallel_power_dot, perpendicular_power_dot, total_power_dot, status)
        real(dp), intent(in) :: gradient(:), parallel_coefficient(:)
        real(dp), intent(in) :: staggered_volumes(:), field(:)
        real(dp), intent(in) :: perpendicular_action(:), canonical_volumes(:)
        real(dp), intent(in) :: gradient_dot(:), parallel_coefficient_dot(:)
        real(dp), intent(in) :: staggered_volumes_dot(:), field_dot(:)
        real(dp), intent(in) :: perpendicular_action_dot(:), canonical_volumes_dot(:)
        real(dp), intent(out) :: parallel_flux_dot(:)
        real(dp), intent(out) :: parallel_power_dot, perpendicular_power_dot
        real(dp), intent(out) :: total_power_dot
        type(fortsparse_status_t), intent(out) :: status

        integer :: support, cell
        real(dp) :: flux_dot, power_dot

        parallel_flux_dot = 0.0_dp
        parallel_power_dot = 0.0_dp
        perpendicular_power_dot = 0.0_dp
        total_power_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI power/flux ledger JVP received incompatible arrays")
        if (.not. valid_value_inputs( &
            gradient, parallel_coefficient, staggered_volumes, field, &
            perpendicular_action, canonical_volumes, parallel_flux_dot)) return
        if (.not. valid_tangent_inputs( &
            gradient_dot, parallel_coefficient_dot, staggered_volumes_dot, &
            field_dot, perpendicular_action_dot, canonical_volumes_dot, &
            size(gradient), size(field))) return

        do support = 1, size(gradient)
            call generated_fci_parallel_flux_power_jvp( &
                gradient(support), parallel_coefficient(support), &
                staggered_volumes(support), gradient_dot(support), &
                parallel_coefficient_dot(support), staggered_volumes_dot(support), &
                flux_dot, power_dot)
            parallel_flux_dot(support) = flux_dot
            parallel_power_dot = parallel_power_dot + power_dot
        end do
        do cell = 1, size(field)
            call generated_fci_perpendicular_power_jvp( &
                field(cell), perpendicular_action(cell), canonical_volumes(cell), &
                field_dot(cell), perpendicular_action_dot(cell), &
                canonical_volumes_dot(cell), power_dot)
            perpendicular_power_dot = perpendicular_power_dot + power_dot
        end do
        total_power_dot = parallel_power_dot + perpendicular_power_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fci_power_flux_ledger_jvp

    subroutine evaluate_fci_power_flux_ledger_vjp( &
            gradient, parallel_coefficient, staggered_volumes, field, &
            perpendicular_action, canonical_volumes, parallel_flux_bar, &
            parallel_power_bar, perpendicular_power_bar, total_power_bar, &
            gradient_bar, parallel_coefficient_bar, staggered_volumes_bar, &
            field_bar, perpendicular_action_bar, canonical_volumes_bar, status)
        real(dp), intent(in) :: gradient(:), parallel_coefficient(:)
        real(dp), intent(in) :: staggered_volumes(:), field(:)
        real(dp), intent(in) :: perpendicular_action(:), canonical_volumes(:)
        real(dp), intent(in) :: parallel_flux_bar(:)
        real(dp), intent(in) :: parallel_power_bar, perpendicular_power_bar
        real(dp), intent(in) :: total_power_bar
        real(dp), intent(out) :: gradient_bar(:), parallel_coefficient_bar(:)
        real(dp), intent(out) :: staggered_volumes_bar(:), field_bar(:)
        real(dp), intent(out) :: perpendicular_action_bar(:), canonical_volumes_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: support, cell
        real(dp) :: local_gradient_bar, local_coefficient_bar
        real(dp) :: local_staggered_volume_bar
        real(dp) :: local_field_bar, local_action_bar, local_canonical_volume_bar
        real(dp) :: effective_parallel_bar, effective_perpendicular_bar

        gradient_bar = 0.0_dp
        parallel_coefficient_bar = 0.0_dp
        staggered_volumes_bar = 0.0_dp
        field_bar = 0.0_dp
        perpendicular_action_bar = 0.0_dp
        canonical_volumes_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI power/flux ledger VJP received incompatible arrays")
        if (.not. valid_value_inputs( &
            gradient, parallel_coefficient, staggered_volumes, field, &
            perpendicular_action, canonical_volumes, parallel_flux_bar)) return
        if (size(parallel_flux_bar) /= size(gradient) .or. &
            .not. ieee_is_finite(parallel_power_bar) .or. &
            .not. ieee_is_finite(perpendicular_power_bar) .or. &
            .not. ieee_is_finite(total_power_bar)) return

        effective_parallel_bar = parallel_power_bar + total_power_bar
        effective_perpendicular_bar = perpendicular_power_bar + total_power_bar
        do support = 1, size(gradient)
            call generated_fci_parallel_flux_power_vjp( &
                gradient(support), parallel_coefficient(support), &
                staggered_volumes(support), parallel_flux_bar(support), &
                effective_parallel_bar, local_gradient_bar, local_coefficient_bar, &
                local_staggered_volume_bar)
            gradient_bar(support) = local_gradient_bar
            parallel_coefficient_bar(support) = local_coefficient_bar
            staggered_volumes_bar(support) = local_staggered_volume_bar
        end do
        do cell = 1, size(field)
            call generated_fci_perpendicular_power_vjp( &
                field(cell), perpendicular_action(cell), canonical_volumes(cell), &
                effective_perpendicular_bar, local_field_bar, local_action_bar, &
                local_canonical_volume_bar)
            field_bar(cell) = local_field_bar
            perpendicular_action_bar(cell) = local_action_bar
            canonical_volumes_bar(cell) = local_canonical_volume_bar
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fci_power_flux_ledger_vjp

    logical function valid_value_inputs( &
            gradient, parallel_coefficient, staggered_volumes, field, &
            perpendicular_action, canonical_volumes, parallel_flux) result(valid)
        real(dp), intent(in) :: gradient(:), parallel_coefficient(:)
        real(dp), intent(in) :: staggered_volumes(:), field(:)
        real(dp), intent(in) :: perpendicular_action(:), canonical_volumes(:)
        real(dp), intent(in) :: parallel_flux(:)

        valid = size(gradient) > 0 .and. size(parallel_coefficient) == size(gradient) &
            .and. size(staggered_volumes) == size(gradient) .and. size(field) > 0 &
            .and. size(perpendicular_action) == size(field) &
            .and. size(canonical_volumes) == size(field) .and. &
            size(parallel_flux) == size(gradient) .and. &
            all(ieee_is_finite(gradient)) .and. &
            all(ieee_is_finite(parallel_coefficient)) .and. &
            all(ieee_is_finite(staggered_volumes)) .and. all(ieee_is_finite(field)) &
            .and. all(ieee_is_finite(perpendicular_action)) .and. &
            all(ieee_is_finite(canonical_volumes)) .and. &
            all(ieee_is_finite(parallel_flux)) .and. &
            all(parallel_coefficient >= 0.0_dp) .and. &
            all(staggered_volumes > 0.0_dp) .and. all(canonical_volumes > 0.0_dp)
    end function valid_value_inputs

    logical function valid_tangent_inputs( &
            gradient_dot, parallel_coefficient_dot, staggered_volumes_dot, &
            field_dot, perpendicular_action_dot, canonical_volumes_dot, &
            gradient_count, field_count) result(valid)
        real(dp), intent(in) :: gradient_dot(:), parallel_coefficient_dot(:)
        real(dp), intent(in) :: staggered_volumes_dot(:), field_dot(:)
        real(dp), intent(in) :: perpendicular_action_dot(:), canonical_volumes_dot(:)
        integer, intent(in) :: gradient_count, field_count

        valid = size(gradient_dot) == gradient_count .and. &
            size(parallel_coefficient_dot) == gradient_count .and. &
            size(staggered_volumes_dot) == gradient_count .and. &
            size(field_dot) == field_count .and. &
            size(perpendicular_action_dot) == field_count .and. &
            size(canonical_volumes_dot) == field_count .and. &
            all(ieee_is_finite(gradient_dot)) .and. &
            all(ieee_is_finite(parallel_coefficient_dot)) .and. &
            all(ieee_is_finite(staggered_volumes_dot)) .and. &
            all(ieee_is_finite(field_dot)) .and. &
            all(ieee_is_finite(perpendicular_action_dot)) .and. &
            all(ieee_is_finite(canonical_volumes_dot))
    end function valid_tangent_inputs

end module fortfem_fci_power_flux_ledger
