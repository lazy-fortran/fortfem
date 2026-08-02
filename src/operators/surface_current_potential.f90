module fortfem_surface_current_potential
    !! Neutral stream-function and harmonic-loop surface-current map.
    !!
    !! For caller-owned tangent gradients g_j(q) and loop representatives
    !! h_l(q), this module evaluates
    !!
    !!   K(q) = s [ n(q) x sum_j g_j(q) a_j
    !!              + sum_l h_l(q) c_l ],
    !!
    !! where s is a fixed orientation sign.  The surface quadrature, cuts,
    !! constitutive law, and physical units remain caller-owned.  The map is
    !! therefore usable for virtual-casing, NESTOR/BIEST-like source traces,
    !! fitted/cut IGA wall patches, or an external resistive-wall model.
    !!
    !! The period ledger is a weighted pairing of each supplied harmonic
    !! representative with K.  It does not normalize the representatives or
    !! assign physical flux units.  A single-valued scalar potential needs one
    !! externally selected gauge; the metadata records that fixed choice.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    real(dp), parameter :: unit_tolerance = 1.0e-12_dp
    real(dp), parameter :: tangent_tolerance = 1.0e-10_dp

    type, public :: surface_current_potential_metadata_t
        integer :: harmonic_count = 0
        integer :: period_count = 0
        integer :: gauge_dof = 0
        integer :: orientation_sign = 1
        logical :: topology_fixed = .false.
        logical :: scalar_potential_single_valued = .true.
    end type surface_current_potential_metadata_t

    public :: initialize_surface_current_potential_metadata
    public :: validate_surface_current_potential_metadata
    public :: assemble_surface_current_potential
    public :: assemble_surface_current_potential_jvp
    public :: assemble_surface_current_potential_vjp

contains

    subroutine initialize_surface_current_potential_metadata( &
            metadata, harmonic_count, period_count, gauge_dof, orientation_sign, &
            topology_fixed, scalar_potential_single_valued, status)
        type(surface_current_potential_metadata_t), intent(out) :: metadata
        integer, intent(in) :: harmonic_count, period_count, gauge_dof
        integer, intent(in) :: orientation_sign
        logical, intent(in) :: topology_fixed, scalar_potential_single_valued
        type(fortsparse_status_t), intent(out) :: status

        metadata%harmonic_count = harmonic_count
        metadata%period_count = period_count
        metadata%gauge_dof = gauge_dof
        metadata%orientation_sign = orientation_sign
        metadata%topology_fixed = topology_fixed
        metadata%scalar_potential_single_valued = scalar_potential_single_valued
        call validate_surface_current_potential_metadata(metadata, status)
    end subroutine initialize_surface_current_potential_metadata

    subroutine validate_surface_current_potential_metadata(metadata, status)
        type(surface_current_potential_metadata_t), intent(in) :: metadata
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface-current potential metadata is incompatible")
        if (metadata%harmonic_count < 0 .or. metadata%period_count < 0 .or. &
            metadata%period_count /= metadata%harmonic_count .or. &
            metadata%gauge_dof < 0 .or. metadata%orientation_sign /= 1 .and. &
            metadata%orientation_sign /= -1 .or. .not. metadata%topology_fixed .or. &
            .not. metadata%scalar_potential_single_valued) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_surface_current_potential_metadata

    subroutine assemble_surface_current_potential( &
            tangent_gradients, scalar_coefficients, harmonic_loop_basis, &
            harmonic_coefficients, normals, surface_weights, metadata, &
            surface_current, integrated_current, period_ledger, status)
        real(dp), intent(in) :: tangent_gradients(:, :, :)
        real(dp), intent(in) :: scalar_coefficients(:)
        real(dp), intent(in) :: harmonic_loop_basis(:, :, :)
        real(dp), intent(in) :: harmonic_coefficients(:)
        real(dp), intent(in) :: normals(:, :), surface_weights(:)
        type(surface_current_potential_metadata_t), intent(in) :: metadata
        real(dp), intent(out) :: surface_current(:, :), integrated_current(:)
        real(dp), intent(out) :: period_ledger(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, scalar_dof, harmonic
        real(dp) :: scalar_gradient(3), current(3)

        surface_current = 0.0_dp
        integrated_current = 0.0_dp
        period_ledger = 0.0_dp
        call validate_potential_inputs( &
            tangent_gradients, scalar_coefficients, harmonic_loop_basis, &
            harmonic_coefficients, normals, surface_weights, metadata, &
            surface_current, integrated_current, period_ledger, status)
        if (status%code /= FORTSPARSE_OK) return

        do quadrature = 1, size(surface_weights)
            scalar_gradient = 0.0_dp
            do scalar_dof = 1, size(scalar_coefficients)
                scalar_gradient = scalar_gradient + &
                    tangent_gradients(quadrature, scalar_dof, :) * &
                    scalar_coefficients(scalar_dof)
            end do
            call cross_product(normals(quadrature, :), scalar_gradient, current)
            do harmonic = 1, metadata%harmonic_count
                current = current + harmonic_loop_basis(quadrature, harmonic, :) * &
                    harmonic_coefficients(harmonic)
            end do
            surface_current(quadrature, :) = metadata%orientation_sign*current
            integrated_current = integrated_current + surface_weights(quadrature)* &
                surface_current(quadrature, :)
            do harmonic = 1, metadata%harmonic_count
                period_ledger(harmonic) = period_ledger(harmonic) + &
                    surface_weights(quadrature)*dot_product( &
                    harmonic_loop_basis(quadrature, harmonic, :), &
                    surface_current(quadrature, :))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_potential

    subroutine assemble_surface_current_potential_jvp( &
            tangent_gradients, scalar_coefficients, harmonic_loop_basis, &
            harmonic_coefficients, normals, surface_weights, metadata, &
            tangent_gradients_dot, scalar_coefficients_dot, harmonic_loop_basis_dot, &
            harmonic_coefficients_dot, normals_dot, surface_weights_dot, &
            surface_current_dot, integrated_current_dot, period_ledger_dot, status)
        real(dp), intent(in) :: tangent_gradients(:, :, :), scalar_coefficients(:)
        real(dp), intent(in) :: harmonic_loop_basis(:, :, :), harmonic_coefficients(:)
        real(dp), intent(in) :: normals(:, :), surface_weights(:)
        type(surface_current_potential_metadata_t), intent(in) :: metadata
        real(dp), intent(in) :: tangent_gradients_dot(:, :, :)
        real(dp), intent(in) :: scalar_coefficients_dot(:)
        real(dp), intent(in) :: harmonic_loop_basis_dot(:, :, :)
        real(dp), intent(in) :: harmonic_coefficients_dot(:)
        real(dp), intent(in) :: normals_dot(:, :), surface_weights_dot(:)
        real(dp), intent(out) :: surface_current_dot(:, :)
        real(dp), intent(out) :: integrated_current_dot(:), period_ledger_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, scalar_dof, harmonic
        real(dp) :: scalar_gradient(3), scalar_gradient_dot(3)
        real(dp) :: current(3), current_dot(3)

        surface_current_dot = 0.0_dp
        integrated_current_dot = 0.0_dp
        period_ledger_dot = 0.0_dp
        call validate_potential_inputs( &
            tangent_gradients, scalar_coefficients, harmonic_loop_basis, &
            harmonic_coefficients, normals, surface_weights, metadata, &
            surface_current_dot, integrated_current_dot, period_ledger_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_direction( &
            tangent_gradients_dot, scalar_coefficients_dot, harmonic_loop_basis_dot, &
            harmonic_coefficients_dot, normals_dot, surface_weights_dot, &
            tangent_gradients, scalar_coefficients, harmonic_loop_basis, &
            harmonic_coefficients, normals, surface_weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current potential JVP has incompatible increments")
            return
        end if

        do quadrature = 1, size(surface_weights)
            scalar_gradient = 0.0_dp
            scalar_gradient_dot = 0.0_dp
            do scalar_dof = 1, size(scalar_coefficients)
                scalar_gradient = scalar_gradient + &
                    tangent_gradients(quadrature, scalar_dof, :) * &
                    scalar_coefficients(scalar_dof)
                scalar_gradient_dot = scalar_gradient_dot + &
                    tangent_gradients_dot(quadrature, scalar_dof, :) * &
                    scalar_coefficients(scalar_dof) + &
                    tangent_gradients(quadrature, scalar_dof, :) * &
                    scalar_coefficients_dot(scalar_dof)
            end do
            call cross_product(normals(quadrature, :), scalar_gradient, current)
            call cross_product(normals_dot(quadrature, :), scalar_gradient, current_dot)
            call cross_product(normals(quadrature, :), scalar_gradient_dot, &
                scalar_gradient)
            current_dot = current_dot + scalar_gradient
            do harmonic = 1, metadata%harmonic_count
                current = current + harmonic_loop_basis(quadrature, harmonic, :) * &
                    harmonic_coefficients(harmonic)
                current_dot = current_dot + harmonic_loop_basis_dot( &
                    quadrature, harmonic, :)*harmonic_coefficients(harmonic) + &
                    harmonic_loop_basis(quadrature, harmonic, :)* &
                    harmonic_coefficients_dot(harmonic)
            end do
            surface_current_dot(quadrature, :) = metadata%orientation_sign*current_dot
            integrated_current_dot = integrated_current_dot + &
                surface_weights_dot(quadrature)*metadata%orientation_sign*current + &
                surface_weights(quadrature)*surface_current_dot(quadrature, :)
            do harmonic = 1, metadata%harmonic_count
                period_ledger_dot(harmonic) = period_ledger_dot(harmonic) + &
                    surface_weights_dot(quadrature)*dot_product( &
                    harmonic_loop_basis(quadrature, harmonic, :), &
                    metadata%orientation_sign*current) + &
                    surface_weights(quadrature)*dot_product( &
                    harmonic_loop_basis_dot(quadrature, harmonic, :), &
                    metadata%orientation_sign*current) + &
                    surface_weights(quadrature)*dot_product( &
                    harmonic_loop_basis(quadrature, harmonic, :), &
                    surface_current_dot(quadrature, :))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_potential_jvp

    subroutine assemble_surface_current_potential_vjp( &
            tangent_gradients, scalar_coefficients, harmonic_loop_basis, &
            harmonic_coefficients, normals, surface_weights, metadata, &
            surface_current, integrated_current, period_ledger, surface_current_bar, &
            integrated_current_bar, period_ledger_bar, tangent_gradients_bar, &
            scalar_coefficients_bar, harmonic_loop_basis_bar, &
            harmonic_coefficients_bar, &
            normals_bar, surface_weights_bar, status)
        real(dp), intent(in) :: tangent_gradients(:, :, :), scalar_coefficients(:)
        real(dp), intent(in) :: harmonic_loop_basis(:, :, :), harmonic_coefficients(:)
        real(dp), intent(in) :: normals(:, :), surface_weights(:)
        type(surface_current_potential_metadata_t), intent(in) :: metadata
        real(dp), intent(in) :: surface_current(:, :), integrated_current(:)
        real(dp), intent(in) :: period_ledger(:), surface_current_bar(:, :)
        real(dp), intent(in) :: integrated_current_bar(:), period_ledger_bar(:)
        real(dp), intent(out) :: tangent_gradients_bar(:, :, :)
        real(dp), intent(out) :: scalar_coefficients_bar(:)
        real(dp), intent(out) :: harmonic_loop_basis_bar(:, :, :)
        real(dp), intent(out) :: harmonic_coefficients_bar(:)
        real(dp), intent(out) :: normals_bar(:, :), surface_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, scalar_dof, harmonic
        real(dp) :: scalar_gradient(3), current(3), total_bar(3), gradient_bar(3)

        tangent_gradients_bar = 0.0_dp
        scalar_coefficients_bar = 0.0_dp
        harmonic_loop_basis_bar = 0.0_dp
        harmonic_coefficients_bar = 0.0_dp
        normals_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        call validate_potential_inputs( &
            tangent_gradients, scalar_coefficients, harmonic_loop_basis, &
            harmonic_coefficients, normals, surface_weights, metadata, &
            surface_current, integrated_current, period_ledger, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(surface_current_bar, 1) /= size(surface_weights) .or. &
            size(surface_current_bar, 2) /= 3 .or. &
            size(integrated_current_bar) /= 3 .or. &
            size(period_ledger_bar) /= metadata%harmonic_count .or. &
            size(tangent_gradients_bar, 1) /= size(tangent_gradients, 1) .or. &
            size(tangent_gradients_bar, 2) /= size(tangent_gradients, 2) .or. &
            size(tangent_gradients_bar, 3) /= 3 .or. &
            size(scalar_coefficients_bar) /= size(scalar_coefficients) .or. &
            size(harmonic_loop_basis_bar, 1) /= size(harmonic_loop_basis, 1) .or. &
            size(harmonic_loop_basis_bar, 2) /= metadata%harmonic_count .or. &
            size(harmonic_loop_basis_bar, 3) /= 3 .or. &
            size(harmonic_coefficients_bar) /= metadata%harmonic_count .or. &
            size(normals_bar, 1) /= size(normals, 1) .or. &
            size(normals_bar, 2) /= 3 .or. &
            size(surface_weights_bar) /= size(surface_weights) .or. &
            .not. all(ieee_is_finite(surface_current_bar)) .or. &
            .not. all(ieee_is_finite(integrated_current_bar)) .or. &
            .not. all(ieee_is_finite(period_ledger_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current potential VJP has incompatible cotangents")
            return
        end if

        do quadrature = 1, size(surface_weights)
            scalar_gradient = 0.0_dp
            do scalar_dof = 1, size(scalar_coefficients)
                scalar_gradient = scalar_gradient + &
                    tangent_gradients(quadrature, scalar_dof, :) * &
                    scalar_coefficients(scalar_dof)
            end do
            current = surface_current(quadrature, :)/metadata%orientation_sign
            total_bar = surface_current_bar(quadrature, :) + &
                surface_weights(quadrature)*integrated_current_bar
            do harmonic = 1, metadata%harmonic_count
                total_bar = total_bar + surface_weights(quadrature)* &
                    period_ledger_bar(harmonic)*harmonic_loop_basis( &
                    quadrature, harmonic, :)
                harmonic_loop_basis_bar(quadrature, harmonic, :) = &
                    harmonic_loop_basis_bar(quadrature, harmonic, :) + &
                    surface_weights(quadrature)*period_ledger_bar(harmonic)* &
                    surface_current(quadrature, :)
            end do
            total_bar = metadata%orientation_sign*total_bar
            call cross_product(total_bar, normals(quadrature, :), gradient_bar)
            call cross_product(scalar_gradient, total_bar, normals_bar(quadrature, :))
            do scalar_dof = 1, size(scalar_coefficients)
                tangent_gradients_bar(quadrature, scalar_dof, :) = &
                    tangent_gradients_bar(quadrature, scalar_dof, :) + &
                    scalar_coefficients(scalar_dof)*gradient_bar
                scalar_coefficients_bar(scalar_dof) = &
                    scalar_coefficients_bar(scalar_dof) + dot_product( &
                    tangent_gradients(quadrature, scalar_dof, :), gradient_bar)
            end do
            do harmonic = 1, metadata%harmonic_count
                harmonic_loop_basis_bar(quadrature, harmonic, :) = &
                    harmonic_loop_basis_bar(quadrature, harmonic, :) + &
                    harmonic_coefficients(harmonic)*total_bar
                harmonic_coefficients_bar(harmonic) = &
                    harmonic_coefficients_bar(harmonic) + dot_product( &
                    harmonic_loop_basis(quadrature, harmonic, :), total_bar)
            end do
            surface_weights_bar(quadrature) = dot_product( &
                integrated_current_bar, surface_current(quadrature, :))
            do harmonic = 1, metadata%harmonic_count
                surface_weights_bar(quadrature) = surface_weights_bar(quadrature) + &
                    period_ledger_bar(harmonic)*dot_product( &
                    harmonic_loop_basis(quadrature, harmonic, :), &
                    surface_current(quadrature, :))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_current_potential_vjp

    subroutine validate_potential_inputs( &
            tangent_gradients, scalar_coefficients, harmonic_loop_basis, &
            harmonic_coefficients, normals, surface_weights, metadata, &
            surface_current, integrated_current, period_ledger, status)
        real(dp), intent(in) :: tangent_gradients(:, :, :), scalar_coefficients(:)
        real(dp), intent(in) :: harmonic_loop_basis(:, :, :), harmonic_coefficients(:)
        real(dp), intent(in) :: normals(:, :), surface_weights(:)
        type(surface_current_potential_metadata_t), intent(in) :: metadata
        real(dp), intent(in) :: surface_current(:, :), integrated_current(:)
        real(dp), intent(in) :: period_ledger(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, scalar_dof, harmonic

        call validate_surface_current_potential_metadata(metadata, status)
        if (status%code /= FORTSPARSE_OK) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface-current potential received incompatible arrays")
        if (size(tangent_gradients, 1) < 1 .or. size(tangent_gradients, 2) < 1 .or. &
            size(tangent_gradients, 3) /= 3 .or. &
            size(normals, 1) /= size(tangent_gradients, 1) .or. &
            size(normals, 2) /= 3 .or. size(surface_weights) /= size(normals, 1) .or. &
            size(scalar_coefficients) /= size(tangent_gradients, 2) .or. &
            size(harmonic_loop_basis, 1) /= size(normals, 1) .or. &
            size(harmonic_loop_basis, 2) /= metadata%harmonic_count .or. &
            size(harmonic_loop_basis, 3) /= 3 .or. &
            size(harmonic_coefficients) /= metadata%harmonic_count .or. &
            size(surface_current, 1) /= size(normals, 1) .or. &
            size(surface_current, 2) /= 3 .or. size(integrated_current) /= 3 .or. &
            size(period_ledger) /= metadata%period_count) return
        if (.not. all(ieee_is_finite(tangent_gradients)) .or. &
            .not. all(ieee_is_finite(scalar_coefficients)) .or. &
            .not. all(ieee_is_finite(harmonic_loop_basis)) .or. &
            .not. all(ieee_is_finite(harmonic_coefficients)) .or. &
            .not. all(ieee_is_finite(normals)) .or. &
            .not. all(ieee_is_finite(surface_weights)) .or. &
            .not. all(ieee_is_finite(surface_current)) .or. &
            .not. all(ieee_is_finite(integrated_current)) .or. &
            .not. all(ieee_is_finite(period_ledger)) .or. &
            any(surface_weights <= 0.0_dp)) return
        do quadrature = 1, size(normals, 1)
            if (abs(dot_product(normals(quadrature, :), &
                normals(quadrature, :)) - 1.0_dp) > &
                unit_tolerance) return
            do scalar_dof = 1, size(scalar_coefficients)
                if (abs(dot_product(normals(quadrature, :), &
                    tangent_gradients(quadrature, scalar_dof, :))) > &
                    tangent_tolerance) return
            end do
            do harmonic = 1, metadata%harmonic_count
                if (abs(dot_product(normals(quadrature, :), &
                    harmonic_loop_basis(quadrature, harmonic, :))) > &
                    tangent_tolerance) return
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_potential_inputs

    logical function valid_direction( &
            tangent_gradients_dot, scalar_coefficients_dot, harmonic_loop_basis_dot, &
            harmonic_coefficients_dot, normals_dot, surface_weights_dot, &
            tangent_gradients, scalar_coefficients, harmonic_loop_basis, &
            harmonic_coefficients, normals, surface_weights) result(valid)
        real(dp), intent(in) :: tangent_gradients_dot(:, :, :), &
            scalar_coefficients_dot(:)
        real(dp), intent(in) :: harmonic_loop_basis_dot(:, :, :), &
            harmonic_coefficients_dot(:)
        real(dp), intent(in) :: normals_dot(:, :), surface_weights_dot(:)
        real(dp), intent(in) :: tangent_gradients(:, :, :), scalar_coefficients(:)
        real(dp), intent(in) :: harmonic_loop_basis(:, :, :), harmonic_coefficients(:)
        real(dp), intent(in) :: normals(:, :), surface_weights(:)

        valid = all(shape(tangent_gradients_dot) == shape(tangent_gradients)) .and. &
            size(scalar_coefficients_dot) == size(scalar_coefficients) .and. &
            all(shape(harmonic_loop_basis_dot) == shape(harmonic_loop_basis)) .and. &
            size(harmonic_coefficients_dot) == size(harmonic_coefficients) .and. &
            all(shape(normals_dot) == shape(normals)) .and. &
            size(surface_weights_dot) == size(surface_weights) .and. &
            all(ieee_is_finite(tangent_gradients_dot)) .and. &
            all(ieee_is_finite(scalar_coefficients_dot)) .and. &
            all(ieee_is_finite(harmonic_loop_basis_dot)) .and. &
            all(ieee_is_finite(harmonic_coefficients_dot)) .and. &
            all(ieee_is_finite(normals_dot)) .and. &
            all(ieee_is_finite(surface_weights_dot))
    end function valid_direction

    pure subroutine cross_product(first, second, result)
        real(dp), intent(in) :: first(3), second(3)
        real(dp), intent(out) :: result(3)

        result = [first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end subroutine cross_product

end module fortfem_surface_current_potential
