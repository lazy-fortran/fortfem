module fortfem_piola_enriched_vector
    !! Batched Piola mapping followed by pointwise vector enrichment.
    !!
    !! The two map kinds are the standard covariant H(curl) and contravariant
    !! H(div) transforms.  Reference values, Jacobians, and activation values
    !! are all caller-owned so the same contract applies to simplex FEEC,
    !! multipatch IGA, and unfitted XFEM/XIGA composition.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det2, det2_jvp, det2_vjp, det3, det3_jvp, &
        det3_vjp, inv2, inv2_jvp, inv2_vjp, inv3, inv3_jvp, inv3_vjp, LINALG_OK
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    integer, parameter, public :: PIOLA_COVARIANT = 1
    integer, parameter, public :: PIOLA_CONTRAVARIANT = 2

    public :: evaluate_piola_enriched_vector_values
    public :: evaluate_piola_enriched_vector_values_jvp
    public :: evaluate_piola_enriched_vector_values_vjp

contains

    subroutine evaluate_piola_enriched_vector_values( &
            map_kind, jacobians, reference_values, activation, physical_values, &
            status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: activation(:, :)
        real(dp), intent(out) :: physical_values(:, :, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: mapped_2(2, size(reference_values, 3))
        real(dp) :: mapped_3(3, size(reference_values, 3))
        integer :: point, basis, local_status

        physical_values = 0.0_dp
        call validate_inputs( &
            map_kind, jacobians, reference_values, activation, physical_values, &
            status)
        if (status%code /= FORTSPARSE_OK) return

        do point = 1, size(reference_values, 2)
            if (size(reference_values, 1) == 2) then
                call map_value_2d( &
                    map_kind, jacobians(:, :, point), reference_values(:, point, :), &
                    mapped_2, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola enrichment failed in the 2D map")
                    return
                end if
                do basis = 1, size(reference_values, 3)
                    physical_values(:, point, basis) = &
                        activation(point, basis)*mapped_2(:, basis)
                end do
            else
                call map_value_3d( &
                    map_kind, jacobians(:, :, point), reference_values(:, point, :), &
                    mapped_3, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola enrichment failed in the 3D map")
                    return
                end if
                do basis = 1, size(reference_values, 3)
                    physical_values(:, point, basis) = &
                        activation(point, basis)*mapped_3(:, basis)
                end do
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_piola_enriched_vector_values

    subroutine evaluate_piola_enriched_vector_values_jvp( &
            map_kind, jacobians, reference_values, activation, jacobians_dot, &
            reference_values_dot, activation_dot, physical_values_dot, status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: activation(:, :)
        real(dp), intent(in) :: jacobians_dot(:, :, :), reference_values_dot(:, :, :)
        real(dp), intent(in) :: activation_dot(:, :)
        real(dp), intent(out) :: physical_values_dot(:, :, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: mapped_2(2, size(reference_values, 3))
        real(dp) :: mapped_dot_2(2, size(reference_values, 3))
        real(dp) :: mapped_3(3, size(reference_values, 3))
        real(dp) :: mapped_dot_3(3, size(reference_values, 3))
        integer :: point, basis, local_status

        physical_values_dot = 0.0_dp
        call validate_inputs( &
            map_kind, jacobians, reference_values, activation, physical_values_dot, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. validate_directions( &
            jacobians, reference_values, activation, jacobians_dot, &
            reference_values_dot, activation_dot, physical_values_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Piola enrichment JVP has incompatible increments")
            return
        end if

        do point = 1, size(reference_values, 2)
            if (size(reference_values, 1) == 2) then
                call map_value_jvp_2d( &
                    map_kind, jacobians(:, :, point), reference_values(:, point, :), &
                    jacobians_dot(:, :, point), reference_values_dot(:, point, :), &
                    mapped_2, mapped_dot_2, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola enrichment JVP failed in the 2D map")
                    return
                end if
                do basis = 1, size(reference_values, 3)
                    physical_values_dot(:, point, basis) = &
                        activation_dot(point, basis)*mapped_2(:, basis) + &
                        activation(point, basis)*mapped_dot_2(:, basis)
                end do
            else
                call map_value_jvp_3d( &
                    map_kind, jacobians(:, :, point), reference_values(:, point, :), &
                    jacobians_dot(:, :, point), reference_values_dot(:, point, :), &
                    mapped_3, mapped_dot_3, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola enrichment JVP failed in the 3D map")
                    return
                end if
                do basis = 1, size(reference_values, 3)
                    physical_values_dot(:, point, basis) = &
                        activation_dot(point, basis)*mapped_3(:, basis) + &
                        activation(point, basis)*mapped_dot_3(:, basis)
                end do
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_piola_enriched_vector_values_jvp

    subroutine evaluate_piola_enriched_vector_values_vjp( &
            map_kind, jacobians, reference_values, activation, physical_values_bar, &
            jacobians_bar, reference_values_bar, activation_bar, status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: activation(:, :), physical_values_bar(:, :, :)
        real(dp), intent(out) :: jacobians_bar(:, :, :), reference_values_bar(:, :, :)
        real(dp), intent(out) :: activation_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: mapped_2(2, size(reference_values, 3))
        real(dp) :: mapped_bar_2(2, size(reference_values, 3))
        real(dp) :: mapped_3(3, size(reference_values, 3))
        real(dp) :: mapped_bar_3(3, size(reference_values, 3))
        real(dp) :: jacobian_bar_2(2, 2), reference_bar_2(2, size(reference_values, 3))
        real(dp) :: jacobian_bar_3(3, 3), reference_bar_3(3, size(reference_values, 3))
        integer :: point, basis, local_status

        jacobians_bar = 0.0_dp
        reference_values_bar = 0.0_dp
        activation_bar = 0.0_dp
        call validate_inputs( &
            map_kind, jacobians, reference_values, activation, physical_values_bar, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. validate_adjoint( &
            jacobians, reference_values, activation, jacobians_bar, &
            reference_values_bar, activation_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Piola enrichment VJP has incompatible cotangents")
            return
        end if

        do point = 1, size(reference_values, 2)
            if (size(reference_values, 1) == 2) then
                call map_value_2d( &
                    map_kind, jacobians(:, :, point), reference_values(:, point, :), &
                    mapped_2, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola enrichment VJP failed in the 2D map")
                    return
                end if
                do basis = 1, size(reference_values, 3)
                    mapped_bar_2(:, basis) = activation(point, basis)* &
                        physical_values_bar(:, point, basis)
                    activation_bar(point, basis) = dot_product( &
                        physical_values_bar(:, point, basis), mapped_2(:, basis))
                end do
                call map_reverse_2d( &
                    map_kind, jacobians(:, :, point), reference_values(:, point, :), &
                    mapped_bar_2, jacobian_bar_2, reference_bar_2, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola enrichment VJP reverse failed in the 2D map")
                    return
                end if
                jacobians_bar(:, :, point) = jacobian_bar_2
                reference_values_bar(:, point, :) = reference_bar_2
            else
                call map_value_3d( &
                    map_kind, jacobians(:, :, point), reference_values(:, point, :), &
                    mapped_3, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola enrichment VJP failed in the 3D map")
                    return
                end if
                do basis = 1, size(reference_values, 3)
                    mapped_bar_3(:, basis) = activation(point, basis)* &
                        physical_values_bar(:, point, basis)
                    activation_bar(point, basis) = dot_product( &
                        physical_values_bar(:, point, basis), mapped_3(:, basis))
                end do
                call map_reverse_3d( &
                    map_kind, jacobians(:, :, point), reference_values(:, point, :), &
                    mapped_bar_3, jacobian_bar_3, reference_bar_3, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola enrichment VJP reverse failed in the 3D map")
                    return
                end if
                jacobians_bar(:, :, point) = jacobian_bar_3
                reference_values_bar(:, point, :) = reference_bar_3
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_piola_enriched_vector_values_vjp

    subroutine map_value_2d(map_kind, jacobian, reference, mapped, status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobian(2, 2), reference(:, :)
        real(dp), intent(out) :: mapped(:, :)
        integer, intent(out) :: status
        real(dp) :: determinant, inverse(2, 2)
        integer :: basis, inverse_status

        mapped = 0.0_dp
        status = LINALG_OK
        determinant = det2(jacobian)
        if (.not. positive_jacobian_2d(jacobian, determinant)) then
            status = 1
            return
        end if
        if (map_kind == PIOLA_COVARIANT) then
            call inv2(jacobian, inverse, inverse_status)
            if (inverse_status /= LINALG_OK) then
                status = inverse_status
                return
            end if
            do basis = 1, size(reference, 2)
                mapped(:, basis) = matmul(transpose(inverse), reference(:, basis))
            end do
        else if (map_kind == PIOLA_CONTRAVARIANT) then
            do basis = 1, size(reference, 2)
                mapped(:, basis) = matmul(jacobian, reference(:, basis))/determinant
            end do
        else
            status = 1
        end if
    end subroutine map_value_2d

    subroutine map_value_3d(map_kind, jacobian, reference, mapped, status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobian(3, 3), reference(:, :)
        real(dp), intent(out) :: mapped(:, :)
        integer, intent(out) :: status
        real(dp) :: determinant, inverse(3, 3)
        integer :: basis, inverse_status

        mapped = 0.0_dp
        status = LINALG_OK
        determinant = det3(jacobian)
        if (.not. positive_jacobian_3d(jacobian, determinant)) then
            status = 1
            return
        end if
        if (map_kind == PIOLA_COVARIANT) then
            call inv3(jacobian, inverse, inverse_status)
            if (inverse_status /= LINALG_OK) then
                status = inverse_status
                return
            end if
            do basis = 1, size(reference, 2)
                mapped(:, basis) = matmul(transpose(inverse), reference(:, basis))
            end do
        else if (map_kind == PIOLA_CONTRAVARIANT) then
            do basis = 1, size(reference, 2)
                mapped(:, basis) = matmul(jacobian, reference(:, basis))/determinant
            end do
        else
            status = 1
        end if
    end subroutine map_value_3d

    subroutine map_value_jvp_2d( &
            map_kind, jacobian, reference, jacobian_dot, reference_dot, mapped, &
            mapped_dot, status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobian(2, 2), reference(:, :), jacobian_dot(2, 2)
        real(dp), intent(in) :: reference_dot(:, :)
        real(dp), intent(out) :: mapped(:, :), mapped_dot(:, :)
        integer, intent(out) :: status
        real(dp) :: determinant, determinant_dot, inverse(2, 2), inverse_dot(2, 2)
        real(dp) :: mapped_raw(2)
        integer :: basis, inverse_status

        mapped = 0.0_dp
        mapped_dot = 0.0_dp
        status = LINALG_OK
        determinant = det2(jacobian)
        if (.not. positive_jacobian_2d(jacobian, determinant)) then
            status = 1
            return
        end if
        if (map_kind == PIOLA_COVARIANT) then
            call inv2_jvp( &
                jacobian, jacobian_dot, inverse, inverse_dot, inverse_status)
            if (inverse_status /= LINALG_OK) then
                status = inverse_status
                return
            end if
            do basis = 1, size(reference, 2)
                mapped(:, basis) = matmul(transpose(inverse), reference(:, basis))
                mapped_dot(:, basis) = matmul(transpose(inverse_dot), &
                    reference(:, basis)) + matmul(transpose(inverse), &
                    reference_dot(:, basis))
            end do
        else if (map_kind == PIOLA_CONTRAVARIANT) then
            call det2_jvp(jacobian, jacobian_dot, determinant_dot)
            do basis = 1, size(reference, 2)
                mapped_raw = matmul(jacobian, reference(:, basis))
                mapped(:, basis) = mapped_raw/determinant
                mapped_dot(:, basis) = (matmul(jacobian_dot, reference(:, basis)) + &
                    matmul(jacobian, reference_dot(:, basis)))/determinant - &
                    mapped_raw*determinant_dot/determinant**2
            end do
        else
            status = 1
        end if
    end subroutine map_value_jvp_2d

    subroutine map_value_jvp_3d( &
            map_kind, jacobian, reference, jacobian_dot, reference_dot, mapped, &
            mapped_dot, status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobian(3, 3), reference(:, :), jacobian_dot(3, 3)
        real(dp), intent(in) :: reference_dot(:, :)
        real(dp), intent(out) :: mapped(:, :), mapped_dot(:, :)
        integer, intent(out) :: status
        real(dp) :: determinant, determinant_dot, inverse(3, 3), inverse_dot(3, 3)
        real(dp) :: mapped_raw(3)
        integer :: basis, inverse_status

        mapped = 0.0_dp
        mapped_dot = 0.0_dp
        status = LINALG_OK
        determinant = det3(jacobian)
        if (.not. positive_jacobian_3d(jacobian, determinant)) then
            status = 1
            return
        end if
        if (map_kind == PIOLA_COVARIANT) then
            call inv3_jvp( &
                jacobian, jacobian_dot, inverse, inverse_dot, inverse_status)
            if (inverse_status /= LINALG_OK) then
                status = inverse_status
                return
            end if
            do basis = 1, size(reference, 2)
                mapped(:, basis) = matmul(transpose(inverse), reference(:, basis))
                mapped_dot(:, basis) = matmul(transpose(inverse_dot), &
                    reference(:, basis)) + matmul(transpose(inverse), &
                    reference_dot(:, basis))
            end do
        else if (map_kind == PIOLA_CONTRAVARIANT) then
            call det3_jvp(jacobian, jacobian_dot, determinant_dot)
            do basis = 1, size(reference, 2)
                mapped_raw = matmul(jacobian, reference(:, basis))
                mapped(:, basis) = mapped_raw/determinant
                mapped_dot(:, basis) = (matmul(jacobian_dot, reference(:, basis)) + &
                    matmul(jacobian, reference_dot(:, basis)))/determinant - &
                    mapped_raw*determinant_dot/determinant**2
            end do
        else
            status = 1
        end if
    end subroutine map_value_jvp_3d

    subroutine map_reverse_2d( &
            map_kind, jacobian, reference, mapped_bar, jacobian_bar, &
            reference_bar, status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobian(2, 2), reference(:, :), mapped_bar(:, :)
        real(dp), intent(out) :: jacobian_bar(2, 2), reference_bar(:, :)
        integer, intent(out) :: status
        real(dp) :: determinant, determinant_bar, inverse(2, 2), inverse_bar(2, 2)
        real(dp) :: determinant_jacobian_bar(2, 2), mapped_raw(2)
        integer :: basis, inverse_status

        jacobian_bar = 0.0_dp
        reference_bar = 0.0_dp
        status = LINALG_OK
        determinant = det2(jacobian)
        if (.not. positive_jacobian_2d(jacobian, determinant)) then
            status = 1
            return
        end if
        if (map_kind == PIOLA_COVARIANT) then
            call inv2(jacobian, inverse, inverse_status)
            if (inverse_status /= LINALG_OK) then
                status = inverse_status
                return
            end if
            inverse_bar = 0.0_dp
            do basis = 1, size(reference, 2)
                reference_bar(:, basis) = matmul(inverse, mapped_bar(:, basis))
                inverse_bar = inverse_bar + spread(reference(:, basis), 2, 2)* &
                    spread(mapped_bar(:, basis), 1, 2)
            end do
            call inv2_vjp( &
                jacobian, inverse_bar, inverse, jacobian_bar, inverse_status)
            if (inverse_status /= LINALG_OK) status = inverse_status
        else if (map_kind == PIOLA_CONTRAVARIANT) then
            determinant_bar = 0.0_dp
            do basis = 1, size(reference, 2)
                mapped_raw = matmul(jacobian, reference(:, basis))
                reference_bar(:, basis) = &
                    matmul(transpose(jacobian), mapped_bar(:, basis))/determinant
                jacobian_bar = jacobian_bar + spread(mapped_bar(:, basis), 2, 2)* &
                    spread(reference(:, basis), 1, 2)/determinant
                determinant_bar = determinant_bar - &
                    dot_product(mapped_bar(:, basis), mapped_raw)/determinant**2
            end do
            call det2_vjp(jacobian, determinant_bar, determinant_jacobian_bar)
            jacobian_bar = jacobian_bar + determinant_jacobian_bar
        else
            status = 1
        end if
    end subroutine map_reverse_2d

    subroutine map_reverse_3d( &
            map_kind, jacobian, reference, mapped_bar, jacobian_bar, &
            reference_bar, status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobian(3, 3), reference(:, :), mapped_bar(:, :)
        real(dp), intent(out) :: jacobian_bar(3, 3), reference_bar(:, :)
        integer, intent(out) :: status
        real(dp) :: determinant, determinant_bar, inverse(3, 3), inverse_bar(3, 3)
        real(dp) :: determinant_jacobian_bar(3, 3), mapped_raw(3)
        integer :: basis, inverse_status

        jacobian_bar = 0.0_dp
        reference_bar = 0.0_dp
        status = LINALG_OK
        determinant = det3(jacobian)
        if (.not. positive_jacobian_3d(jacobian, determinant)) then
            status = 1
            return
        end if
        if (map_kind == PIOLA_COVARIANT) then
            call inv3(jacobian, inverse, inverse_status)
            if (inverse_status /= LINALG_OK) then
                status = inverse_status
                return
            end if
            inverse_bar = 0.0_dp
            do basis = 1, size(reference, 2)
                reference_bar(:, basis) = matmul(inverse, mapped_bar(:, basis))
                inverse_bar = inverse_bar + spread(reference(:, basis), 2, 3)* &
                    spread(mapped_bar(:, basis), 1, 3)
            end do
            call inv3_vjp( &
                jacobian, inverse_bar, inverse, jacobian_bar, inverse_status)
            if (inverse_status /= LINALG_OK) status = inverse_status
        else if (map_kind == PIOLA_CONTRAVARIANT) then
            determinant_bar = 0.0_dp
            do basis = 1, size(reference, 2)
                mapped_raw = matmul(jacobian, reference(:, basis))
                reference_bar(:, basis) = &
                    matmul(transpose(jacobian), mapped_bar(:, basis))/determinant
                jacobian_bar = jacobian_bar + spread(mapped_bar(:, basis), 2, 3)* &
                    spread(reference(:, basis), 1, 3)/determinant
                determinant_bar = determinant_bar - &
                    dot_product(mapped_bar(:, basis), mapped_raw)/determinant**2
            end do
            call det3_vjp(jacobian, determinant_bar, determinant_jacobian_bar)
            jacobian_bar = jacobian_bar + determinant_jacobian_bar
        else
            status = 1
        end if
    end subroutine map_reverse_3d

    subroutine validate_inputs( &
            map_kind, jacobians, reference_values, activation, target, status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: activation(:, :), target(:, :, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: dimension, points, basis, point

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Piola enrichment has incompatible arrays")
        if (map_kind /= PIOLA_COVARIANT .and. &
            map_kind /= PIOLA_CONTRAVARIANT) return
        dimension = size(reference_values, 1)
        points = size(reference_values, 2)
        basis = size(reference_values, 3)
        if (dimension < 2 .or. dimension > 3) return
        if (points < 1 .or. basis < 1) return
        if (size(jacobians, 1) /= dimension) return
        if (size(jacobians, 2) /= dimension) return
        if (size(jacobians, 3) /= points) return
        if (size(activation, 1) /= points) return
        if (size(activation, 2) /= basis) return
        if (size(target, 1) /= dimension) return
        if (size(target, 2) /= points) return
        if (size(target, 3) /= basis) return
        if (any(.not. ieee_is_finite(jacobians))) return
        if (any(.not. ieee_is_finite(reference_values))) return
        if (any(.not. ieee_is_finite(activation))) return
        if (any(.not. ieee_is_finite(target))) return
        do point = 1, points
            if (dimension == 2) then
                if (.not. positive_jacobian_2d( &
                    jacobians(:, :, point), det2(jacobians(:, :, point)))) return
            else
                if (.not. positive_jacobian_3d( &
                    jacobians(:, :, point), det3(jacobians(:, :, point)))) return
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_inputs

    logical function validate_directions( &
            jacobians, reference_values, activation, jacobians_dot, &
            reference_values_dot, activation_dot, target) result(valid)
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: activation(:, :), jacobians_dot(:, :, :)
        real(dp), intent(in) :: reference_values_dot(:, :, :), activation_dot(:, :)
        real(dp), intent(in) :: target(:, :, :)

        valid = size(jacobians_dot, 1) == size(jacobians, 1) .and. &
            size(jacobians_dot, 2) == size(jacobians, 2) .and. &
            size(jacobians_dot, 3) == size(jacobians, 3) .and. &
            size(reference_values_dot, 1) == size(reference_values, 1) .and. &
            size(reference_values_dot, 2) == size(reference_values, 2) .and. &
            size(reference_values_dot, 3) == size(reference_values, 3) .and. &
            size(activation_dot, 1) == size(activation, 1) .and. &
            size(activation_dot, 2) == size(activation, 2) .and. &
            all(ieee_is_finite(jacobians_dot)) .and. &
            all(ieee_is_finite(reference_values_dot)) .and. &
            all(ieee_is_finite(activation_dot)) .and. &
            size(target, 1) == size(reference_values, 1)
    end function validate_directions

    logical function validate_adjoint( &
            jacobians, reference_values, activation, jacobians_bar, &
            reference_values_bar, activation_bar) result(valid)
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: activation(:, :)
        real(dp), intent(in) :: jacobians_bar(:, :, :), reference_values_bar(:, :, :)
        real(dp), intent(in) :: activation_bar(:, :)

        valid = size(jacobians_bar, 1) == size(jacobians, 1) .and. &
            size(jacobians_bar, 2) == size(jacobians, 2) .and. &
            size(jacobians_bar, 3) == size(jacobians, 3) .and. &
            size(reference_values_bar, 1) == size(reference_values, 1) .and. &
            size(reference_values_bar, 2) == size(reference_values, 2) .and. &
            size(reference_values_bar, 3) == size(reference_values, 3) .and. &
            size(activation_bar, 1) == size(activation, 1) .and. &
            size(activation_bar, 2) == size(activation, 2)
    end function validate_adjoint

    pure logical function positive_jacobian_2d(jacobian, determinant) result(valid)
        real(dp), intent(in) :: jacobian(2, 2), determinant
        real(dp) :: scale, tolerance

        scale = max(1.0_dp, maxval(abs(jacobian)))
        tolerance = 64.0_dp*epsilon(1.0_dp)*scale**2
        valid = determinant > tolerance
    end function positive_jacobian_2d

    pure logical function positive_jacobian_3d(jacobian, determinant) result(valid)
        real(dp), intent(in) :: jacobian(3, 3), determinant
        real(dp) :: scale, tolerance

        scale = max(1.0_dp, maxval(abs(jacobian)))
        tolerance = 64.0_dp*epsilon(1.0_dp)*scale**3
        valid = determinant > tolerance
    end function positive_jacobian_3d

end module fortfem_piola_enriched_vector
