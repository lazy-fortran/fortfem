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
    public :: evaluate_piola_enriched_vector_differential_3d
    public :: evaluate_piola_enriched_vector_differential_3d_jvp
    public :: evaluate_piola_enriched_vector_differential_3d_vjp

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

    subroutine evaluate_piola_enriched_vector_differential_3d( &
            map_kind, jacobians, reference_values, reference_curls, &
            reference_divergence, activation, activation_gradient, physical_values, &
            curl_values, divergence, status)
        !! Evaluate affine Piola H(curl)/H(div) enrichment differentials.
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: reference_curls(:, :, :), reference_divergence(:, :)
        real(dp), intent(in) :: activation(:, :), activation_gradient(:, :, :)
        real(dp), intent(out) :: physical_values(:, :, :), curl_values(:, :, :)
        real(dp), intent(out) :: divergence(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: mapped(3, size(reference_values, 3))
        real(dp) :: mapped_curl(3, size(reference_values, 3))
        real(dp) :: mapped_divergence
        integer :: point, basis, local_status

        physical_values = 0.0_dp
        curl_values = 0.0_dp
        divergence = 0.0_dp
        call validate_differential_inputs( &
            map_kind, jacobians, reference_values, reference_curls, reference_divergence, &
            activation, activation_gradient, physical_values, curl_values, divergence, status)
        if (status%code /= FORTSPARSE_OK) return
        do point = 1, size(reference_values, 2)
            call map_value_3d( &
                map_kind, jacobians(:, :, point), reference_values(:, point, :), &
                mapped, local_status)
            if (local_status /= LINALG_OK) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "Piola differential value map failed")
                return
            end if
            do basis = 1, size(reference_values, 3)
                physical_values(:, point, basis) = activation(point, basis)*mapped(:, basis)
            end do
            if (map_kind == PIOLA_COVARIANT) then
                call map_value_3d( &
                    PIOLA_CONTRAVARIANT, jacobians(:, :, point), reference_curls(:, point, :), &
                    mapped_curl, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola differential curl map failed")
                    return
                end if
                do basis = 1, size(reference_values, 3)
                    curl_values(:, point, basis) = activation(point, basis)* &
                        mapped_curl(:, basis) + cross3( &
                        activation_gradient(:, point, basis), mapped(:, basis))
                end do
            else
                do basis = 1, size(reference_values, 3)
                    mapped_divergence = reference_divergence(point, basis)/ &
                        det3(jacobians(:, :, point))
                    divergence(point, basis) = activation(point, basis)*mapped_divergence + &
                        dot_product(activation_gradient(:, point, basis), mapped(:, basis))
                end do
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_piola_enriched_vector_differential_3d

    subroutine evaluate_piola_enriched_vector_differential_3d_jvp( &
            map_kind, jacobians, reference_values, reference_curls, reference_divergence, &
            activation, activation_gradient, jacobians_dot, reference_values_dot, &
            reference_curls_dot, reference_divergence_dot, activation_dot, &
            activation_gradient_dot, physical_values_dot, curl_values_dot, divergence_dot, &
            status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: reference_curls(:, :, :), reference_divergence(:, :)
        real(dp), intent(in) :: activation(:, :), activation_gradient(:, :, :)
        real(dp), intent(in) :: jacobians_dot(:, :, :), reference_values_dot(:, :, :)
        real(dp), intent(in) :: reference_curls_dot(:, :, :), reference_divergence_dot(:, :)
        real(dp), intent(in) :: activation_dot(:, :), activation_gradient_dot(:, :, :)
        real(dp), intent(out) :: physical_values_dot(:, :, :), curl_values_dot(:, :, :)
        real(dp), intent(out) :: divergence_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: mapped(3, size(reference_values, 3))
        real(dp) :: mapped_dot(3, size(reference_values, 3))
        real(dp) :: mapped_curl(3, size(reference_values, 3))
        real(dp) :: mapped_curl_dot(3, size(reference_values, 3))
        real(dp) :: determinant, determinant_dot, mapped_divergence, mapped_divergence_dot
        integer :: point, basis, local_status

        physical_values_dot = 0.0_dp
        curl_values_dot = 0.0_dp
        divergence_dot = 0.0_dp
        call validate_differential_inputs( &
            map_kind, jacobians, reference_values, reference_curls, reference_divergence, &
            activation, activation_gradient, physical_values_dot, curl_values_dot, &
            divergence_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. validate_differential_direction( &
            jacobians, reference_values, reference_curls, reference_divergence, activation, &
            activation_gradient, jacobians_dot, reference_values_dot, reference_curls_dot, &
            reference_divergence_dot, activation_dot, activation_gradient_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Piola differential JVP has incompatible increments")
            return
        end if
        do point = 1, size(reference_values, 2)
            call map_value_jvp_3d( &
                map_kind, jacobians(:, :, point), reference_values(:, point, :), &
                jacobians_dot(:, :, point), reference_values_dot(:, point, :), mapped, &
                mapped_dot, local_status)
            if (local_status /= LINALG_OK) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "Piola differential JVP value map failed")
                return
            end if
            if (map_kind == PIOLA_COVARIANT) then
                call map_value_jvp_3d( &
                    PIOLA_CONTRAVARIANT, jacobians(:, :, point), reference_curls(:, point, :), &
                    jacobians_dot(:, :, point), reference_curls_dot(:, point, :), mapped_curl, &
                    mapped_curl_dot, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola differential JVP curl map failed")
                    return
                end if
            else
                determinant = det3(jacobians(:, :, point))
                call det3_jvp(jacobians(:, :, point), jacobians_dot(:, :, point), determinant_dot)
            end if
            do basis = 1, size(reference_values, 3)
                physical_values_dot(:, point, basis) = activation_dot(point, basis)* &
                    mapped(:, basis) + activation(point, basis)*mapped_dot(:, basis)
                if (map_kind == PIOLA_COVARIANT) then
                    curl_values_dot(:, point, basis) = activation_dot(point, basis)* &
                        mapped_curl(:, basis) + activation(point, basis)*mapped_curl_dot(:, basis) + &
                        cross3(activation_gradient_dot(:, point, basis), mapped(:, basis)) + &
                        cross3(activation_gradient(:, point, basis), mapped_dot(:, basis))
                else
                    mapped_divergence = reference_divergence(point, basis)/determinant
                    mapped_divergence_dot = reference_divergence_dot(point, basis)/determinant - &
                        reference_divergence(point, basis)*determinant_dot/determinant**2
                    divergence_dot(point, basis) = activation_dot(point, basis)*mapped_divergence + &
                        activation(point, basis)*mapped_divergence_dot + &
                        dot_product(activation_gradient_dot(:, point, basis), mapped(:, basis)) + &
                        dot_product(activation_gradient(:, point, basis), mapped_dot(:, basis))
                end if
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_piola_enriched_vector_differential_3d_jvp

    subroutine evaluate_piola_enriched_vector_differential_3d_vjp( &
            map_kind, jacobians, reference_values, reference_curls, reference_divergence, &
            activation, activation_gradient, physical_values_bar, curl_values_bar, &
            divergence_bar, jacobians_bar, reference_values_bar, reference_curls_bar, &
            reference_divergence_bar, activation_bar, activation_gradient_bar, status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: reference_curls(:, :, :), reference_divergence(:, :)
        real(dp), intent(in) :: activation(:, :), activation_gradient(:, :, :)
        real(dp), intent(in) :: physical_values_bar(:, :, :), curl_values_bar(:, :, :)
        real(dp), intent(in) :: divergence_bar(:, :)
        real(dp), intent(out) :: jacobians_bar(:, :, :), reference_values_bar(:, :, :)
        real(dp), intent(out) :: reference_curls_bar(:, :, :), reference_divergence_bar(:, :)
        real(dp), intent(out) :: activation_bar(:, :), activation_gradient_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: mapped(3, size(reference_values, 3))
        real(dp) :: mapped_curl(3, size(reference_values, 3))
        real(dp) :: mapped_bar(3, size(reference_values, 3))
        real(dp) :: mapped_curl_bar(3, size(reference_values, 3))
        real(dp) :: jacobian_bar_local(3, 3), reference_bar_local(3, size(reference_values, 3))
        real(dp) :: jacobian_bar_curl(3, 3), reference_curl_bar_local(3, size(reference_values, 3))
        real(dp) :: determinant, determinant_bar, determinant_jacobian_bar(3, 3)
        real(dp) :: mapped_divergence, mapped_divergence_bar
        integer :: point, basis, local_status

        jacobians_bar = 0.0_dp
        reference_values_bar = 0.0_dp
        reference_curls_bar = 0.0_dp
        reference_divergence_bar = 0.0_dp
        activation_bar = 0.0_dp
        activation_gradient_bar = 0.0_dp
        call validate_differential_inputs( &
            map_kind, jacobians, reference_values, reference_curls, reference_divergence, &
            activation, activation_gradient, physical_values_bar, curl_values_bar, &
            divergence_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. validate_differential_adjoint( &
            jacobians, reference_values, reference_curls, reference_divergence, activation, &
            activation_gradient, jacobians_bar, reference_values_bar, reference_curls_bar, &
            reference_divergence_bar, activation_bar, activation_gradient_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Piola differential VJP has incompatible cotangents")
            return
        end if
        do point = 1, size(reference_values, 2)
            call map_value_3d( &
                map_kind, jacobians(:, :, point), reference_values(:, point, :), mapped, local_status)
            if (local_status /= LINALG_OK) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "Piola differential VJP value map failed")
                return
            end if
            mapped_bar = 0.0_dp
            mapped_curl_bar = 0.0_dp
            do basis = 1, size(reference_values, 3)
                mapped_bar(:, basis) = activation(point, basis)* &
                    physical_values_bar(:, point, basis)
                activation_bar(point, basis) = dot_product( &
                    physical_values_bar(:, point, basis), mapped(:, basis))
            end do
            if (map_kind == PIOLA_COVARIANT) then
                call map_value_3d( &
                    PIOLA_CONTRAVARIANT, jacobians(:, :, point), reference_curls(:, point, :), &
                    mapped_curl, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola differential VJP curl map failed")
                    return
                end if
                do basis = 1, size(reference_values, 3)
                    mapped_bar(:, basis) = mapped_bar(:, basis) + &
                        cross3(curl_values_bar(:, point, basis), activation_gradient(:, point, basis))
                    activation_bar(point, basis) = activation_bar(point, basis) + &
                        dot_product(curl_values_bar(:, point, basis), mapped_curl(:, basis))
                    activation_gradient_bar(:, point, basis) = &
                        cross3(mapped(:, basis), curl_values_bar(:, point, basis))
                    mapped_curl_bar(:, basis) = activation(point, basis)* &
                        curl_values_bar(:, point, basis)
                end do
                call map_reverse_3d( &
                    map_kind, jacobians(:, :, point), reference_values(:, point, :), mapped_bar, &
                    jacobian_bar_local, reference_bar_local, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola differential VJP value reverse failed")
                    return
                end if
                jacobians_bar(:, :, point) = jacobians_bar(:, :, point) + jacobian_bar_local
                reference_values_bar(:, point, :) = reference_bar_local
                call map_reverse_3d( &
                    PIOLA_CONTRAVARIANT, jacobians(:, :, point), reference_curls(:, point, :), &
                    mapped_curl_bar, jacobian_bar_curl, reference_curl_bar_local, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola differential VJP curl reverse failed")
                    return
                end if
                jacobians_bar(:, :, point) = jacobians_bar(:, :, point) + jacobian_bar_curl
                reference_curls_bar(:, point, :) = reference_curl_bar_local
            else
                determinant = det3(jacobians(:, :, point))
                do basis = 1, size(reference_values, 3)
                    mapped_divergence = reference_divergence(point, basis)/determinant
                    mapped_divergence_bar = activation(point, basis)*divergence_bar(point, basis)
                    activation_bar(point, basis) = activation_bar(point, basis) + &
                        divergence_bar(point, basis)*mapped_divergence
                    mapped_bar(:, basis) = mapped_bar(:, basis) + &
                        divergence_bar(point, basis)*activation_gradient(:, point, basis)
                    activation_gradient_bar(:, point, basis) = &
                        divergence_bar(point, basis)*mapped(:, basis)
                    reference_divergence_bar(point, basis) = mapped_divergence_bar/determinant
                    determinant_bar = -mapped_divergence_bar*reference_divergence(point, basis)/determinant**2
                    call det3_vjp(jacobians(:, :, point), determinant_bar, determinant_jacobian_bar)
                    jacobians_bar(:, :, point) = jacobians_bar(:, :, point) + determinant_jacobian_bar
                end do
                call map_reverse_3d( &
                    map_kind, jacobians(:, :, point), reference_values(:, point, :), mapped_bar, &
                    jacobian_bar_local, reference_bar_local, local_status)
                if (local_status /= LINALG_OK) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        "Piola differential VJP value reverse failed")
                    return
                end if
                jacobians_bar(:, :, point) = jacobians_bar(:, :, point) + jacobian_bar_local
                reference_values_bar(:, point, :) = reference_bar_local
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_piola_enriched_vector_differential_3d_vjp

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

    subroutine validate_differential_inputs( &
            map_kind, jacobians, reference_values, reference_curls, reference_divergence, &
            activation, activation_gradient, values, curl_values, divergence, status)
        integer, intent(in) :: map_kind
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: reference_curls(:, :, :), reference_divergence(:, :)
        real(dp), intent(in) :: activation(:, :), activation_gradient(:, :, :)
        real(dp), intent(in) :: values(:, :, :), curl_values(:, :, :), divergence(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: points, basis, point

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Piola differential has incompatible arrays")
        if (map_kind /= PIOLA_COVARIANT .and. map_kind /= PIOLA_CONTRAVARIANT) return
        if (size(reference_values, 1) /= 3) return
        points = size(reference_values, 2)
        basis = size(reference_values, 3)
        if (points < 1 .or. basis < 1) return
        if (size(jacobians, 1) /= 3 .or. size(jacobians, 2) /= 3 .or. &
            size(jacobians, 3) /= points .or. size(reference_curls, 1) /= 3 .or. &
            size(reference_curls, 2) /= points .or. size(reference_curls, 3) /= basis .or. &
            size(reference_divergence, 1) /= points .or. &
            size(reference_divergence, 2) /= basis .or. size(activation, 1) /= points .or. &
            size(activation, 2) /= basis .or. size(activation_gradient, 1) /= 3 .or. &
            size(activation_gradient, 2) /= points .or. size(activation_gradient, 3) /= basis .or. &
            size(values, 1) /= 3 .or. size(values, 2) /= points .or. &
            size(values, 3) /= basis .or. size(curl_values, 1) /= 3 .or. &
            size(curl_values, 2) /= points .or. size(curl_values, 3) /= basis .or. &
            size(divergence, 1) /= points .or. size(divergence, 2) /= basis) return
        if (.not. all(ieee_is_finite(jacobians)) .or. &
            .not. all(ieee_is_finite(reference_values)) .or. &
            .not. all(ieee_is_finite(reference_curls)) .or. &
            .not. all(ieee_is_finite(reference_divergence)) .or. &
            .not. all(ieee_is_finite(activation)) .or. &
            .not. all(ieee_is_finite(activation_gradient)) .or. &
            .not. all(ieee_is_finite(values)) .or. .not. all(ieee_is_finite(curl_values)) .or. &
            .not. all(ieee_is_finite(divergence))) return
        do point = 1, points
            if (.not. positive_jacobian_3d( &
                jacobians(:, :, point), det3(jacobians(:, :, point)))) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_differential_inputs

    logical function validate_differential_direction( &
            jacobians, reference_values, reference_curls, reference_divergence, activation, &
            activation_gradient, jacobians_dot, reference_values_dot, reference_curls_dot, &
            reference_divergence_dot, activation_dot, activation_gradient_dot) result(valid)
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: reference_curls(:, :, :), reference_divergence(:, :)
        real(dp), intent(in) :: activation(:, :), activation_gradient(:, :, :)
        real(dp), intent(in) :: jacobians_dot(:, :, :), reference_values_dot(:, :, :)
        real(dp), intent(in) :: reference_curls_dot(:, :, :), reference_divergence_dot(:, :)
        real(dp), intent(in) :: activation_dot(:, :), activation_gradient_dot(:, :, :)

        valid = size(jacobians_dot, 1) == size(jacobians, 1) .and. &
            size(jacobians_dot, 2) == size(jacobians, 2) .and. &
            size(jacobians_dot, 3) == size(jacobians, 3) .and. &
            size(reference_values_dot, 1) == size(reference_values, 1) .and. &
            size(reference_values_dot, 2) == size(reference_values, 2) .and. &
            size(reference_values_dot, 3) == size(reference_values, 3) .and. &
            size(reference_curls_dot, 1) == size(reference_curls, 1) .and. &
            size(reference_curls_dot, 2) == size(reference_curls, 2) .and. &
            size(reference_curls_dot, 3) == size(reference_curls, 3) .and. &
            size(reference_divergence_dot, 1) == size(reference_divergence, 1) .and. &
            size(reference_divergence_dot, 2) == size(reference_divergence, 2) .and. &
            size(activation_dot, 1) == size(activation, 1) .and. &
            size(activation_dot, 2) == size(activation, 2) .and. &
            size(activation_gradient_dot, 1) == size(activation_gradient, 1) .and. &
            size(activation_gradient_dot, 2) == size(activation_gradient, 2) .and. &
            size(activation_gradient_dot, 3) == size(activation_gradient, 3) .and. &
            all(ieee_is_finite(jacobians_dot)) .and. all(ieee_is_finite(reference_values_dot)) .and. &
            all(ieee_is_finite(reference_curls_dot)) .and. all(ieee_is_finite(reference_divergence_dot)) .and. &
            all(ieee_is_finite(activation_dot)) .and. all(ieee_is_finite(activation_gradient_dot))
    end function validate_differential_direction

    logical function validate_differential_adjoint( &
            jacobians, reference_values, reference_curls, reference_divergence, activation, &
            activation_gradient, jacobians_bar, reference_values_bar, reference_curls_bar, &
            reference_divergence_bar, activation_bar, activation_gradient_bar) result(valid)
        real(dp), intent(in) :: jacobians(:, :, :), reference_values(:, :, :)
        real(dp), intent(in) :: reference_curls(:, :, :), reference_divergence(:, :)
        real(dp), intent(in) :: activation(:, :), activation_gradient(:, :, :)
        real(dp), intent(in) :: jacobians_bar(:, :, :), reference_values_bar(:, :, :)
        real(dp), intent(in) :: reference_curls_bar(:, :, :), reference_divergence_bar(:, :)
        real(dp), intent(in) :: activation_bar(:, :), activation_gradient_bar(:, :, :)

        valid = size(jacobians_bar, 1) == size(jacobians, 1) .and. &
            size(jacobians_bar, 2) == size(jacobians, 2) .and. &
            size(jacobians_bar, 3) == size(jacobians, 3) .and. &
            size(reference_values_bar, 1) == size(reference_values, 1) .and. &
            size(reference_values_bar, 2) == size(reference_values, 2) .and. &
            size(reference_values_bar, 3) == size(reference_values, 3) .and. &
            size(reference_curls_bar, 1) == size(reference_curls, 1) .and. &
            size(reference_curls_bar, 2) == size(reference_curls, 2) .and. &
            size(reference_curls_bar, 3) == size(reference_curls, 3) .and. &
            size(reference_divergence_bar, 1) == size(reference_divergence, 1) .and. &
            size(reference_divergence_bar, 2) == size(reference_divergence, 2) .and. &
            size(activation_bar, 1) == size(activation, 1) .and. &
            size(activation_bar, 2) == size(activation, 2) .and. &
            size(activation_gradient_bar, 1) == size(activation_gradient, 1) .and. &
            size(activation_gradient_bar, 2) == size(activation_gradient, 2) .and. &
            size(activation_gradient_bar, 3) == size(activation_gradient, 3)
    end function validate_differential_adjoint

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

    pure function cross3(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3)

        product = [left(2)*right(3) - left(3)*right(2), &
            left(3)*right(1) - left(1)*right(3), &
            left(1)*right(2) - left(2)*right(1)]
    end function cross3

end module fortfem_piola_enriched_vector
