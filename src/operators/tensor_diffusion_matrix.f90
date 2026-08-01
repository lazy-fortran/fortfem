module fortfem_tensor_diffusion_matrix
    !! Tensor-weighted scalar diffusion contraction at fixed quadrature points.
    !!
    !! For basis gradients g_i and a caller-owned constitutive tensor K, the
    !! returned dense block is
    !!
    !!   A_ij = sum_q w_q g_i(q)^T K(q) g_j(q).
    !!
    !! The block is neutral: K may be a strongly anisotropic conductivity,
    !! resistivity, elastic tangent, or another symmetric or nonsymmetric
    !! tensor.  Geometry and constitutive clients own the quadrature data.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_tensor_diffusion_matrix
    public :: assemble_tensor_diffusion_matrix_jvp
    public :: assemble_tensor_diffusion_matrix_vjp

contains

    subroutine assemble_tensor_diffusion_matrix( &
            basis_gradients, tensor, weights, matrix, status)
        real(dp), intent(in) :: basis_gradients(:, :, :), tensor(:, :, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, first_basis, second_basis, row, column

        matrix = 0.0_dp
        if (.not. validate_inputs( &
            basis_gradients, tensor, weights, matrix, status)) return
        do quadrature = 1, size(weights)
            do first_basis = 1, size(matrix, 1)
                do second_basis = 1, size(matrix, 2)
                    do row = 1, 2
                        do column = 1, 2
                            matrix(first_basis, second_basis) = &
                                matrix(first_basis, second_basis) + &
                                weights(quadrature)*basis_gradients( &
                                quadrature, first_basis, row)*tensor( &
                                quadrature, row, column)*basis_gradients( &
                                quadrature, second_basis, column)
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_tensor_diffusion_matrix

    subroutine assemble_tensor_diffusion_matrix_jvp( &
            basis_gradients, tensor, weights, basis_gradients_dot, tensor_dot, &
            weights_dot, matrix_dot, status)
        real(dp), intent(in) :: basis_gradients(:, :, :), tensor(:, :, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: basis_gradients_dot(:, :, :), tensor_dot(:, :, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: matrix_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, first_basis, second_basis, row, column

        matrix_dot = 0.0_dp
        if (.not. validate_inputs( &
            basis_gradients, tensor, weights, matrix_dot, status)) return
        if (.not. validate_direction( &
            basis_gradients_dot, tensor_dot, weights_dot, &
            size(weights), size(matrix_dot, 1), size(matrix_dot, 2))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "tensor diffusion matrix JVP has incompatible increments")
            return
        end if
        do quadrature = 1, size(weights)
            do first_basis = 1, size(matrix_dot, 1)
                do second_basis = 1, size(matrix_dot, 2)
                    do row = 1, 2
                        do column = 1, 2
                            matrix_dot(first_basis, second_basis) = &
                                matrix_dot(first_basis, second_basis) + &
                                weights_dot(quadrature)*basis_gradients( &
                                quadrature, first_basis, row)*tensor( &
                                quadrature, row, column)*basis_gradients( &
                                quadrature, second_basis, column) + &
                                weights(quadrature)*basis_gradients_dot( &
                                quadrature, first_basis, row)*tensor( &
                                quadrature, row, column)*basis_gradients( &
                                quadrature, second_basis, column) + &
                                weights(quadrature)*basis_gradients( &
                                quadrature, first_basis, row)*tensor_dot( &
                                quadrature, row, column)*basis_gradients( &
                                quadrature, second_basis, column) + &
                                weights(quadrature)*basis_gradients( &
                                quadrature, first_basis, row)*tensor( &
                                quadrature, row, column)*basis_gradients_dot( &
                                quadrature, second_basis, column)
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_tensor_diffusion_matrix_jvp

    subroutine assemble_tensor_diffusion_matrix_vjp( &
            basis_gradients, tensor, weights, matrix_bar, basis_gradients_bar, &
            tensor_bar, weights_bar, status)
        real(dp), intent(in) :: basis_gradients(:, :, :), tensor(:, :, :)
        real(dp), intent(in) :: weights(:), matrix_bar(:, :)
        real(dp), intent(out) :: basis_gradients_bar(:, :, :), tensor_bar(:, :, :)
        real(dp), intent(out) :: weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, first_basis, second_basis, row, column
        real(dp) :: cotangent

        basis_gradients_bar = 0.0_dp
        tensor_bar = 0.0_dp
        weights_bar = 0.0_dp
        if (.not. validate_inputs( &
            basis_gradients, tensor, weights, matrix_bar, status)) return
        if (.not. validate_adjoint( &
            basis_gradients_bar, tensor_bar, weights_bar, matrix_bar, &
            size(weights), size(matrix_bar, 1), size(matrix_bar, 2))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "tensor diffusion matrix VJP has incompatible cotangents")
            return
        end if
        do quadrature = 1, size(weights)
            do first_basis = 1, size(matrix_bar, 1)
                do second_basis = 1, size(matrix_bar, 2)
                    cotangent = matrix_bar(first_basis, second_basis)
                    do row = 1, 2
                        do column = 1, 2
                            basis_gradients_bar(quadrature, first_basis, row) = &
                                basis_gradients_bar(quadrature, first_basis, row) + &
                                cotangent*weights(quadrature)*tensor( &
                                quadrature, row, column)*basis_gradients( &
                                quadrature, second_basis, column)
                            basis_gradients_bar(quadrature, second_basis, column) = &
                                basis_gradients_bar(quadrature, second_basis, column) + &
                                cotangent*weights(quadrature)*basis_gradients( &
                                quadrature, first_basis, row)*tensor( &
                                quadrature, row, column)
                            tensor_bar(quadrature, row, column) = &
                                tensor_bar(quadrature, row, column) + cotangent*&
                                weights(quadrature)*basis_gradients( &
                                quadrature, first_basis, row)*basis_gradients( &
                                quadrature, second_basis, column)
                            weights_bar(quadrature) = weights_bar(quadrature) + &
                                cotangent*basis_gradients(quadrature, first_basis, row)*&
                                tensor(quadrature, row, column)*basis_gradients( &
                                quadrature, second_basis, column)
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_tensor_diffusion_matrix_vjp

    logical function validate_inputs( &
            basis_gradients, tensor, weights, matrix, status) result(valid)
        real(dp), intent(in) :: basis_gradients(:, :, :), tensor(:, :, :)
        real(dp), intent(in) :: weights(:), matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: quadrature_count, basis_count

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "tensor diffusion matrix has incompatible arrays")
        quadrature_count = size(basis_gradients, 1)
        basis_count = size(basis_gradients, 2)
        if (quadrature_count < 1 .or. basis_count < 1 .or. &
            size(basis_gradients, 3) /= 2 .or. &
            size(tensor, 1) /= quadrature_count .or. &
            size(tensor, 2) /= 2 .or. size(tensor, 3) /= 2 .or. &
            size(weights) /= quadrature_count .or. &
            size(matrix, 1) /= basis_count .or. size(matrix, 2) /= basis_count) return
        if (any(.not. ieee_is_finite(basis_gradients)) .or. &
            any(.not. ieee_is_finite(tensor)) .or. &
            any(.not. ieee_is_finite(weights)) .or. any(weights <= 0.0_dp) .or. &
            any(.not. ieee_is_finite(matrix))) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_inputs

    logical function validate_direction( &
            basis_gradients_dot, tensor_dot, weights_dot, quadrature_count, &
            basis_count, matrix_count) result(valid)
        real(dp), intent(in) :: basis_gradients_dot(:, :, :), tensor_dot(:, :, :)
        real(dp), intent(in) :: weights_dot(:)
        integer, intent(in) :: quadrature_count, basis_count, matrix_count

        valid = size(basis_gradients_dot, 1) == quadrature_count .and. &
            size(basis_gradients_dot, 2) == basis_count .and. &
            size(basis_gradients_dot, 3) == 2 .and. &
            size(tensor_dot, 1) == quadrature_count .and. &
            size(tensor_dot, 2) == 2 .and. size(tensor_dot, 3) == 2 .and. &
            size(weights_dot) == quadrature_count .and. matrix_count == basis_count &
            .and. all(ieee_is_finite(basis_gradients_dot)) .and. &
            all(ieee_is_finite(tensor_dot)) .and. all(ieee_is_finite(weights_dot))
    end function validate_direction

    logical function validate_adjoint( &
            basis_gradients_bar, tensor_bar, weights_bar, matrix_bar, &
            quadrature_count, basis_count, matrix_count) result(valid)
        real(dp), intent(in) :: basis_gradients_bar(:, :, :), tensor_bar(:, :, :)
        real(dp), intent(in) :: weights_bar(:), matrix_bar(:, :)
        integer, intent(in) :: quadrature_count, basis_count, matrix_count

        valid = size(basis_gradients_bar, 1) == quadrature_count .and. &
            size(basis_gradients_bar, 2) == basis_count .and. &
            size(basis_gradients_bar, 3) == 2 .and. &
            size(tensor_bar, 1) == quadrature_count .and. &
            size(tensor_bar, 2) == 2 .and. size(tensor_bar, 3) == 2 .and. &
            size(weights_bar) == quadrature_count .and. &
            size(matrix_bar, 1) == basis_count .and. &
            size(matrix_bar, 2) == matrix_count .and. &
            all(ieee_is_finite(matrix_bar))
    end function validate_adjoint

end module fortfem_tensor_diffusion_matrix
