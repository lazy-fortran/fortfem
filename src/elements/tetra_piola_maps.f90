module fortfem_tetra_piola_maps
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3, det3_jvp, det3_vjp, inv3, inv3_jvp, &
        inv3_vjp
    implicit none

    private

    public :: map_tetra_nedelec_covariant
    public :: map_tetra_nedelec_covariant_jvp
    public :: map_tetra_nedelec_covariant_vjp
    public :: map_tetra_rt_contravariant

contains

    pure subroutine map_tetra_nedelec_covariant( &
            jacobian, reference_values, reference_curls, physical_values, &
            physical_curls, status)
        real(dp), intent(in) :: jacobian(3, 3)
        real(dp), intent(in) :: reference_values(:, :)
        real(dp), intent(in) :: reference_curls(:, :)
        real(dp), intent(out) :: physical_values(:, :)
        real(dp), intent(out) :: physical_curls(:, :)
        integer, intent(out) :: status

        real(dp) :: determinant, inverse_jacobian(3, 3), tolerance
        integer :: basis, dof_count, inverse_status

        physical_values = 0.0_dp
        physical_curls = 0.0_dp
        status = 1
        dof_count = size(reference_values, 2)
        if (size(reference_values, 1) /= 3) return
        if (size(reference_curls, 1) /= 3) return
        if (size(reference_curls, 2) /= dof_count) return
        if (size(physical_values, 1) /= 3) return
        if (size(physical_values, 2) /= dof_count) return
        if (size(physical_curls, 1) /= 3) return
        if (size(physical_curls, 2) /= dof_count) return

        determinant = det3(jacobian)
        tolerance = 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, maxval(abs(jacobian))**3)
        if (determinant <= tolerance) return
        call inv3(jacobian, inverse_jacobian, inverse_status)
        if (inverse_status /= 0) return

        do basis = 1, dof_count
            physical_values(:, basis) = matmul( &
                transpose(inverse_jacobian), reference_values(:, basis))
            physical_curls(:, basis) = &
                matmul(jacobian, reference_curls(:, basis)) / determinant
        end do
        status = 0
    end subroutine map_tetra_nedelec_covariant

    pure subroutine map_tetra_nedelec_covariant_jvp( &
            jacobian, reference_values, reference_curls, jacobian_dot, &
            reference_values_dot, reference_curls_dot, physical_values_dot, &
            physical_curls_dot, status)
        real(dp), intent(in) :: jacobian(3, 3), jacobian_dot(3, 3)
        real(dp), intent(in) :: reference_values(:, :), reference_curls(:, :)
        real(dp), intent(in) :: reference_values_dot(:, :)
        real(dp), intent(in) :: reference_curls_dot(:, :)
        real(dp), intent(out) :: physical_values_dot(:, :)
        real(dp), intent(out) :: physical_curls_dot(:, :)
        integer, intent(out) :: status

        real(dp) :: determinant, determinant_dot
        real(dp) :: inverse_jacobian(3, 3), inverse_jacobian_dot(3, 3)
        real(dp) :: mapped_curls(3, size(reference_curls, 2)), tolerance
        integer :: inverse_status

        physical_values_dot = 0.0_dp
        physical_curls_dot = 0.0_dp
        call validate_nedelec_shapes( &
            reference_values, reference_curls, physical_values_dot, &
            physical_curls_dot, status)
        if (status /= 0) return
        if (any(shape(reference_values_dot) /= shape(reference_values)) .or. &
            any(shape(reference_curls_dot) /= shape(reference_curls))) then
            status = 1
            return
        end if
        determinant = det3(jacobian)
        tolerance = 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)
        if (determinant <= tolerance) return
        call det3_jvp(jacobian, jacobian_dot, determinant_dot)
        call inv3_jvp( &
            jacobian, jacobian_dot, inverse_jacobian, inverse_jacobian_dot, &
            inverse_status)
        if (inverse_status /= 0) return
        physical_values_dot = &
            matmul(transpose(inverse_jacobian_dot), reference_values) + &
            matmul(transpose(inverse_jacobian), reference_values_dot)
        mapped_curls = matmul(jacobian, reference_curls)
        physical_curls_dot = ( &
            matmul(jacobian_dot, reference_curls) + &
            matmul(jacobian, reference_curls_dot))/determinant - &
            mapped_curls*determinant_dot/determinant**2
        status = 0
    end subroutine map_tetra_nedelec_covariant_jvp

    pure subroutine map_tetra_nedelec_covariant_vjp( &
            jacobian, reference_values, reference_curls, physical_values_bar, &
            physical_curls_bar, jacobian_bar, reference_values_bar, &
            reference_curls_bar, status)
        real(dp), intent(in) :: jacobian(3, 3)
        real(dp), intent(in) :: reference_values(:, :), reference_curls(:, :)
        real(dp), intent(in) :: physical_values_bar(:, :)
        real(dp), intent(in) :: physical_curls_bar(:, :)
        real(dp), intent(out) :: jacobian_bar(3, 3)
        real(dp), intent(out) :: reference_values_bar(:, :)
        real(dp), intent(out) :: reference_curls_bar(:, :)
        integer, intent(out) :: status

        real(dp) :: determinant, determinant_bar
        real(dp) :: determinant_jacobian_bar(3, 3)
        real(dp) :: inverse_jacobian(3, 3), inverse_jacobian_bar(3, 3)
        real(dp) :: inverse_jacobian_jacobian_bar(3, 3)
        real(dp) :: mapped_curls(3, size(reference_curls, 2))
        real(dp) :: mapped_curls_bar(3, size(reference_curls, 2)), tolerance
        integer :: inverse_status

        jacobian_bar = 0.0_dp
        reference_values_bar = 0.0_dp
        reference_curls_bar = 0.0_dp
        call validate_nedelec_shapes( &
            reference_values, reference_curls, physical_values_bar, &
            physical_curls_bar, status)
        if (status /= 0) return
        if (any(shape(reference_values_bar) /= shape(reference_values)) .or. &
            any(shape(reference_curls_bar) /= shape(reference_curls))) then
            status = 1
            return
        end if
        determinant = det3(jacobian)
        tolerance = 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)
        if (determinant <= tolerance) return

        call inv3(jacobian, inverse_jacobian, inverse_status)
        if (inverse_status /= 0) return
        reference_values_bar = &
            matmul(inverse_jacobian, physical_values_bar)
        inverse_jacobian_bar = &
            matmul(reference_values, transpose(physical_values_bar))
        call inv3_vjp( &
            jacobian, inverse_jacobian_bar, inverse_jacobian, &
            inverse_jacobian_jacobian_bar, inverse_status)
        if (inverse_status /= 0) return

        mapped_curls = matmul(jacobian, reference_curls)
        mapped_curls_bar = physical_curls_bar/determinant
        reference_curls_bar = &
            matmul(transpose(jacobian), mapped_curls_bar)
        jacobian_bar = inverse_jacobian_jacobian_bar + &
            matmul(mapped_curls_bar, transpose(reference_curls))
        determinant_bar = &
            -sum(physical_curls_bar*mapped_curls)/determinant**2
        call det3_vjp(jacobian, determinant_bar, determinant_jacobian_bar)
        jacobian_bar = jacobian_bar + determinant_jacobian_bar
        status = 0
    end subroutine map_tetra_nedelec_covariant_vjp

    pure subroutine validate_nedelec_shapes( &
            reference_values, reference_curls, physical_values, &
            physical_curls, status)
        real(dp), intent(in) :: reference_values(:, :), reference_curls(:, :)
        real(dp), intent(in) :: physical_values(:, :), physical_curls(:, :)
        integer, intent(out) :: status

        status = 1
        if (size(reference_values, 1) /= 3) return
        if (size(reference_curls, 1) /= 3) return
        if (size(reference_curls, 2) /= size(reference_values, 2)) return
        if (any(shape(physical_values) /= shape(reference_values))) return
        if (any(shape(physical_curls) /= shape(reference_curls))) return
        status = 0
    end subroutine validate_nedelec_shapes

    pure subroutine map_tetra_rt_contravariant( &
            jacobian, reference_values, reference_divergences, &
            physical_values, physical_divergences, status)
        real(dp), intent(in) :: jacobian(3, 3)
        real(dp), intent(in) :: reference_values(:, :)
        real(dp), intent(in) :: reference_divergences(:)
        real(dp), intent(out) :: physical_values(:, :)
        real(dp), intent(out) :: physical_divergences(:)
        integer, intent(out) :: status

        real(dp) :: determinant, tolerance
        integer :: basis, dof_count

        physical_values = 0.0_dp
        physical_divergences = 0.0_dp
        status = 1
        dof_count = size(reference_values, 2)
        if (size(reference_values, 1) /= 3) return
        if (size(reference_divergences) /= dof_count) return
        if (size(physical_values, 1) /= 3) return
        if (size(physical_values, 2) /= dof_count) return
        if (size(physical_divergences) /= dof_count) return

        determinant = det3(jacobian)
        tolerance = 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)
        if (determinant <= tolerance) return

        do basis = 1, dof_count
            physical_values(:, basis) = &
                matmul(jacobian, reference_values(:, basis))/determinant
        end do
        physical_divergences = reference_divergences/determinant
        status = 0
    end subroutine map_tetra_rt_contravariant

end module fortfem_tetra_piola_maps
