module fortfem_triangle_piola_maps
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: &
        det2, det2_jvp, det2_vjp, inv2_jvp, inv2_vjp
    implicit none

    private

    public :: map_triangle_nedelec_covariant
    public :: map_triangle_nedelec_covariant_jvp
    public :: map_triangle_nedelec_covariant_vjp
    public :: map_triangle_rt_contravariant
    public :: map_triangle_rt_contravariant_jvp
    public :: map_triangle_rt_contravariant_vjp

contains

    pure subroutine map_triangle_nedelec_covariant( &
            jacobian, reference_values, reference_curls, physical_values, &
            physical_curls, status)
        real(dp), intent(in) :: jacobian(2, 2)
        real(dp), intent(in) :: reference_values(:, :), reference_curls(:)
        real(dp), intent(out) :: physical_values(:, :), physical_curls(:)
        integer, intent(out) :: status

        real(dp) :: determinant, tolerance
        integer :: basis_dof, dof_count

        physical_values = 0.0_dp
        physical_curls = 0.0_dp
        status = 1
        dof_count = size(reference_values, 2)
        if (size(reference_values, 1) /= 2) return
        if (size(reference_curls) /= dof_count) return
        if (size(physical_values, 1) /= 2 .or. &
            size(physical_values, 2) /= dof_count) return
        if (size(physical_curls) /= dof_count) return

        determinant = jacobian(1, 1) * jacobian(2, 2) - &
            jacobian(1, 2) * jacobian(2, 1)
        tolerance = 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, maxval(abs(jacobian))**2)
        if (determinant <= tolerance) return

        do basis_dof = 1, dof_count
            physical_values(1, basis_dof) = ( &
                jacobian(2, 2) * reference_values(1, basis_dof) - &
                jacobian(2, 1) * reference_values(2, basis_dof)) / determinant
            physical_values(2, basis_dof) = ( &
                -jacobian(1, 2) * reference_values(1, basis_dof) + &
                jacobian(1, 1) * reference_values(2, basis_dof)) / determinant
        end do
        physical_curls = reference_curls / determinant
        status = 0
    end subroutine map_triangle_nedelec_covariant

    pure subroutine map_triangle_nedelec_covariant_jvp( &
            jacobian, reference_values, reference_curls, jacobian_dot, &
            reference_values_dot, reference_curls_dot, physical_values_dot, &
            physical_curls_dot, status)
        real(dp), intent(in) :: jacobian(2, 2), jacobian_dot(2, 2)
        real(dp), intent(in) :: reference_values(:, :), reference_curls(:)
        real(dp), intent(in) :: reference_values_dot(:, :)
        real(dp), intent(in) :: reference_curls_dot(:)
        real(dp), intent(out) :: physical_values_dot(:, :)
        real(dp), intent(out) :: physical_curls_dot(:)
        integer, intent(out) :: status

        real(dp) :: determinant, determinant_dot
        real(dp) :: inverse(2, 2), inverse_dot(2, 2)
        integer :: basis_dof, dof_count, local_status

        physical_values_dot = 0.0_dp
        physical_curls_dot = 0.0_dp
        call validate_products( &
            reference_values, reference_curls, reference_values_dot, &
            reference_curls_dot, physical_values_dot, physical_curls_dot, &
            status)
        if (status /= 0) return
        call inv2_jvp( &
            jacobian, jacobian_dot, inverse, inverse_dot, local_status)
        if (local_status /= 0) then
            status = 1
            return
        end if
        determinant = det2(jacobian)
        call det2_jvp(jacobian, jacobian_dot, determinant_dot)
        if (.not. valid_determinant(jacobian, determinant)) then
            status = 1
            return
        end if
        dof_count = size(reference_values, 2)
        do basis_dof = 1, dof_count
            physical_values_dot(:, basis_dof) = matmul( &
                transpose(inverse_dot), reference_values(:, basis_dof)) + &
                matmul( &
                transpose(inverse), reference_values_dot(:, basis_dof))
        end do
        physical_curls_dot = reference_curls_dot/determinant - &
            reference_curls*determinant_dot/determinant**2
        status = 0
    end subroutine map_triangle_nedelec_covariant_jvp

    pure subroutine map_triangle_nedelec_covariant_vjp( &
            jacobian, reference_values, reference_curls, physical_values_bar, &
            physical_curls_bar, jacobian_bar, reference_values_bar, &
            reference_curls_bar, status)
        real(dp), intent(in) :: jacobian(2, 2)
        real(dp), intent(in) :: reference_values(:, :), reference_curls(:)
        real(dp), intent(in) :: physical_values_bar(:, :)
        real(dp), intent(in) :: physical_curls_bar(:)
        real(dp), intent(out) :: jacobian_bar(2, 2)
        real(dp), intent(out) :: reference_values_bar(:, :)
        real(dp), intent(out) :: reference_curls_bar(:)
        integer, intent(out) :: status

        real(dp) :: determinant, determinant_bar, determinant_jacobian_bar(2, 2)
        real(dp) :: inverse(2, 2), inverse_bar(2, 2)
        integer :: basis_dof, dof_count, local_status

        jacobian_bar = 0.0_dp
        reference_values_bar = 0.0_dp
        reference_curls_bar = 0.0_dp
        call validate_products( &
            reference_values, reference_curls, physical_values_bar, &
            physical_curls_bar, reference_values_bar, reference_curls_bar, &
            status)
        if (status /= 0) return
        determinant = det2(jacobian)
        if (.not. valid_determinant(jacobian, determinant)) then
            status = 1
            return
        end if
        inverse_bar = 0.0_dp
        dof_count = size(reference_values, 2)
        do basis_dof = 1, dof_count
            inverse_bar = inverse_bar + spread( &
                reference_values(:, basis_dof), 2, 2)*spread( &
                physical_values_bar(:, basis_dof), 1, 2)
        end do
        call inv2_vjp( &
            jacobian, inverse_bar, inverse, jacobian_bar, local_status)
        if (local_status /= 0) then
            status = 1
            return
        end if
        do basis_dof = 1, dof_count
            reference_values_bar(:, basis_dof) = &
                matmul(inverse, physical_values_bar(:, basis_dof))
        end do
        reference_curls_bar = physical_curls_bar/determinant
        determinant_bar = -dot_product( &
            physical_curls_bar, reference_curls)/determinant**2
        call det2_vjp( &
            jacobian, determinant_bar, determinant_jacobian_bar)
        jacobian_bar = jacobian_bar + determinant_jacobian_bar
        status = 0
    end subroutine map_triangle_nedelec_covariant_vjp

    pure subroutine map_triangle_rt_contravariant( &
            jacobian, reference_values, reference_divergences, &
            physical_values, physical_divergences, status)
        real(dp), intent(in) :: jacobian(2, 2)
        real(dp), intent(in) :: reference_values(:, :)
        real(dp), intent(in) :: reference_divergences(:)
        real(dp), intent(out) :: physical_values(:, :)
        real(dp), intent(out) :: physical_divergences(:)
        integer, intent(out) :: status

        real(dp) :: determinant, tolerance
        integer :: basis_dof, dof_count

        physical_values = 0.0_dp
        physical_divergences = 0.0_dp
        status = 1
        dof_count = size(reference_values, 2)
        if (size(reference_values, 1) /= 2) return
        if (size(reference_divergences) /= dof_count) return
        if (size(physical_values, 1) /= 2 .or. &
            size(physical_values, 2) /= dof_count) return
        if (size(physical_divergences) /= dof_count) return

        determinant = jacobian(1, 1) * jacobian(2, 2) - &
            jacobian(1, 2) * jacobian(2, 1)
        tolerance = 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, maxval(abs(jacobian))**2)
        if (determinant <= tolerance) return

        do basis_dof = 1, dof_count
            physical_values(1, basis_dof) = ( &
                jacobian(1, 1) * reference_values(1, basis_dof) + &
                jacobian(1, 2) * reference_values(2, basis_dof)) / determinant
            physical_values(2, basis_dof) = ( &
                jacobian(2, 1) * reference_values(1, basis_dof) + &
                jacobian(2, 2) * reference_values(2, basis_dof)) / determinant
        end do
        physical_divergences = reference_divergences / determinant
        status = 0
    end subroutine map_triangle_rt_contravariant

    pure subroutine map_triangle_rt_contravariant_jvp( &
            jacobian, reference_values, reference_divergences, jacobian_dot, &
            reference_values_dot, reference_divergences_dot, &
            physical_values_dot, physical_divergences_dot, status)
        real(dp), intent(in) :: jacobian(2, 2), jacobian_dot(2, 2)
        real(dp), intent(in) :: reference_values(:, :)
        real(dp), intent(in) :: reference_divergences(:)
        real(dp), intent(in) :: reference_values_dot(:, :)
        real(dp), intent(in) :: reference_divergences_dot(:)
        real(dp), intent(out) :: physical_values_dot(:, :)
        real(dp), intent(out) :: physical_divergences_dot(:)
        integer, intent(out) :: status

        real(dp) :: determinant, determinant_dot
        integer :: basis_dof, dof_count

        physical_values_dot = 0.0_dp
        physical_divergences_dot = 0.0_dp
        call validate_products( &
            reference_values, reference_divergences, reference_values_dot, &
            reference_divergences_dot, physical_values_dot, &
            physical_divergences_dot, status)
        if (status /= 0) return
        determinant = det2(jacobian)
        call det2_jvp(jacobian, jacobian_dot, determinant_dot)
        if (.not. valid_determinant(jacobian, determinant)) then
            status = 1
            return
        end if
        dof_count = size(reference_values, 2)
        do basis_dof = 1, dof_count
            physical_values_dot(:, basis_dof) = ( &
                matmul(jacobian_dot, reference_values(:, basis_dof)) + &
                matmul(jacobian, reference_values_dot(:, basis_dof)))/ &
                determinant - matmul( &
                jacobian, reference_values(:, basis_dof))*determinant_dot/ &
                determinant**2
        end do
        physical_divergences_dot = reference_divergences_dot/determinant - &
            reference_divergences*determinant_dot/determinant**2
        status = 0
    end subroutine map_triangle_rt_contravariant_jvp

    pure subroutine map_triangle_rt_contravariant_vjp( &
            jacobian, reference_values, reference_divergences, &
            physical_values_bar, physical_divergences_bar, jacobian_bar, &
            reference_values_bar, reference_divergences_bar, status)
        real(dp), intent(in) :: jacobian(2, 2)
        real(dp), intent(in) :: reference_values(:, :)
        real(dp), intent(in) :: reference_divergences(:)
        real(dp), intent(in) :: physical_values_bar(:, :)
        real(dp), intent(in) :: physical_divergences_bar(:)
        real(dp), intent(out) :: jacobian_bar(2, 2)
        real(dp), intent(out) :: reference_values_bar(:, :)
        real(dp), intent(out) :: reference_divergences_bar(:)
        integer, intent(out) :: status

        real(dp) :: determinant, determinant_bar, determinant_jacobian_bar(2, 2)
        real(dp) :: physical_value(2)
        integer :: basis_dof, dof_count

        jacobian_bar = 0.0_dp
        reference_values_bar = 0.0_dp
        reference_divergences_bar = 0.0_dp
        call validate_products( &
            reference_values, reference_divergences, physical_values_bar, &
            physical_divergences_bar, reference_values_bar, &
            reference_divergences_bar, status)
        if (status /= 0) return
        determinant = det2(jacobian)
        if (.not. valid_determinant(jacobian, determinant)) then
            status = 1
            return
        end if
        determinant_bar = 0.0_dp
        dof_count = size(reference_values, 2)
        do basis_dof = 1, dof_count
            physical_value = &
                matmul(jacobian, reference_values(:, basis_dof))/determinant
            jacobian_bar = jacobian_bar + spread( &
                physical_values_bar(:, basis_dof), 2, 2)*spread( &
                reference_values(:, basis_dof), 1, 2)/determinant
            reference_values_bar(:, basis_dof) = matmul( &
                transpose(jacobian), physical_values_bar(:, basis_dof))/ &
                determinant
            determinant_bar = determinant_bar - dot_product( &
                physical_values_bar(:, basis_dof), physical_value)/determinant
        end do
        reference_divergences_bar = physical_divergences_bar/determinant
        determinant_bar = determinant_bar - dot_product( &
            physical_divergences_bar, reference_divergences)/determinant**2
        call det2_vjp( &
            jacobian, determinant_bar, determinant_jacobian_bar)
        jacobian_bar = jacobian_bar + determinant_jacobian_bar
        status = 0
    end subroutine map_triangle_rt_contravariant_vjp

    pure subroutine validate_products( &
            reference_values, reference_scalars, reference_values_product, &
            reference_scalars_product, physical_values_product, &
            physical_scalars_product, status)
        real(dp), intent(in) :: reference_values(:, :), reference_scalars(:)
        real(dp), intent(in) :: reference_values_product(:, :)
        real(dp), intent(in) :: reference_scalars_product(:)
        real(dp), intent(in) :: physical_values_product(:, :)
        real(dp), intent(in) :: physical_scalars_product(:)
        integer, intent(out) :: status

        integer :: dof_count

        status = 1
        dof_count = size(reference_values, 2)
        if (size(reference_values, 1) /= 2) return
        if (size(reference_scalars) /= dof_count) return
        if (any(shape(reference_values_product) /= shape(reference_values))) &
            return
        if (size(reference_scalars_product) /= dof_count) return
        if (any(shape(physical_values_product) /= shape(reference_values))) &
            return
        if (size(physical_scalars_product) /= dof_count) return
        status = 0
    end subroutine validate_products

    pure logical function valid_determinant(jacobian, determinant) result(valid)
        real(dp), intent(in) :: jacobian(2, 2), determinant

        valid = determinant > 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**2)
    end function valid_determinant

end module fortfem_triangle_piola_maps
