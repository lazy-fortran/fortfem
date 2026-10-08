module fortfem_triangle_piola_maps
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det2
    use fortfem_generated_triangle_covariant_vjp, only: &
        generated_triangle_covariant_vjp
    use fortfem_generated_triangle_contravariant_vjp, only: &
        generated_triangle_contravariant_vjp
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
        real(dp) :: determinant
        integer :: basis_dof

        physical_values = 0.0_dp
        physical_curls = 0.0_dp
        call validate_products(reference_values, reference_curls, &
            reference_values, reference_curls, physical_values, &
            physical_curls, status)
        if (status /= 0) return
        status = 1
        determinant = det2(jacobian)
        if (.not. valid_determinant(jacobian, determinant)) return
        do basis_dof = 1, size(reference_values, 2)
            include "../generated/fortfem_triangle_covariant_primal.inc"
        end do
        status = 0
    end subroutine map_triangle_nedelec_covariant

    pure subroutine map_triangle_nedelec_covariant_jvp( &
            jacobian, reference_values, reference_curls, jacobian_dot, &
            reference_values_dot, reference_curls_dot, physical_values_dot, &
            physical_curls_dot, status)
        real(dp), intent(in) :: jacobian(2, 2), jacobian_dot(2, 2)
        real(dp), intent(in) :: reference_values(:, :), reference_curls(:)
        real(dp), intent(in) :: reference_values_dot(:, :), reference_curls_dot(:)
        real(dp), intent(out) :: physical_values_dot(:, :), physical_curls_dot(:)
        integer, intent(out) :: status
        real(dp) :: determinant
        integer :: basis_dof

        physical_values_dot = 0.0_dp
        physical_curls_dot = 0.0_dp
        call validate_products(reference_values, reference_curls, &
            reference_values_dot, reference_curls_dot, physical_values_dot, &
            physical_curls_dot, status)
        if (status /= 0) return
        status = 1
        determinant = det2(jacobian)
        if (.not. valid_determinant(jacobian, determinant)) return
        do basis_dof = 1, size(reference_values, 2)
            include "../generated/fortfem_triangle_covariant_jvp.inc"
        end do
        status = 0
    end subroutine map_triangle_nedelec_covariant_jvp

    pure subroutine map_triangle_nedelec_covariant_vjp( &
            jacobian, reference_values, reference_curls, physical_values_bar, &
            physical_curls_bar, jacobian_bar, reference_values_bar, &
            reference_curls_bar, status)
        real(dp), intent(in) :: jacobian(2, 2)
        real(dp), intent(in) :: reference_values(:, :), reference_curls(:)
        real(dp), intent(in) :: physical_values_bar(:, :), physical_curls_bar(:)
        real(dp), intent(out) :: jacobian_bar(2, 2)
        real(dp), intent(out) :: reference_values_bar(:, :), reference_curls_bar(:)
        integer, intent(out) :: status
        real(dp) :: determinant, product(7)
        integer :: basis_dof

        jacobian_bar = 0.0_dp
        reference_values_bar = 0.0_dp
        reference_curls_bar = 0.0_dp
        call validate_products(reference_values, reference_curls, &
            physical_values_bar, physical_curls_bar, reference_values_bar, &
            reference_curls_bar, status)
        if (status /= 0) return
        status = 1
        determinant = det2(jacobian)
        if (.not. valid_determinant(jacobian, determinant)) return
        do basis_dof = 1, size(reference_values, 2)
            call generated_triangle_covariant_vjp( &
                jacobian(1, 1), jacobian(2, 1), jacobian(1, 2), jacobian(2, 2), &
                reference_values(1, basis_dof), reference_values(2, basis_dof), &
                reference_curls(basis_dof), &
                physical_values_bar(1, basis_dof), &
                physical_values_bar(2, basis_dof), physical_curls_bar(basis_dof), &
                product)
            jacobian_bar = jacobian_bar + reshape(product(:4), [2, 2])
            reference_values_bar(:, basis_dof) = product(5:6)
            reference_curls_bar(basis_dof) = product(7)
        end do
        status = 0
    end subroutine map_triangle_nedelec_covariant_vjp

    pure subroutine map_triangle_rt_contravariant( &
            jacobian, reference_values, reference_divergences, physical_values, &
            physical_divergences, status)
        real(dp), intent(in) :: jacobian(2, 2)
        real(dp), intent(in) :: reference_values(:, :), reference_divergences(:)
        real(dp), intent(out) :: physical_values(:, :), physical_divergences(:)
        integer, intent(out) :: status
        real(dp) :: determinant
        integer :: basis_dof

        physical_values = 0.0_dp
        physical_divergences = 0.0_dp
        call validate_products(reference_values, reference_divergences, &
            reference_values, reference_divergences, physical_values, &
            physical_divergences, status)
        if (status /= 0) return
        status = 1
        determinant = det2(jacobian)
        if (.not. valid_determinant(jacobian, determinant)) return
        do basis_dof = 1, size(reference_values, 2)
            include "../generated/fortfem_triangle_contravariant_primal.inc"
        end do
        status = 0
    end subroutine map_triangle_rt_contravariant

    pure subroutine map_triangle_rt_contravariant_jvp( &
            jacobian, reference_values, reference_divergences, jacobian_dot, &
            reference_values_dot, reference_divergences_dot, physical_values_dot, &
            physical_divergences_dot, status)
        real(dp), intent(in) :: jacobian(2, 2), jacobian_dot(2, 2)
        real(dp), intent(in) :: reference_values(:, :), reference_divergences(:)
        real(dp), intent(in) :: reference_values_dot(:, :), reference_divergences_dot(:)
        real(dp), intent(out) :: physical_values_dot(:, :), physical_divergences_dot(:)
        integer, intent(out) :: status
        real(dp) :: determinant
        integer :: basis_dof

        physical_values_dot = 0.0_dp
        physical_divergences_dot = 0.0_dp
        call validate_products(reference_values, reference_divergences, &
            reference_values_dot, reference_divergences_dot, physical_values_dot, &
            physical_divergences_dot, status)
        if (status /= 0) return
        status = 1
        determinant = det2(jacobian)
        if (.not. valid_determinant(jacobian, determinant)) return
        do basis_dof = 1, size(reference_values, 2)
            include "../generated/fortfem_triangle_contravariant_jvp.inc"
        end do
        status = 0
    end subroutine map_triangle_rt_contravariant_jvp

    pure subroutine map_triangle_rt_contravariant_vjp( &
            jacobian, reference_values, reference_divergences, physical_values_bar, &
            physical_divergences_bar, jacobian_bar, reference_values_bar, &
            reference_divergences_bar, status)
        real(dp), intent(in) :: jacobian(2, 2)
        real(dp), intent(in) :: reference_values(:, :), reference_divergences(:)
        real(dp), intent(in) :: physical_values_bar(:, :), physical_divergences_bar(:)
        real(dp), intent(out) :: jacobian_bar(2, 2)
        real(dp), intent(out) :: reference_values_bar(:, :), reference_divergences_bar(:)
        integer, intent(out) :: status
        real(dp) :: determinant, product(7)
        integer :: basis_dof

        jacobian_bar = 0.0_dp
        reference_values_bar = 0.0_dp
        reference_divergences_bar = 0.0_dp
        call validate_products(reference_values, reference_divergences, &
            physical_values_bar, physical_divergences_bar, reference_values_bar, &
            reference_divergences_bar, status)
        if (status /= 0) return
        status = 1
        determinant = det2(jacobian)
        if (.not. valid_determinant(jacobian, determinant)) return
        do basis_dof = 1, size(reference_values, 2)
            call generated_triangle_contravariant_vjp( &
                jacobian(1, 1), jacobian(2, 1), jacobian(1, 2), jacobian(2, 2), &
                reference_values(1, basis_dof), reference_values(2, basis_dof), &
                reference_divergences(basis_dof), &
                physical_values_bar(1, basis_dof), &
                physical_values_bar(2, basis_dof), physical_divergences_bar(basis_dof), &
                product)
            jacobian_bar = jacobian_bar + reshape(product(:4), [2, 2])
            reference_values_bar(:, basis_dof) = product(5:6)
            reference_divergences_bar(basis_dof) = product(7)
        end do
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
