module fortfem_triangle_piola_maps
    use fortfem_kinds, only: dp
    implicit none

    private

    public :: map_triangle_nedelec_covariant
    public :: map_triangle_rt_contravariant

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

end module fortfem_triangle_piola_maps
