module fortfem_tetra_piola_maps
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3, inv3
    implicit none

    private

    public :: map_tetra_nedelec_covariant

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

end module fortfem_tetra_piola_maps
