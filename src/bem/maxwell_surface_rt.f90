module fortfem_maxwell_surface_rt
    !! Arbitrary-order divergence-conforming currents on parametric panels.
    use fortfem_kinds, only: dp
    use fortfem_triangle_rt_arbitrary_order, only: &
        evaluate_triangle_raviart_thomas, triangle_rt_basis_t, &
        triangle_rt_dof_count
    implicit none
    private

    public :: evaluate_maxwell_surface_rt_basis

contains

    subroutine evaluate_maxwell_surface_rt_basis( &
            basis, xi, eta, tangent_xi, tangent_eta, surface_jacobian, &
            values, surface_divergences, status)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta, tangent_xi(3), tangent_eta(3)
        real(dp), intent(in) :: surface_jacobian
        real(dp), intent(out) :: values(:, :), surface_divergences(:)
        integer, intent(out) :: status

        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :)
        integer :: dof, dof_count

        status = 1
        values = 0.0_dp
        surface_divergences = 0.0_dp
        dof_count = triangle_rt_dof_count(basis)
        if (dof_count < 1 .or. surface_jacobian <= 0.0_dp) return
        if (size(values, 1) /= 3 .or. size(values, 2) /= dof_count) return
        if (size(surface_divergences) /= dof_count) return
        allocate( &
            reference_values(2, dof_count), &
            reference_divergences(dof_count))
        call evaluate_triangle_raviart_thomas( &
            basis, xi, eta, reference_values, reference_divergences, status)
        if (status /= 0) return
        do dof = 1, dof_count
            values(:, dof) = ( &
                tangent_xi*reference_values(1, dof) + &
                tangent_eta*reference_values(2, dof))/surface_jacobian
        end do
        surface_divergences = reference_divergences/surface_jacobian
        status = 0
    end subroutine evaluate_maxwell_surface_rt_basis

end module fortfem_maxwell_surface_rt
