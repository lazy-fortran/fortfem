module fortfem_assembly_rt_arbitrary_order_2d
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_piola_maps, only: map_triangle_rt_contravariant
    use fortfem_triangle_rt_arbitrary_order, only: &
        evaluate_triangle_raviart_thomas, initialize_triangle_raviart_thomas, &
        triangle_rt_basis_t, triangle_rt_dof_count
    implicit none

    private

    public :: assemble_triangle_rt_div_mass_element

contains

    subroutine assemble_triangle_rt_div_mass_element( &
            vertices, degree, quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(2, 3)
        integer, intent(in) :: degree, quadrature_degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        type(triangle_rt_basis_t) :: basis
        real(dp), allocatable :: eta(:), physical_divergences(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :), weights(:), xi(:)
        real(dp) :: determinant, jacobian(2, 2), physical_weight
        integer :: basis_dof_count, column, point, row

        status = 1
        if (degree < 0 .or. quadrature_degree < 0) return

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        determinant = jacobian(1, 1) * jacobian(2, 2) - &
            jacobian(1, 2) * jacobian(2, 1)
        if (determinant <= 64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, maxval(abs(jacobian))**2)) return

        call initialize_triangle_raviart_thomas(degree, basis, status)
        if (status /= 0) return
        basis_dof_count = triangle_rt_dof_count(basis)
        allocate(matrix(basis_dof_count, basis_dof_count))
        allocate(reference_values(2, basis_dof_count))
        allocate(reference_divergences(basis_dof_count))
        allocate(physical_values(2, basis_dof_count))
        allocate(physical_divergences(basis_dof_count))
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return

        matrix = 0.0_dp
        do point = 1, size(weights)
            call evaluate_triangle_raviart_thomas( &
                basis, xi(point), eta(point), reference_values, &
                reference_divergences, status)
            if (status /= 0) return
            call map_triangle_rt_contravariant( &
                jacobian, reference_values, reference_divergences, &
                physical_values, physical_divergences, status)
            if (status /= 0) return
            physical_weight = determinant * weights(point)
            do column = 1, basis_dof_count
                do row = 1, basis_dof_count
                    matrix(row, column) = matrix(row, column) + &
                        physical_weight * ( &
                        physical_divergences(row) * &
                        physical_divergences(column) + &
                        dot_product( &
                        physical_values(:, row), physical_values(:, column)))
                end do
            end do
        end do
        status = 0
    end subroutine assemble_triangle_rt_div_mass_element

end module fortfem_assembly_rt_arbitrary_order_2d
