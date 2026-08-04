subroutine fortfem_h1_geometry(jacobian_11, jacobian_21, jacobian_12, &
                               jacobian_22, gradient_row_1, gradient_row_2, &
                               gradient_column_1, gradient_column_2, &
                               basis_row, basis_column, &
                               stiffness_coefficient, mass_coefficient, &
                               quadrature_weight, contribution)
    !! One H1 element-matrix entry on a mapped cell.
    !!
    !! The primal for the fortad path. It states the same computation the
    !! Enzyme fixture states as an objective over packed `geometry` and `data`
    !! arrays, and the same one the fortsym generator states symbolically, but
    !! with the arguments named - which is what lets the three be compared
    !! argument for argument rather than only in aggregate.
    !!
    !! Only the four Jacobian entries are active. The reference gradients, the
    !! basis values, the coefficients and the quadrature weight select the
    !! entry being assembled; they are data, not geometry.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: jacobian_11, jacobian_21, jacobian_12, jacobian_22
    real(dp), intent(in) :: gradient_row_1, gradient_row_2
    real(dp), intent(in) :: gradient_column_1, gradient_column_2
    real(dp), intent(in) :: basis_row, basis_column
    real(dp), intent(in) :: stiffness_coefficient, mass_coefficient
    real(dp), intent(in) :: quadrature_weight
    real(dp), intent(out) :: contribution
    real(dp) :: determinant, row_1, row_2, column_1, column_2, stiffness, mass

    determinant = jacobian_11*jacobian_22 - jacobian_12*jacobian_21
    row_1 = (jacobian_22*gradient_row_1 - jacobian_21*gradient_row_2)/determinant
    row_2 = (-jacobian_12*gradient_row_1 + jacobian_11*gradient_row_2)/determinant
    column_1 = (jacobian_22*gradient_column_1 - &
                jacobian_21*gradient_column_2)/determinant
    column_2 = (-jacobian_12*gradient_column_1 + &
                jacobian_11*gradient_column_2)/determinant
    stiffness = stiffness_coefficient*(row_1*column_1 + row_2*column_2)
    mass = mass_coefficient*basis_row*basis_column
    contribution = quadrature_weight*determinant*(stiffness + mass)
end subroutine fortfem_h1_geometry
