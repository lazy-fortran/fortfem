program test_tetra_rt_piola_map
    use check, only: check_condition, check_summary
    use fortfem_feec, only: map_tetra_rt_contravariant
    use fortfem_kinds, only: dp
    implicit none

    real(dp) :: jacobian(3, 3), physical_divergences(3)
    real(dp) :: physical_normal(3), physical_tangent_one(3)
    real(dp) :: physical_tangent_two(3), physical_values(3, 3)
    real(dp) :: reference_divergences(3), reference_normal(3)
    real(dp) :: reference_tangent_one(3), reference_tangent_two(3)
    real(dp) :: reference_values(3, 3)
    integer :: basis, status
    logical :: all_passed

    all_passed = .true.
    jacobian(:, 1) = [2.0_dp, -0.4_dp, 0.3_dp]
    jacobian(:, 2) = [0.2_dp, 1.5_dp, -0.1_dp]
    jacobian(:, 3) = [-0.3_dp, 0.4_dp, 1.2_dp]
    reference_values(:, 1) = [1.0_dp, -2.0_dp, 0.5_dp]
    reference_values(:, 2) = [0.25_dp, 0.75_dp, -1.0_dp]
    reference_values(:, 3) = [-1.2_dp, 0.4_dp, 2.0_dp]
    reference_divergences = [-1.0_dp, 2.5_dp, 0.125_dp]
    reference_tangent_one = [0.6_dp, -0.2_dp, 0.3_dp]
    reference_tangent_two = [-0.1_dp, 0.4_dp, 0.7_dp]
    reference_normal = cross_product( &
        reference_tangent_one, reference_tangent_two)
    physical_tangent_one = matmul(jacobian, reference_tangent_one)
    physical_tangent_two = matmul(jacobian, reference_tangent_two)
    physical_normal = cross_product( &
        physical_tangent_one, physical_tangent_two)

    call map_tetra_rt_contravariant( &
        jacobian, reference_values, reference_divergences, physical_values, &
        physical_divergences, status)
    call record_condition(status == 0, &
        "Tetrahedral contravariant Piola map accepts a skew tetrahedron")
    do basis = 1, 3
        call record_condition(abs( &
            dot_product(physical_values(:, basis), physical_normal) - &
            dot_product(reference_values(:, basis), reference_normal)) < &
            3.0e-14_dp, &
            "Tetrahedral contravariant Piola map preserves normal flux")
    end do
    call record_condition(maxval(abs(physical_divergences - &
        reference_divergences/determinant(jacobian))) < 3.0e-14_dp, &
        "Tetrahedral contravariant Piola map scales divergence exactly")

    jacobian(:, 3) = jacobian(:, 1) + jacobian(:, 2)
    call map_tetra_rt_contravariant( &
        jacobian, reference_values, reference_divergences, physical_values, &
        physical_divergences, status)
    call record_condition(status /= 0, &
        "Tetrahedral contravariant Piola map rejects a singular transform")

    call check_summary("Tetrahedral RT contravariant Piola map")
    if (.not. all_passed) error stop 1

contains

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    pure function determinant(matrix) result(value)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp) :: value

        value = matrix(1, 1)*( &
            matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) - &
            matrix(1, 2)*( &
            matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
            matrix(1, 3)*( &
            matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
    end function determinant

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_rt_piola_map
