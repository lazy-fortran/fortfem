program test_arbitrary_order_piola_maps
    use check, only: check_condition, check_summary
    use fortfem_api, only: map_triangle_nedelec_covariant, &
        map_triangle_rt_contravariant
    use fortfem_kinds, only: dp
    implicit none

    real(dp) :: jacobian(2, 2), physical_curls(3), physical_divergences(3)
    real(dp) :: physical_edge(2), physical_normal(2)
    real(dp) :: physical_values(2, 3), reference_curls(3)
    real(dp) :: reference_divergences(3), reference_edge(2)
    real(dp) :: reference_normal(2), reference_values(2, 3)
    integer :: basis_dof, status
    logical :: all_passed

    all_passed = .true.
    jacobian(1, :) = [2.0_dp, 0.3_dp]
    jacobian(2, :) = [-0.4_dp, 1.5_dp]
    reference_values(:, 1) = [1.0_dp, -2.0_dp]
    reference_values(:, 2) = [0.25_dp, 0.75_dp]
    reference_values(:, 3) = [-1.2_dp, 0.4_dp]
    reference_curls = [2.0_dp, -0.5_dp, 3.25_dp]
    reference_divergences = [-1.0_dp, 2.5_dp, 0.125_dp]
    reference_edge = [0.6_dp, -0.2_dp]
    reference_normal = [reference_edge(2), -reference_edge(1)]
    physical_edge = matmul(jacobian, reference_edge)
    physical_normal = [physical_edge(2), -physical_edge(1)]

    call map_triangle_nedelec_covariant( &
        jacobian, reference_values, reference_curls, physical_values, &
        physical_curls, status)
    do basis_dof = 1, 3
        call record_condition(abs( &
            dot_product(physical_values(:, basis_dof), physical_edge) - &
            dot_product(reference_values(:, basis_dof), reference_edge)) < &
            2.0e-14_dp, &
            "Covariant Piola map preserves an arbitrary tangential moment")
    end do
    call record_condition(maxval(abs(physical_curls - &
        reference_curls / determinant(jacobian))) < 2.0e-14_dp, &
        "Covariant Piola map applies exact curl scaling")

    call map_triangle_rt_contravariant( &
        jacobian, reference_values, reference_divergences, physical_values, &
        physical_divergences, status)
    do basis_dof = 1, 3
        call record_condition(abs( &
            dot_product(physical_values(:, basis_dof), physical_normal) - &
            dot_product(reference_values(:, basis_dof), reference_normal)) < &
            2.0e-14_dp, &
            "Contravariant Piola map preserves an arbitrary normal moment")
    end do
    call record_condition(maxval(abs(physical_divergences - &
        reference_divergences / determinant(jacobian))) < 2.0e-14_dp, &
        "Contravariant Piola map applies exact divergence scaling")

    jacobian(2, :) = 2.0_dp * jacobian(1, :)
    call map_triangle_nedelec_covariant( &
        jacobian, reference_values, reference_curls, physical_values, &
        physical_curls, status)
    call record_condition(status /= 0, &
        "Piola map rejects a singular affine transformation")

    call check_summary("Arbitrary-order affine Piola maps")
    if (.not. all_passed) error stop 1

contains

    pure function determinant(matrix) result(value)
        real(dp), intent(in) :: matrix(2, 2)
        real(dp) :: value

        value = matrix(1, 1) * matrix(2, 2) - &
            matrix(1, 2) * matrix(2, 1)
    end function determinant

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_arbitrary_order_piola_maps
