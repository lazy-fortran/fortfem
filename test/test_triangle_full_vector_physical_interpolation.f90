program test_triangle_full_vector_physical_interpolation
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_triangle_bdm_interpolant, &
        evaluate_triangle_nedelec_second_kind_interpolant, &
        initialize_triangle_bdm, initialize_triangle_nedelec_second_kind, &
        interpolate_triangle_bdm, interpolate_triangle_nedelec_second_kind, &
        triangle_bdm_basis_t, triangle_nedelec_second_kind_t
    use fortfem_kinds, only: dp
    implicit none

    type(triangle_bdm_basis_t) :: bdm_basis
    type(triangle_nedelec_second_kind_t) :: nedelec_basis
    real(dp), allocatable :: bdm_dofs(:), nedelec_dofs(:)
    real(dp) :: bdm_divergence, bdm_value(2), exact_divergence
    real(dp) :: exact_value(2), jacobian(2, 2), nedelec_curl
    real(dp) :: nedelec_value(2), physical_point(2), vertices(2, 3)
    real(dp) :: xi, eta
    integer :: active_degree, degree, status
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.2_dp, -0.1_dp]
    vertices(:, 2) = [1.7_dp, 0.3_dp]
    vertices(:, 3) = [-0.2_dp, 1.4_dp]
    jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
    jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
    xi = 0.27_dp
    eta = 0.19_dp
    physical_point = vertices(:, 1) + matmul(jacobian, [xi, eta])

    do degree = 1, 4
        active_degree = degree
        call manufactured_field( &
            physical_point(1), physical_point(2), exact_value)
        exact_divergence = real(degree, dp) * ( &
            physical_point(1)**(degree - 1) + &
            physical_point(2)**(degree - 1))

        call initialize_triangle_nedelec_second_kind( &
            degree, nedelec_basis, status)
        call interpolate_triangle_nedelec_second_kind( &
            vertices, degree, 2 * degree + 4, manufactured_field, &
            nedelec_dofs, status)
        call evaluate_triangle_nedelec_second_kind_interpolant( &
            vertices, nedelec_basis, nedelec_dofs, xi, eta, nedelec_value, &
            nedelec_curl, status)
        call record_condition(maxval(abs(nedelec_value - exact_value)) < &
            2.0e-9_dp, &
            "Physical second-kind interpolation reproduces [P_k]^2")
        call record_condition(abs(nedelec_curl + 2.0_dp) < 3.0e-9_dp, &
            "Physical second-kind interpolation reproduces exact curl")

        call initialize_triangle_bdm(degree, bdm_basis, status)
        call interpolate_triangle_bdm( &
            vertices, degree, 2 * degree + 4, manufactured_field, &
            bdm_dofs, status)
        call evaluate_triangle_bdm_interpolant( &
            vertices, bdm_basis, bdm_dofs, xi, eta, bdm_value, &
            bdm_divergence, status)
        call record_condition(maxval(abs(bdm_value - exact_value)) < &
            2.0e-9_dp, "Physical BDM interpolation reproduces [P_k]^2")
        call record_condition(abs(bdm_divergence - exact_divergence) < &
            3.0e-9_dp, &
            "Physical BDM interpolation reproduces exact divergence")

        deallocate(bdm_dofs, nedelec_dofs)
    end do

    call interpolate_triangle_bdm( &
        vertices, 0, 2, manufactured_field, bdm_dofs, status)
    call record_condition(status /= 0, &
        "Physical BDM interpolation rejects degree zero")

    call check_summary("Full polynomial vector physical interpolation")
    if (.not. all_passed) error stop 1

contains

    subroutine manufactured_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [x**active_degree + y, y**active_degree - x]
    end subroutine manufactured_field

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_full_vector_physical_interpolation
