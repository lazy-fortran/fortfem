program test_tetra_nedelec_affine_map
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_tetra_nedelec_first_order, &
        map_tetra_nedelec_covariant
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: inv3
    implicit none

    real(dp) :: inverse_jacobian(3, 3), jacobian(3, 3)
    real(dp) :: physical_curls(3, 6), physical_values(3, 6)
    real(dp) :: reference_curls(3, 6), reference_values(3, 6)
    real(dp) :: reference_vertices(3, 4), midpoint(3), tangent(3)
    integer :: edge, edge_vertices(2, 6), info, other_edge, status
    logical :: all_passed

    all_passed = .true.
    jacobian(:, 1) = [1.7_dp, 0.2_dp, -0.1_dp]
    jacobian(:, 2) = [0.3_dp, 1.4_dp, 0.2_dp]
    jacobian(:, 3) = [-0.2_dp, 0.4_dp, 1.2_dp]
    call inv3(jacobian, inverse_jacobian, info)
    reference_vertices = 0.0_dp
    reference_vertices(1, 2) = 1.0_dp
    reference_vertices(2, 3) = 1.0_dp
    reference_vertices(3, 4) = 1.0_dp
    edge_vertices(:, 1) = [1, 2]
    edge_vertices(:, 2) = [1, 3]
    edge_vertices(:, 3) = [1, 4]
    edge_vertices(:, 4) = [2, 3]
    edge_vertices(:, 5) = [2, 4]
    edge_vertices(:, 6) = [3, 4]

    do edge = 1, 6
        midpoint = 0.5_dp * ( &
            reference_vertices(:, edge_vertices(1, edge)) + &
            reference_vertices(:, edge_vertices(2, edge)))
        tangent = matmul(jacobian, &
            reference_vertices(:, edge_vertices(2, edge)) - &
            reference_vertices(:, edge_vertices(1, edge)))
        call evaluate_tetra_nedelec_first_order( &
            midpoint, reference_values, reference_curls, status)
        call map_tetra_nedelec_covariant( &
            jacobian, reference_values, reference_curls, physical_values, &
            physical_curls, status)
        do other_edge = 1, 6
            call record_condition(abs( &
                dot_product(physical_values(:, other_edge), tangent) - &
                merge(1.0_dp, 0.0_dp, other_edge == edge)) < 3.0e-14_dp, &
                "Covariant tetrahedral map preserves oriented edge moments")
        end do
    end do
    call check_physical_curls(jacobian, inverse_jacobian)

    call check_summary("Affine tetrahedral Nedelec map")
    if (.not. all_passed) error stop 1

contains

    subroutine check_physical_curls(jacobian, inverse_jacobian)
        real(dp), intent(in) :: jacobian(3, 3), inverse_jacobian(3, 3)

        real(dp), parameter :: step = 1.0e-6_dp
        real(dp) :: center(3), center_physical(3), numerical_curl(3)
        real(dp) :: shifted_physical(3), values_minus(3, 6)
        real(dp) :: values_plus(3, 6)
        integer :: basis, direction, local_status

        center = [0.2_dp, 0.2_dp, 0.2_dp]
        center_physical = matmul(jacobian, center)
        call mapped_values( &
            center_physical, jacobian, inverse_jacobian, physical_values, &
            physical_curls, local_status)
        do basis = 1, 6
            numerical_curl = 0.0_dp
            do direction = 1, 3
                shifted_physical = center_physical
                shifted_physical(direction) = &
                    shifted_physical(direction) + step
                call mapped_values( &
                    shifted_physical, jacobian, inverse_jacobian, &
                    values_plus, reference_curls, local_status)
                shifted_physical(direction) = &
                    center_physical(direction) - step
                call mapped_values( &
                    shifted_physical, jacobian, inverse_jacobian, &
                    values_minus, reference_curls, local_status)
                call accumulate_curl_difference( &
                    direction, basis, values_plus, values_minus, step, &
                    numerical_curl)
            end do
            call record_condition(maxval(abs( &
                numerical_curl - physical_curls(:, basis))) < 3.0e-9_dp, &
                "Mapped curls agree with physical finite differences")
        end do
    end subroutine check_physical_curls

    subroutine mapped_values( &
            physical_point, jacobian, inverse_jacobian, values, curls, status)
        real(dp), intent(in) :: physical_point(3), jacobian(3, 3)
        real(dp), intent(in) :: inverse_jacobian(3, 3)
        real(dp), intent(out) :: values(3, 6), curls(3, 6)
        integer, intent(out) :: status

        real(dp) :: point(3), reference_curls(3, 6)
        real(dp) :: reference_values(3, 6)

        point = matmul(inverse_jacobian, physical_point)
        call evaluate_tetra_nedelec_first_order( &
            point, reference_values, reference_curls, status)
        call map_tetra_nedelec_covariant( &
            jacobian, reference_values, reference_curls, values, curls, status)
    end subroutine mapped_values

    subroutine accumulate_curl_difference( &
            direction, basis, plus, minus, step, curl)
        integer, intent(in) :: direction, basis
        real(dp), intent(in) :: plus(3, 6), minus(3, 6), step
        real(dp), intent(inout) :: curl(3)

        real(dp) :: derivative(3)

        derivative = (plus(:, basis) - minus(:, basis)) / (2.0_dp * step)
        select case (direction)
        case (1)
            curl(2) = curl(2) - derivative(3)
            curl(3) = curl(3) + derivative(2)
        case (2)
            curl(1) = curl(1) + derivative(3)
            curl(3) = curl(3) - derivative(1)
        case (3)
            curl(1) = curl(1) - derivative(2)
            curl(2) = curl(2) + derivative(1)
        end select
    end subroutine accumulate_curl_difference

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_affine_map
