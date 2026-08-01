program test_paper_magnetic_box_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_magnetic_box_3d
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: reference = 0.07367135328_dp
    type(fortsparse_status_t) :: status
    real(dp) :: coarse_field(3), exact_field(3), fine_field(3)
    real(dp) :: p_fields(3, 2:4)
    real(dp) :: field_point(3)
    real(dp) :: coarse_value, fine_value, p_values(2:4)
    integer :: coarse_dofs, degree, fine_dofs, p_dofs
    logical :: all_passed

    all_passed = .true.
    field_point = [0.23_dp, 0.41_dp, 1.47_dp]
    exact_field = membrane_magnetic_field(field_point)
    call solve_magnetic_box_3d( &
        6, 4, coarse_value, coarse_dofs, status, field_point, coarse_field)
    call record_condition(status%code == 0, &
        "Coarse magnetic-paper box solve succeeds")
    call solve_magnetic_box_3d( &
        12, 8, fine_value, fine_dofs, status, field_point, fine_field)
    call record_condition(status%code == 0, &
        "Paper-resolution magnetic box solve succeeds")
    call record_condition(abs(fine_value - reference) / reference < 0.05_dp, &
        "Magnetic box centre value matches the membrane series")
    call record_condition(abs(fine_value - reference) < &
        abs(coarse_value - reference), &
        "Magnetic box centre value improves under refinement")
    call record_condition(fine_dofs > coarse_dofs, &
        "Magnetic box refinement increases edge degrees of freedom")
    call record_condition(norm2(fine_field - exact_field) / &
        norm2(exact_field) < 0.15_dp, &
        "Magnetic box curl reaches the differentiated membrane series")
    call record_condition(norm2(fine_field - exact_field) < &
        norm2(coarse_field - exact_field), &
        "Magnetic box curl improves under refinement")
    do degree = 2, 4
        call solve_magnetic_box_3d( &
            2, 1, p_values(degree), p_dofs, status, field_point, &
            p_fields(:, degree), order=degree)
        call record_condition(status%code == 0, &
            "Higher-order magnetic-paper box solve succeeds")
    end do
    write (*, '(a,3(es12.4,1x))') &
        "Magnetic box p=2:4 centre values ", p_values
    call record_condition(all(abs(p_values(3:4) - reference) < &
        abs(p_values(2:3) - reference)), &
        "Magnetic box centre value improves with polynomial order")
    call record_condition(abs(p_values(4) - reference)/reference < 0.05_dp, &
        "Order-four magnetic box reaches the membrane series")
    call record_condition( &
        norm2(p_fields(:, 4) - exact_field) < &
        norm2(p_fields(:, 2) - exact_field), &
        "Magnetic box curl improves with polynomial order")
    call record_condition( &
        norm2(p_fields(:, 4) - exact_field)/norm2(exact_field) < 0.05_dp, &
        "Order-four magnetic box curl reaches the differentiated series")

    call check_summary("Paper magnetic three-dimensional box")
    if (.not. all_passed) error stop 1

contains

    pure function membrane_magnetic_field(point) result(field)
        real(dp), intent(in) :: point(3)
        real(dp) :: field(3)

        real(dp), parameter :: pi = acos(-1.0_dp)
        integer, parameter :: maximum_mode = 399
        real(dp) :: denominator
        integer :: m, n

        field = 0.0_dp
        do n = 1, maximum_mode, 2
            do m = 1, maximum_mode, 2
                denominator = real(m * m + n * n, dp)
                field(1) = field(1) + &
                    sin(real(m, dp) * pi * point(1)) * &
                    cos(real(n, dp) * pi * point(2)) / &
                    (real(m, dp) * denominator)
                field(2) = field(2) - &
                    cos(real(m, dp) * pi * point(1)) * &
                    sin(real(n, dp) * pi * point(2)) / &
                    (real(n, dp) * denominator)
            end do
        end do
        field = 16.0_dp * field / pi**3
    end function membrane_magnetic_field

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_paper_magnetic_box_3d
