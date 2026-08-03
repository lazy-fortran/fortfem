program test_cartesian_pml_geometry
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: build_cartesian_pml_element_stretch, &
        build_cartesian_pml_element_stretch_jvp, &
        build_cartesian_pml_element_stretch_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(2, 9) = reshape([ &
        0.2_dp, 0.2_dp, 0.4_dp, 0.2_dp, 0.2_dp, 0.4_dp, &
        1.2_dp, 0.2_dp, 1.4_dp, 0.2_dp, 1.2_dp, 0.4_dp, &
        1.2_dp, 1.2_dp, 1.4_dp, 1.2_dp, 1.2_dp, 1.4_dp], [2, 9])
    integer, parameter :: cells(3, 3) = reshape([ &
        1, 2, 3, 4, 5, 6, 7, 8, 9], [3, 3])
    real(dp), parameter :: physical_min(2) = [0.0_dp, 0.0_dp]
    real(dp), parameter :: physical_max(2) = [1.0_dp, 1.0_dp]
    real(dp), parameter :: outer_min(2) = [-0.5_dp, -0.5_dp]
    real(dp), parameter :: outer_max(2) = [1.5_dp, 1.5_dp]
    complex(dp), allocatable :: stretch(:, :), stretch_bar(:, :)
    complex(dp), allocatable :: stretch_dot(:, :)
    complex(dp), allocatable :: stretch_minus(:, :), stretch_plus(:, :)
    real(dp) :: expected_sigma, forward_pairing, reverse_pairing
    real(dp) :: outer_max_bar(2), outer_max_dot(2)
    real(dp) :: outer_min_bar(2), outer_min_dot(2)
    real(dp) :: physical_max_bar(2), physical_max_dot(2)
    real(dp) :: physical_min_bar(2), physical_min_dot(2)
    real(dp) :: sigma_max_bar, sigma_max_dot
    real(dp) :: vertices_bar(2, 9), vertices_dot(2, 9)
    real(dp) :: wave_number_bar, wave_number_dot
    real(dp), parameter :: difference_step = 1.0e-6_dp
    integer :: i, status
    integer :: status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    call build_cartesian_pml_element_stretch( &
        vertices, cells, physical_min, physical_max, outer_min, outer_max, &
        2.0_dp, 4.0_dp, 2, stretch, status)
    call record_condition(status == 0, "Cartesian PML geometry is accepted")
    call record_condition(maxval(abs(stretch(:, 1) - &
        cmplx(1.0_dp, 0.0_dp, dp))) < 1.0e-14_dp, &
        "physical-domain element retains unit stretch")
    expected_sigma = 4.0_dp*((4.0_dp/15.0_dp)/0.5_dp)**2
    call record_condition(abs(stretch(1, 2) - &
        cmplx(1.0_dp, expected_sigma/2.0_dp, dp)) < 1.0e-14_dp .and. &
        abs(stretch(2, 2) - cmplx(1.0_dp, 0.0_dp, dp)) < 1.0e-14_dp, &
        "face-layer element receives the polynomial normal stretch")
    call record_condition(abs(stretch(1, 3) - &
        cmplx(1.0_dp, expected_sigma/2.0_dp, dp)) < 1.0e-14_dp .and. &
        abs(stretch(2, 3) - &
        cmplx(1.0_dp, expected_sigma/2.0_dp, dp)) < 1.0e-14_dp, &
        "corner-layer element receives independent coordinate stretches")

    vertices_dot = reshape([ &
        (0.01_dp*sin(real(i, dp)), i=1, 18)], [2, 9])
    physical_min_dot = [0.02_dp, -0.01_dp]
    physical_max_dot = [-0.03_dp, 0.01_dp]
    outer_min_dot = [-0.01_dp, 0.02_dp]
    outer_max_dot = [0.02_dp, -0.02_dp]
    wave_number_dot = 0.04_dp
    sigma_max_dot = -0.05_dp
    call build_cartesian_pml_element_stretch_jvp( &
        vertices, cells, physical_min, physical_max, outer_min, outer_max, &
        2.0_dp, 4.0_dp, 2, vertices_dot, physical_min_dot, &
        physical_max_dot, outer_min_dot, outer_max_dot, wave_number_dot, &
        sigma_max_dot, stretch_dot, status)
    call build_cartesian_pml_element_stretch( &
        vertices - difference_step*vertices_dot, cells, &
        physical_min - difference_step*physical_min_dot, &
        physical_max - difference_step*physical_max_dot, &
        outer_min - difference_step*outer_min_dot, &
        outer_max - difference_step*outer_max_dot, &
        2.0_dp - difference_step*wave_number_dot, &
        4.0_dp - difference_step*sigma_max_dot, 2, stretch_minus, status_minus)
    call build_cartesian_pml_element_stretch( &
        vertices + difference_step*vertices_dot, cells, &
        physical_min + difference_step*physical_min_dot, &
        physical_max + difference_step*physical_max_dot, &
        outer_min + difference_step*outer_min_dot, &
        outer_max + difference_step*outer_max_dot, &
        2.0_dp + difference_step*wave_number_dot, &
        4.0_dp + difference_step*sigma_max_dot, 2, stretch_plus, status_plus)
    call record_condition(status == 0 .and. status_minus == 0 .and. &
        status_plus == 0, "Cartesian PML geometry JVP accepts valid inputs")
    call record_condition(maxval(abs(stretch_dot - &
        (stretch_plus - stretch_minus)/(2.0_dp*difference_step))) < 2.0e-9_dp, &
        "Cartesian PML geometry JVP matches complete central difference")

    allocate (stretch_bar(2, 3))
    stretch_bar = reshape([ &
        cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.1_dp, 0.5_dp, dp), cmplx(-0.2_dp, 0.3_dp, dp), &
        cmplx(0.4_dp, -0.2_dp, dp), cmplx(0.3_dp, 0.1_dp, dp)], [2, 3])
    call build_cartesian_pml_element_stretch_vjp( &
        vertices, cells, physical_min, physical_max, outer_min, outer_max, &
        2.0_dp, 4.0_dp, 2, stretch_bar, vertices_bar, physical_min_bar, &
        physical_max_bar, outer_min_bar, outer_max_bar, wave_number_bar, &
        sigma_max_bar, status)
    forward_pairing = real(sum(conjg(stretch_bar)*stretch_dot), dp)
    reverse_pairing = sum(vertices_bar*vertices_dot) + &
        dot_product(physical_min_bar, physical_min_dot) + &
        dot_product(physical_max_bar, physical_max_dot) + &
        dot_product(outer_min_bar, outer_min_dot) + &
        dot_product(outer_max_bar, outer_max_dot) + &
        wave_number_bar*wave_number_dot + sigma_max_bar*sigma_max_dot
    call record_condition(abs(forward_pairing - reverse_pairing) < 2.0e-12_dp, &
        "Cartesian PML geometry JVP and VJP satisfy the dot identity")
    call check_summary("Cartesian PML geometry")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_cartesian_pml_geometry
