program test_cartesian_pml_chain_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: build_cartesian_pml_element_stretch, &
        build_cartesian_pml_element_stretch_jvp, &
        build_cartesian_pml_element_stretch_vjp, &
        cartesian_scalar_helmholtz_pml_coefficients, &
        cartesian_scalar_helmholtz_pml_coefficients_jvp, &
        cartesian_scalar_helmholtz_pml_coefficients_vjp
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
    real(dp), parameter :: difference_step = 1.0e-6_dp
    complex(dp), allocatable :: stretch(:, :), stretch_bar(:, :)
    complex(dp), allocatable :: stretch_dot(:, :)
    complex(dp) :: coefficient_stretch(3), coefficient_stretch_dot(3)
    complex(dp) :: gradient_bar(3, 3), gradient_dot(3), local_stretch_bar(3)
    complex(dp) :: mass_bar(3), mass_dot
    real(dp) :: forward_product, objective_minus, objective_plus
    real(dp) :: outer_max_bar(2), outer_max_dot(2)
    real(dp) :: outer_min_bar(2), outer_min_dot(2)
    real(dp) :: physical_max_bar(2), physical_max_dot(2)
    real(dp) :: physical_min_bar(2), physical_min_dot(2)
    real(dp) :: reverse_product, sigma_max_bar, sigma_max_dot
    real(dp) :: vertices_bar(2, 9), vertices_dot(2, 9)
    real(dp) :: wave_number_bar, wave_number_dot
    integer :: cell, i, status
    logical :: all_passed

    all_passed = .true.
    vertices_dot = reshape([ &
        (0.01_dp*cos(0.3_dp*real(i, dp)), i=1, 18)], [2, 9])
    physical_min_dot = [0.02_dp, -0.01_dp]
    physical_max_dot = [-0.03_dp, 0.01_dp]
    outer_min_dot = [-0.01_dp, 0.02_dp]
    outer_max_dot = [0.02_dp, -0.02_dp]
    wave_number_dot = 0.04_dp
    sigma_max_dot = -0.05_dp
    gradient_bar = reshape([ &
        (cmplx(0.07_dp*real(i, dp), -0.03_dp*real(i + 1, dp), dp), &
        i=1, 9)], [3, 3])
    mass_bar = [ &
        cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.1_dp, 0.2_dp, dp)]

    call build_cartesian_pml_element_stretch( &
        vertices, cells, physical_min, physical_max, outer_min, outer_max, &
        2.0_dp, 4.0_dp, 2, stretch, status)
    call record_condition(status == 0, "PML chain primal geometry succeeds")
    call build_cartesian_pml_element_stretch_jvp( &
        vertices, cells, physical_min, physical_max, outer_min, outer_max, &
        2.0_dp, 4.0_dp, 2, vertices_dot, physical_min_dot, &
        physical_max_dot, outer_min_dot, outer_max_dot, wave_number_dot, &
        sigma_max_dot, stretch_dot, status)
    call record_condition(status == 0, "PML chain geometry JVP succeeds")

    forward_product = 0.0_dp
    allocate (stretch_bar(2, size(cells, 2)), source=cmplx(0.0_dp, 0.0_dp, dp))
    do cell = 1, size(cells, 2)
        coefficient_stretch = [ &
            stretch(:, cell), cmplx(1.0_dp, 0.0_dp, dp)]
        coefficient_stretch_dot = [ &
            stretch_dot(:, cell), cmplx(0.0_dp, 0.0_dp, dp)]
        call cartesian_scalar_helmholtz_pml_coefficients_jvp( &
            coefficient_stretch, coefficient_stretch_dot, gradient_dot, &
            mass_dot, status)
        if (status /= 0) exit
        forward_product = forward_product + &
            real(sum(conjg(gradient_bar(:, cell))*gradient_dot) + &
            conjg(mass_bar(cell))*mass_dot, dp)
        call cartesian_scalar_helmholtz_pml_coefficients_vjp( &
            coefficient_stretch, gradient_bar(:, cell), mass_bar(cell), &
            local_stretch_bar, status)
        if (status /= 0) exit
        stretch_bar(:, cell) = local_stretch_bar(1:2)
    end do
    call record_condition(status == 0, "PML coefficient products compose")

    call evaluate_objective( &
        vertices - difference_step*vertices_dot, &
        physical_min - difference_step*physical_min_dot, &
        physical_max - difference_step*physical_max_dot, &
        outer_min - difference_step*outer_min_dot, &
        outer_max - difference_step*outer_max_dot, &
        2.0_dp - difference_step*wave_number_dot, &
        4.0_dp - difference_step*sigma_max_dot, objective_minus, status)
    call evaluate_objective( &
        vertices + difference_step*vertices_dot, &
        physical_min + difference_step*physical_min_dot, &
        physical_max + difference_step*physical_max_dot, &
        outer_min + difference_step*outer_min_dot, &
        outer_max + difference_step*outer_max_dot, &
        2.0_dp + difference_step*wave_number_dot, &
        4.0_dp + difference_step*sigma_max_dot, objective_plus, status)
    call record_condition(abs(forward_product - &
        (objective_plus - objective_minus)/(2.0_dp*difference_step)) < &
        2.0e-9_dp, "Complete PML JVP matches primal central difference")

    call build_cartesian_pml_element_stretch_vjp( &
        vertices, cells, physical_min, physical_max, outer_min, outer_max, &
        2.0_dp, 4.0_dp, 2, stretch_bar, vertices_bar, physical_min_bar, &
        physical_max_bar, outer_min_bar, outer_max_bar, wave_number_bar, &
        sigma_max_bar, status)
    reverse_product = sum(vertices_bar*vertices_dot) + &
        dot_product(physical_min_bar, physical_min_dot) + &
        dot_product(physical_max_bar, physical_max_dot) + &
        dot_product(outer_min_bar, outer_min_dot) + &
        dot_product(outer_max_bar, outer_max_dot) + &
        wave_number_bar*wave_number_dot + sigma_max_bar*sigma_max_dot
    call record_condition(abs(reverse_product - forward_product) < 2.0e-12_dp, &
        "Complete PML reverse sweep satisfies the dot identity")

    call check_summary("Cartesian PML geometry-to-coefficient derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine evaluate_objective( &
            active_vertices, active_physical_min, active_physical_max, &
            active_outer_min, active_outer_max, wave_number, sigma_max, &
            objective, operator_status)
        real(dp), intent(in) :: active_vertices(:, :)
        real(dp), intent(in) :: active_physical_min(:), active_physical_max(:)
        real(dp), intent(in) :: active_outer_min(:), active_outer_max(:)
        real(dp), intent(in) :: wave_number, sigma_max
        real(dp), intent(out) :: objective
        integer, intent(out) :: operator_status

        complex(dp), allocatable :: active_stretch(:, :)
        complex(dp) :: gradient(3), local_stretch(3), mass
        integer :: active_cell

        objective = 0.0_dp
        call build_cartesian_pml_element_stretch( &
            active_vertices, cells, active_physical_min, active_physical_max, &
            active_outer_min, active_outer_max, wave_number, sigma_max, 2, &
            active_stretch, operator_status)
        if (operator_status /= 0) return
        do active_cell = 1, size(cells, 2)
            local_stretch = [ &
                active_stretch(:, active_cell), cmplx(1.0_dp, 0.0_dp, dp)]
            call cartesian_scalar_helmholtz_pml_coefficients( &
                local_stretch, gradient, mass, operator_status)
            if (operator_status /= 0) return
            objective = objective + &
                real(sum(conjg(gradient_bar(:, active_cell))*gradient) + &
                conjg(mass_bar(active_cell))*mass, dp)
        end do
    end subroutine evaluate_objective

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_cartesian_pml_chain_ad
