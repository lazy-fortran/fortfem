program test_tetra_lagrange_curvilinear_pml_state_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        solve_tetra_lagrange_curvilinear_pml, &
        solve_tetra_lagrange_curvilinear_pml_jvp, &
        solve_tetra_lagrange_curvilinear_pml_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: difference_step = 1.0e-6_dp
    real(dp), parameter :: vertices(3, 4) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    integer, parameter :: tetrahedra(4, 1) = reshape([1, 2, 3, 4], [4, 1])
    integer, parameter :: degree = 1
    real(dp), parameter :: wave_number = 1.7_dp
    real(dp), parameter :: wave_number_dot = -0.13_dp
    logical, parameter :: constrained(4) = [.true., .true., .true., .false.]
    complex(dp), parameter :: constrained_values(4) = [ &
        cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.1_dp, 0.3_dp, dp), &
        cmplx(0.15_dp, 0.05_dp, dp), cmplx(0.0_dp, 0.0_dp, dp)]
    complex(dp), parameter :: load(4) = [ &
        cmplx(0.3_dp, 0.2_dp, dp), cmplx(-0.1_dp, 0.4_dp, dp), &
        cmplx(0.2_dp, -0.3_dp, dp), cmplx(0.8_dp, 0.1_dp, dp)]
    real(dp) :: vertices_dot(3, 4), wave_number_bar
    real(dp) :: vertices_bar(3, 4)
    complex(dp) :: stretch(3, 3, 1), stretch_dot(3, 3, 1)
    complex(dp) :: stretch_bar(3, 3, 1)
    complex(dp) :: load_dot(4), constrained_values_dot(4)
    complex(dp) :: load_bar(4), constrained_values_bar(4)
    complex(dp), allocatable :: solution(:), solution_dot(:)
    complex(dp), allocatable :: solution_minus(:), solution_plus(:)
    integer :: i, status, status_minus, status_plus
    logical :: all_passed
    real(dp) :: forward_pairing, reverse_pairing

    all_passed = .true.
    vertices_dot = reshape([ &
        0.03_dp, -0.02_dp, 0.04_dp, -0.01_dp, 0.05_dp, 0.02_dp, &
        -0.04_dp, 0.01_dp, 0.02_dp, 0.03_dp, -0.05_dp, 0.01_dp], [3, 4])
    load_dot = [ &
        cmplx(0.02_dp, -0.01_dp, dp), cmplx(-0.03_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(-0.02_dp, 0.03_dp, dp)]
    constrained_values_dot = [ &
        cmplx(-0.01_dp, 0.02_dp, dp), cmplx(0.03_dp, -0.02_dp, dp), &
        cmplx(0.02_dp, 0.01_dp, dp), cmplx(0.0_dp, 0.0_dp, dp)]
    stretch = cmplx(0.0_dp, 0.0_dp, dp)
    stretch(1, 1, 1) = cmplx(1.0_dp, 0.08_dp, dp)
    stretch(2, 2, 1) = cmplx(1.0_dp, 0.13_dp, dp)
    stretch(3, 3, 1) = cmplx(1.0_dp, 0.05_dp, dp)
    stretch(1, 2, 1) = cmplx(0.03_dp, -0.02_dp, dp)
    stretch(2, 1, 1) = cmplx(-0.01_dp, 0.04_dp, dp)
    stretch_dot(:, :, 1) = reshape([ &
        cmplx(0.02_dp, -0.01_dp, dp), cmplx(-0.03_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(0.03_dp, 0.01_dp, dp), &
        cmplx(-0.02_dp, 0.03_dp, dp), cmplx(0.01_dp, -0.04_dp, dp), &
        cmplx(0.02_dp, 0.01_dp, dp), cmplx(-0.01_dp, 0.02_dp, dp), &
        cmplx(0.03_dp, -0.02_dp, dp)], [3, 3])

    call solve_tetra_lagrange_curvilinear_pml( &
        vertices, tetrahedra, degree, stretch, wave_number, load, constrained, &
        constrained_values, solution, status)
    call record_condition(status == 0, "curvilinear scalar PML state solves")
    do i = 1, 4
        if (constrained(i)) call record_condition( &
            abs(solution(i) - constrained_values(i)) < 1.0e-14_dp, &
            "curvilinear scalar PML state preserves constrained values")
    end do

    call solve_tetra_lagrange_curvilinear_pml_jvp( &
        vertices, tetrahedra, degree, stretch, wave_number, load, constrained, &
        constrained_values, vertices_dot, stretch_dot, wave_number_dot, &
        load_dot, constrained_values_dot, solution_dot, status)
    call solve_tetra_lagrange_curvilinear_pml( &
        vertices - difference_step*vertices_dot, tetrahedra, degree, &
        stretch - difference_step*stretch_dot, &
        wave_number - difference_step*wave_number_dot, &
        load - difference_step*load_dot, constrained, &
        constrained_values - difference_step*constrained_values_dot, &
        solution_minus, status_minus)
    call solve_tetra_lagrange_curvilinear_pml( &
        vertices + difference_step*vertices_dot, tetrahedra, degree, &
        stretch + difference_step*stretch_dot, &
        wave_number + difference_step*wave_number_dot, &
        load + difference_step*load_dot, constrained, &
        constrained_values + difference_step*constrained_values_dot, &
        solution_plus, status_plus)
    call record_condition(status == 0 .and. status_minus == 0 .and. &
        status_plus == 0, "curvilinear scalar PML state JVP inputs are accepted")
    call record_condition(maxval(abs(solution_dot - &
        (solution_plus - solution_minus)/(2.0_dp*difference_step))) < 4.0e-8_dp, &
        "curvilinear scalar PML state JVP matches re-solves")

    call solve_tetra_lagrange_curvilinear_pml_vjp( &
        vertices, tetrahedra, degree, stretch, wave_number, load, constrained, &
        constrained_values, solution, solution_dot, vertices_bar, stretch_bar, &
        wave_number_bar, load_bar, constrained_values_bar, status)
    forward_pairing = real(sum(conjg(solution_dot)*solution_dot), dp)
    reverse_pairing = real(sum(conjg(load_bar)*load_dot) + &
        sum(conjg(constrained_values_bar)*constrained_values_dot) + &
        sum(conjg(stretch_bar)*stretch_dot), dp) + &
        sum(vertices_bar*vertices_dot) + wave_number_bar*wave_number_dot
    call record_condition(abs(forward_pairing - reverse_pairing) < 3.0e-10_dp, &
        "curvilinear scalar PML state VJP satisfies the complex dot identity")

    call check_summary("Tetrahedral scalar curvilinear PML state derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_lagrange_curvilinear_pml_state_ad
