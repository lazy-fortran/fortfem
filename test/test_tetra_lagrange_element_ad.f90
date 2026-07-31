program test_tetra_lagrange_element_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: assemble_tetra_lagrange_stiffness_element, &
        assemble_tetra_lagrange_stiffness_element_jvp, &
        assemble_tetra_lagrange_stiffness_element_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: degree = 2, quadrature_degree = 6
    integer, parameter :: dof_count = (degree + 1)*(degree + 2)*(degree + 3)/6
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 4), vertices_dot(3, 4), vertices_bar(3, 4)
    real(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    real(dp), allocatable :: minus(:, :), plus(:, :)
    real(dp) :: matrix_bar(dof_count, dof_count)
    real(dp) :: stiffness_coefficient, mass_coefficient
    real(dp) :: stiffness_coefficient_dot, mass_coefficient_dot
    real(dp) :: stiffness_coefficient_bar, mass_coefficient_bar
    real(dp) :: lhs, relative_error, rhs
    integer :: column, failures, row, status

    failures = 0
    vertices = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, &
        1.3_dp, -0.1_dp, 0.25_dp, &
        0.25_dp, 0.9_dp, 0.35_dp, &
        0.0_dp, 0.1_dp, 1.4_dp], [3, 4])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, &
        -0.03_dp, 0.01_dp, 0.02_dp, &
        0.025_dp, -0.015_dp, 0.005_dp, &
        -0.01_dp, 0.03_dp, -0.02_dp], [3, 4])
    stiffness_coefficient = 1.7_dp
    mass_coefficient = 0.8_dp
    stiffness_coefficient_dot = 0.13_dp
    mass_coefficient_dot = -0.09_dp
    do column = 1, dof_count
        do row = 1, dof_count
            matrix_bar(row, column) = &
                0.002_dp*row - 0.0015_dp*column
        end do
    end do

    call assemble_tetra_lagrange_stiffness_element( &
        vertices, degree, quadrature_degree, matrix, status, &
        stiffness_coefficient, mass_coefficient)
    call check(status == 0, "H1 element primal assembly succeeds")
    call assemble_tetra_lagrange_stiffness_element_jvp( &
        vertices, degree, quadrature_degree, stiffness_coefficient, &
        mass_coefficient, vertices_dot, stiffness_coefficient_dot, &
        mass_coefficient_dot, matrix_dot, status)
    call check(status == 0, "H1 element JVP succeeds")
    call assemble_tetra_lagrange_stiffness_element( &
        vertices + step*vertices_dot, degree, quadrature_degree, plus, status, &
        stiffness_coefficient + step*stiffness_coefficient_dot, &
        mass_coefficient + step*mass_coefficient_dot)
    call assemble_tetra_lagrange_stiffness_element( &
        vertices - step*vertices_dot, degree, quadrature_degree, minus, status, &
        stiffness_coefficient - step*stiffness_coefficient_dot, &
        mass_coefficient - step*mass_coefficient_dot)
    relative_error = maxval(abs( &
        matrix_dot - (plus - minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot)))
    call check(relative_error < 5.0e-8_dp, &
        "H1 element JVP matches independent complete reassembly")

    call assemble_tetra_lagrange_stiffness_element_vjp( &
        vertices, degree, quadrature_degree, stiffness_coefficient, &
        mass_coefficient, matrix_bar, vertices_bar, &
        stiffness_coefficient_bar, mass_coefficient_bar, status)
    call check(status == 0, "H1 element VJP succeeds")
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(vertices_bar*vertices_dot) + &
        stiffness_coefficient_bar*stiffness_coefficient_dot + &
        mass_coefficient_bar*mass_coefficient_dot
    call check(abs(lhs - rhs) < 2.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "H1 element products satisfy the complete adjoint identity")

    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_tetra_lagrange_element_ad
