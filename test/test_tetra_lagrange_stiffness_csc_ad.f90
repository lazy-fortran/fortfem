program test_tetra_lagrange_stiffness_csc_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_feec, only: assemble_tetra_lagrange_stiffness_csc, &
        assemble_tetra_lagrange_stiffness_csc_jvp, &
        assemble_tetra_lagrange_stiffness_csc_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: degree = 2, quadrature_degree = 6
    integer, parameter :: tetrahedra(4, 2) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    real(dp), allocatable :: matrix_values_bar(:)
    type(csc_t) :: matrix, matrix_dot, minus, plus
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: stiffness_coefficient, mass_coefficient
    real(dp) :: stiffness_coefficient_dot, mass_coefficient_dot
    real(dp) :: stiffness_coefficient_bar, mass_coefficient_bar
    real(dp) :: lhs, relative_error, rhs
    integer :: entry, failures

    failures = 0
    vertices = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp], [3, 5])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, &
        -0.03_dp, 0.01_dp, 0.02_dp, &
        0.025_dp, -0.015_dp, 0.005_dp, &
        -0.01_dp, 0.03_dp, -0.02_dp, &
        0.02_dp, -0.01_dp, 0.025_dp], [3, 5])
    stiffness_coefficient = 1.7_dp
    mass_coefficient = 0.8_dp
    stiffness_coefficient_dot = 0.13_dp
    mass_coefficient_dot = -0.09_dp

    call assemble_tetra_lagrange_stiffness_csc( &
        vertices, tetrahedra, degree, quadrature_degree, matrix, &
        sparse_status, stiffness_coefficient, mass_coefficient)
    call check(sparse_status%code == 0, "global H1 primal assembly succeeds")
    call assemble_tetra_lagrange_stiffness_csc_jvp( &
        vertices, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient, mass_coefficient, vertices_dot, &
        stiffness_coefficient_dot, mass_coefficient_dot, matrix_dot, &
        sparse_status)
    call check(sparse_status%code == 0, "global H1 JVP succeeds")
    call assemble_tetra_lagrange_stiffness_csc( &
        vertices + step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        plus, sparse_status, &
        stiffness_coefficient + step*stiffness_coefficient_dot, &
        mass_coefficient + step*mass_coefficient_dot)
    call assemble_tetra_lagrange_stiffness_csc( &
        vertices - step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        minus, sparse_status, &
        stiffness_coefficient - step*stiffness_coefficient_dot, &
        mass_coefficient - step*mass_coefficient_dot)
    call check( &
        matrix_dot%nnz == matrix%nnz .and. plus%nnz == matrix%nnz .and. &
        minus%nnz == matrix%nnz .and. &
        all(matrix_dot%col_ptr == matrix%col_ptr) .and. &
        all(matrix_dot%row_idx == matrix%row_idx), &
        "global H1 JVP preserves the merged CSC pattern")
    relative_error = maxval(abs( &
        matrix_dot%val - (plus%val - minus%val)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(matrix_dot%val)))
    call check(relative_error < 5.0e-8_dp, &
        "global H1 JVP matches independent complete mesh reassembly")

    allocate(matrix_values_bar(matrix%nnz))
    do entry = 1, matrix%nnz
        matrix_values_bar(entry) = 0.003_dp*entry - 0.00001_dp*entry**2
    end do
    call assemble_tetra_lagrange_stiffness_csc_vjp( &
        vertices, tetrahedra, degree, quadrature_degree, &
        stiffness_coefficient, mass_coefficient, matrix_values_bar, &
        vertices_bar, stiffness_coefficient_bar, mass_coefficient_bar, &
        sparse_status)
    call check(sparse_status%code == 0, "global H1 VJP succeeds")
    lhs = dot_product(matrix_values_bar, matrix_dot%val)
    rhs = sum(vertices_bar*vertices_dot) + &
        stiffness_coefficient_bar*stiffness_coefficient_dot + &
        mass_coefficient_bar*mass_coefficient_dot
    call check(abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "global H1 products accumulate the shared-vertex adjoint")

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

end program test_tetra_lagrange_stiffness_csc_ad
