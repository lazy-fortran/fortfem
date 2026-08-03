program test_helmholtz_bem_hypersingular
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_helmholtz_hypersingular_linear
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), parameter :: diagonal_reference = &
        (0.17401430520778820_dp, 0.17904385560359232_dp)
    complex(dp), parameter :: off_diagonal_reference = &
        (-0.29713338712646714_dp, -0.29896381168979896_dp)
    complex(dp) :: expected(2, 2), matrix(4, 4)
    integer :: panel_nodes(2, 4), status
    real(dp) :: panel_end(2, 4), panel_start(2, 4)
    logical :: all_passed

    all_passed = .true.
    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    panel_nodes(:, 1) = [1, 2]
    expected(1, 1) = diagonal_reference
    expected(2, 2) = diagonal_reference
    expected(1, 2) = off_diagonal_reference
    expected(2, 1) = off_diagonal_reference

    call assemble_helmholtz_hypersingular_linear( &
        panel_start(:, 1:1), panel_end(:, 1:1), panel_nodes(:, 1:1), &
        2, 1.0_dp, 48, matrix(1:2, 1:2), status)
    call record_condition(status == 0 .and. &
        maxval(abs(matrix(1:2, 1:2) - expected)) < 8.0e-11_dp, &
        "Helmholtz hypersingular self panel matches the regularized SciPy form")

    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    panel_start(:, 2) = [1.0_dp, 0.0_dp]
    panel_end(:, 2) = [1.0_dp, 1.0_dp]
    panel_start(:, 3) = [1.0_dp, 1.0_dp]
    panel_end(:, 3) = [0.0_dp, 1.0_dp]
    panel_start(:, 4) = [0.0_dp, 1.0_dp]
    panel_end(:, 4) = [0.0_dp, 0.0_dp]
    panel_nodes(:, 1) = [1, 2]
    panel_nodes(:, 2) = [2, 3]
    panel_nodes(:, 3) = [3, 4]
    panel_nodes(:, 4) = [4, 1]
    call assemble_helmholtz_hypersingular_linear( &
        panel_start, panel_end, panel_nodes, 4, 1.0_dp, 24, matrix, status)
    call record_condition(status == 0 .and. &
        maxval(abs(matrix - transpose(matrix))) < 3.0e-13_dp, &
        "Helmholtz hypersingular Galerkin matrix is complex symmetric")

    call check_summary("Helmholtz BEM hypersingular")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_helmholtz_bem_hypersingular
