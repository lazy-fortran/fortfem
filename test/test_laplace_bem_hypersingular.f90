program test_laplace_bem_hypersingular
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_laplace_hypersingular_linear
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: panel_nodes(2, 4)
    real(dp) :: expected, matrix(4, 4), panel_end(2, 4), panel_start(2, 4)
    real(dp) :: alternating(4)
    integer :: i, status
    logical :: all_passed

    all_passed = .true.

    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [2.0_dp, 0.0_dp]
    panel_nodes(:, 1) = [1, 2]
    call assemble_laplace_hypersingular_linear( &
        panel_start(:, 1:1), panel_end(:, 1:1), panel_nodes(:, 1:1), &
        2, 24, matrix(1:2, 1:2), status)
    expected = (1.5_dp - log(2.0_dp)) / (2.0_dp * pi)
    call record_condition(status == 0 .and. &
        maxval(abs(matrix(1:2, 1:2) - reshape( &
        [expected, -expected, -expected, expected], [2, 2]))) < 2.0e-14_dp, &
        "Hypersingular linear form matches an exact one-panel integral")

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
    call assemble_laplace_hypersingular_linear( &
        panel_start, panel_end, panel_nodes, 4, 24, matrix, status)
    call record_condition(status == 0 .and. &
        maxval(abs(matrix - transpose(matrix))) < 2.0e-14_dp, &
        "Hypersingular Galerkin matrix is symmetric")
    call record_condition(maxval(abs(matmul(matrix, [(1.0_dp, i=1, 4)]))) < &
        2.0e-14_dp, "Hypersingular operator annihilates a constant trace")
    alternating = [1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp]
    call record_condition(dot_product( &
        alternating, matmul(matrix, alternating)) > 0.0_dp, &
        "Hypersingular energy is positive for a nonconstant trace")

    panel_nodes(2, 1) = 5
    call assemble_laplace_hypersingular_linear( &
        panel_start(:, 1:1), panel_end(:, 1:1), panel_nodes(:, 1:1), &
        2, 24, matrix(1:2, 1:2), status)
    call record_condition(status /= 0, &
        "Hypersingular assembly rejects an invalid node index")

    call check_summary("Laplace BEM hypersingular")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_laplace_bem_hypersingular
