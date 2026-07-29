program test_laplace_bem_mixed_double_layer
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_laplace_double_layer_constant, &
        assemble_laplace_double_layer_mixed_linear
    use fortfem_kinds, only: dp
    implicit none

    real(dp) :: constant_matrix(4, 4), mixed_matrix(4, 8)
    real(dp) :: continuous_matrix(4, 4)
    real(dp) :: panel_end(2, 4), panel_start(2, 4)
    integer :: continuous_nodes(2, 4), discontinuous_nodes(2, 4)
    integer :: panel, status
    logical :: all_passed

    all_passed = .true.
    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    panel_start(:, 2) = [1.0_dp, 0.0_dp]
    panel_end(:, 2) = [1.0_dp, 1.0_dp]
    panel_start(:, 3) = [1.0_dp, 1.0_dp]
    panel_end(:, 3) = [0.0_dp, 1.0_dp]
    panel_start(:, 4) = [0.0_dp, 1.0_dp]
    panel_end(:, 4) = [0.0_dp, 0.0_dp]
    discontinuous_nodes = reshape([1, 2, 3, 4, 5, 6, 7, 8], [2, 4])

    call assemble_laplace_double_layer_constant( &
        panel_start, panel_end, 24, constant_matrix, status)
    call assemble_laplace_double_layer_mixed_linear( &
        panel_start, panel_end, discontinuous_nodes, 8, 24, mixed_matrix, status)
    do panel = 1, 4
        call record_condition(maxval(abs( &
            mixed_matrix(:, 2 * panel - 1) + mixed_matrix(:, 2 * panel) - &
            constant_matrix(:, panel))) < 3.0e-14_dp, &
            "Linear source partition of unity recovers each constant panel")
    end do

    continuous_nodes(:, 1) = [1, 2]
    continuous_nodes(:, 2) = [2, 3]
    continuous_nodes(:, 3) = [3, 4]
    continuous_nodes(:, 4) = [4, 1]
    call assemble_laplace_double_layer_mixed_linear( &
        panel_start, panel_end, continuous_nodes, 4, 24, &
        continuous_matrix, status)
    do panel = 1, 4
        call record_condition( &
            abs(sum(continuous_matrix(panel, :)) + 0.5_dp) < 3.0e-14_dp, &
            "Closed-boundary mixed double layer maps one to minus one half")
    end do

    continuous_nodes(1, 1) = 0
    call assemble_laplace_double_layer_mixed_linear( &
        panel_start, panel_end, continuous_nodes, 4, 24, &
        continuous_matrix, status)
    call record_condition(status /= 0, "Invalid source node is rejected")

    call check_summary("Laplace BEM mixed double-layer")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_laplace_bem_mixed_double_layer
