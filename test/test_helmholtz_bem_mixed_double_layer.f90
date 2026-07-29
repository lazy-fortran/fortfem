program test_helmholtz_bem_mixed_double_layer
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_helmholtz_double_layer_constant, &
        assemble_helmholtz_double_layer_mixed_linear, &
        assemble_laplace_double_layer_mixed_linear
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: panel_start(2, 3) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.2_dp, 0.9_dp], [2, 3])
    real(dp), parameter :: panel_end(2, 3) = reshape([ &
        1.0_dp, 0.0_dp, 0.2_dp, 0.9_dp, 0.0_dp, 0.0_dp], [2, 3])
    integer, parameter :: panel_nodes(2, 3) = reshape([1, 2, 2, 3, 3, 1], [2, 3])
    complex(dp) :: constant_matrix(3, 3), mixed_matrix(3, 3)
    complex(dp) :: low_frequency(3, 3)
    real(dp) :: laplace_matrix(3, 3)
    integer :: source, status
    logical :: all_passed

    all_passed = .true.
    call assemble_helmholtz_double_layer_constant( &
        panel_start, panel_end, 1.7_dp, 28, constant_matrix, status)
    call assemble_helmholtz_double_layer_mixed_linear( &
        panel_start, panel_end, panel_nodes, 3, 1.7_dp, 28, &
        mixed_matrix, status)
    call record_condition(status == 0, "Mixed Helmholtz double layer assembles")
    call record_condition(maxval(abs(matmul(mixed_matrix, &
        [(cmplx(1.0_dp, 0.0_dp, dp), source=1, 3)]) - &
        matmul(constant_matrix, &
        [(cmplx(1.0_dp, 0.0_dp, dp), source=1, 3)]))) < 2.0e-11_dp, &
        "Linear partition of unity reproduces the constant operator")

    call assemble_laplace_double_layer_mixed_linear( &
        panel_start, panel_end, panel_nodes, 3, 28, laplace_matrix, status)
    call assemble_helmholtz_double_layer_mixed_linear( &
        panel_start, panel_end, panel_nodes, 3, 1.0e-5_dp, 28, &
        low_frequency, status)
    call record_condition(maxval(abs(low_frequency - laplace_matrix)) < &
        2.0e-9_dp, "Mixed Helmholtz operator has the Laplace low-frequency limit")

    call check_summary("Mixed Helmholtz BEM double-layer")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record_condition
end program test_helmholtz_bem_mixed_double_layer
