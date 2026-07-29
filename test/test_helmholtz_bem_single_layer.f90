program test_helmholtz_bem_single_layer
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_helmholtz_single_layer_constant
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), parameter :: self_reference = &
        (0.24623816383735733_dp, 0.23983991217241335_dp)
    complex(dp), parameter :: separated_reference = &
        (-0.003630378324870141_dp, -0.07066865243283912_dp)
    complex(dp), parameter :: adjacent_reference = &
        (-0.02690275227956196_dp, 0.16774273359470834_dp)
    complex(dp) :: matrix(2, 2)
    real(dp) :: panel_end(2, 2), panel_start(2, 2)
    integer :: status
    logical :: all_passed

    all_passed = .true.

    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    call assemble_helmholtz_single_layer_constant( &
        panel_start(:, 1:1), panel_end(:, 1:1), 1.0_dp, 48, &
        matrix(1:1, 1:1), status)
    call record_condition(status == 0 .and. &
        abs(matrix(1, 1) - self_reference) < 5.0e-11_dp, &
        "Helmholtz self panel matches SciPy 1.18 quadrature")

    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    panel_start(:, 2) = [2.0_dp, 0.0_dp]
    panel_end(:, 2) = [3.0_dp, 0.0_dp]
    call assemble_helmholtz_single_layer_constant( &
        panel_start, panel_end, 2.0_dp, 48, matrix, status)
    call record_condition(status == 0 .and. &
        abs(matrix(1, 2) - separated_reference) < 5.0e-11_dp, &
        "Regular Helmholtz panel interaction matches SciPy 1.18 quadrature")
    call record_condition(abs(matrix(2, 1) - matrix(1, 2)) < 1.0e-14_dp, &
        "Helmholtz single-layer Galerkin matrix is complex symmetric")

    panel_start(:, 1) = [1.0_dp, 0.0_dp]
    panel_end(:, 1) = [0.0_dp, 0.0_dp]
    panel_start(:, 2) = [0.0_dp, 0.0_dp]
    panel_end(:, 2) = [0.0_dp, 1.0_dp]
    call assemble_helmholtz_single_layer_constant( &
        panel_start, panel_end, 1.5_dp, 48, matrix, status)
    call record_condition(status == 0 .and. &
        abs(matrix(1, 2) - adjacent_reference) < 2.0e-9_dp, &
        "Endpoint-singular Helmholtz interaction matches SciPy 1.18 quadrature")

    call assemble_helmholtz_single_layer_constant( &
        panel_start, panel_end, 0.0_dp, 48, matrix, status)
    call record_condition(status /= 0, &
        "Helmholtz single-layer assembly rejects zero wavenumber")

    call check_summary("Helmholtz BEM single-layer")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_helmholtz_bem_single_layer
