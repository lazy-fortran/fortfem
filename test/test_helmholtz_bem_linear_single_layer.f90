program test_helmholtz_bem_linear_single_layer
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_helmholtz_single_layer_linear
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), parameter :: diagonal_reference = &
        (0.07222385862956911_dp, 0.06079605656882104_dp)
    complex(dp), parameter :: off_diagonal_reference = &
        (0.050895223289109816_dp, 0.05912389951738563_dp)
    complex(dp), parameter :: self_15_diagonal = &
        (0.05465153714629806_dp, 0.058754210190784244_dp)
    complex(dp), parameter :: self_15_off_diagonal = &
        (0.03265906095922309_dp, 0.05516462422252709_dp)
    complex(dp), parameter :: cross_11 = &
        (-0.009813119182056268_dp, 0.04169961664572386_dp)
    complex(dp), parameter :: cross_12 = &
        (-0.020471521792898992_dp, 0.03266982040411914_dp)
    complex(dp), parameter :: cross_21 = &
        (0.013195007877415087_dp, 0.051673679899141456_dp)
    complex(dp), parameter :: cross_22 = &
        (-0.009813119182024717_dp, 0.04169961664572387_dp)
    complex(dp) :: expected(3, 3), matrix(3, 3)
    integer :: panel_nodes(2, 2), status
    real(dp) :: panel_end(2, 2), panel_start(2, 2)
    logical :: all_passed

    all_passed = .true.
    panel_start(:, 1) = [0.0_dp, 0.0_dp]
    panel_end(:, 1) = [1.0_dp, 0.0_dp]
    panel_nodes(:, 1) = [1, 2]
    expected(1, 1) = diagonal_reference
    expected(2, 2) = diagonal_reference
    expected(1, 2) = off_diagonal_reference
    expected(2, 1) = off_diagonal_reference

    call assemble_helmholtz_single_layer_linear( &
        panel_start(:, 1:1), panel_end(:, 1:1), panel_nodes(:, 1:1), &
        2, 1.0_dp, 48, matrix(1:2, 1:2), status)
    call record_condition(status == 0 .and. &
        maxval(abs(matrix(1:2, 1:2) - expected(1:2, 1:2))) < 5.0e-11_dp, &
        "Linear Helmholtz self panel matches SciPy 1.18 weighted quadrature")
    call record_condition(maxval(abs( &
        matrix(1:2, 1:2) - transpose(matrix(1:2, 1:2)))) < &
        1.0e-14_dp, "Linear Helmholtz single layer is complex symmetric")

    panel_start(:, 1) = [1.0_dp, 0.0_dp]
    panel_end(:, 1) = [0.0_dp, 0.0_dp]
    panel_start(:, 2) = [0.0_dp, 0.0_dp]
    panel_end(:, 2) = [0.0_dp, 1.0_dp]
    panel_nodes(:, 1) = [1, 2]
    panel_nodes(:, 2) = [2, 3]
    expected = (0.0_dp, 0.0_dp)
    expected(1, 1) = self_15_diagonal
    expected(3, 3) = self_15_diagonal
    expected(1, 2) = self_15_off_diagonal + cross_11
    expected(2, 1) = expected(1, 2)
    expected(1, 3) = cross_12
    expected(3, 1) = cross_12
    expected(2, 2) = 2.0_dp * (self_15_diagonal + cross_21)
    expected(2, 3) = self_15_off_diagonal + cross_22
    expected(3, 2) = expected(2, 3)
    call assemble_helmholtz_single_layer_linear( &
        panel_start, panel_end, panel_nodes, 3, 1.5_dp, 48, matrix, status)
    call record_condition(status == 0 .and. &
        maxval(abs(matrix - expected)) < 2.0e-9_dp, &
        "Linear adjacent-panel moments match SciPy 1.18 weighted quadrature")

    panel_nodes(2, 1) = 4
    call assemble_helmholtz_single_layer_linear( &
        panel_start, panel_end, panel_nodes, 3, 1.0_dp, 48, matrix, status)
    call record_condition(status /= 0, &
        "Linear Helmholtz single layer rejects an invalid node index")

    call check_summary("Linear Helmholtz BEM single-layer")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_helmholtz_bem_linear_single_layer
