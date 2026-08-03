program test_maxwell_cfie_regularized_3d
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_maxwell_regularized_cfie_rwg_3d
    use fortfem_feec, only: assemble_maxwell_rwg_rbc_pairing
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: efie(:, :), manual_product(:, :)
    complex(dp), allocatable :: matrix(:, :), mfie(:, :), product(:, :)
    complex(dp), allocatable :: regularizer(:, :)
    complex(dp), allocatable :: scaled_efie(:, :), scaled_matrix(:, :)
    complex(dp), allocatable :: scaled_mfie(:, :), scaled_product(:, :)
    complex(dp), allocatable :: scaled_regularizer(:, :)
    real(dp), allocatable :: mass(:, :)
    real(dp) :: vertices(3, 4)
    integer :: index, status, triangles(3, 4)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    triangles(:, 1) = [1, 3, 2]
    triangles(:, 2) = [1, 2, 4]
    triangles(:, 3) = [1, 4, 3]
    triangles(:, 4) = [2, 3, 4]

    call assemble_maxwell_regularized_cfie_rwg_3d( &
        vertices, triangles, 1.1_dp, 1.6_dp, 2, 1.0e-5_dp, 1, matrix, &
        efie, mfie, regularizer, product, status)
    call assemble_maxwell_rwg_rbc_pairing( &
        vertices, triangles, 3, mass, status)
    allocate(manual_product(6, 6))
    manual_product = cmplx(0.0_dp, 0.0_dp, dp)
    do index = 1, 6
        manual_product = manual_product + outer_product( &
            regularizer(:, index), efie(index, :))/mass(index, index)
    end do
    call record_condition(status == 0 .and. maxval(abs( &
        mass - diagonal_matrix([(mass(index, index), index=1, 6)]))) < &
        2.0e-14_dp .and. maxval(abs(product - manual_product)) < 3.0e-13_dp, &
        "CFIE product uses the stable tetrahedral RWG-RBC mass inverse")
    call record_condition(maxval(abs(matrix - (mfie - product))) < 3.0e-14_dp, &
        "regularized CFIE combines compatible MFIE and electric product forms")

    call assemble_maxwell_regularized_cfie_rwg_3d( &
        2.0_dp*vertices, triangles, 0.55_dp, 1.6_dp, 2, 1.0e-5_dp, 1, &
        scaled_matrix, scaled_efie, scaled_mfie, scaled_regularizer, &
        scaled_product, status)
    call record_condition(maxval(abs(scaled_matrix - 2.0_dp*matrix)) < &
        2.0e-10_dp, &
        "regularized CFIE obeys exact reciprocal wave-geometry scaling")

    call check_summary("Three-dimensional regularized Maxwell CFIE")
    if (.not. all_passed) error stop 1

contains

    pure function outer_product(first, second) result(matrix_result)
        complex(dp), intent(in) :: first(:), second(:)
        complex(dp) :: matrix_result(size(first), size(second))
        integer :: column

        do column = 1, size(second)
            matrix_result(:, column) = first*second(column)
        end do
    end function outer_product

    pure function diagonal_matrix(diagonal) result(matrix_result)
        real(dp), intent(in) :: diagonal(:)
        real(dp) :: matrix_result(size(diagonal), size(diagonal))
        integer :: local_index

        matrix_result = 0.0_dp
        do local_index = 1, size(diagonal)
            matrix_result(local_index, local_index) = diagonal(local_index)
        end do
    end function diagonal_matrix

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_cfie_regularized_3d
