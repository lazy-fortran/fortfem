program test_maxwell_bc_surface
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_maxwell_rwg_rbc_pairing, build_maxwell_bc_transformation
    use fortfem_kinds, only: dp
    implicit none

    real(dp), allocatable :: pairing(:, :), refined_vertices(:, :)
    real(dp), allocatable :: transformation(:, :)
    real(dp) :: column_norms(6), reversed_norms(6), vertices(3, 4)
    integer, allocatable :: refined_triangles(:, :)
    integer :: column, status, triangles(3, 4)
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

    call build_maxwell_bc_transformation( &
        vertices, triangles, refined_vertices, refined_triangles, &
        transformation, status)
    call record_condition(status == 0 .and. &
        size(transformation, 1) == 72 .and. size(transformation, 2) == 6, &
        "closed tetrahedral surface builds six BC basis functions")
    do column = 1, 6
        call record_condition( &
            count(abs(transformation(:, column)) > 2.0e-14_dp) == 20, &
            "each valence-three BC function has the Bempp reference support")
        column_norms(column) = sum(transformation(:, column)**2)
    end do
    call sort_real(column_norms)
    call record_condition(maxval(abs(column_norms - [ &
        10.666666666666666_dp, 10.666666666666666_dp, &
        10.666666666666666_dp, 14.133333333333333_dp, &
        14.133333333333333_dp, 14.133333333333333_dp])) < 2.0e-13_dp, &
        "BC coefficients reproduce the independent Bempp tetrahedron fixture")
    triangles = triangles([1, 3, 2], :)
    call build_maxwell_bc_transformation( &
        vertices, triangles, refined_vertices, refined_triangles, &
        transformation, status)
    do column = 1, 6
        reversed_norms(column) = sum(transformation(:, column)**2)
    end do
    call sort_real(reversed_norms)
    call record_condition(status == 0 .and. &
        maxval(abs(reversed_norms - column_norms)) < 2.0e-13_dp, &
        "BC coefficient norm is invariant under global normal reversal")
    triangles = triangles([1, 3, 2], :)
    call assemble_maxwell_rwg_rbc_pairing( &
        vertices, triangles, 4, pairing, status)
    column_norms = [(pairing(column, column), column=1, 6)]
    call sort_real(column_norms)
    call record_condition(status == 0 .and. maxval(abs( &
        pairing - diagonal_matrix([(pairing(column, column), column=1, 6)]))) < &
        2.0e-14_dp, &
        "RWG-RBC dual pairing is diagonal on the reference tetrahedron")
    call record_condition(maxval(abs(column_norms - [ &
        0.4814814814814815_dp, 0.4814814814814815_dp, &
        0.4814814814814815_dp, 0.6809176411426006_dp, &
        0.6809176411426006_dp, 0.6809176411426006_dp])) < 2.0e-13_dp, &
        "RWG-RBC pairing reproduces the independent Bempp Gram fixture")

    call check_summary("Maxwell Buffa-Christiansen surface space")
    if (.not. all_passed) error stop 1

contains

    pure function diagonal_matrix(diagonal) result(matrix)
        real(dp), intent(in) :: diagonal(:)
        real(dp) :: matrix(size(diagonal), size(diagonal))
        integer :: index

        matrix = 0.0_dp
        do index = 1, size(diagonal)
            matrix(index, index) = diagonal(index)
        end do
    end function diagonal_matrix

    subroutine sort_real(values)
        real(dp), intent(inout) :: values(:)
        real(dp) :: temporary
        integer :: first, second

        do first = 1, size(values) - 1
            do second = first + 1, size(values)
                if (values(second) >= values(first)) cycle
                temporary = values(first)
                values(first) = values(second)
                values(second) = temporary
            end do
        end do
    end subroutine sort_real

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_bc_surface
