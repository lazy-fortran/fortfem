program test_structured_tetra_box_mesh
    use check, only: check_condition, check_summary
    use fortfem_core, only: generate_structured_tetra_box_mesh
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    implicit none

    real(dp), parameter :: bounds(3, 2) = reshape([ &
        -1.0_dp, 2.0_dp, 4.0_dp, &
        1.0_dp, 5.0_dp, 10.0_dp], [3, 2])
    integer, parameter :: counts(3) = [2, 1, 3]
    integer, allocatable :: tetrahedra(:, :)
    real(dp), allocatable :: vertices(:, :)
    real(dp) :: jacobian(3, 3), volume
    integer :: cell, status
    logical :: all_passed

    all_passed = .true.
    call generate_structured_tetra_box_mesh( &
        bounds, counts, vertices, tetrahedra, status)
    call record_condition(status == 0, "structured tetra box generation succeeds")
    call record_condition( &
        size(vertices, 2) == product(counts + 1), &
        "structured tetra box has the tensor-product vertex count")
    call record_condition( &
        size(tetrahedra, 2) == 6*product(counts), &
        "structured tetra box has six tetrahedra per brick")
    call record_condition( &
        all(abs(minval(vertices, dim=2) - bounds(:, 1)) < 1.0e-14_dp) .and. &
        all(abs(maxval(vertices, dim=2) - bounds(:, 2)) < 1.0e-14_dp), &
        "structured tetra box reaches every requested bound")

    volume = 0.0_dp
    do cell = 1, size(tetrahedra, 2)
        jacobian(:, 1) = &
            vertices(:, tetrahedra(2, cell)) - &
            vertices(:, tetrahedra(1, cell))
        jacobian(:, 2) = &
            vertices(:, tetrahedra(3, cell)) - &
            vertices(:, tetrahedra(1, cell))
        jacobian(:, 3) = &
            vertices(:, tetrahedra(4, cell)) - &
            vertices(:, tetrahedra(1, cell))
        call record_condition( &
            det3(jacobian) > 0.0_dp, &
            "structured tetra box orients every tetrahedron positively")
        volume = volume + det3(jacobian)/6.0_dp
    end do
    call record_condition( &
        abs(volume - product(bounds(:, 2) - bounds(:, 1))) < 1.0e-12_dp, &
        "structured tetra box exactly partitions the requested volume")
    call check_summary("structured tetrahedral box mesh")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_structured_tetra_box_mesh
