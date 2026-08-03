program test_laplace_bem_adaptive_3d
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_laplace_single_layer_p0_adaptive_3d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), allocatable :: matrix(:, :), scaled_matrix(:, :)
    real(dp) :: vertices(3, 4)
    integer :: status, triangles(3, 2)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [1.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 1.0_dp, 0.0_dp]
    triangles(:, 1) = [1, 2, 3]
    triangles(:, 2) = [1, 3, 4]
    call assemble_laplace_single_layer_p0_adaptive_3d( &
        vertices, triangles, 1.0e-3_dp, 3, matrix, status)
    call record_condition(status == 0 .and. &
        maxval(abs(matrix - transpose(matrix))) < 2.0e-14_dp .and. &
        all(matrix > 0.0_dp), &
        "Adaptive 3D single layer resolves a shared-edge panel pair")

    call assemble_laplace_single_layer_p0_adaptive_3d( &
        2.0_dp*vertices, triangles, 1.0e-3_dp, 3, scaled_matrix, status)
    call record_condition(maxval(abs( &
        scaled_matrix - 8.0_dp*matrix)) < 2.0e-8_dp, &
        "Adaptive 3D single layer obeys exact geometric cubic scaling")

    call assemble_laplace_single_layer_p0_adaptive_3d( &
        vertices, triangles, -1.0_dp, 8, scaled_matrix, status)
    call record_condition(status /= 0, &
        "Adaptive 3D single layer rejects a negative tolerance")
    call check_summary("Adaptive three-dimensional Laplace BEM")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_laplace_bem_adaptive_3d
