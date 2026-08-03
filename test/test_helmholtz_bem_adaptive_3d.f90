program test_helmholtz_bem_adaptive_3d
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_helmholtz_single_layer_p0_adaptive_3d, &
        assemble_laplace_single_layer_p0_adaptive_3d
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: helmholtz(:, :), scaled(:, :)
    real(dp), allocatable :: laplace(:, :)
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
        vertices, triangles, 1.0e-3_dp, 3, laplace, status)
    call assemble_helmholtz_single_layer_p0_adaptive_3d( &
        vertices, triangles, 0.0_dp, 6, 1.0e-3_dp, 3, helmholtz, status)
    call record_condition(status == 0 .and. &
        maxval(abs(helmholtz - cmplx(laplace, 0.0_dp, dp))) < 2.0e-14_dp, &
        "Adaptive Helmholtz single layer has the exact Laplace limit")

    call assemble_helmholtz_single_layer_p0_adaptive_3d( &
        vertices, triangles, 0.7_dp, 6, 1.0e-3_dp, 3, helmholtz, status)
    call record_condition(status == 0 .and. &
        maxval(abs(helmholtz - transpose(helmholtz))) < 2.0e-14_dp, &
        "Adaptive Helmholtz single layer is complex symmetric")

    call assemble_helmholtz_single_layer_p0_adaptive_3d( &
        2.0_dp*vertices, triangles, 0.35_dp, 6, 1.0e-3_dp, 3, scaled, &
        status)
    call record_condition(maxval(abs(scaled - 8.0_dp*helmholtz)) < 2.0e-8_dp, &
        "Adaptive Helmholtz single layer obeys exact wave-scaled cubic scaling")

    call assemble_helmholtz_single_layer_p0_adaptive_3d( &
        vertices, triangles, 0.7_dp, 6, -1.0_dp, 3, scaled, status)
    call record_condition(status /= 0, &
        "Adaptive Helmholtz single layer rejects a negative tolerance")
    call check_summary("Adaptive three-dimensional Helmholtz BEM")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_helmholtz_bem_adaptive_3d
