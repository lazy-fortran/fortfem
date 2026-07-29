program test_helmholtz_symmetric_coupling
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_helmholtz_symmetric_coupling_p1_p0, &
        assemble_laplace_symmetric_coupling_p1_p0
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: vertices(2, 3) = reshape([ &
        0.0_dp, 0.0_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.25_dp], [2, 3])
    integer, parameter :: triangles(3, 1) = reshape([1, 2, 3], [3, 1])
    integer, parameter :: panel_nodes(2, 3) = reshape([1, 2, 2, 3, 3, 1], [2, 3])
    real(dp) :: laplace(6, 6), panel_end(2, 3), panel_start(2, 3)
    real(dp) :: panel_length(3)
    complex(dp) :: helmholtz(6, 6), low_frequency(6, 6)
    complex(dp) :: zero_mean_flux(3)
    integer :: panel, status
    logical :: all_passed

    all_passed = .true.
    do panel = 1, 3
        panel_start(:, panel) = vertices(:, panel_nodes(1, panel))
        panel_end(:, panel) = vertices(:, panel_nodes(2, panel))
        panel_length(panel) = norm2(panel_end(:, panel) - panel_start(:, panel))
    end do

    call assemble_helmholtz_symmetric_coupling_p1_p0( &
        vertices, triangles, panel_start, panel_end, panel_nodes, &
        1.2_dp, 28, helmholtz, status)
    call record_condition(status == 0, "Helmholtz Costabel coupling assembles")
    call record_condition(maxval(abs(helmholtz(1:3, 4:6) + &
        transpose(helmholtz(4:6, 1:3)))) < 3.0e-14_dp, &
        "Helmholtz transmission cross blocks have opposite transpose signs")
    call record_condition(maxval(abs(helmholtz(4:6, 4:6) - &
        transpose(helmholtz(4:6, 4:6)))) < 3.0e-14_dp, &
        "Helmholtz single-layer block is complex symmetric")

    call assemble_laplace_symmetric_coupling_p1_p0( &
        vertices, triangles, panel_start, panel_end, panel_nodes, 28, &
        laplace, status)
    call assemble_helmholtz_symmetric_coupling_p1_p0( &
        vertices, triangles, panel_start, panel_end, panel_nodes, &
        1.0e-5_dp, 28, low_frequency, status)
    call record_condition(maxval(abs(low_frequency(1:3, 1:6) - &
        laplace(1:3, 1:6))) < 3.0e-9_dp .and. &
        maxval(abs(low_frequency(4:6, 1:3) - &
        laplace(4:6, 1:3))) < 3.0e-9_dp, &
        "Helmholtz stiffness, W, and K have the Laplace low-frequency limit")
    zero_mean_flux = [cmplx(panel_length(2), 0.0_dp, dp), &
        cmplx(-panel_length(1), 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp)]
    call record_condition(maxval(abs(matmul(low_frequency(4:6, 4:6) - &
        laplace(4:6, 4:6), zero_mean_flux))) < 3.0e-9_dp, &
        "Helmholtz V tends to Laplace V on zero-total-flux traces")

    call check_summary("Helmholtz symmetric FEM-BEM coupling")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record_condition
end program test_helmholtz_symmetric_coupling
