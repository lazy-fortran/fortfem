program test_surface_delta_load
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_surface_delta_load
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: trace_basis(3, 2) = reshape([ &
        1.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 1.0_dp, 0.5_dp], [3, 2])
    real(dp), parameter :: surface_weights(3) = [2.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: surface_source(3) = [3.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: expected_load(2) = [11.0_dp, 9.0_dp]
    real(dp) :: load(2)
    real(dp), parameter :: bad_weights(2) = [1.0_dp, 2.0_dp]
    type(fortsparse_status_t) :: status

    call assemble_surface_delta_load( &
        trace_basis, surface_weights, surface_source, load, status)
    call check_condition(status%code == 0, &
        "surface delta load accepts a positive trace quadrature")
    call check_condition(maxval(abs(load - expected_load)) < 1.0e-14_dp, &
        "surface delta load matches the explicit trace-transpose oracle")

    call assemble_surface_delta_load( &
        trace_basis, bad_weights, surface_source, load, status)
    call check_condition(status%code /= 0, &
        "surface delta load rejects incompatible quadrature sizes")
    call check_summary("surface delta load")
end program test_surface_delta_load
