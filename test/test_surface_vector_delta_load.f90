program test_surface_vector_delta_load
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_surface_vector_delta_load
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: tangential_trace_basis(2, 2, 3) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], [2, 2, 3])
    real(dp), parameter :: surface_weights(2) = [2.0_dp, 1.0_dp]
    real(dp), parameter :: surface_current(2, 3) = reshape([ &
        3.0_dp, 4.0_dp, 5.0_dp, 1.0_dp, -2.0_dp, 2.0_dp], [2, 3])
    real(dp), parameter :: expected_load(2) = [7.0_dp, 1.0_dp]
    real(dp) :: load(2)
    real(dp), parameter :: bad_weights(1) = [1.0_dp]
    type(fortsparse_status_t) :: status

    call assemble_surface_vector_delta_load( &
        tangential_trace_basis, surface_weights, surface_current, load, status)
    call check_condition(status%code == 0, &
        "surface vector delta load accepts tangential trace data")
    call check_condition(maxval(abs(load - expected_load)) < 1.0e-14_dp, &
        "surface vector delta load matches the current-sheet pairing oracle")

    call assemble_surface_vector_delta_load( &
        tangential_trace_basis, bad_weights, surface_current, load, status)
    call check_condition(status%code /= 0, &
        "surface vector delta load rejects incompatible quadrature sizes")
    call check_summary("surface vector delta load")
end program test_surface_vector_delta_load
