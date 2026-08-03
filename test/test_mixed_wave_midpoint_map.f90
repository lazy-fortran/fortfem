program test_mixed_wave_midpoint_map
    !! Independent oracle for the reusable mixed-wave Cayley map.
    !!
    !! The expected map is assembled from the two block matrices in the
    !! midpoint equations, rather than by probing the state update.  This
    !! keeps the map contract useful to structure diagnostics and clients
    !! that need a matrix representation for small modal systems.
    use check, only: check_condition, check_summary
    use fortfem_time, only: advance_mixed_wave_midpoint, &
        assemble_mixed_wave_midpoint_map
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_SINGULAR
    implicit none

    real(dp), parameter :: mass_q(1, 1) = reshape([2.0_dp], [1, 1])
    real(dp), parameter :: mass_v(1, 1) = reshape([3.0_dp], [1, 1])
    real(dp), parameter :: coupling(1, 1) = reshape([1.25_dp], [1, 1])
    real(dp), parameter :: time_step = 0.4_dp
    real(dp) :: map(2, 2), map_small(1, 1), expected(2, 2), determinant
    real(dp) :: q(1), v(1), mapped_state(2)
    real(dp) :: singular_q(1, 1), singular_v(1, 1), singular_coupling(1, 1)
    real(dp) :: denominator, half_coupling
    type(fortsparse_status_t) :: status

    call assemble_mixed_wave_midpoint_map( &
        mass_q, mass_v, coupling, time_step, map, status)
    denominator = mass_q(1, 1)*mass_v(1, 1) + &
        (0.5_dp*time_step*coupling(1, 1))**2
    half_coupling = 0.5_dp*time_step*coupling(1, 1)
    expected(1, 1) = (mass_q(1, 1)*mass_v(1, 1) - half_coupling**2)/denominator
    expected(1, 2) = -mass_v(1, 1)*2.0_dp*half_coupling/denominator
    expected(2, 1) = mass_q(1, 1)*2.0_dp*half_coupling/denominator
    expected(2, 2) = (mass_q(1, 1)*mass_v(1, 1) - half_coupling**2)/denominator
    determinant = map(1, 1)*map(2, 2) - map(1, 2)*map(2, 1)

    call check_condition(status%code == 0 .and. &
        maxval(abs(map - expected)) < 2.0e-13_dp, &
        "mixed-wave midpoint map matches independent Cayley oracle")
    call check_condition(abs(determinant - 1.0_dp) < 2.0e-13_dp, &
        "mixed-wave midpoint map has unit determinant")

    q = [0.7_dp]
    v = [-0.2_dp]
    mapped_state = matmul(map, [q(1), v(1)])
    call advance_mixed_wave_midpoint( &
        mass_q, mass_v, coupling, time_step, q, v, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs([q(1), v(1)] - mapped_state)) < 2.0e-13_dp, &
        "mixed-wave midpoint update agrees with its explicit map")

    call assemble_mixed_wave_midpoint_map( &
        mass_q, mass_v, coupling, time_step, map_small, status)
    call check_condition(status%code == FORTSPARSE_INVALID_MATRIX .and. &
        maxval(abs(map_small)) == 0.0_dp, &
        "mixed-wave midpoint map rejects an output with the wrong shape")

    singular_q = 0.0_dp
    singular_v = 0.0_dp
    singular_coupling = 0.0_dp
    call assemble_mixed_wave_midpoint_map( &
        singular_q, singular_v, singular_coupling, time_step, map, status)
    call check_condition(status%code == FORTSPARSE_SINGULAR, &
        "mixed-wave midpoint map reports a singular block")

    call check_summary("Mixed-wave midpoint map")
end program test_mixed_wave_midpoint_map
