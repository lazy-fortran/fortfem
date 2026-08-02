program test_collocation_transform_contract
    use check, only: check_condition, check_summary
    use fortfem_collocation_grid, only: &
        COLLOCATION_GRID_CONCENTRIC, COLLOCATION_GRID_LINEAR, &
        COLLOCATION_GRID_QUADRATURE, collocation_grid_chunk_bounds, &
        collocation_grid_flat_index, collocation_grid_metadata, &
        collocation_grid_point_count, collocation_grid_t, &
        collocation_grid_unflatten_index, initialize_collocation_grid, &
        validate_collocation_grid
    use fortfem_direct_fourier_transform, only: &
        direct_fourier_adjoint, direct_fourier_forward, &
        direct_fourier_plan_chunk_bounds, direct_fourier_plan_metadata, &
        direct_fourier_plan_mode_count, direct_fourier_plan_sample_count, &
        direct_fourier_plan_t, initialize_direct_fourier_plan, &
        validate_direct_fourier_plan
    use fortfem_kinds, only: dp, pi
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: n_first = 3, n_second = 4, n_point = n_first*n_second
    real(dp), parameter :: first_nodes(n_first) = [0.0_dp, 0.5_dp, 1.0_dp]
    real(dp), parameter :: second_nodes(n_second) = [0.0_dp, pi/3.0_dp, &
        2.0_dp*pi/3.0_dp, 2.0_dp*pi]
    real(dp), parameter :: first_weights(n_first) = [0.25_dp, 0.5_dp, 0.25_dp]
    real(dp), parameter :: second_weights(n_second) = [0.5_dp, 1.0_dp, 1.0_dp, 0.5_dp]
    integer, parameter :: endpoint_map(n_second) = [1, 2, 3, 1]
    integer, parameter :: axis_map(n_first) = [1, 1, 3]
    logical :: all_passed
    integer :: first, second, flat, mapped_first, mapped_second
    integer :: grid_kind, chunk_size, chunk, first_point, last_point
    real(dp) :: expected_weight, round_trip_error
    real(dp), allocatable :: copied_first(:), copied_second(:), copied_weights(:)
    integer, allocatable :: copied_endpoint(:), copied_axis(:), copied_start(:), copied_end(:)
    type(collocation_grid_t) :: grid, bad_grid, quadrature, concentric
    type(direct_fourier_plan_t) :: plan, bad_plan
    type(fortsparse_status_t) :: status
    real(dp), parameter :: angles(4) = [0.0_dp, 0.5_dp*pi, pi, 1.5_dp*pi]
    integer, parameter :: modes(4) = [0, 1, 2, 3]
    complex(dp), parameter :: samples(4) = [ &
        cmplx(1.0_dp, 0.5_dp, dp), cmplx(-2.0_dp, 1.0_dp, dp), &
        cmplx(0.25_dp, -0.25_dp, dp), cmplx(3.0_dp, -1.5_dp, dp)]
    complex(dp) :: coefficients(4), recovered(4)
    complex(dp) :: expected_coefficient

    all_passed = .true.
    call initialize_collocation_grid(grid, COLLOCATION_GRID_LINEAR, &
        first_nodes, second_nodes, first_weights, second_weights, endpoint_map, &
        axis_map, 5, status)
    call record(status%code == 0 .and. validate_collocation_grid(grid, status), &
        "linear collocation grid accepts positive tensor-product metadata")
    call record(collocation_grid_point_count(grid) == n_point, &
        "grid point count is the tensor-product cardinality")

    call collocation_grid_metadata(grid, grid_kind, copied_first, copied_second, &
        copied_weights, copied_endpoint, copied_axis, chunk_size, copied_start, &
        copied_end, status)
    call record(status%code == 0 .and. grid_kind == COLLOCATION_GRID_LINEAR .and. &
        chunk_size == 5 .and. all(copied_first == first_nodes) .and. &
        all(copied_second == second_nodes) .and. all(copied_endpoint == endpoint_map) .and. &
        all(copied_axis == axis_map), &
        "grid metadata exports independent coordinate and endpoint copies")
    call record(all(copied_weights > 0.0_dp), &
        "tensor-product quadrature weights are strictly positive")

    do second = 1, n_second
        do first = 1, n_first
            flat = collocation_grid_flat_index(grid, first, second, status)
            call collocation_grid_unflatten_index(grid, flat, mapped_first, mapped_second, status)
            call record(flat == first + (second - 1)*n_first .and. &
                mapped_first == first .and. mapped_second == second, &
                "flatten and unflatten maps are deterministic and inverse")
            expected_weight = first_weights(first)*second_weights(second)
            call record(abs(copied_weights(flat) - expected_weight) < 1.0e-14_dp, &
                "flattened weight agrees with the independent tensor-product oracle")
        end do
    end do
    do chunk = 1, size(copied_start)
        call collocation_grid_chunk_bounds(grid, chunk, first_point, last_point, status)
        call record(status%code == 0 .and. first_point == copied_start(chunk) .and. &
            last_point == copied_end(chunk) .and. last_point - first_point + 1 <= chunk_size, &
            "grid chunk metadata is contiguous and bounded")
    end do

    call initialize_collocation_grid(quadrature, COLLOCATION_GRID_QUADRATURE, &
        first_nodes, second_nodes, first_weights, second_weights, chunk_size=4, status=status)
    call record(status%code == 0, "quadrature grid kind shares the neutral contract")
    call initialize_collocation_grid(concentric, COLLOCATION_GRID_CONCENTRIC, &
        [0.0_dp, 0.2_dp, 0.8_dp], second_nodes, first_weights, second_weights, &
        chunk_size=4, status=status)
    call record(status%code == 0, "concentric grid kind accepts nonnegative radii")

    call initialize_collocation_grid(bad_grid, COLLOCATION_GRID_LINEAR, first_nodes, &
        second_nodes, first_weights(1:2), second_weights, endpoint_map, axis_map, 5, status)
    call record(status%code /= 0, &
        "grid initialization rejects incompatible weight shapes")
    call initialize_collocation_grid(bad_grid, COLLOCATION_GRID_LINEAR, first_nodes, &
        second_nodes, [-1.0_dp, 0.5_dp, 0.25_dp], second_weights, endpoint_map, axis_map, &
        5, status)
    call record(status%code /= 0, "grid initialization rejects nonpositive weights")
    call initialize_collocation_grid(bad_grid, COLLOCATION_GRID_LINEAR, first_nodes, &
        second_nodes, first_weights, second_weights, [1, 2, 4, 1], axis_map, 5, status)
    call record(status%code /= 0, "grid initialization rejects a non-idempotent endpoint map")
    call initialize_collocation_grid(bad_grid, COLLOCATION_GRID_LINEAR, first_nodes, &
        second_nodes, first_weights, second_weights, endpoint_map, axis_map, n_point + 1, status)
    call record(status%code /= 0, "grid initialization bounds chunk size by point count")

    call initialize_direct_fourier_plan(plan, angles, modes, 3, status=status)
    call record(status%code == 0 .and. validate_direct_fourier_plan(plan, status) .and. &
        direct_fourier_plan_sample_count(plan) == 4 .and. &
        direct_fourier_plan_mode_count(plan) == 4, &
        "uniform square Fourier plan accepts bounded direct-transform metadata")
    call direct_fourier_forward(plan, samples, coefficients, status)
    call record(status%code == 0, "direct Fourier forward action accepts a finite sample vector")
    do first = 1, 4
        expected_coefficient = cmplx(0.0_dp, 0.0_dp, dp)
        do second = 1, 4
            expected_coefficient = expected_coefficient + samples(second)*cmplx( &
                cos(real(modes(first), dp)*angles(second)), &
                -sin(real(modes(first), dp)*angles(second)), dp)/2.0_dp
        end do
        call record(abs(coefficients(first) - expected_coefficient) < 1.0e-13_dp, &
            "forward action agrees with an independently evaluated Fourier sum")
    end do
    call direct_fourier_adjoint(plan, coefficients, recovered, status)
    round_trip_error = maxval(abs(recovered - samples))
    call record(status%code == 0 .and. round_trip_error < 1.0e-13_dp, &
        "normalized uniform direct Fourier forward/adjoint round trip is deterministic")

    call direct_fourier_plan_chunk_bounds(plan, 1, first_point, last_point, status)
    call record(status%code == 0 .and. first_point == 1 .and. last_point == 3, &
        "Fourier chunk metadata reports the requested bounded first chunk")
    call direct_fourier_plan_metadata(plan, copied_first, copied_axis, copied_weights, &
        chunk_size, copied_start, copied_end, status)
    call record(status%code == 0 .and. chunk_size == 3 .and. all(copied_start == [1, 4]) .and. &
        all(copied_end == [3, 4]), "Fourier plan metadata exports contiguous chunk bounds")
    call direct_fourier_forward(plan, samples(1:3), coefficients, status)
    call record(status%code /= 0, "Fourier forward rejects an invalid sample shape")
    call initialize_direct_fourier_plan(bad_plan, angles, modes, 5, status=status)
    call record(status%code /= 0, "Fourier plan rejects a chunk larger than the sample set")

    call check_summary("collocation and direct Fourier contract")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_collocation_transform_contract
