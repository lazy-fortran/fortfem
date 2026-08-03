program test_fourier_collocation_facade
    !! Downstream-style behavioral oracle for collocation metadata exported by
    !! the canonical Fourier facade.  No compatibility umbrella is imported.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: &
        COLLOCATION_GRID_CONCENTRIC, COLLOCATION_GRID_LINEAR, &
        COLLOCATION_GRID_QUADRATURE, collocation_grid_chunk_bounds, &
        collocation_grid_flat_index, collocation_grid_metadata, &
        collocation_grid_point_count, collocation_grid_t, &
        collocation_grid_unflatten_index, initialize_collocation_grid, &
        validate_collocation_grid
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    real(dp), parameter :: first_nodes(2) = [-1.0_dp, 1.0_dp]
    real(dp), parameter :: second_nodes(3) = [0.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: first_weights(2) = [2.0_dp, 3.0_dp]
    real(dp), parameter :: second_weights(3) = [5.0_dp, 7.0_dp, 11.0_dp]
    type(collocation_grid_t) :: grid
    type(fortsparse_status_t) :: status
    real(dp), allocatable :: copied_first(:), copied_second(:), weights(:)
    integer, allocatable :: endpoint_map(:), axis_map(:), chunk_start(:), chunk_end(:)
    integer :: grid_kind, chunk_size, flat, first_index, second_index
    integer :: first_point, last_point

    call initialize_collocation_grid( &
        grid, COLLOCATION_GRID_LINEAR, first_nodes, second_nodes, &
        first_weights, second_weights, chunk_size=4, status=status)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        validate_collocation_grid(grid, status), &
        "canonical Fourier facade initializes valid collocation metadata")
    call check_condition(all([COLLOCATION_GRID_LINEAR, &
        COLLOCATION_GRID_QUADRATURE, COLLOCATION_GRID_CONCENTRIC] == [1, 2, 3]), &
        "canonical Fourier facade preserves the public grid-kind protocol")
    call check_condition(collocation_grid_point_count(grid) == 6, &
        "canonical Fourier facade reports tensor-product cardinality")

    flat = collocation_grid_flat_index(grid, 2, 3, status)
    call collocation_grid_unflatten_index( &
        grid, flat, first_index, second_index, status)
    call check_condition(flat == 6 .and. first_index == 2 .and. second_index == 3, &
        "facade flattening follows the independent column-major index oracle")

    call collocation_grid_chunk_bounds(grid, 2, first_point, last_point, status)
    call check_condition(first_point == 5 .and. last_point == 6, &
        "facade chunk bounds agree with the requested four-point partition")

    call collocation_grid_metadata( &
        grid, grid_kind, copied_first, copied_second, weights, endpoint_map, &
        axis_map, chunk_size, chunk_start, chunk_end, status)
    call check_condition(grid_kind == COLLOCATION_GRID_LINEAR .and. &
        chunk_size == 4 .and. all(weights == [10.0_dp, 15.0_dp, &
        14.0_dp, 21.0_dp, 22.0_dp, 33.0_dp]), &
        "facade metadata agrees with independent tensor-product weights")

    call check_summary("canonical Fourier collocation facade")

end program test_fourier_collocation_facade
