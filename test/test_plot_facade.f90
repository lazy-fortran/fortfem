program test_plot_facade
    !! Independent smoke test for the canonical plotting/sampling facade.
    !!
    !! The oracle is the affine field u=x+y.  P1 interpolation must reproduce
    !! that field at every point in the sampled unit-square grid exactly up to
    !! roundoff.  The mesh-path contract is checked independently as well;
    !! plotting itself is intentionally not called, so this test writes no
    !! media files.
    use check, only: check_condition, check_summary
    use fortfem_core, only: mesh_t, unit_square_mesh, function_space_t
    use fortfem_feec, only: function, function_space, function_t
    use fortfem_kinds, only: dp
    use fortfem_plot, only: build_mesh_edge_path, compute_scalar_plot_grid
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none

    integer, parameter :: nx = 3, ny = 3
    real(dp), parameter :: tolerance = 128.0_dp*epsilon(1.0_dp)
    type(mesh_t) :: mesh
    type(function_space_t) :: space
    type(function_t) :: field
    real(dp) :: x_grid(nx + 1), y_grid(ny + 1), z_grid(nx, ny)
    real(dp), allocatable :: x_path(:), y_path(:)
    real(dp) :: expected, error
    integer :: i, j, finite_count
    logical :: all_passed

    all_passed = .true.
    mesh = unit_square_mesh(2)
    space = function_space(mesh, "Lagrange", 1)
    field = function(space)
    field%values = mesh%data%vertices(1, :) + mesh%data%vertices(2, :)

    call compute_scalar_plot_grid(field, nx, ny, x_grid, y_grid, z_grid)
    call record_condition(abs(x_grid(1)) <= tolerance .and. &
        abs(x_grid(nx + 1) - 1.0_dp) <= tolerance, &
        "scalar plot grid spans the mesh in x")
    call record_condition(abs(y_grid(1)) <= tolerance .and. &
        abs(y_grid(ny + 1) - 1.0_dp) <= tolerance, &
        "scalar plot grid spans the mesh in y")

    error = 0.0_dp
    do i = 1, nx
        do j = 1, ny
            expected = x_grid(i) + y_grid(j)
            error = max(error, abs(z_grid(i, j) - expected))
        end do
    end do
    call record_condition(error <= tolerance, &
        "scalar plot sampling reproduces an affine P1 field")

    call build_mesh_edge_path(mesh%data%vertices, mesh%data%edges, x_path, &
        y_path)
    finite_count = count(.not. ieee_is_nan(x_path) .and. &
        .not. ieee_is_nan(y_path))
    call record_condition(size(x_path) == 3*mesh%data%n_edges .and. &
        size(y_path) == size(x_path), &
        "mesh plot path reserves one segment per edge")
    call record_condition(finite_count == 2*mesh%data%n_edges, &
        "mesh plot path separates edges with NaN sentinels")

    call check_summary("Canonical plot facade")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_plot_facade
