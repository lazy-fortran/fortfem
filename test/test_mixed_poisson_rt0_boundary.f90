program test_mixed_poisson_rt0_boundary
    use check, only: check_condition, check_summary
    use fortfem_core, only: mesh_t, rectangle_mesh
    use fortfem_feec, only: solve_mixed_poisson_rt0
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    type(fortsparse_status_t) :: status
    type(mesh_t) :: mesh
    real(dp), allocatable :: boundary_pressure(:), cell_pressure(:)
    real(dp), allocatable :: cell_source(:), flux_dofs(:)
    real(dp) :: centroid(2), edge_vector(2), expected_flux
    real(dp) :: flux_error, pressure_error
    integer :: degree_of_freedom, edge, triangle
    logical :: all_passed

    all_passed = .true.
    mesh = rectangle_mesh(5, 5, [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
    allocate(cell_source(mesh%data%n_triangles))
    allocate(boundary_pressure(mesh%data%n_edges))
    cell_source = 0.0_dp
    boundary_pressure = 0.0_dp
    do edge = 1, mesh%data%n_edges
        if (mesh%data%is_boundary_edge(edge)) then
            boundary_pressure(edge) = 0.5_dp * sum( &
                mesh%data%vertices(1, mesh%data%edges(:, edge)) + &
                mesh%data%vertices(2, mesh%data%edges(:, edge)))
        end if
    end do

    call solve_mixed_poisson_rt0( &
        mesh%data, cell_source, flux_dofs, cell_pressure, status, &
        boundary_pressure)
    call record_condition(status%code == 0, &
        "RT0-DG0 mixed solve accepts nonzero boundary pressure")

    flux_error = 0.0_dp
    do edge = 1, mesh%data%n_edges
        degree_of_freedom = mesh%data%edge_to_dof(edge) + 1
        edge_vector = mesh%data%vertices(:, mesh%data%edges(2, edge)) - &
            mesh%data%vertices(:, mesh%data%edges(1, edge))
        expected_flux = edge_vector(1) - edge_vector(2)
        flux_error = max(flux_error, abs( &
            flux_dofs(degree_of_freedom) - expected_flux))
    end do
    call record_condition(flux_error < 2.0e-13_dp, &
        "Affine pressure produces its exact constant Darcy flux")

    pressure_error = 0.0_dp
    do triangle = 1, mesh%data%n_triangles
        centroid = sum(mesh%data%vertices(:, &
            mesh%data%triangles(:, triangle)), dim=2) / 3.0_dp
        pressure_error = max(pressure_error, abs( &
            cell_pressure(triangle) - sum(centroid)))
    end do
    call record_condition(pressure_error < 2.0e-13_dp, &
        "Affine pressure produces its exact DG0 cell averages")

    call check_summary("RT0-DG0 mixed Poisson boundary data")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_mixed_poisson_rt0_boundary
