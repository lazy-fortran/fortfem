program test_mixed_poisson_rt0_conservation
    use check, only: check_condition, check_summary
    use fortfem_core, only: mesh_t, rectangle_mesh
    use fortfem_feec, only: solve_mixed_poisson_rt0
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    type(fortsparse_status_t) :: status
    type(mesh_t) :: mesh
    real(dp), allocatable :: cell_pressure(:), flux_dofs(:)
    real(dp) :: cell_source(2), determinant, local_balance
    real(dp) :: vertex_a(2), vertex_b(2), vertex_c(2)
    integer :: edge, edge_dofs(3), edge_orientations(3), triangle
    logical :: all_passed

    all_passed = .true.
    mesh = rectangle_mesh(2, 2, [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
    cell_source = [2.0_dp, -1.0_dp]
    call solve_mixed_poisson_rt0( &
        mesh%data, cell_source, flux_dofs, cell_pressure, status)

    call record_condition(status%code == 0, &
        "RT0-DG0 mixed Poisson system solves")
    call record_condition(allocated(flux_dofs), &
        "Mixed Poisson solve returns flux degrees of freedom")
    call record_condition(allocated(cell_pressure), &
        "Mixed Poisson solve returns cell pressures")
    if (allocated(flux_dofs) .and. allocated(cell_pressure)) then
        call record_condition(size(flux_dofs) == mesh%data%n_edges, &
            "Mixed Poisson flux uses one global moment per edge")
        call record_condition( &
            size(cell_pressure) == mesh%data%n_triangles, &
            "Mixed Poisson pressure uses one value per cell")
    end if

    do triangle = 1, mesh%data%n_triangles
        call mesh%data%get_triangle_edge_dofs( &
            triangle, edge_dofs, edge_orientations)
        local_balance = 0.0_dp
        do edge = 1, 3
            local_balance = local_balance + &
                real(edge_orientations(edge), dp) * &
                flux_dofs(edge_dofs(edge) + 1)
        end do
        vertex_a = mesh%data%vertices(:, mesh%data%triangles(1, triangle))
        vertex_b = mesh%data%vertices(:, mesh%data%triangles(2, triangle))
        vertex_c = mesh%data%vertices(:, mesh%data%triangles(3, triangle))
        determinant = &
            (vertex_b(1) - vertex_a(1)) * &
            (vertex_c(2) - vertex_a(2)) - &
            (vertex_b(2) - vertex_a(2)) * &
            (vertex_c(1) - vertex_a(1))
        call record_condition(abs( &
            local_balance - 0.5_dp * determinant * cell_source(triangle)) < &
            3.0e-13_dp, "Mixed flux conserves its cell source exactly")
    end do

    call check_summary("RT0-DG0 mixed Poisson conservation")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_mixed_poisson_rt0_conservation
