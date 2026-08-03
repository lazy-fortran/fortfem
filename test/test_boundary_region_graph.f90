program test_boundary_region_graph
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        boundary_region_graph_components, boundary_region_graph_cycle_basis, &
        boundary_region_graph_incidence, boundary_region_graph_interface_samples, &
        boundary_region_graph_interface_metadata, &
        boundary_region_graph_t, initialize_boundary_region_graph, &
        validate_boundary_region_graph
    use fortfem_kinds, only: dp
    implicit none

    type(boundary_region_graph_t) :: graph
    integer :: plus_region(2), minus_region(2), interface_genus(2)
    logical :: exterior_interface(2)
    integer :: sample_offsets(3), incidence(3, 2), cycle_count, component_count, status
    integer, allocatable :: incidence_alloc(:, :), components(:), cycle_basis(:, :)
    real(dp) :: sample_points(3, 3), sample_normals(3, 3), sample_weights(3)
    real(dp), allocatable :: points(:, :), normals(:, :), weights(:)
    integer, allocatable :: genus_copy(:)
    logical, allocatable :: exterior_copy(:)
    logical :: all_passed

    all_passed = .true.
    plus_region = [1, 2]
    minus_region = [2, 3]
    interface_genus = [0, 1]
    exterior_interface = [.false., .true.]
    sample_offsets = [1, 3, 4]
    sample_points = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
    sample_normals = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
    sample_weights = [0.5_dp, 0.75_dp, 1.25_dp]

    call initialize_boundary_region_graph( &
        graph, 3, plus_region, minus_region, interface_genus, exterior_interface, &
        sample_offsets, sample_points, sample_normals, sample_weights, status)
    call record_condition(status == 0, &
        "boundary region graph accepts oriented interface samples")
    call validate_boundary_region_graph(graph, status)
    call record_condition(status == 0, "boundary region graph validates fixed metadata")
    call boundary_region_graph_incidence(graph, incidence_alloc, status)
    incidence = reshape([1, -1, 0, 0, 1, -1], [3, 2])
    call record_condition(status == 0 .and. all(incidence_alloc == incidence), &
        "boundary region graph preserves oriented region incidence")
    call boundary_region_graph_components(graph, components, component_count, status)
    call record_condition(status == 0 .and. component_count == 1 .and. &
        all(components == 1), "boundary region graph delegates connectivity")
    call boundary_region_graph_cycle_basis(graph, cycle_basis, cycle_count, status)
    call record_condition(status == 0 .and. cycle_count == 0 .and. &
        size(cycle_basis, 1) == 2, "boundary region graph exposes cycle metadata")
    call boundary_region_graph_interface_samples(graph, 1, points, normals, weights, status)
    call record_condition(status == 0 .and. size(points, 2) == 2 .and. &
        maxval(abs(points - sample_points(:, 1:2))) < 2.0e-14_dp .and. &
        maxval(abs(normals - sample_normals(:, 1:2))) < 2.0e-14_dp .and. &
        maxval(abs(weights - sample_weights(1:2))) < 2.0e-14_dp, &
        "boundary region graph extracts per-interface sample ranges")
    call boundary_region_graph_interface_samples(graph, 2, points, normals, weights, status)
    call record_condition(status == 0 .and. size(points, 2) == 1 .and. &
        maxval(abs(points(:, 1) - sample_points(:, 3))) < 2.0e-14_dp, &
        "boundary region graph extracts exterior interface samples")
    call boundary_region_graph_interface_metadata( &
        graph, genus_copy, exterior_copy, status)
    call record_condition(status == 0 .and. all(genus_copy == interface_genus) .and. &
        all(exterior_copy .eqv. exterior_interface), &
        "boundary region graph exposes interface topology metadata")

    sample_offsets = [1, 4, 5]
    call initialize_boundary_region_graph( &
        graph, 3, plus_region, minus_region, interface_genus, exterior_interface, &
        sample_offsets, sample_points, sample_normals, sample_weights, status)
    call record_condition(status /= 0, "boundary region graph rejects inconsistent offsets")

    call check_summary("boundary region graph")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_boundary_region_graph
