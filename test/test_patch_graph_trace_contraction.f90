program test_patch_graph_trace_contraction
    use check, only: check_condition, check_summary
    use fortfem_boundary_region_graph, only: &
        boundary_region_graph_cycle_basis, boundary_region_graph_incidence, &
        boundary_region_graph_interface_samples, boundary_region_graph_t, &
        initialize_boundary_region_graph
    use fortfem_kinds, only: dp
    use fortfem_patch_graph_trace_contraction, only: &
        assemble_patch_graph_trace_contraction, &
        assemble_patch_graph_trace_contraction_jvp, &
        assemble_patch_graph_trace_contraction_vjp
    implicit none

    type(boundary_region_graph_t) :: graph
    integer :: plus_region(3), minus_region(3), interface_genus(3)
    logical :: exterior_interface(3)
    integer :: sample_offsets(4), status, cycle_count, edge_count
    integer, allocatable :: cycle_basis(:, :), incidence(:, :)
    real(dp) :: sample_points(3, 5), sample_normals(3, 5), sample_weights(5)
    real(dp) :: trace_values(2, 5), trace_dot(2, 5), trace_bar(2, 2)
    real(dp), allocatable :: contraction(:, :), contraction_dot(:, :)
    real(dp), allocatable :: contraction_plus(:, :), contraction_minus(:, :)
    real(dp), allocatable :: trace_values_bar(:, :), expected(:, :)
    real(dp), allocatable :: points(:, :), normals(:, :), weights(:)
    real(dp) :: epsilon, lhs, rhs
    integer :: interface_id, cycle_id, component, offset, sample_count
    logical :: all_passed

    all_passed = .true.
    plus_region = [1, 2, 1]
    minus_region = [2, 1, 1]
    interface_genus = [0, 0, 1]
    exterior_interface = [.false., .false., .false.]
    sample_offsets = [1, 3, 4, 6]
    sample_points = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 1.0_dp], [3, 5])
    sample_normals = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp], [3, 5])
    sample_weights = [0.5_dp, 0.75_dp, 1.25_dp, 0.8_dp, 1.1_dp]
    trace_values = reshape([ &
        1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp, -1.5_dp, 0.25_dp, &
        2.0_dp, 1.0_dp, -0.75_dp, 0.5_dp], [2, 5])
    trace_dot = reshape([ &
        -0.4_dp, 0.9_dp, 0.25_dp, -0.3_dp, 0.8_dp, 0.2_dp, &
        0.7_dp, -0.5_dp, 0.6_dp, -0.1_dp], [2, 5])
    trace_bar = reshape([1.2_dp, -0.3_dp, 0.8_dp, 0.4_dp], [2, 2])

    call initialize_boundary_region_graph( &
        graph, 2, plus_region, minus_region, interface_genus, &
        exterior_interface, sample_offsets, sample_points, sample_normals, &
        sample_weights, status)
    call record_condition(status == 0, &
        "patch graph with periodic self-identification initializes")
    call boundary_region_graph_cycle_basis( &
        graph, cycle_basis, cycle_count, status)
    edge_count = size(cycle_basis, 1)
    call boundary_region_graph_incidence(graph, incidence, status)
    call record_condition(status == 0 .and. edge_count == 3 .and. &
        cycle_count == 2 .and. maxval(abs(matmul( &
        real(incidence, dp), real(cycle_basis, dp)))) < 2.0e-14_dp, &
        "cycle basis closes the oriented arbitrary patch graph")

    call assemble_patch_graph_trace_contraction( &
        graph, trace_values, contraction, status)
    call record_condition(status == 0 .and. size(contraction, 1) == 2 .and. &
        size(contraction, 2) == cycle_count, &
        "trace contraction returns one value per component and cycle")
    call manual_contraction( &
        graph, cycle_basis, trace_values, expected, status)
    call record_condition(status == 0 .and. maxval(abs(contraction - expected)) < &
        2.0e-14_dp, "trace contraction matches independent weighted ledger")

    call assemble_patch_graph_trace_contraction_jvp( &
        graph, trace_dot, contraction_dot, status)
    call manual_contraction( &
        graph, cycle_basis, trace_dot, expected, status)
    call record_condition(status == 0 .and. maxval(abs(contraction_dot - expected)) < &
        2.0e-14_dp, "trace contraction JVP matches independent ledger")

    epsilon = 3.0e-7_dp
    call assemble_patch_graph_trace_contraction( &
        graph, trace_values + epsilon * trace_dot, contraction_plus, status)
    call assemble_patch_graph_trace_contraction( &
        graph, trace_values - epsilon * trace_dot, contraction_minus, status)
    call record_condition(maxval(abs( &
        (contraction_plus - contraction_minus) / (2.0_dp * epsilon) - &
        contraction_dot)) < 2.0e-9_dp, &
        "trace contraction JVP agrees with a centered directional check")

    call assemble_patch_graph_trace_contraction_vjp( &
        graph, trace_bar, trace_values_bar, status)
    rhs = sum(trace_values_bar * trace_dot)
    call assemble_patch_graph_trace_contraction_jvp( &
        graph, trace_dot, contraction_dot, status)
    lhs = sum(trace_bar * contraction_dot)
    call record_condition(status == 0 .and. abs(lhs - rhs) < 2.0e-14_dp, &
        "trace contraction VJP satisfies the dot-product oracle")

    call assemble_patch_graph_trace_contraction( &
        graph, trace_values(:, :4), contraction, status)
    call record_condition(status /= 0 .and. size(contraction, 2) == 0, &
        "trace contraction rejects a mismatched global sample ledger")

    call check_summary("patch graph trace contraction")
    if (.not. all_passed) error stop 1

contains

    subroutine manual_contraction(graph, cycle_basis, values, result, status)
        type(boundary_region_graph_t), intent(in) :: graph
        integer, intent(in) :: cycle_basis(:, :)
        real(dp), intent(in) :: values(:, :)
        real(dp), allocatable, intent(out) :: result(:, :)
        integer, intent(out) :: status
        integer :: interface_id, cycle_id, component, offset, sample_count
        real(dp), allocatable :: points(:, :), normals(:, :), weights(:)

        status = 1
        allocate(result(size(values, 1), size(cycle_basis, 2)))
        result = 0.0_dp
        offset = 0
        do interface_id = 1, size(cycle_basis, 1)
            call boundary_region_graph_interface_samples( &
                graph, interface_id, points, normals, weights, status)
            if (status /= 0) return
            sample_count = size(weights)
            if (offset + sample_count > size(values, 2)) return
            do cycle_id = 1, size(cycle_basis, 2)
                do component = 1, size(values, 1)
                    result(component, cycle_id) = result(component, cycle_id) + &
                        real(cycle_basis(interface_id, cycle_id), dp) * &
                        sum(weights * values( &
                        component, offset + 1:offset + sample_count))
                end do
            end do
            offset = offset + sample_count
        end do
        status = merge(0, 1, offset == size(values, 2))
    end subroutine manual_contraction

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_patch_graph_trace_contraction
