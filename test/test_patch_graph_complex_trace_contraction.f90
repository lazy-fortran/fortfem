program test_patch_graph_complex_trace_contraction
    use check, only: check_condition, check_summary
    use fortfem_boundary_region_graph, only: &
        boundary_region_graph_cycle_basis, &
        boundary_region_graph_interface_samples, boundary_region_graph_t, &
        initialize_boundary_region_graph
    use fortfem_kinds, only: dp
    use fortfem_patch_graph_complex_trace_contraction, only: &
        assemble_patch_graph_complex_trace_contraction, &
        assemble_patch_graph_complex_trace_contraction_jvp, &
        assemble_patch_graph_complex_trace_contraction_vjp
    implicit none

    type(boundary_region_graph_t) :: graph, malformed_graph
    integer :: plus_region(3), minus_region(3), interface_genus(3)
    logical :: exterior_interface(3)
    integer :: sample_offsets(4), status, cycle_count
    integer, allocatable :: cycle_basis(:, :)
    real(dp) :: sample_points(3, 5), sample_normals(3, 5)
    real(dp) :: sample_weights(5)
    complex(dp) :: trace_values(2, 5), trace_dot(2, 5), trace_bar(2, 2)
    complex(dp), allocatable :: contraction(:, :), contraction_dot(:, :)
    complex(dp), allocatable :: contraction_minus(:, :)
    complex(dp), allocatable :: contraction_plus(:, :), expected(:, :)
    complex(dp), allocatable :: trace_values_bar(:, :)
    real(dp) :: epsilon, lhs, rhs
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
        cmplx(1.0_dp, 0.2_dp, dp), cmplx(-2.0_dp, 0.7_dp, dp), &
        cmplx(0.5_dp, -0.4_dp, dp), cmplx(3.0_dp, 0.1_dp, dp), &
        cmplx(-1.5_dp, 0.6_dp, dp), cmplx(0.25_dp, -0.9_dp, dp), &
        cmplx(2.0_dp, 0.3_dp, dp), cmplx(1.0_dp, -0.8_dp, dp), &
        cmplx(-0.75_dp, 0.5_dp, dp), cmplx(0.5_dp, -0.2_dp, dp)], [2, 5])
    trace_dot = reshape([ &
        cmplx(-0.4_dp, 0.6_dp, dp), cmplx(0.9_dp, -0.2_dp, dp), &
        cmplx(0.25_dp, 0.3_dp, dp), cmplx(-0.3_dp, 0.8_dp, dp), &
        cmplx(0.8_dp, -0.5_dp, dp), cmplx(0.2_dp, 0.7_dp, dp), &
        cmplx(0.7_dp, -0.1_dp, dp), cmplx(-0.5_dp, 0.4_dp, dp), &
        cmplx(0.6_dp, 0.9_dp, dp), cmplx(-0.1_dp, -0.6_dp, dp)], [2, 5])
    trace_bar = reshape([ &
        cmplx(1.2_dp, -0.5_dp, dp), cmplx(-0.3_dp, 0.9_dp, dp), &
        cmplx(0.8_dp, 0.2_dp, dp), cmplx(0.4_dp, -0.7_dp, dp)], [2, 2])

    call initialize_boundary_region_graph( &
        graph, 2, plus_region, minus_region, interface_genus, &
        exterior_interface, sample_offsets, sample_points, sample_normals, &
        sample_weights, status)
    call record_condition(status == 0, &
        "complex contraction patch graph initializes")
    call boundary_region_graph_cycle_basis( &
        graph, cycle_basis, cycle_count, status)

    call assemble_patch_graph_complex_trace_contraction( &
        graph, trace_values, contraction, status)
    call record_condition(status == 0, &
        "complex patch-graph trace contraction succeeds")
    call manual_contraction( &
        graph, cycle_basis, trace_values, expected, status)
    call record_condition(status == 0, &
        "independent complex cycle-weight ledger succeeds")
    call record_condition(size(contraction, 2) == cycle_count, &
        "complex contraction returns one column per independent cycle")
    call record_condition(maxval(abs(contraction - expected)) < 2.0e-14_dp, &
        "complex contraction matches an independent cycle-weight ledger")

    call assemble_patch_graph_complex_trace_contraction_jvp( &
        graph, trace_dot, contraction_dot, status)
    call record_condition(status == 0, &
        "complex patch-graph trace contraction JVP succeeds")
    call manual_contraction(graph, cycle_basis, trace_dot, expected, status)
    call record_condition(status == 0 .and. &
        maxval(abs(contraction_dot - expected)) < 2.0e-14_dp, &
        "complex contraction JVP matches the independent ledger")

    epsilon = 3.0e-7_dp
    call assemble_patch_graph_complex_trace_contraction( &
        graph, trace_values + epsilon*trace_dot, contraction_plus, status)
    call assemble_patch_graph_complex_trace_contraction( &
        graph, trace_values - epsilon*trace_dot, contraction_minus, status)
    call record_condition(maxval(abs( &
        (contraction_plus - contraction_minus)/(2.0_dp*epsilon) - &
        contraction_dot)) < 2.0e-9_dp, &
        "complex contraction JVP agrees with a centered directional oracle")

    call assemble_patch_graph_complex_trace_contraction_vjp( &
        graph, trace_bar, trace_values_bar, status)
    lhs = real(sum(conjg(trace_bar)*contraction_dot), dp)
    rhs = real(sum(conjg(trace_values_bar)*trace_dot), dp)
    call record_condition(status == 0 .and. abs(lhs - rhs) < 2.0e-14_dp, &
        "complex contraction VJP satisfies the real-part adjoint oracle")

    call assemble_patch_graph_complex_trace_contraction( &
        graph, trace_values(:, :4), contraction, status)
    call record_condition(status /= 0 .and. size(contraction, 2) == 0, &
        "complex contraction rejects a mismatched sample ledger")
    call assemble_patch_graph_complex_trace_contraction( &
        malformed_graph, trace_values, contraction, status)
    call record_condition(status /= 0 .and. size(contraction) == 0, &
        "complex contraction rejects malformed graph metadata")

    call check_summary("complex patch graph trace contraction")
    if (.not. all_passed) error stop 1

contains

    subroutine manual_contraction(graph, cycle_basis, values, result, status)
        type(boundary_region_graph_t), intent(in) :: graph
        integer, intent(in) :: cycle_basis(:, :)
        complex(dp), intent(in) :: values(:, :)
        complex(dp), allocatable, intent(out) :: result(:, :)
        integer, intent(out) :: status
        integer :: interface_id, cycle_id, component, offset, sample_count
        real(dp), allocatable :: points(:, :), normals(:, :), weights(:)

        status = 1
        allocate(result(size(values, 1), size(cycle_basis, 2)))
        result = cmplx(0.0_dp, 0.0_dp, dp)
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
                        real(cycle_basis(interface_id, cycle_id), dp)* &
                        sum(weights*values( &
                        component, offset + 1:offset + sample_count))
                end do
            end do
            offset = offset + sample_count
        end do
        if (offset /= size(values, 2)) return
        status = 0
    end subroutine manual_contraction

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_patch_graph_complex_trace_contraction
