module fortfem_patch_graph_trace_contraction
    !! Weighted cycle contractions of traces on an oriented patch graph.
    !!
    !! Interfaces are stored as contiguous sample blocks by
    !! `boundary_region_graph_t`.  This module contracts a component-wise
    !! trace against each graph cycle, preserving interface orientation from
    !! the integer cycle basis.  The operation is deliberately agnostic to
    !! the field type: a component can be a scalar, vector, or modal value.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_boundary_region_graph, only: &
        boundary_region_graph_cycle_basis, boundary_region_graph_interface_samples, &
        boundary_region_graph_t, validate_boundary_region_graph
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: assemble_patch_graph_trace_contraction
    public :: assemble_patch_graph_trace_contraction_jvp
    public :: assemble_patch_graph_trace_contraction_vjp

contains

    subroutine assemble_patch_graph_trace_contraction( &
            graph, trace_values, cycle_contraction, status)
        type(boundary_region_graph_t), intent(in) :: graph
        real(dp), intent(in) :: trace_values(:, :)
        real(dp), allocatable, intent(out) :: cycle_contraction(:, :)
        integer, intent(out) :: status
        integer, allocatable :: cycle_basis(:,:), edge_offsets(:)
        real(dp), allocatable :: all_weights(:)
        integer :: component_count, cycle_count, sample_count

        if (allocated(cycle_contraction)) deallocate(cycle_contraction)
        call prepare_ledger( &
            graph, cycle_basis, edge_offsets, all_weights, status)
        if (status /= 0) then
            allocate(cycle_contraction(0, 0))
            return
        end if
        component_count = size(trace_values, 1)
        cycle_count = size(cycle_basis, 2)
        sample_count = size(all_weights)
        if (.not. valid_trace_shape(trace_values, component_count, sample_count)) then
            allocate(cycle_contraction(component_count, 0))
            status = 2
            return
        end if
        allocate(cycle_contraction(component_count, cycle_count))
        call apply_contraction( &
            cycle_basis, edge_offsets, all_weights, trace_values, &
            cycle_contraction)
        status = 0
    end subroutine assemble_patch_graph_trace_contraction

    subroutine assemble_patch_graph_trace_contraction_jvp( &
            graph, trace_values_dot, cycle_contraction_dot, status)
        type(boundary_region_graph_t), intent(in) :: graph
        real(dp), intent(in) :: trace_values_dot(:, :)
        real(dp), allocatable, intent(out) :: cycle_contraction_dot(:, :)
        integer, intent(out) :: status
        integer, allocatable :: cycle_basis(:,:), edge_offsets(:)
        real(dp), allocatable :: all_weights(:)
        integer :: component_count, cycle_count, sample_count

        if (allocated(cycle_contraction_dot)) deallocate(cycle_contraction_dot)
        call prepare_ledger( &
            graph, cycle_basis, edge_offsets, all_weights, status)
        if (status /= 0) then
            allocate(cycle_contraction_dot(0, 0))
            return
        end if
        component_count = size(trace_values_dot, 1)
        cycle_count = size(cycle_basis, 2)
        sample_count = size(all_weights)
        if (.not. valid_trace_shape( &
                trace_values_dot, component_count, sample_count)) then
            allocate(cycle_contraction_dot(component_count, 0))
            status = 2
            return
        end if
        allocate(cycle_contraction_dot(component_count, cycle_count))
        call apply_contraction( &
            cycle_basis, edge_offsets, all_weights, trace_values_dot, &
            cycle_contraction_dot)
        status = 0
    end subroutine assemble_patch_graph_trace_contraction_jvp

    subroutine assemble_patch_graph_trace_contraction_vjp( &
            graph, cycle_contraction_bar, trace_values_bar, status)
        type(boundary_region_graph_t), intent(in) :: graph
        real(dp), intent(in) :: cycle_contraction_bar(:, :)
        real(dp), allocatable, intent(out) :: trace_values_bar(:, :)
        integer, intent(out) :: status
        integer, allocatable :: cycle_basis(:,:), edge_offsets(:)
        real(dp), allocatable :: all_weights(:)
        integer :: component_count, cycle_count, sample_count

        if (allocated(trace_values_bar)) deallocate(trace_values_bar)
        call prepare_ledger( &
            graph, cycle_basis, edge_offsets, all_weights, status)
        if (status /= 0) then
            allocate(trace_values_bar(0, 0))
            return
        end if
        component_count = size(cycle_contraction_bar, 1)
        cycle_count = size(cycle_basis, 2)
        sample_count = size(all_weights)
        if (.not. valid_cycle_shape( &
                cycle_contraction_bar, component_count, cycle_count)) then
            allocate(trace_values_bar(0, 0))
            status = 2
            return
        end if
        allocate(trace_values_bar(component_count, sample_count))
        call apply_contraction_vjp( &
            cycle_basis, edge_offsets, all_weights, cycle_contraction_bar, &
            trace_values_bar)
        status = 0
    end subroutine assemble_patch_graph_trace_contraction_vjp

    subroutine prepare_ledger( &
            graph, cycle_basis, edge_offsets, all_weights, status)
        type(boundary_region_graph_t), intent(in) :: graph
        integer, allocatable, intent(out) :: cycle_basis(:, :)
        integer, allocatable, intent(out) :: edge_offsets(:)
        real(dp), allocatable, intent(out) :: all_weights(:)
        integer, intent(out) :: status
        real(dp), allocatable :: points(:, :), normals(:, :), weights(:)
        integer :: edge_count, edge, sample_count, cycle_count, total_samples

        if (allocated(cycle_basis)) deallocate(cycle_basis)
        if (allocated(edge_offsets)) deallocate(edge_offsets)
        if (allocated(all_weights)) deallocate(all_weights)
        call validate_boundary_region_graph(graph, status)
        if (status /= 0) then
            allocate(cycle_basis(0, 0), edge_offsets(0), all_weights(0))
            return
        end if
        call boundary_region_graph_cycle_basis( &
            graph, cycle_basis, cycle_count, status)
        if (status /= 0) then
            allocate(edge_offsets(0), all_weights(0))
            return
        end if
        edge_count = size(cycle_basis, 1)
        allocate(edge_offsets(edge_count + 1))
        edge_offsets(1) = 1
        total_samples = 0
        do edge = 1, edge_count
            call boundary_region_graph_interface_samples( &
                graph, edge, points, normals, weights, status)
            if (status /= 0) then
                deallocate(edge_offsets)
                allocate(edge_offsets(0), all_weights(0))
                return
            end if
            sample_count = size(weights)
            total_samples = total_samples + sample_count
            edge_offsets(edge + 1) = total_samples + 1
        end do
        allocate(all_weights(total_samples))
        total_samples = 0
        do edge = 1, edge_count
            call boundary_region_graph_interface_samples( &
                graph, edge, points, normals, weights, status)
            if (status /= 0) then
                deallocate(edge_offsets, all_weights)
                allocate(edge_offsets(0), all_weights(0))
                return
            end if
            sample_count = size(weights)
            if (sample_count > 0) all_weights( &
                total_samples + 1:total_samples + sample_count) = weights
            total_samples = total_samples + sample_count
        end do
        status = 0
    end subroutine prepare_ledger

    subroutine apply_contraction( &
            cycle_basis, edge_offsets, all_weights, values, contraction)
        integer, intent(in) :: cycle_basis(:, :), edge_offsets(:)
        real(dp), intent(in) :: all_weights(:), values(:, :)
        real(dp), intent(out) :: contraction(:, :)
        integer :: edge, cycle, component, first_sample, last_sample

        contraction = 0.0_dp
        do edge = 1, size(cycle_basis, 1)
            first_sample = edge_offsets(edge)
            last_sample = edge_offsets(edge + 1) - 1
            if (last_sample < first_sample) cycle
            do cycle = 1, size(cycle_basis, 2)
                do component = 1, size(values, 1)
                    contraction(component, cycle) = &
                        contraction(component, cycle) + &
                        real(cycle_basis(edge, cycle), dp) * sum( &
                        all_weights(first_sample:last_sample) * &
                        values(component, first_sample:last_sample))
                end do
            end do
        end do
    end subroutine apply_contraction

    subroutine apply_contraction_vjp( &
            cycle_basis, edge_offsets, all_weights, contraction_bar, values_bar)
        integer, intent(in) :: cycle_basis(:, :), edge_offsets(:)
        real(dp), intent(in) :: all_weights(:), contraction_bar(:, :)
        real(dp), intent(out) :: values_bar(:, :)
        integer :: edge, cycle, component, sample, first_sample, last_sample

        values_bar = 0.0_dp
        do edge = 1, size(cycle_basis, 1)
            first_sample = edge_offsets(edge)
            last_sample = edge_offsets(edge + 1) - 1
            if (last_sample < first_sample) cycle
            do sample = first_sample, last_sample
                do cycle = 1, size(cycle_basis, 2)
                    do component = 1, size(contraction_bar, 1)
                        values_bar(component, sample) = &
                            values_bar(component, sample) + &
                            real(cycle_basis(edge, cycle), dp) * &
                            all_weights(sample) * contraction_bar(component, cycle)
                    end do
                end do
            end do
        end do
    end subroutine apply_contraction_vjp

    logical function valid_trace_shape(values, component_count, sample_count) &
            result(valid)
        real(dp), intent(in) :: values(:, :)
        integer, intent(in) :: component_count, sample_count

        valid = component_count > 0 .and. size(values, 2) == sample_count .and. &
            all(ieee_is_finite(values))
    end function valid_trace_shape

    logical function valid_cycle_shape(values, component_count, cycle_count) &
            result(valid)
        real(dp), intent(in) :: values(:, :)
        integer, intent(in) :: component_count, cycle_count

        valid = component_count > 0 .and. size(values, 2) == cycle_count .and. &
            all(ieee_is_finite(values))
    end function valid_cycle_shape

end module fortfem_patch_graph_trace_contraction
