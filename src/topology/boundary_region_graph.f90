module fortfem_boundary_region_graph
    !! Neutral region/interface graph with caller-owned boundary samples.
    !!
    !! The topology delegates to `fortfem_region_interface_graph`.  This
    !! wrapper adds fixed interface metadata and a contiguous physical sampler:
    !! sample_offsets(i):sample_offsets(i+1)-1 belongs to interface i.  Points,
    !! normals, and weights remain application-owned geometry data; no mesh,
    !! coil, equilibrium, or file-format convention is selected here.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_region_interface_graph, only: &
        initialize_region_interface_graph, region_interface_graph_components, &
        region_interface_graph_cycle_basis, region_interface_graph_incidence, &
        region_interface_graph_t, validate_region_interface_graph
    implicit none
    private

    type, public :: boundary_region_graph_t
        private
        type(region_interface_graph_t) :: topology
        integer, allocatable :: interface_genus(:)
        logical, allocatable :: exterior_interface(:)
        integer, allocatable :: sample_offsets(:)
        real(dp), allocatable :: sample_points(:, :)
        real(dp), allocatable :: sample_normals(:, :)
        real(dp), allocatable :: sample_weights(:)
    end type boundary_region_graph_t

    public :: initialize_boundary_region_graph
    public :: validate_boundary_region_graph
    public :: boundary_region_graph_incidence
    public :: boundary_region_graph_components
    public :: boundary_region_graph_cycle_basis
    public :: boundary_region_graph_interface_samples
    public :: boundary_region_graph_interface_metadata

contains

    subroutine initialize_boundary_region_graph( &
            graph, region_count, plus_region, minus_region, interface_genus, &
            exterior_interface, sample_offsets, sample_points, sample_normals, &
            sample_weights, status)
        type(boundary_region_graph_t), intent(inout) :: graph
        integer, intent(in) :: region_count
        integer, intent(in) :: plus_region(:), minus_region(:), interface_genus(:)
        logical, intent(in) :: exterior_interface(:)
        integer, intent(in) :: sample_offsets(:)
        real(dp), intent(in) :: sample_points(:, :), sample_normals(:, :)
        real(dp), intent(in) :: sample_weights(:)
        integer, intent(out) :: status

        call clear_boundary_region_graph(graph)
        status = 1
        call initialize_region_interface_graph( &
            graph%topology, region_count, plus_region, minus_region, status)
        if (status /= 0) return
        if (.not. valid_geometry( &
            plus_region, interface_genus, exterior_interface, sample_offsets, &
            sample_points, sample_normals, sample_weights)) then
            status = 1
            return
        end if
        allocate(graph%interface_genus(size(interface_genus)), &
            graph%exterior_interface(size(exterior_interface)), &
            graph%sample_offsets(size(sample_offsets)), &
            graph%sample_points(3, size(sample_points, 2)), &
            graph%sample_normals(3, size(sample_normals, 2)), &
            graph%sample_weights(size(sample_weights)))
        graph%interface_genus = interface_genus
        graph%exterior_interface = exterior_interface
        graph%sample_offsets = sample_offsets
        graph%sample_points = sample_points
        graph%sample_normals = sample_normals
        graph%sample_weights = sample_weights
        status = 0
    end subroutine initialize_boundary_region_graph

    subroutine validate_boundary_region_graph(graph, status)
        type(boundary_region_graph_t), intent(in) :: graph
        integer, intent(out) :: status

        status = 1
        call validate_region_interface_graph(graph%topology, status)
        if (status /= 0) return
        if (.not. allocated(graph%interface_genus) .or. &
            .not. allocated(graph%exterior_interface) .or. &
            .not. allocated(graph%sample_offsets) .or. &
            .not. allocated(graph%sample_points) .or. &
            .not. allocated(graph%sample_normals) .or. &
            .not. allocated(graph%sample_weights)) return
        if (.not. valid_geometry_from_graph(graph)) return
        status = 0
    end subroutine validate_boundary_region_graph

    subroutine boundary_region_graph_incidence(graph, incidence, status)
        type(boundary_region_graph_t), intent(in) :: graph
        integer, allocatable, intent(out) :: incidence(:, :)
        integer, intent(out) :: status

        if (allocated(incidence)) deallocate(incidence)
        call validate_boundary_region_graph(graph, status)
        if (status /= 0) then
            allocate(incidence(0, 0))
            return
        end if
        call region_interface_graph_incidence(graph%topology, incidence, status)
    end subroutine boundary_region_graph_incidence

    subroutine boundary_region_graph_components( &
            graph, components, component_count, status)
        type(boundary_region_graph_t), intent(in) :: graph
        integer, allocatable, intent(out) :: components(:)
        integer, intent(out) :: component_count
        integer, intent(out) :: status

        if (allocated(components)) deallocate(components)
        component_count = 0
        call validate_boundary_region_graph(graph, status)
        if (status /= 0) then
            allocate(components(0))
            return
        end if
        call region_interface_graph_components( &
            graph%topology, components, component_count, status)
    end subroutine boundary_region_graph_components

    subroutine boundary_region_graph_cycle_basis( &
            graph, cycle_basis, cycle_count, status)
        type(boundary_region_graph_t), intent(in) :: graph
        integer, allocatable, intent(out) :: cycle_basis(:, :)
        integer, intent(out) :: cycle_count
        integer, intent(out) :: status

        if (allocated(cycle_basis)) deallocate(cycle_basis)
        cycle_count = 0
        call validate_boundary_region_graph(graph, status)
        if (status /= 0) then
            allocate(cycle_basis(0, 0))
            return
        end if
        call region_interface_graph_cycle_basis( &
            graph%topology, cycle_basis, cycle_count, status)
    end subroutine boundary_region_graph_cycle_basis

    subroutine boundary_region_graph_interface_samples( &
            graph, interface_id, points, normals, weights, status)
        type(boundary_region_graph_t), intent(in) :: graph
        integer, intent(in) :: interface_id
        real(dp), allocatable, intent(out) :: points(:, :), normals(:, :)
        real(dp), allocatable, intent(out) :: weights(:)
        integer, intent(out) :: status
        integer :: first_sample, last_sample, sample_count

        if (allocated(points)) deallocate(points)
        if (allocated(normals)) deallocate(normals)
        if (allocated(weights)) deallocate(weights)
        call validate_boundary_region_graph(graph, status)
        if (status /= 0 .or. interface_id < 1 .or. &
            interface_id >= size(graph%sample_offsets)) then
            allocate(points(3, 0), normals(3, 0), weights(0))
            status = 1
            return
        end if
        first_sample = graph%sample_offsets(interface_id)
        last_sample = graph%sample_offsets(interface_id + 1) - 1
        sample_count = max(0, last_sample - first_sample + 1)
        allocate(points(3, sample_count), normals(3, sample_count), weights(sample_count))
        if (sample_count > 0) then
            points = graph%sample_points(:, first_sample:last_sample)
            normals = graph%sample_normals(:, first_sample:last_sample)
            weights = graph%sample_weights(first_sample:last_sample)
        end if
        status = 0
    end subroutine boundary_region_graph_interface_samples

    subroutine boundary_region_graph_interface_metadata( &
            graph, interface_genus, exterior_interface, status)
        type(boundary_region_graph_t), intent(in) :: graph
        integer, allocatable, intent(out) :: interface_genus(:)
        logical, allocatable, intent(out) :: exterior_interface(:)
        integer, intent(out) :: status

        if (allocated(interface_genus)) deallocate(interface_genus)
        if (allocated(exterior_interface)) deallocate(exterior_interface)
        call validate_boundary_region_graph(graph, status)
        if (status /= 0) then
            allocate(interface_genus(0), exterior_interface(0))
            return
        end if
        allocate(interface_genus(size(graph%interface_genus)), &
            exterior_interface(size(graph%exterior_interface)))
        interface_genus = graph%interface_genus
        exterior_interface = graph%exterior_interface
        status = 0
    end subroutine boundary_region_graph_interface_metadata

    subroutine clear_boundary_region_graph(graph)
        type(boundary_region_graph_t), intent(inout) :: graph

        if (allocated(graph%interface_genus)) deallocate(graph%interface_genus)
        if (allocated(graph%exterior_interface)) deallocate(graph%exterior_interface)
        if (allocated(graph%sample_offsets)) deallocate(graph%sample_offsets)
        if (allocated(graph%sample_points)) deallocate(graph%sample_points)
        if (allocated(graph%sample_normals)) deallocate(graph%sample_normals)
        if (allocated(graph%sample_weights)) deallocate(graph%sample_weights)
    end subroutine clear_boundary_region_graph

    logical function valid_geometry( &
            plus_region, interface_genus, exterior_interface, sample_offsets, &
            sample_points, sample_normals, sample_weights) result(valid)
        integer, intent(in) :: plus_region(:), interface_genus(:)
        logical, intent(in) :: exterior_interface(:)
        integer, intent(in) :: sample_offsets(:)
        real(dp), intent(in) :: sample_points(:, :), sample_normals(:, :)
        real(dp), intent(in) :: sample_weights(:)
        integer :: interface_count, sample_count

        interface_count = size(plus_region)
        sample_count = size(sample_weights)
        valid = size(interface_genus) == interface_count .and. &
            size(exterior_interface) == interface_count .and. &
            size(sample_offsets) == interface_count + 1 .and. &
            size(sample_points, 1) == 3 .and. size(sample_normals, 1) == 3 .and. &
            size(sample_points, 2) == sample_count .and. &
            size(sample_normals, 2) == sample_count .and. &
            sample_offsets(1) == 1 .and. sample_offsets(interface_count + 1) == &
            sample_count + 1 .and. all(interface_genus >= 0) .and. &
            all(sample_offsets(1:) >= 1) .and. all(sample_offsets(1:) <= sample_count + 1)
        if (.not. valid) return
        if (interface_count == 0) then
            valid = finite_real(sample_points) .and. finite_real(sample_normals) .and. &
                finite_vector(sample_weights) .and. all(sample_weights > 0.0_dp)
            return
        end if
        valid = all(sample_offsets(2:) >= sample_offsets(:interface_count)) .and. &
            finite_real(sample_points) .and. finite_real(sample_normals) .and. &
            finite_vector(sample_weights) .and. all(sample_weights > 0.0_dp)
        if (.not. valid) return
        if (sample_count > 0) then
            valid = all(sqrt(sum(sample_normals*sample_normals, dim=1)) > tiny(1.0_dp))
        end if
    end function valid_geometry

    logical function valid_geometry_from_graph(graph) result(valid)
        type(boundary_region_graph_t), intent(in) :: graph
        integer :: status

        call validate_region_interface_graph(graph%topology, status)
        if (status /= 0) then
            valid = .false.
            return
        end if
        valid = valid_geometry( &
            graph%interface_genus, graph%interface_genus, graph%exterior_interface, &
            graph%sample_offsets, graph%sample_points, graph%sample_normals, &
            graph%sample_weights)
    end function valid_geometry_from_graph

    pure logical function finite_real(values) result(valid)
        real(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(values))
    end function finite_real

    pure logical function finite_vector(values) result(valid)
        real(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(values))
    end function finite_vector

end module fortfem_boundary_region_graph
