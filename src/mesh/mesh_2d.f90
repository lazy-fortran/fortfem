module fortfem_mesh_2d
    use fortfem_kinds
    implicit none
    private
    
    type, public :: mesh_2d_t
        ! Vertices
        integer :: n_vertices = 0
        real(dp), allocatable :: vertices(:,:)  ! (2, n_vertices) - x,y coordinates
        
        ! Triangles  
        integer :: n_triangles = 0
        integer, allocatable :: triangles(:,:)  ! (3, n_triangles) - vertex indices
        
        ! Edges
        integer :: n_edges = 0
        integer, allocatable :: edges(:,:)      ! (2, n_edges) - vertex indices
        integer, allocatable :: edge_to_triangles(:,:)  ! (2, n_edges) - triangle indices
        
        ! Edge DOF numbering for H(curl) elements
        integer :: n_edge_dofs = 0
        integer, allocatable :: edge_to_dof(:)    ! (n_edges) - maps edge to DOF number
        integer, allocatable :: dof_to_edge(:)    ! (n_edge_dofs) - maps DOF to edge number
        integer :: n_interior_dofs = 0
        integer :: n_boundary_dofs = 0
        
        ! Boundary
        integer :: n_boundary_edges = 0
        integer, allocatable :: boundary_edges(:)     ! Indices of boundary edges
        logical, allocatable :: is_boundary_vertex(:) ! Flag for boundary vertices
        
        ! Connectivity
        type(sparse_int_list_t), allocatable :: vertex_to_triangles(:)
        type(sparse_int_list_t), allocatable :: vertex_to_vertices(:)
        
    contains
        procedure :: create_rectangular
        procedure :: create_unit_disk
        procedure :: create_from_boundary
        procedure :: build_connectivity
        procedure :: build_edge_connectivity
        procedure :: get_edge_vertices
        procedure :: get_edge_triangles
        procedure :: get_edge_length_tangent
        procedure :: is_boundary_edge
        procedure :: build_edge_dof_numbering
        procedure :: get_n_edge_dofs
        procedure :: get_last_interior_dof
        procedure :: get_first_boundary_dof
        procedure :: get_triangle_edge_dofs
        procedure :: find_boundary
        procedure :: build_edge_to_triangle_mapping
        procedure :: compute_areas
        procedure :: save_to_file
        procedure :: load_from_file
        procedure :: destroy
    end type mesh_2d_t
    
    ! Helper type for sparse connectivity
    type :: sparse_int_list_t
        integer :: n = 0
        integer, allocatable :: items(:)
    end type sparse_int_list_t
    
contains

    subroutine create_rectangular(this, nx, ny, x_min, x_max, y_min, y_max)
        class(mesh_2d_t), intent(out) :: this
        integer, intent(in) :: nx, ny  ! Number of vertices in x,y directions
        real(dp), intent(in) :: x_min, x_max, y_min, y_max
        
        integer :: i, j, k, v1, v2, v3, v4
        real(dp) :: dx, dy
        
        ! Total vertices and triangles
        this%n_vertices = nx * ny
        this%n_triangles = 2 * (nx-1) * (ny-1)
        
        ! Allocate arrays
        allocate(this%vertices(2, this%n_vertices))
        allocate(this%triangles(3, this%n_triangles))
        
        ! Grid spacing
        dx = (x_max - x_min) / real(nx - 1, dp)
        dy = (y_max - y_min) / real(ny - 1, dp)
        
        ! Create vertices
        k = 0
        do j = 1, ny
            do i = 1, nx
                k = k + 1
                this%vertices(1, k) = x_min + (i-1) * dx
                this%vertices(2, k) = y_min + (j-1) * dy
            end do
        end do
        
        ! Create triangles (each rectangle split into 2 triangles)
        k = 0
        do j = 1, ny-1
            do i = 1, nx-1
                ! Vertices of rectangle
                v1 = (j-1)*nx + i       ! bottom-left
                v2 = (j-1)*nx + i + 1   ! bottom-right
                v3 = j*nx + i + 1       ! top-right
                v4 = j*nx + i           ! top-left
                
                ! Lower triangle
                k = k + 1
                this%triangles(1, k) = v1
                this%triangles(2, k) = v2
                this%triangles(3, k) = v3
                
                ! Upper triangle
                k = k + 1
                this%triangles(1, k) = v1
                this%triangles(2, k) = v3
                this%triangles(3, k) = v4
            end do
        end do
        
    end subroutine create_rectangular
    
    subroutine build_connectivity(this)
        class(mesh_2d_t), intent(inout) :: this
        
        integer :: t, v, i, j
        integer, allocatable :: temp_list(:)
        
        ! Allocate connectivity arrays
        if (allocated(this%vertex_to_triangles)) deallocate(this%vertex_to_triangles)
        if (allocated(this%vertex_to_vertices)) deallocate(this%vertex_to_vertices)
        
        allocate(this%vertex_to_triangles(this%n_vertices))
        allocate(this%vertex_to_vertices(this%n_vertices))
        allocate(temp_list(20))  ! Temporary list for building connectivity
        
        ! Build vertex-to-triangle connectivity
        do v = 1, this%n_vertices
            this%vertex_to_triangles(v)%n = 0
            
            ! Count triangles containing this vertex
            do t = 1, this%n_triangles
                do i = 1, 3
                    if (this%triangles(i, t) == v) then
                        this%vertex_to_triangles(v)%n = this%vertex_to_triangles(v)%n + 1
                        temp_list(this%vertex_to_triangles(v)%n) = t
                        exit
                    end if
                end do
            end do
            
            ! Store triangle list
            if (this%vertex_to_triangles(v)%n > 0) then
                allocate(this%vertex_to_triangles(v)%items(this%vertex_to_triangles(v)%n))
                this%vertex_to_triangles(v)%items = temp_list(1:this%vertex_to_triangles(v)%n)
            end if
        end do
        
        ! Build edges from triangles
        call build_edges_from_triangles(this)
        
        deallocate(temp_list)
        
    end subroutine build_connectivity
    
    subroutine build_edge_connectivity(this)
        class(mesh_2d_t), intent(inout) :: this
        
        ! Build edge connectivity for H(curl) finite elements
        ! This calls the existing build_connectivity but ensures edges are built
        call this%build_connectivity()
        
        ! Ensure boundary edges are identified
        call this%find_boundary()
    end subroutine build_edge_connectivity
    
    subroutine get_edge_vertices(this, edge_index, vertices)
        class(mesh_2d_t), intent(in) :: this
        integer, intent(in) :: edge_index
        integer, intent(out) :: vertices(2)
        
        if (edge_index < 1 .or. edge_index > this%n_edges) then
            error stop "Edge index out of bounds"
        end if
        
        if (.not. allocated(this%edges)) then
            error stop "Edge connectivity not built. Call build_edge_connectivity first"
        end if
        
        vertices(1) = this%edges(1, edge_index)
        vertices(2) = this%edges(2, edge_index)
    end subroutine get_edge_vertices
    
    function is_boundary_edge(this, edge_index) result(is_boundary)
        class(mesh_2d_t), intent(in) :: this
        integer, intent(in) :: edge_index
        logical :: is_boundary
        
        integer :: i
        
        if (edge_index < 1 .or. edge_index > this%n_edges) then
            error stop "Edge index out of bounds"
        end if
        
        if (.not. allocated(this%boundary_edges)) then
            error stop "Boundary edges not identified. Call build_edge_connectivity first"
        end if
        
        ! Check if edge is in boundary_edges array
        is_boundary = .false.
        do i = 1, this%n_boundary_edges
            if (this%boundary_edges(i) == edge_index) then
                is_boundary = .true.
                return
            end if
        end do
    end function is_boundary_edge
    
    subroutine get_edge_triangles(this, edge_index, triangles)
        class(mesh_2d_t), intent(in) :: this
        integer, intent(in) :: edge_index
        integer, allocatable, intent(out) :: triangles(:)
        
        integer :: t, i, j, v1, v2, count
        integer, allocatable :: temp_triangles(:)
        
        if (edge_index < 1 .or. edge_index > this%n_edges) then
            error stop "Edge index out of bounds"
        end if
        
        if (.not. allocated(this%edges)) then
            error stop "Edge connectivity not built. Call build_edge_connectivity first"
        end if
        
        ! Get edge vertices
        v1 = this%edges(1, edge_index)
        v2 = this%edges(2, edge_index)
        
        ! Allocate temporary array for worst case
        allocate(temp_triangles(2))
        count = 0
        
        ! Find triangles containing this edge
        do t = 1, this%n_triangles
            do i = 1, 3
                j = mod(i, 3) + 1
                if ((this%triangles(i, t) == v1 .and. this%triangles(j, t) == v2) .or. &
                    (this%triangles(i, t) == v2 .and. this%triangles(j, t) == v1)) then
                    count = count + 1
                    if (count > 2) then
                        error stop "More than 2 triangles share an edge - invalid mesh"
                    end if
                    temp_triangles(count) = t
                    exit
                end if
            end do
        end do
        
        ! Copy to output array
        allocate(triangles(count))
        triangles = temp_triangles(1:count)
        
        deallocate(temp_triangles)
    end subroutine get_edge_triangles
    
    subroutine get_edge_length_tangent(this, edge_index, length, tangent)
        class(mesh_2d_t), intent(in) :: this
        integer, intent(in) :: edge_index
        real(dp), intent(out) :: length, tangent(2)
        
        integer :: v1, v2
        real(dp) :: dx, dy
        
        if (edge_index < 1 .or. edge_index > this%n_edges) then
            error stop "Edge index out of bounds"
        end if
        
        if (.not. allocated(this%edges)) then
            error stop "Edge connectivity not built. Call build_edge_connectivity first"
        end if
        
        ! Get edge vertices
        v1 = this%edges(1, edge_index)
        v2 = this%edges(2, edge_index)
        
        ! Compute edge vector
        dx = this%vertices(1, v2) - this%vertices(1, v1)
        dy = this%vertices(2, v2) - this%vertices(2, v1)
        
        ! Compute edge length
        length = sqrt(dx**2 + dy**2)
        
        ! Compute unit tangent vector
        if (length > 0.0_dp) then
            tangent(1) = dx / length
            tangent(2) = dy / length
        else
            ! Zero-length edge (degenerate)
            tangent(1) = 0.0_dp
            tangent(2) = 0.0_dp
        end if
    end subroutine get_edge_length_tangent
    
    subroutine build_edge_dof_numbering(this)
        class(mesh_2d_t), intent(inout) :: this
        
        integer :: i, dof_count
        
        if (.not. allocated(this%edges)) then
            error stop "Edge connectivity not built. Call build_edge_connectivity first"
        end if
        
        ! RT0/Nédélec lowest order: one DOF per edge
        this%n_edge_dofs = this%n_edges
        
        ! Allocate DOF mapping arrays
        allocate(this%edge_to_dof(this%n_edges))
        allocate(this%dof_to_edge(this%n_edge_dofs))
        
        ! Number interior DOFs first (0-based indexing)
        dof_count = 0
        do i = 1, this%n_edges
            if (.not. this%is_boundary_edge(i)) then
                this%edge_to_dof(i) = dof_count
                this%dof_to_edge(dof_count + 1) = i
                dof_count = dof_count + 1
            end if
        end do
        
        this%n_interior_dofs = dof_count
        
        ! Number boundary DOFs next
        do i = 1, this%n_edges
            if (this%is_boundary_edge(i)) then
                this%edge_to_dof(i) = dof_count
                this%dof_to_edge(dof_count + 1) = i
                dof_count = dof_count + 1
            end if
        end do
        
        this%n_boundary_dofs = dof_count - this%n_interior_dofs
    end subroutine build_edge_dof_numbering
    
    function get_n_edge_dofs(this) result(n_dofs)
        class(mesh_2d_t), intent(in) :: this
        integer :: n_dofs
        
        n_dofs = this%n_edge_dofs
    end function get_n_edge_dofs
    
    function get_last_interior_dof(this) result(last_dof)
        class(mesh_2d_t), intent(in) :: this
        integer :: last_dof
        
        last_dof = this%n_interior_dofs - 1
    end function get_last_interior_dof
    
    function get_first_boundary_dof(this) result(first_dof)
        class(mesh_2d_t), intent(in) :: this
        integer :: first_dof
        
        first_dof = this%n_interior_dofs
    end function get_first_boundary_dof
    
    subroutine get_triangle_edge_dofs(this, triangle_index, edge_dofs)
        class(mesh_2d_t), intent(in) :: this
        integer, intent(in) :: triangle_index
        integer, intent(out) :: edge_dofs(3)
        
        integer :: i, j, v1, v2, edge_index
        
        if (triangle_index < 1 .or. triangle_index > this%n_triangles) then
            error stop "Triangle index out of bounds"
        end if
        
        if (.not. allocated(this%edge_to_dof)) then
            error stop "Edge DOF numbering not built. Call build_edge_dof_numbering first"
        end if
        
        ! Find the 3 edges of this triangle
        do i = 1, 3
            j = mod(i, 3) + 1
            v1 = min(this%triangles(i, triangle_index), this%triangles(j, triangle_index))
            v2 = max(this%triangles(i, triangle_index), this%triangles(j, triangle_index))
            
            ! Find this edge in the edge list
            edge_index = 0
            do edge_index = 1, this%n_edges
                if (this%edges(1, edge_index) == v1 .and. this%edges(2, edge_index) == v2) then
                    exit
                end if
            end do
            
            if (edge_index > this%n_edges) then
                error stop "Edge not found in triangle"
            end if
            
            edge_dofs(i) = this%edge_to_dof(edge_index)
        end do
    end subroutine get_triangle_edge_dofs
    
    subroutine build_edges_from_triangles(this)
        class(mesh_2d_t), intent(inout) :: this
        
        integer :: t, i, j, v1, v2, e, found
        integer :: max_edges
        integer, allocatable :: temp_edges(:,:)
        
        ! Maximum possible edges
        max_edges = 3 * this%n_triangles
        allocate(temp_edges(2, max_edges))
        
        this%n_edges = 0
        
        ! Extract unique edges from triangles
        do t = 1, this%n_triangles
            do i = 1, 3
                j = mod(i, 3) + 1
                v1 = min(this%triangles(i, t), this%triangles(j, t))
                v2 = max(this%triangles(i, t), this%triangles(j, t))
                
                ! Check if edge already exists
                found = 0
                do e = 1, this%n_edges
                    if (temp_edges(1, e) == v1 .and. temp_edges(2, e) == v2) then
                        found = e
                        exit
                    end if
                end do
                
                if (found == 0) then
                    this%n_edges = this%n_edges + 1
                    temp_edges(1, this%n_edges) = v1
                    temp_edges(2, this%n_edges) = v2
                end if
            end do
        end do
        
        ! Copy to final array
        allocate(this%edges(2, this%n_edges))
        this%edges = temp_edges(:, 1:this%n_edges)
        
        deallocate(temp_edges)
        
        ! Build edge-to-triangle mapping
        call this%build_edge_to_triangle_mapping()
        
    end subroutine build_edges_from_triangles
    
    subroutine build_edge_to_triangle_mapping(this)
        class(mesh_2d_t), intent(inout) :: this
        
        integer :: t, i, j, v1, v2, e
        integer, allocatable :: edge_triangle_count(:)
        
        if (.not. allocated(this%edges)) then
            error stop "Edges must be built before edge-to-triangle mapping"
        end if
        
        ! Allocate edge-to-triangle mapping (2 triangles max per edge)
        allocate(this%edge_to_triangles(2, this%n_edges))
        allocate(edge_triangle_count(this%n_edges))
        
        ! Initialize
        this%edge_to_triangles = 0
        edge_triangle_count = 0
        
        ! For each triangle, find its edges and record the triangle-edge relationship
        do t = 1, this%n_triangles
            do i = 1, 3
                j = mod(i, 3) + 1
                v1 = min(this%triangles(i, t), this%triangles(j, t))
                v2 = max(this%triangles(i, t), this%triangles(j, t))
                
                ! Find the edge with these vertices
                do e = 1, this%n_edges
                    if (this%edges(1, e) == v1 .and. this%edges(2, e) == v2) then
                        edge_triangle_count(e) = edge_triangle_count(e) + 1
                        if (edge_triangle_count(e) <= 2) then
                            this%edge_to_triangles(edge_triangle_count(e), e) = t
                        else
                            error stop "Edge belongs to more than 2 triangles - invalid mesh"
                        end if
                        exit
                    end if
                end do
            end do
        end do
        
        deallocate(edge_triangle_count)
        
    end subroutine build_edge_to_triangle_mapping
    
    subroutine find_boundary(this)
        class(mesh_2d_t), intent(inout) :: this
        
        integer :: e, t, i, j, v1, v2, count
        integer, allocatable :: temp_boundary(:)
        
        if (.not. allocated(this%edges)) then
            call build_edges_from_triangles(this)
        end if
        
        allocate(temp_boundary(this%n_edges))
        allocate(this%is_boundary_vertex(this%n_vertices))
        this%is_boundary_vertex = .false.
        
        this%n_boundary_edges = 0
        
        ! Find edges that belong to only one triangle
        do e = 1, this%n_edges
            v1 = this%edges(1, e)
            v2 = this%edges(2, e)
            count = 0
            
            ! Count triangles containing this edge
            do t = 1, this%n_triangles
                do i = 1, 3
                    j = mod(i, 3) + 1
                    if ((this%triangles(i, t) == v1 .and. this%triangles(j, t) == v2) .or. &
                        (this%triangles(i, t) == v2 .and. this%triangles(j, t) == v1)) then
                        count = count + 1
                        exit
                    end if
                end do
            end do
            
            ! Boundary edge belongs to exactly one triangle
            if (count == 1) then
                this%n_boundary_edges = this%n_boundary_edges + 1
                temp_boundary(this%n_boundary_edges) = e
                this%is_boundary_vertex(v1) = .true.
                this%is_boundary_vertex(v2) = .true.
            end if
        end do
        
        ! Store boundary edges
        allocate(this%boundary_edges(this%n_boundary_edges))
        this%boundary_edges = temp_boundary(1:this%n_boundary_edges)
        
        deallocate(temp_boundary)
        
    end subroutine find_boundary
    
    function compute_areas(this) result(areas)
        class(mesh_2d_t), intent(in) :: this
        real(dp), allocatable :: areas(:)
        
        integer :: t
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        allocate(areas(this%n_triangles))
        
        do t = 1, this%n_triangles
            x1 = this%vertices(1, this%triangles(1, t))
            y1 = this%vertices(2, this%triangles(1, t))
            x2 = this%vertices(1, this%triangles(2, t))
            y2 = this%vertices(2, this%triangles(2, t))
            x3 = this%vertices(1, this%triangles(3, t))
            y3 = this%vertices(2, this%triangles(3, t))
            
            ! Area = 0.5 * |det([[x1-x3, x2-x3], [y1-y3, y2-y3]])|
            areas(t) = 0.5_dp * abs((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3))
        end do
        
    end function compute_areas
    
    subroutine save_to_file(this, filename)
        class(mesh_2d_t), intent(in) :: this
        character(len=*), intent(in) :: filename
        
        integer :: unit, i
        
        open(newunit=unit, file=filename, status='replace', action='write')
        
        ! Write header
        write(unit, '(a)') '# FortFEM 2D Mesh'
        write(unit, '(a,i0)') '# Vertices: ', this%n_vertices
        write(unit, '(a,i0)') '# Triangles: ', this%n_triangles
        
        ! Write vertices
        write(unit, '(a)') 'VERTICES'
        write(unit, '(i0)') this%n_vertices
        do i = 1, this%n_vertices
            write(unit, '(2(es23.16,1x))') this%vertices(:, i)
        end do
        
        ! Write triangles
        write(unit, '(a)') 'TRIANGLES'
        write(unit, '(i0)') this%n_triangles
        do i = 1, this%n_triangles
            write(unit, '(3(i0,1x))') this%triangles(:, i)
        end do
        
        close(unit)
        
    end subroutine save_to_file
    
    subroutine load_from_file(this, filename)
        class(mesh_2d_t), intent(out) :: this
        character(len=*), intent(in) :: filename
        
        integer :: unit, i
        character(len=100) :: line
        
        open(newunit=unit, file=filename, status='old', action='read')
        
        ! Skip header
        do
            read(unit, '(a)') line
            if (line == 'VERTICES') exit
        end do
        
        ! Read vertices
        read(unit, *) this%n_vertices
        allocate(this%vertices(2, this%n_vertices))
        do i = 1, this%n_vertices
            read(unit, *) this%vertices(:, i)
        end do
        
        ! Read triangles
        read(unit, '(a)') line  ! Should be 'TRIANGLES'
        read(unit, *) this%n_triangles
        allocate(this%triangles(3, this%n_triangles))
        do i = 1, this%n_triangles
            read(unit, *) this%triangles(:, i)
        end do
        
        close(unit)
        
    end subroutine load_from_file
    
    subroutine destroy(this)
        class(mesh_2d_t), intent(inout) :: this
        
        integer :: i
        
        if (allocated(this%vertices)) deallocate(this%vertices)
        if (allocated(this%triangles)) deallocate(this%triangles)
        if (allocated(this%edges)) deallocate(this%edges)
        if (allocated(this%edge_to_triangles)) deallocate(this%edge_to_triangles)
        if (allocated(this%boundary_edges)) deallocate(this%boundary_edges)
        if (allocated(this%is_boundary_vertex)) deallocate(this%is_boundary_vertex)
        if (allocated(this%edge_to_dof)) deallocate(this%edge_to_dof)
        if (allocated(this%dof_to_edge)) deallocate(this%dof_to_edge)
        
        if (allocated(this%vertex_to_triangles)) then
            do i = 1, size(this%vertex_to_triangles)
                if (allocated(this%vertex_to_triangles(i)%items)) then
                    deallocate(this%vertex_to_triangles(i)%items)
                end if
            end do
            deallocate(this%vertex_to_triangles)
        end if
        
        if (allocated(this%vertex_to_vertices)) then
            do i = 1, size(this%vertex_to_vertices)
                if (allocated(this%vertex_to_vertices(i)%items)) then
                    deallocate(this%vertex_to_vertices(i)%items)
                end if
            end do
            deallocate(this%vertex_to_vertices)
        end if
        
        this%n_vertices = 0
        this%n_triangles = 0
        this%n_edges = 0
        this%n_boundary_edges = 0
        this%n_edge_dofs = 0
        this%n_interior_dofs = 0
        this%n_boundary_dofs = 0
        
    end subroutine destroy

    ! Clean mesh generation stubs
    subroutine create_unit_disk(this, max_element_size)
        class(mesh_2d_t), intent(out) :: this
        real(dp), intent(in) :: max_element_size
        
        write(*,*) "STUB: create_unit_disk with resolution", max_element_size
        call this%create_rectangular(10, 10, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp)
    end subroutine create_unit_disk

    subroutine create_from_boundary(this, boundary, resolution)
        class(mesh_2d_t), intent(out) :: this
        class(*), intent(in) :: boundary
        real(dp), intent(in) :: resolution
        
        write(*,*) "STUB: create_from_boundary with resolution", resolution
        call this%create_rectangular(10, 10, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    end subroutine create_from_boundary

end module fortfem_mesh_2d