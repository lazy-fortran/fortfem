module triangle_io
    !> Triangle mesh file format I/O routines.
    !
    !  This module provides routines to read and write mesh files in the
    !  Triangle format developed by Jonathan Shewchuk.
    !
    !  Triangle file format reference:
    !    https://www.cs.cmu.edu/~quake/triangle.html
    !
    !  File formats supported:
    !    .node - Vertex coordinates
    !    .ele  - Triangle connectivity
    !    .poly - Planar Straight Line Graph (PSLG) for input
    !
    !  The .node file format:
    !    First line: <# of vertices> <dimension (2)> <# of attributes> <# of
    !                boundary markers (0 or 1)>
    !    Remaining lines: <vertex #> <x> <y> [attributes] [boundary marker]
    !
    !  The .ele file format:
    !    First line: <# of triangles> <nodes per triangle (3)> <# of attributes>
    !    Remaining lines: <triangle #> <node 1> <node 2> <node 3> [attributes]
    !
    !  The .poly file format:
    !    First line: <# of vertices> <dimension (2)> <# of attributes> <# of
    !                boundary markers (0 or 1)>
    !    Following lines: vertex list (same as .node) OR 0 if vertices in
    !                     separate .node file
    !    One line: <# of segments> <# of boundary markers (0 or 1)>
    !    Following lines: <segment #> <endpoint 1> <endpoint 2> [boundary marker]
    !    One line: <# of holes>
    !    Following lines: <hole #> <x> <y>
    !
    !  Note: Triangle uses 1-based indexing by default (-z for 0-based).
    !
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: read_triangle_node_file
    public :: read_triangle_ele_file
    public :: read_triangle_poly_file
    public :: write_triangle_poly_file
    public :: read_triangle_mesh
    public :: ensure_triangle_available

contains

    subroutine read_triangle_node_file(filename, vertices, n_vertices, stat)
        !> Read vertices from Triangle .node file format.
        !
        !  Arguments:
        !    filename   - Path to .node file
        !    vertices   - Output array (2, n_vertices) of x,y coordinates
        !    n_vertices - Number of vertices read
        !    stat       - Status: 0 = success, nonzero = error
        !
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: vertices(:,:)
        integer, intent(out) :: n_vertices
        integer, intent(out) :: stat

        integer :: unit_num, i, idx, dim, n_attr, n_bm
        real(dp) :: x, y
        logical :: file_exists
        character(len=256) :: line

        stat = 0
        n_vertices = 0

        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            stat = 1
            return
        end if

        open(newunit=unit_num, file=filename, status="old", action="read",    &
             iostat=stat)
        if (stat /= 0) return

        ! Skip comment lines starting with #
        do
            read(unit_num, '(A)', iostat=stat) line
            if (stat /= 0) then
                close(unit_num)
                return
            end if
            line = adjustl(line)
            if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
        end do

        ! Parse header: n_vertices dimension n_attributes n_boundary_markers
        read(line, *, iostat=stat) n_vertices, dim, n_attr, n_bm
        if (stat /= 0) then
            close(unit_num)
            return
        end if

        if (dim /= 2) then
            stat = 2  ! Only 2D supported
            close(unit_num)
            return
        end if

        allocate(vertices(2, n_vertices))

        do i = 1, n_vertices
            ! Skip comment lines
            do
                read(unit_num, '(A)', iostat=stat) line
                if (stat /= 0) then
                    close(unit_num)
                    return
                end if
                line = adjustl(line)
                if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
            end do

            read(line, *, iostat=stat) idx, x, y
            if (stat /= 0) then
                close(unit_num)
                return
            end if
            vertices(1, i) = x
            vertices(2, i) = y
        end do

        close(unit_num)
        stat = 0
    end subroutine read_triangle_node_file

    subroutine read_triangle_ele_file(filename, triangles, n_triangles, stat)
        !> Read triangles from Triangle .ele file format.
        !
        !  Arguments:
        !    filename    - Path to .ele file
        !    triangles   - Output array (3, n_triangles) of vertex indices
        !    n_triangles - Number of triangles read
        !    stat        - Status: 0 = success, nonzero = error
        !
        character(len=*), intent(in) :: filename
        integer, allocatable, intent(out) :: triangles(:,:)
        integer, intent(out) :: n_triangles
        integer, intent(out) :: stat

        integer :: unit_num, i, idx, nodes_per_tri, n_attr
        integer :: v1, v2, v3
        logical :: file_exists
        character(len=256) :: line

        stat = 0
        n_triangles = 0

        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            stat = 1
            return
        end if

        open(newunit=unit_num, file=filename, status="old", action="read",    &
             iostat=stat)
        if (stat /= 0) return

        ! Skip comment lines
        do
            read(unit_num, '(A)', iostat=stat) line
            if (stat /= 0) then
                close(unit_num)
                return
            end if
            line = adjustl(line)
            if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
        end do

        ! Parse header: n_triangles nodes_per_triangle n_attributes
        read(line, *, iostat=stat) n_triangles, nodes_per_tri, n_attr
        if (stat /= 0) then
            close(unit_num)
            return
        end if

        if (nodes_per_tri /= 3) then
            stat = 2  ! Only linear triangles supported
            close(unit_num)
            return
        end if

        allocate(triangles(3, n_triangles))

        do i = 1, n_triangles
            ! Skip comment lines
            do
                read(unit_num, '(A)', iostat=stat) line
                if (stat /= 0) then
                    close(unit_num)
                    return
                end if
                line = adjustl(line)
                if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
            end do

            read(line, *, iostat=stat) idx, v1, v2, v3
            if (stat /= 0) then
                close(unit_num)
                return
            end if
            triangles(1, i) = v1
            triangles(2, i) = v2
            triangles(3, i) = v3
        end do

        close(unit_num)
        stat = 0
    end subroutine read_triangle_ele_file

    subroutine read_triangle_poly_file(filename, vertices, segments,          &
                                       n_vertices, n_segments, stat)
        !> Read PSLG from Triangle .poly file format.
        !
        !  Arguments:
        !    filename   - Path to .poly file
        !    vertices   - Output array (2, n_vertices) of x,y coordinates
        !    segments   - Output array (2, n_segments) of endpoint indices
        !    n_vertices - Number of vertices read
        !    n_segments - Number of segments read
        !    stat       - Status: 0 = success, nonzero = error
        !
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: vertices(:,:)
        integer, allocatable, intent(out) :: segments(:,:)
        integer, intent(out) :: n_vertices
        integer, intent(out) :: n_segments
        integer, intent(out) :: stat

        integer :: unit_num, i, idx, dim, n_attr, n_bm
        real(dp) :: x, y
        integer :: v1, v2
        logical :: file_exists
        character(len=256) :: line

        stat = 0
        n_vertices = 0
        n_segments = 0

        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            stat = 1
            return
        end if

        open(newunit=unit_num, file=filename, status="old", action="read",    &
             iostat=stat)
        if (stat /= 0) return

        ! Skip comment lines and read vertex header
        do
            read(unit_num, '(A)', iostat=stat) line
            if (stat /= 0) then
                close(unit_num)
                return
            end if
            line = adjustl(line)
            if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
        end do

        read(line, *, iostat=stat) n_vertices, dim, n_attr, n_bm
        if (stat /= 0) then
            close(unit_num)
            return
        end if

        ! Read vertices if present (n_vertices > 0)
        if (n_vertices > 0) then
            allocate(vertices(2, n_vertices))
            do i = 1, n_vertices
                do
                    read(unit_num, '(A)', iostat=stat) line
                    if (stat /= 0) then
                        close(unit_num)
                        return
                    end if
                    line = adjustl(line)
                    if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
                end do
                read(line, *, iostat=stat) idx, x, y
                if (stat /= 0) then
                    close(unit_num)
                    return
                end if
                vertices(1, i) = x
                vertices(2, i) = y
            end do
        end if

        ! Read segment header
        do
            read(unit_num, '(A)', iostat=stat) line
            if (stat /= 0) then
                close(unit_num)
                return
            end if
            line = adjustl(line)
            if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
        end do

        read(line, *, iostat=stat) n_segments, n_bm
        if (stat /= 0) then
            close(unit_num)
            return
        end if

        if (n_segments > 0) then
            allocate(segments(2, n_segments))
            do i = 1, n_segments
                do
                    read(unit_num, '(A)', iostat=stat) line
                    if (stat /= 0) then
                        close(unit_num)
                        return
                    end if
                    line = adjustl(line)
                    if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
                end do
                read(line, *, iostat=stat) idx, v1, v2
                if (stat /= 0) then
                    close(unit_num)
                    return
                end if
                segments(1, i) = v1
                segments(2, i) = v2
            end do
        end if

        close(unit_num)
        stat = 0
    end subroutine read_triangle_poly_file

    subroutine write_triangle_poly_file(filename, vertices, segments,         &
                                        n_vertices, n_segments, stat)
        !> Write PSLG to Triangle .poly file format.
        !
        !  Arguments:
        !    filename   - Path to .poly file
        !    vertices   - Array (2, n_vertices) of x,y coordinates
        !    segments   - Array (2, n_segments) of endpoint indices (1-based)
        !    n_vertices - Number of vertices
        !    n_segments - Number of segments
        !    stat       - Status: 0 = success, nonzero = error
        !
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: vertices(:,:)
        integer, intent(in) :: segments(:,:)
        integer, intent(in) :: n_vertices
        integer, intent(in) :: n_segments
        integer, intent(out) :: stat

        integer :: unit_num, i

        open(newunit=unit_num, file=filename, status="replace", action="write",&
             iostat=stat)
        if (stat /= 0) return

        ! Write vertex header: n_vertices dimension n_attributes n_bm
        write(unit_num, '(I0,A)') n_vertices, " 2 0 0"

        ! Write vertices
        do i = 1, n_vertices
            write(unit_num, '(I0,A,ES23.15,A,ES23.15)')                        &
                i, " ", vertices(1, i), " ", vertices(2, i)
        end do

        ! Write segment header: n_segments n_boundary_markers
        write(unit_num, '(I0,A)') n_segments, " 0"

        ! Write segments
        do i = 1, n_segments
            write(unit_num, '(I0,A,I0,A,I0)')                                  &
                i, " ", segments(1, i), " ", segments(2, i)
        end do

        ! Write holes header (no holes)
        write(unit_num, '(A)') "0"

        close(unit_num)
        stat = 0
    end subroutine write_triangle_poly_file

    subroutine read_triangle_mesh(basename, vertices, triangles,              &
                                  n_vertices, n_triangles, stat)
        !> Read complete mesh from Triangle output files (.node and .ele).
        !
        !  Triangle outputs meshes with a numeric suffix (e.g., .1.node, .1.ele).
        !  This routine reads both files and returns the mesh data.
        !
        !  Arguments:
        !    basename    - Base filename (e.g., "/tmp/mesh" for mesh.1.node)
        !    vertices    - Output array (2, n_vertices) of coordinates
        !    triangles   - Output array (3, n_triangles) of vertex indices
        !    n_vertices  - Number of vertices
        !    n_triangles - Number of triangles
        !    stat        - Status: 0 = success, nonzero = error
        !
        character(len=*), intent(in) :: basename
        real(dp), allocatable, intent(out) :: vertices(:,:)
        integer, allocatable, intent(out) :: triangles(:,:)
        integer, intent(out) :: n_vertices
        integer, intent(out) :: n_triangles
        integer, intent(out) :: stat

        character(len=512) :: node_file, ele_file

        node_file = trim(basename) // ".1.node"
        ele_file = trim(basename) // ".1.ele"

        call read_triangle_node_file(node_file, vertices, n_vertices, stat)
        if (stat /= 0) return

        call read_triangle_ele_file(ele_file, triangles, n_triangles, stat)
    end subroutine read_triangle_mesh

    subroutine ensure_triangle_available(triangle_path, stat)
        !> Ensure Triangle executable is available, downloading if necessary.
        !
        !  This routine checks if Triangle is available at the specified path.
        !  If not, it downloads from netlib and compiles it.
        !
        !  Arguments:
        !    triangle_path - Path where Triangle executable should be
        !    stat          - Status: 0 = success, nonzero = error
        !
        !  Note: Requires wget and gcc to be available in PATH.
        !
        character(len=*), intent(in) :: triangle_path
        integer, intent(out) :: stat

        character(len=1024) :: cmd
        logical :: file_exists

        inquire(file=triangle_path, exist=file_exists)
        if (file_exists) then
            stat = 0
            return
        end if

        ! Download Triangle source from netlib
        write(*, '(A)') "   Downloading Triangle from netlib.org..."
        cmd = "wget -q http://www.netlib.org/voronoi/triangle.zip "           &
              // "-O /tmp/triangle_download.zip 2>/dev/null"
        call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
        if (stat /= 0) then
            write(*, '(A)') "   Warning: Failed to download Triangle"
            return
        end if

        ! Extract
        cmd = "cd /tmp && unzip -o triangle_download.zip >/dev/null 2>&1"
        call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
        if (stat /= 0) then
            write(*, '(A)') "   Warning: Failed to extract Triangle"
            return
        end if

        ! Compile (using gnu89 for K&R C compatibility)
        write(*, '(A)') "   Compiling Triangle..."
        cmd = "cd /tmp && gcc -std=gnu89 -O2 -DLINUX -o " // trim(triangle_path) &
              // " triangle.c -lm 2>/dev/null"
        call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
        if (stat /= 0) then
            write(*, '(A)') "   Warning: Failed to compile Triangle"
            return
        end if

        inquire(file=triangle_path, exist=file_exists)
        if (.not. file_exists) then
            stat = 1
            return
        end if

        stat = 0
    end subroutine ensure_triangle_available

end module triangle_io
