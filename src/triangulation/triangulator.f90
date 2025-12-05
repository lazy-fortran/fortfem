module triangulator
    ! Clean interface between FortFEM mesh API and triangulation backends
    use fortfem_kinds, only: dp
    use triangulation_fortran, only: triangulation_result_t, triangulate_fortran, &
                                     cleanup_triangulation
    use triangle_io, only: write_triangle_poly_file, ensure_triangle_available, &
        read_triangle_mesh
    implicit none
    
    private
    public :: triangulate_boundary, triangulate_points
    
contains

    subroutine triangulate_boundary(boundary_points, segments, mesh_points, &
                                   mesh_triangles, n_points, n_triangles)
        !> Triangulate a boundary defined by points and segments using
        !  constrained Delaunay triangulation.
        !
        !  Preferred backend is the Triangle library (external binary),
        !  with a fallback to the internal Fortran CDT implementation
        !  if Triangle is not available.
        real(dp), intent(in) :: boundary_points(:,:)    ! (2, n_boundary_points)
        integer, intent(in) :: segments(:,:)            ! (2, n_segments)
        real(dp), allocatable, intent(out) :: mesh_points(:,:)     ! (2, n_points)
        integer, allocatable, intent(out) :: mesh_triangles(:,:)   ! (3, n_triangles)
        integer, intent(out) :: n_points, n_triangles
        
        type(triangulation_result_t) :: result
        character(len=*), parameter :: TRIANGLE_PATH = "/tmp/triangle_bin"
        character(len=*), parameter :: POLY_BASE = "/tmp/fortfem_boundary"
        character(len=512) :: cmd
        integer :: stat
        logical :: triangle_ok

        ! Try Triangle external library first
        call ensure_triangle_available(TRIANGLE_PATH, stat)
        triangle_ok = (stat == 0)

        if (triangle_ok) then
            ! Write PSLG to .poly file
            call write_triangle_poly_file(POLY_BASE // ".poly",              &
                                          boundary_points, segments,         &
                                          size(boundary_points, 2),          &
                                          size(segments, 2), stat)

            if (stat == 0) then
                cmd = trim(TRIANGLE_PATH) // " -pQ " // POLY_BASE // ".poly"
                call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
            end if

            if (stat == 0) then
                call read_triangle_mesh(POLY_BASE, mesh_points,              &
                                        mesh_triangles, n_points,           &
                                        n_triangles, stat)
                if (stat == 0 .and. n_triangles > 0) then
                    return
                end if
            end if
        end if

        ! Fallback: internal Fortran CDT implementation
        call triangulate_fortran(boundary_points, segments, result)

        n_points = result%npoints
        n_triangles = result%ntriangles

        allocate(mesh_points(2, n_points))
        allocate(mesh_triangles(3, n_triangles))

        mesh_points = result%points(:, 1:n_points)
        mesh_triangles = result%triangles(:, 1:n_triangles)

        call cleanup_triangulation(result)
        
    end subroutine triangulate_boundary

    subroutine triangulate_points(input_points, mesh_points, mesh_triangles, &
                                 n_points, n_triangles)
        !> Triangulate a set of points (unconstrained Delaunay)
        real(dp), intent(in) :: input_points(:,:)       ! (2, n_input_points)
        real(dp), allocatable, intent(out) :: mesh_points(:,:)     ! (2, n_points)
        integer, allocatable, intent(out) :: mesh_triangles(:,:)   ! (3, n_triangles)
        integer, intent(out) :: n_points, n_triangles
        
        type(triangulation_result_t) :: result
        integer, allocatable :: empty_segments(:,:)
        
        allocate(empty_segments(2, 0))
        
        call triangulate_fortran(input_points, empty_segments, result)
        
        n_points = result%npoints
        n_triangles = result%ntriangles
        
        allocate(mesh_points(2, n_points))
        allocate(mesh_triangles(3, n_triangles))
        
        mesh_points = result%points(:, 1:n_points)
        mesh_triangles = result%triangles(:, 1:n_triangles)
        
        call cleanup_triangulation(result)
        deallocate(empty_segments)
        
    end subroutine triangulate_points

end module triangulator
