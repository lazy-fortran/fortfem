module triangulator
    ! Clean interface between FortFEM mesh API and Delaunay triangulation backend
    use fortfem_kinds, only: dp
    use triangulation_fortran
    implicit none
    
    private
    public :: triangulate_boundary, triangulate_points
    
contains

    subroutine triangulate_boundary(boundary_points, segments, mesh_points, &
                                   mesh_triangles, n_points, n_triangles)
        !> Triangulate a boundary defined by points and segments
        real(dp), intent(in) :: boundary_points(:,:)    ! (2, n_boundary_points)
        integer, intent(in) :: segments(:,:)            ! (2, n_segments)
        real(dp), allocatable, intent(out) :: mesh_points(:,:)     ! (2, n_points)
        integer, allocatable, intent(out) :: mesh_triangles(:,:)   ! (3, n_triangles)
        integer, intent(out) :: n_points, n_triangles
        
        type(triangulation_result_t) :: result
        
        ! Call constrained Delaunay triangulation
        call triangulate_fortran(boundary_points, segments, result)
        
        ! Extract results
        n_points = result%npoints
        n_triangles = result%ntriangles
        
        ! Allocate output arrays
        allocate(mesh_points(2, n_points))
        allocate(mesh_triangles(3, n_triangles))
        
        ! Copy results
        mesh_points = result%points(:, 1:n_points)
        mesh_triangles = result%triangles(:, 1:n_triangles)
        
        ! Cleanup
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
        
        ! No segments for unconstrained triangulation
        allocate(empty_segments(2, 0))
        
        ! Call triangulation
        call triangulate_fortran(input_points, empty_segments, result)
        
        ! Extract results
        n_points = result%npoints
        n_triangles = result%ntriangles
        
        ! Allocate output arrays
        allocate(mesh_points(2, n_points))
        allocate(mesh_triangles(3, n_triangles))
        
        ! Copy results
        mesh_points = result%points(:, 1:n_points)
        mesh_triangles = result%triangles(:, 1:n_triangles)
        
        ! Cleanup
        call cleanup_triangulation(result)
        deallocate(empty_segments)
        
    end subroutine triangulate_points

end module triangulator