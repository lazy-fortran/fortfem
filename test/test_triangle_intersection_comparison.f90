program test_triangle_intersection_comparison
    !> Compare FortFEM CDT with Triangle on a PSLG with crossing segments.
    !
    !  Geometry:
    !    - Outer square: (0,0)-(2,0)-(2,2)-(0,2)
    !    - Horizontal segment through center: (0,1)-(2,1)
    !    - Vertical segment through center:   (1,0)-(1,2)
    !
    !  Both FortFEM and Triangle should:
    !    - Insert a Steiner vertex near (1,1),
    !    - Split the crossing segments into a plus-shaped graph,
    !    - Produce triangulations with no edge crossing the constraints.
    !
    use fortfem_kinds, only : dp
    use fortfem_api,  only : mesh_t, mesh_from_arrays, mesh_from_triangle_files, plot
    use triangulation_fortran, only : triangulation_result_t, triangulate_fortran
    use triangle_io, only : write_triangle_poly_file, ensure_triangle_available
    use check, only : check_condition, check_summary
    implicit none

    real(dp), parameter :: pts(2,8) = reshape([ &
        0.0_dp, 0.0_dp, &  ! 1: outer square
        2.0_dp, 0.0_dp, &  ! 2
        2.0_dp, 2.0_dp, &  ! 3
        0.0_dp, 2.0_dp, &  ! 4
        0.0_dp, 1.0_dp, &  ! 5: mid left
        2.0_dp, 1.0_dp, &  ! 6: mid right
        1.0_dp, 0.0_dp, &  ! 7: mid bottom
        1.0_dp, 2.0_dp ], [2, 8])  ! 8: mid top

    integer, parameter :: segs(2,6) = reshape([ &
        1, 2, &  ! outer square
        2, 3, &
        3, 4, &
        4, 1, &
        5, 6, &  ! horizontal
        7, 8 ], [2, 6])   ! vertical

    type(triangulation_result_t) :: fortfem_result
    type(mesh_t) :: fortfem_mesh, triangle_mesh
    logical :: triangle_ok

    write(*,*) "=== Triangle vs FortFEM Segment Intersection Test ==="
    write(*,*) ""

    call run_fortfem_intersection(pts, segs, fortfem_result, fortfem_mesh)
    call check_condition(fortfem_result%ntriangles > 0,                       &
        "FortFEM intersection CDT produces triangles")
    call check_condition(has_vertex_near(fortfem_result%points,              &
        1.0_dp, 1.0_dp, 1.0e-6_dp),                                           &
        "FortFEM intersection CDT has Steiner vertex near (1,1)")

    call run_triangle_intersection(pts, segs, triangle_mesh, triangle_ok)
    if (triangle_ok) then
        call check_condition(triangle_mesh%data%n_triangles > 0,              &
            "Triangle intersection CDT produces triangles")
        call check_condition(has_vertex_near(triangle_mesh%data%vertices,     &
            1.0_dp, 1.0_dp, 1.0e-6_dp),                                      &
            "Triangle intersection CDT has Steiner vertex near (1,1)")
    else
        write(*,*) "   Triangle intersection run skipped (triangle not available)"
    end if

    write(*,*) ""
    write(*,*) "=== Writing intersection comparison plots ==="
    call plot(fortfem_mesh, filename="build/intersection_fortfem.png",        &
              title="FortFEM CDT with Segment Intersection")
    write(*,*) "   Saved: build/intersection_fortfem.png"

    if (triangle_ok) then
        call plot(triangle_mesh, filename="build/intersection_triangle.png",  &
                  title="Triangle CDT with Segment Intersection")
        write(*,*) "   Saved: build/intersection_triangle.png"
    end if

    call check_summary("Triangle Intersection Comparison")

contains

    subroutine run_fortfem_intersection(points, segments, result, mesh)
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: segments(:,:)
        type(triangulation_result_t), intent(out) :: result
        type(mesh_t), intent(out) :: mesh

        integer :: status

        call triangulate_fortran(points, segments, result, status)
        if (status /= 0) then
            write(*,*) "Warning: triangulate_fortran (intersection) status =", &
                status
        end if

        mesh = mesh_from_arrays(result%points, result%triangles)
    end subroutine run_fortfem_intersection

    subroutine run_triangle_intersection(points, segments, mesh, success)
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: segments(:,:)
        type(mesh_t), intent(out) :: mesh
        logical, intent(out) :: success

        character(len=*), parameter :: TRIANGLE_PATH = "/tmp/triangle_bin"
        character(len=*), parameter :: POLY_FILE = "/tmp/intersection_comparison"
        integer :: stat, nverts, nsegs

        success = .false.

        call ensure_triangle_available(TRIANGLE_PATH, stat)
        if (stat /= 0) return

        nverts = size(points, 2)
        nsegs  = size(segments, 2)

        call write_triangle_poly_file(POLY_FILE // ".poly", points, segments, &
                                      nverts, nsegs, stat)
        if (stat /= 0) then
            write(*,*) "   Warning: failed to write .poly for intersection test"
            return
        end if

        call run_triangle_binary(TRIANGLE_PATH, POLY_FILE, stat)
        if (stat /= 0) then
            write(*,*) "   Warning: Triangle execution failed for intersection test"
            return
        end if

        mesh = mesh_from_triangle_files(POLY_FILE)
        if (mesh%data%n_triangles == 0) then
            write(*,*) "   Warning: failed to read Triangle intersection mesh"
            return
        end if

        success = .true.
    end subroutine run_triangle_intersection

    subroutine run_triangle_binary(bin_path, poly_basename, exit_stat)
        character(len=*), intent(in) :: bin_path
        character(len=*), intent(in) :: poly_basename
        integer, intent(out) :: exit_stat

        character(len=512) :: cmd

        cmd = trim(bin_path) // " -pQ " // trim(poly_basename) // ".poly"
        call execute_command_line(trim(cmd), wait=.true., exitstat=exit_stat)
    end subroutine run_triangle_binary

    logical function has_vertex_near(points, x, y, tol) result(found)
        real(dp), intent(in) :: points(:,:)
        real(dp), intent(in) :: x, y, tol

        integer :: i, n
        real(dp) :: dx, dy, dist2

        n = size(points, 2)
        found = .false.

        do i = 1, n
            dx = points(1, i) - x
            dy = points(2, i) - y
            dist2 = dx*dx + dy*dy
            if (dist2 <= tol*tol) then
                found = .true.
                return
            end if
        end do
    end function has_vertex_near

end program test_triangle_intersection_comparison
