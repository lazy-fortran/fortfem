program test_triangle_hole_comparison
    !> Compare FortFEM CDT-with-holes against Triangle on a square-with-hole domain.
    !
    !  This test:
    !    1. Builds a square-with-hole PSLG (outer and inner squares)
    !    2. Runs FortFEM constrained Delaunay with a hole seed
    !    3. Runs Triangle on the same PSLG with a hole seed
    !    4. Writes comparison plots for visual inspection
    !    5. Checks that both meshes have no triangles whose centroid lies inside the hole
    !
    use fortfem_kinds, only : dp
    use fortfem_api,  only : mesh_t, mesh_from_arrays, mesh_from_triangle_files, plot
    use triangulation_fortran, only : triangulation_result_t,                 &
                                      triangulate_with_hole_fortran
    use triangle_io, only : write_triangle_poly_file, ensure_triangle_available
    use check, only : check_condition, check_summary
    implicit none

    real(dp), parameter :: points(2,8) = reshape([ &
        0.0_dp, 0.0_dp, &  ! outer square
        2.0_dp, 0.0_dp, &
        2.0_dp, 2.0_dp, &
        0.0_dp, 2.0_dp, &
        0.5_dp, 0.5_dp, &  ! inner square (hole)
        1.5_dp, 0.5_dp, &
        1.5_dp, 1.5_dp, &
        0.5_dp, 1.5_dp ], [2, 8])

    integer, parameter :: segments(2,8) = reshape([ &
        1, 2, &  ! outer segments
        2, 3, &
        3, 4, &
        4, 1, &
        5, 6, &  ! inner segments
        6, 7, &
        7, 8, &
        8, 5 ], [2, 8])

    real(dp), parameter :: hole_point(2) = [1.0_dp, 1.0_dp]

    type(triangulation_result_t) :: fortfem_result
    type(mesh_t) :: fortfem_mesh, triangle_mesh
    logical :: triangle_ok

    write(*,*) "=== Triangle vs FortFEM Hole Comparison Test ==="
    write(*,*) ""

    call run_fortfem_with_hole(points, segments, hole_point, fortfem_result,  &
        fortfem_mesh)
    call check_condition(fortfem_result%ntriangles > 0,                       &
        "FortFEM-with-hole produces triangles")
    call check_condition(all_outside_hole(fortfem_result%points,             &
        fortfem_result%triangles),                                            &
        "FortFEM triangles avoid interior of hole")

    call run_triangle_with_hole(points, segments, hole_point, triangle_mesh,  &
        triangle_ok)
    if (triangle_ok) then
        call check_condition(triangle_mesh%data%n_triangles > 0,              &
            "Triangle-with-hole produces triangles")
        call check_condition(all_outside_hole(triangle_mesh%data%vertices,    &
            triangle_mesh%data%triangles),                                    &
            "Triangle triangles avoid interior of hole")
    else
        write(*,*) "   Triangle-with-hole run skipped (triangle not available)"
    end if

    write(*,*) ""
    write(*,*) "=== Writing comparison plots ==="
    call plot(fortfem_mesh, filename="build/hole_fortfem.png",                &
              title="FortFEM CDT with Hole")
    write(*,*) "   Saved: build/hole_fortfem.png"

    if (triangle_ok) then
        call plot(triangle_mesh, filename="build/hole_triangle.png",          &
                  title="Triangle CDT with Hole")
        write(*,*) "   Saved: build/hole_triangle.png"
    end if

    call check_summary("Triangle Hole Comparison")

contains

    subroutine run_fortfem_with_hole(pts, segs, hpt, result, mesh)
        real(dp), intent(in) :: pts(:,:)
        integer, intent(in) :: segs(:,:)
        real(dp), intent(in) :: hpt(:)
        type(triangulation_result_t), intent(out) :: result
        type(mesh_t), intent(out) :: mesh

        integer :: status
        real(dp) :: hole_pts(2,1)

        hole_pts(:,1) = hpt(:)
        call triangulate_with_hole_fortran(pts, segs, hole_pts, result, status)

        if (status /= 0) then
            write(*,*) "Warning: triangulate_with_hole_fortran returned ",    &
                status
        end if

        mesh = mesh_from_arrays(result%points, result%triangles)
    end subroutine run_fortfem_with_hole

    subroutine run_triangle_with_hole(pts, segs, hpt, mesh, success)
        real(dp), intent(in) :: pts(:,:)
        integer, intent(in) :: segs(:,:)
        real(dp), intent(in) :: hpt(:)
        type(mesh_t), intent(out) :: mesh
        logical, intent(out) :: success

        character(len=*), parameter :: TRIANGLE_PATH = "/tmp/triangle_bin"
        character(len=*), parameter :: POLY_FILE = "/tmp/square_hole_comparison"
        real(dp) :: hole_points(2,1)
        integer :: n, i, stat

        success = .false.

        call ensure_triangle_available(TRIANGLE_PATH, stat)
        if (stat /= 0) return

        n = size(pts, 2)
        hole_points(:,1) = hpt(:)

        call write_triangle_poly_file(POLY_FILE // ".poly", pts, segs,        &
                                      n, size(segs,2), stat, hole_points)
        if (stat /= 0) then
            write(*,*) "   Warning: Failed to write .poly file with hole"
            return
        end if

        call run_triangle_binary(TRIANGLE_PATH, POLY_FILE, stat)
        if (stat /= 0) then
            write(*,*) "   Warning: Triangle execution with hole failed"
            return
        end if

        mesh = mesh_from_triangle_files(POLY_FILE)
        if (mesh%data%n_triangles == 0) then
            write(*,*) "   Warning: Failed to read Triangle hole mesh"
            return
        end if

        success = .true.
    end subroutine run_triangle_with_hole

    subroutine run_triangle_binary(bin_path, poly_basename, exit_stat)
        character(len=*), intent(in) :: bin_path
        character(len=*), intent(in) :: poly_basename
        integer, intent(out) :: exit_stat

        character(len=512) :: cmd

        ! -p: read PSLG from .poly, -Q: quiet, no quality refinement (-q)
        cmd = trim(bin_path) // " -pQ " // trim(poly_basename) // ".poly"
        call execute_command_line(trim(cmd), wait=.true., exitstat=exit_stat)
    end subroutine run_triangle_binary

    logical function all_outside_hole(pts, tris) result(ok)
        real(dp), intent(in) :: pts(:,:)
        integer, intent(in) :: tris(:,:)

        integer :: i
        real(dp) :: cx, cy

        ok = .true.

        do i = 1, size(tris, 2)
            call triangle_centroid(pts, tris(:,i), cx, cy)

            ! Hole is inner square with corners (0.5,0.5) to (1.5,1.5)
            if (cx > 0.5_dp .and. cx < 1.5_dp .and.                            &
                cy > 0.5_dp .and. cy < 1.5_dp) then
                ok = .false.
                return
            end if
        end do
    end function all_outside_hole

    subroutine triangle_centroid(pts, tri, cx, cy)
        real(dp), intent(in) :: pts(:,:)
        integer, intent(in) :: tri(3)
        real(dp), intent(out) :: cx, cy

        cx = (pts(1, tri(1)) + pts(1, tri(2)) + pts(1, tri(3))) / 3.0_dp
        cy = (pts(2, tri(1)) + pts(2, tri(2)) + pts(2, tri(3))) / 3.0_dp
    end subroutine triangle_centroid

end program test_triangle_hole_comparison
