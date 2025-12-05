program test_lshape_hole_refinement_comparison
    !> Compare FortFEM vs Triangle on an L-shaped domain with a hole,
    !> including one level of uniform mesh refinement, and write plots.
    use fortfem_kinds, only : dp
    use fortfem_api,  only : mesh_t, boundary_t, l_shape_boundary,          &
                             mesh_from_boundary, mesh_from_arrays,          &
                             mesh_from_triangle_files, refine_uniform,      &
                             plot
    use triangulation_fortran, only : triangulation_result_t,               &
                                      triangulate_with_hole_fortran
    use triangle_io, only : write_triangle_poly_file, ensure_triangle_available
    use check, only : check_condition, check_summary
    implicit none

    type(boundary_t) :: boundary
    real(dp), allocatable :: points(:,:)
    integer, allocatable :: segments(:,:)
    real(dp) :: hole_point(2)

    type(triangulation_result_t) :: fortfem_result
    type(mesh_t) :: fortfem_mesh, fortfem_refined
    type(mesh_t) :: triangle_mesh, triangle_refined
    logical :: triangle_ok

    write(*,*) "=== L-shape-with-hole Refinement: FortFEM vs Triangle ==="
    write(*,*) ""

    call build_lshape_with_hole_pslg(boundary, points, segments, hole_point)

    call run_fortfem_with_hole(points, segments, hole_point,                 &
        fortfem_result, fortfem_mesh)
    call check_condition(fortfem_result%ntriangles > 0,                      &
        "FortFEM L-shape-with-hole produces triangles")

    fortfem_refined = refine_uniform(fortfem_mesh)

    call run_triangle_with_hole(points, segments, hole_point,                &
        triangle_mesh, triangle_ok)
    if (triangle_ok) then
        triangle_refined = refine_uniform(triangle_mesh)
    else
        write(*,*) "   Triangle run skipped (triangle not available)"
    end if

    write(*,*) ""
    write(*,*) "=== Writing comparison plots (base meshes) ==="
    call plot(fortfem_mesh,                                                  &
              filename="build/lshape_hole_fortfem.png",                      &
              title="FortFEM L-shape with Hole")
    write(*,*) "   Saved: build/lshape_hole_fortfem.png"

    if (triangle_ok) then
        call plot(triangle_mesh,                                             &
                  filename="build/lshape_hole_triangle.png",                 &
                  title="Triangle L-shape with Hole")
        write(*,*) "   Saved: build/lshape_hole_triangle.png"
    end if

    write(*,*) ""
    write(*,*) "=== Writing comparison plots (refined meshes) ==="
    call plot(fortfem_refined,                                               &
              filename="build/lshape_hole_fortfem_refined.png",              &
              title="FortFEM L-shape with Hole (refined)")
    write(*,*) "   Saved: build/lshape_hole_fortfem_refined.png"

    if (triangle_ok) then
        call plot(triangle_refined,                                          &
                  filename="build/lshape_hole_triangle_refined.png",         &
                  title="Triangle L-shape with Hole (refined)")
        write(*,*) "   Saved: build/lshape_hole_triangle_refined.png"
    end if

    call check_condition(fortfem_refined%data%n_triangles >                  &
                         fortfem_mesh%data%n_triangles,                      &
        "FortFEM refined mesh has more triangles")

    if (triangle_ok) then
        call check_condition(triangle_refined%data%n_triangles >             &
                             triangle_mesh%data%n_triangles,                 &
            "Triangle refined mesh has more triangles")
    end if

    call check_summary("L-shape Hole Refinement Comparison")

contains

    subroutine build_lshape_with_hole_pslg(boundary, pts, segs, hpt)
        type(boundary_t), intent(out) :: boundary
        real(dp), allocatable, intent(out) :: pts(:,:)
        integer, allocatable, intent(out) :: segs(:,:)
        real(dp), intent(out) :: hpt(2)

        integer :: n_outer, n_points, i, idx

        boundary = l_shape_boundary(1.0_dp, 24)

        n_outer = boundary%n_points
        n_points = n_outer

        allocate(pts(2, n_points))
        pts(:, :) = boundary%points(:, :)

        allocate(segs(2, n_outer))
        do i = 1, n_outer
            segs(1, i) = i
            if (i < n_outer) then
                segs(2, i) = i + 1
            else
                segs(2, i) = 1
            end if
        end do

        hpt(1) = 1.5_dp
        hpt(2) = 0.5_dp
    end subroutine build_lshape_with_hole_pslg

    subroutine run_fortfem_with_hole(pts, segs, hpt, result, mesh)
        real(dp), intent(in) :: pts(:,:)
        integer, intent(in) :: segs(:,:)
        real(dp), intent(in) :: hpt(2)
        type(triangulation_result_t), intent(out) :: result
        type(mesh_t), intent(out) :: mesh

        integer :: status
        real(dp) :: hole_pts(2,1)

        hole_pts(:,1) = hpt(:)
        call triangulate_with_hole_fortran(pts, segs, hole_pts,              &
                                           result, status)

        if (status /= 0) then
            write(*,*) "Warning: triangulate_with_hole_fortran returned ",   &
                status
        end if

        mesh = mesh_from_arrays(result%points, result%triangles)
    end subroutine run_fortfem_with_hole

    subroutine run_triangle_with_hole(pts, segs, hpt, mesh, success)
        real(dp), intent(in) :: pts(:,:)
        integer, intent(in) :: segs(:,:)
        real(dp), intent(in) :: hpt(2)
        type(mesh_t), intent(out) :: mesh
        logical, intent(out) :: success

        character(len=*), parameter :: TRIANGLE_PATH = "/tmp/triangle_bin"
        character(len=*), parameter :: POLY_FILE =                               &
            "build/lshape_hole_comparison"
        real(dp) :: hole_points(2,1)
        integer :: n, stat

        success = .false.

        call ensure_triangle_available(TRIANGLE_PATH, stat)
        if (stat /= 0) return

        n = size(pts, 2)
        hole_points(:,1) = hpt(:)

        call write_triangle_poly_file(POLY_FILE // ".poly", pts, segs,       &
                                      n, size(segs,2), stat, hole_points)
        if (stat /= 0) then
            write(*,*) "   Warning: Failed to write .poly file for L-shape hole"
            return
        end if

        call run_triangle_binary(TRIANGLE_PATH, POLY_FILE, stat)
        if (stat /= 0) then
            write(*,*) "   Warning: Triangle execution for L-shape hole failed"
            return
        end if

        mesh = mesh_from_triangle_files(POLY_FILE)
        if (mesh%data%n_triangles == 0) then
            write(*,*) "   Warning: Failed to read Triangle L-shape hole mesh"
            return
        end if

        success = .true.
    end subroutine run_triangle_with_hole

    subroutine run_triangle_binary(bin_path, poly_basename, exit_stat)
        character(len=*), intent(in) :: bin_path
        character(len=*), intent(in) :: poly_basename
        integer, intent(out) :: exit_stat

        character(len=512) :: cmd

        cmd = trim(bin_path) // " -pQ " // trim(poly_basename) // ".poly"
        call execute_command_line(trim(cmd), wait=.true., exitstat=exit_stat)
    end subroutine run_triangle_binary

end program test_lshape_hole_refinement_comparison

