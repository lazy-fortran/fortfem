program test_triangle_quality_comparison
    !> Compare FortFEM Delaunay refinement against Triangle on a simple PSLG.
    !
    !  Geometry:
    !    - Outer square: (0,0)-(2,0)-(2,2)-(0,2)
    !    - One interior point near the lower edge to create skinny triangles
    !
    !  This test:
    !    1. Builds a simple PSLG (scaled unit square).
    !    2. Runs plain FortFEM CDT and measures minimum triangle angle.
    !    3. Runs FortFEM quality CDT via triangulate_with_quality_fortran
    !       (currently a quality-aware wrapper around triangulate_fortran).
    !    4. Runs Triangle with quality options on the same PSLG.
    !    5. Writes comparison plots for visual inspection.
    !
    use fortfem_kinds, only : dp
    use fortfem_api,  only : mesh_t, mesh_from_arrays, mesh_from_triangle_files, &
                             plot
    use triangulation_fortran, only : triangulation_result_t,                   &
                                      triangulate_fortran,                      &
                                      triangulate_with_quality_fortran
    use triangle_io, only : write_triangle_poly_file, ensure_triangle_available
    use check, only : check_condition, check_summary
    implicit none

    real(dp), parameter :: points(2,5) = reshape([ &
        0.0_dp, 0.0_dp, &  ! 1: outer square
        2.0_dp, 0.0_dp, &  ! 2
        2.0_dp, 2.0_dp, &  ! 3
        0.0_dp, 2.0_dp, &  ! 4
        0.5_dp, 0.15_dp ], [2, 5])  ! 5: interior point near lower edge

    integer, parameter :: segments(2,4) = reshape([ &
        1, 2, &  ! outer boundary
        2, 3, &
        3, 4, &
        4, 1 ], [2, 4])

    real(dp), parameter :: target_min_angle = 20.0_dp

    type(triangulation_result_t) :: base_result, quality_result
    type(mesh_t) :: base_mesh, quality_mesh, triangle_mesh
    real(dp) :: base_min_angle, quality_min_angle, triangle_min_angle
    integer :: status_quality
    logical :: triangle_ok

    write(*,*) "=== Triangle vs FortFEM Quality Comparison Test ==="
    write(*,*) ""

    call run_fortfem_base(points, segments, base_result, base_mesh,           &
        base_min_angle)
    call check_condition(base_result%ntriangles > 0,                          &
        "Base CDT: has triangles")

    call run_fortfem_quality(points, segments, target_min_angle,             &
        quality_result, quality_mesh, quality_min_angle, status_quality)

    call check_condition(status_quality == 0,                                 &
        "Quality CDT: status indicates acceptable mesh")
    call check_condition(quality_result%ntriangles >= base_result%ntriangles, &
        "Quality CDT: refinement does not reduce triangle count")
    call check_condition(quality_min_angle >= base_min_angle,                &
        "Quality CDT: minimum angle improves")
    call check_condition(quality_min_angle >= 0.9_dp * target_min_angle,     &
        "Quality CDT: minimum angle near target")

    call run_triangle_quality(points, segments, target_min_angle,            &
        triangle_mesh, triangle_min_angle, triangle_ok)

    if (triangle_ok) then
        call check_condition(triangle_min_angle >= 0.9_dp * target_min_angle,&
            "Triangle quality CDT: minimum angle near target")
        call check_condition(quality_min_angle >= 0.8_dp * triangle_min_angle,&
            "FortFEM quality CDT: minimum angle comparable to Triangle")
    else
        write(*,*) "   Triangle quality run skipped (triangle not available)"
    end if

    write(*,'(A,F8.3)') "   Base CDT min angle (deg):       ", base_min_angle
    write(*,'(A,F8.3)') "   FortFEM quality CDT min angle:  ", quality_min_angle
    if (triangle_ok) then
        write(*,'(A,F8.3)') "   Triangle quality CDT min angle:", triangle_min_angle
    end if

    write(*,*) ""
    write(*,*) "=== Writing quality comparison plots ==="
    call plot(base_mesh, filename="build/quality_base_fortfem.png",          &
              title="FortFEM Base CDT")
    write(*,*) "   Saved: build/quality_base_fortfem.png"

    call plot(quality_mesh, filename="build/quality_refined_fortfem.png",    &
              title="FortFEM Quality CDT")
    write(*,*) "   Saved: build/quality_refined_fortfem.png"

    if (triangle_ok) then
        call plot(triangle_mesh, filename="build/quality_triangle.png",      &
                  title="Triangle Quality CDT")
        write(*,*) "   Saved: build/quality_triangle.png"
    end if

    call check_summary("Triangle Quality Comparison")

contains

    subroutine run_fortfem_base(pts, segs, result, mesh, min_angle_deg)
        real(dp), intent(in) :: pts(:,:)
        integer, intent(in) :: segs(:,:)
        type(triangulation_result_t), intent(out) :: result
        type(mesh_t), intent(out) :: mesh
        real(dp), intent(out) :: min_angle_deg

        integer :: status_local

        call triangulate_fortran(pts, segs, result, status_local)
        mesh = mesh_from_arrays(result%points, result%triangles)
        min_angle_deg = compute_min_result_angle(result)
    end subroutine run_fortfem_base

    subroutine run_fortfem_quality(pts, segs, min_angle_target, result,      &
                                   mesh, min_angle_deg, status)
        real(dp), intent(in) :: pts(:,:)
        integer, intent(in) :: segs(:,:)
        real(dp), intent(in) :: min_angle_target
        type(triangulation_result_t), intent(out) :: result
        type(mesh_t), intent(out) :: mesh
        real(dp), intent(out) :: min_angle_deg
        integer, intent(out) :: status

        call triangulate_with_quality_fortran(pts, segs, min_angle_target,   &
                                             result, status)
        mesh = mesh_from_arrays(result%points, result%triangles)
        min_angle_deg = compute_min_result_angle(result)
    end subroutine run_fortfem_quality

    subroutine run_triangle_quality(pts, segs, min_angle_target, mesh,       &
                                    min_angle_deg, success)
        real(dp), intent(in) :: pts(:,:)
        integer, intent(in) :: segs(:,:)
        real(dp), intent(in) :: min_angle_target
        type(mesh_t), intent(out) :: mesh
        real(dp), intent(out) :: min_angle_deg
        logical, intent(out) :: success

        character(len=*), parameter :: TRIANGLE_PATH = "/tmp/triangle_bin"
        character(len=*), parameter :: POLY_FILE = "/tmp/quality_comparison"
        integer :: nverts, nsegs, stat
        real(dp) :: min_angle_measured

        success = .false.
        min_angle_deg = 0.0_dp

        call ensure_triangle_available(TRIANGLE_PATH, stat)
        if (stat /= 0) return

        nverts = size(pts, 2)
        nsegs  = size(segs, 2)

        call write_triangle_poly_file(POLY_FILE // ".poly", pts, segs,       &
                                      nverts, nsegs, stat)
        if (stat /= 0) then
            write(*,*) "   Warning: failed to write .poly for quality test"
            return
        end if

        call run_triangle_binary_quality(TRIANGLE_PATH, POLY_FILE,           &
                                         min_angle_target, stat)
        if (stat /= 0) then
            write(*,*) "   Warning: Triangle execution failed for quality test"
            return
        end if

        mesh = mesh_from_triangle_files(POLY_FILE)
        if (mesh%data%n_triangles == 0) then
            write(*,*) "   Warning: failed to read Triangle quality mesh"
            return
        end if

        min_angle_measured = compute_min_mesh_angle(mesh)
        min_angle_deg = min_angle_measured

        success = .true.
    end subroutine run_triangle_quality

    subroutine run_triangle_binary_quality(bin_path, poly_basename,          &
                                           min_angle_deg, exit_stat)
        character(len=*), intent(in) :: bin_path
        character(len=*), intent(in) :: poly_basename
        real(dp), intent(in) :: min_angle_deg
        integer, intent(out) :: exit_stat

        character(len=512) :: cmd

        ! -p: PSLG input, -qX: quality with minimum angle X degrees, -Q: quiet
        write(cmd, '(A,F6.2,A)') " -pq", min_angle_deg, " "
        cmd = trim(bin_path) // trim(cmd) // " " // trim(poly_basename) // ".poly"
        call execute_command_line(trim(cmd), wait=.true., exitstat=exit_stat)
    end subroutine run_triangle_binary_quality

    real(dp) function compute_min_result_angle(result) result(min_angle_deg)
        type(triangulation_result_t), intent(in) :: result

        integer :: i
        real(dp) :: ax, ay, bx, by, cx, cy
        real(dp) :: angle1, angle2, angle3, tri_min
        real(dp), parameter :: pi = acos(-1.0_dp)

        min_angle_deg = 180.0_dp

        do i = 1, result%ntriangles
            ax = result%points(1, result%triangles(1, i))
            ay = result%points(2, result%triangles(1, i))
            bx = result%points(1, result%triangles(2, i))
            by = result%points(2, result%triangles(2, i))
            cx = result%points(1, result%triangles(3, i))
            cy = result%points(2, result%triangles(3, i))

            call triangle_angles_from_coords(ax, ay, bx, by, cx, cy,         &
                                            angle1, angle2, angle3, pi)

            tri_min = min(angle1, min(angle2, angle3))
            if (tri_min < min_angle_deg) min_angle_deg = tri_min
        end do
    end function compute_min_result_angle

    real(dp) function compute_min_mesh_angle(mesh) result(min_angle_deg)
        type(mesh_t), intent(in) :: mesh

        integer :: i, v1, v2, v3
        real(dp) :: ax, ay, bx, by, cx, cy
        real(dp) :: angle1, angle2, angle3, tri_min
        real(dp), parameter :: pi = acos(-1.0_dp)

        min_angle_deg = 180.0_dp

        do i = 1, mesh%data%n_triangles
            v1 = mesh%data%triangles(1, i)
            v2 = mesh%data%triangles(2, i)
            v3 = mesh%data%triangles(3, i)

            ax = mesh%data%vertices(1, v1)
            ay = mesh%data%vertices(2, v1)
            bx = mesh%data%vertices(1, v2)
            by = mesh%data%vertices(2, v2)
            cx = mesh%data%vertices(1, v3)
            cy = mesh%data%vertices(2, v3)

            call triangle_angles_from_coords(ax, ay, bx, by, cx, cy,         &
                                            angle1, angle2, angle3, pi)

            tri_min = min(angle1, min(angle2, angle3))
            if (tri_min < min_angle_deg) min_angle_deg = tri_min
        end do
    end function compute_min_mesh_angle

    subroutine triangle_angles_from_coords(ax, ay, bx, by, cx, cy,           &
                                           angle1, angle2, angle3, pi)
        real(dp), intent(in) :: ax, ay, bx, by, cx, cy, pi
        real(dp), intent(out) :: angle1, angle2, angle3

        real(dp) :: a, b, c
        real(dp) :: cos_a, cos_b, cos_c

        a = sqrt((bx - cx)*(bx - cx) + (by - cy)*(by - cy))
        b = sqrt((ax - cx)*(ax - cx) + (ay - cy)*(ay - cy))
        c = sqrt((ax - bx)*(ax - bx) + (ay - by)*(ay - by))

        if (a <= 0.0_dp .or. b <= 0.0_dp .or. c <= 0.0_dp) then
            angle1 = 0.0_dp
            angle2 = 0.0_dp
            angle3 = 0.0_dp
            return
        end if

        cos_a = (b*b + c*c - a*a) / (2.0_dp * b * c)
        cos_b = (a*a + c*c - b*b) / (2.0_dp * a * c)
        cos_c = (a*a + b*b - c*c) / (2.0_dp * a * b)

        cos_a = max(-1.0_dp, min(1.0_dp, cos_a))
        cos_b = max(-1.0_dp, min(1.0_dp, cos_b))
        cos_c = max(-1.0_dp, min(1.0_dp, cos_c))

        angle1 = acos(cos_a) * 180.0_dp / pi
        angle2 = acos(cos_b) * 180.0_dp / pi
        angle3 = acos(cos_c) * 180.0_dp / pi
    end subroutine triangle_angles_from_coords

end program test_triangle_quality_comparison
