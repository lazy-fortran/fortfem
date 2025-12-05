program test_triangle_comparison
    !> Compare FortFEM CDT with Triangle library on L-shape domain.
    !
    !  This test:
    !    1. Generates an L-shape mesh with FortFEM
    !    2. Downloads and compiles Triangle (if needed)
    !    3. Runs Triangle on the same boundary
    !    4. Plots both meshes for visual comparison
    !
    use fortfem_kinds, only: dp
    use fortfem_api, only: mesh_t, boundary_t, l_shape_boundary,              &
        mesh_from_boundary, mesh_from_triangle_files, plot
    use triangle_io, only: write_triangle_poly_file, ensure_triangle_available
    use check, only: check_condition, check_summary
    implicit none

    type(boundary_t) :: boundary
    type(mesh_t) :: fortfem_mesh, triangle_mesh
    logical :: triangle_ok
    integer :: stat

    write(*,*) "=== Triangle vs FortFEM Comparison Test ==="
    write(*,*) ""

    boundary = l_shape_boundary(1.0_dp, 24)

    write(*,*) "1. Generating FortFEM mesh..."
    fortfem_mesh = mesh_from_boundary(boundary, resolution=0.08_dp)
    write(*,'(A,I0,A,I0,A)') "   FortFEM: ", fortfem_mesh%data%n_vertices,    &
        " vertices, ", fortfem_mesh%data%n_triangles, " triangles"

    write(*,*) ""
    write(*,*) "2. Running Triangle library..."
    call run_triangle_on_boundary(boundary, triangle_mesh, triangle_ok)

    if (triangle_ok) then
        write(*,'(A,I0,A,I0,A)') "   Triangle: ",                             &
            triangle_mesh%data%n_vertices, " vertices, ",                      &
            triangle_mesh%data%n_triangles, " triangles"
    else
        write(*,*) "   Triangle: FAILED (not available or error)"
    end if

    write(*,*) ""
    write(*,*) "3. Generating comparison plots..."

    call plot(fortfem_mesh, filename="build/comparison_fortfem.png",          &
              title="FortFEM CDT")
    write(*,*) "   Saved: build/comparison_fortfem.png"

    if (triangle_ok) then
        call plot(triangle_mesh, filename="build/comparison_triangle.png",    &
                  title="Triangle Library")
        write(*,*) "   Saved: build/comparison_triangle.png"
    end if

    call check_condition(fortfem_mesh%data%n_triangles > 0,                   &
        "FortFEM produces triangles")

    if (triangle_ok) then
        call check_condition(triangle_mesh%data%n_triangles > 0,              &
            "Triangle produces triangles")
    end if

    call check_summary("Triangle Comparison")

contains

    subroutine run_triangle_on_boundary(boundary, mesh, success)
        !> Run Triangle on a boundary and return the resulting mesh.
        type(boundary_t), intent(in) :: boundary
        type(mesh_t), intent(out) :: mesh
        logical, intent(out) :: success

        character(len=*), parameter :: TRIANGLE_PATH = "/tmp/triangle_bin"
        character(len=*), parameter :: POLY_FILE = "/tmp/lshape_comparison"
        character(len=512) :: cmd
        integer :: stat, i, n
        integer, allocatable :: segments(:,:)

        success = .false.

        ! Ensure Triangle is available
        call ensure_triangle_available(TRIANGLE_PATH, stat)
        if (stat /= 0) return

        ! Create segments array from boundary (closed polygon)
        n = boundary%n_points
        allocate(segments(2, n))
        do i = 1, n
            segments(1, i) = i
            if (i < n) then
                segments(2, i) = i + 1
            else
                segments(2, i) = 1
            end if
        end do

        ! Write .poly file for Triangle
        call write_triangle_poly_file(POLY_FILE // ".poly", boundary%points,  &
                                      segments, n, n, stat)
        if (stat /= 0) then
            write(*,*) "   Warning: Failed to write .poly file"
            return
        end if

        ! Run Triangle with ONLY CDT (no quality refinement) for fair comparison
        ! -p: read PSLG from .poly file
        ! -Q: quiet mode
        ! Note: Without -q, Triangle does pure CDT like our implementation
        cmd = TRIANGLE_PATH // " -pQ " // POLY_FILE // ".poly"
        call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
        if (stat /= 0) then
            write(*,*) "   Warning: Triangle execution failed"
            return
        end if

        ! Load the resulting mesh
        mesh = mesh_from_triangle_files(POLY_FILE)
        if (mesh%data%n_triangles == 0) then
            write(*,*) "   Warning: Failed to read Triangle output"
            return
        end if

        success = .true.
    end subroutine run_triangle_on_boundary

end program test_triangle_comparison
