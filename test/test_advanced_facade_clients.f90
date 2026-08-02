program test_advanced_facade_clients
    !! Compile-and-value oracle for the advanced canonical capability surfaces.
    !!
    !! The gallery programs exercise the full numerical paths.  This small
    !! client keeps the facade contract independently testable without tying
    !! the release gate to a large toroidal solve or to plot generation.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        apply_maxwell_trace_to_flux_map, &
        assemble_laplace_single_layer_p0_3d
    use fortfem_core, only: &
        generate_sphere_surface_mesh, generate_torus_surface_mesh, &
        invert_tetra_affine_map
    use fortfem_feec, only: build_bspline_polar_feec_2d_operators
    use fortfem_kinds, only: dp
    implicit none

    call check_geometry_facade()
    call check_polar_feec_facade()
    call check_boundary_facade()
    call check_summary("advanced canonical facade clients")

contains

    subroutine check_geometry_facade()
        integer, allocatable :: sphere_triangles(:, :), torus_triangles(:, :)
        real(dp), allocatable :: sphere_vertices(:, :), torus_vertices(:, :)
        real(dp) :: tetra_vertices(3, 4), point(3), reference(3)
        integer :: status

        tetra_vertices = reshape([ &
            0.0_dp, 0.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
        point = [0.2_dp, 0.3_dp, 0.1_dp]
        call invert_tetra_affine_map( &
            tetra_vertices, point, reference, status)
        call check_condition(status == 0 .and. &
            maxval(abs(reference - point)) < 2.0e-14_dp, &
            "core facade inverts an affine tetrahedron analytically")

        call generate_sphere_surface_mesh(1.0_dp, 0, sphere_vertices, &
            sphere_triangles)
        call generate_torus_surface_mesh(2.0_dp, 0.5_dp, 3, 4, &
            torus_vertices, torus_triangles)
        call check_condition(size(sphere_vertices, 2) == 6 .and. &
            size(sphere_triangles, 2) == 8 .and. &
            size(torus_vertices, 2) == 12 .and. &
            size(torus_triangles, 2) == 24, &
            "core facade exposes sphere and torus mesh constructors")
    end subroutine check_geometry_facade

    subroutine check_polar_feec_facade()
        real(dp), allocatable :: gradient(:, :), curl(:, :)
        integer :: status

        call build_bspline_polar_feec_2d_operators( &
            6, 5, gradient, curl, status)
        call check_condition(status == 0 .and. size(gradient, 1) == 38 .and. &
            size(curl, 1) == 18 .and. &
            maxval(abs(matmul(curl, gradient))) < 2.0e-14_dp, &
            "FEEC facade preserves the polar spline curl-gradient complex")
    end subroutine check_polar_feec_facade

    subroutine check_boundary_facade()
        real(dp), allocatable :: vertices(:, :), matrix(:, :)
        integer, allocatable :: triangles(:, :)
        complex(dp) :: map(2, 2), trace(2), flux(2)
        integer :: status

        call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
        call assemble_laplace_single_layer_p0_3d( &
            vertices, triangles, 4, matrix, status)
        call check_condition(status == 0 .and. &
            size(matrix, 1) == size(triangles, 2) .and. &
            maxval(abs(matrix - transpose(matrix))) < 2.0e-13_dp, &
            "boundary facade assembles a reciprocal 3-D Laplace layer")

        map = reshape([ &
            cmplx(2.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 1.0_dp, dp), &
            cmplx(-1.0_dp, 0.0_dp, dp), cmplx(0.5_dp, 0.0_dp, dp)], [2, 2])
        trace = [cmplx(1.0_dp, 0.0_dp, dp), cmplx(2.0_dp, -1.0_dp, dp)]
        call apply_maxwell_trace_to_flux_map(map, trace, flux, status)
        call check_condition(status == 0 .and. &
            maxval(abs(flux - matmul(map, trace))) < 2.0e-14_dp, &
            "boundary facade applies a Maxwell trace-to-flux map")
    end subroutine check_boundary_facade

end program test_advanced_facade_clients
