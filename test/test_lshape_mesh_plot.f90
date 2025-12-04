program test_lshape_mesh_plot
    use fortfem_kinds, only: dp
    use fortfem_api, only: mesh_t, boundary_t, l_shape_boundary,          &
        mesh_from_boundary, plot
    use check, only: check_condition, check_summary
    implicit none

    type(boundary_t) :: boundary
    type(mesh_t) :: mesh
    integer :: n_boundary_vertices

    write(*,*) "Testing L-shape boundary Delaunay mesh and plotting..."

    boundary = l_shape_boundary(1.0_dp, 24)
    mesh = mesh_from_boundary(boundary, resolution=0.08_dp)

    n_boundary_vertices = count(mesh%data%is_boundary_vertex)

    call check_condition(mesh%data%n_vertices > 0, &
        "L-shape mesh: positive vertex count")
    call check_condition(mesh%data%n_triangles > 0, &
        "L-shape mesh: positive triangle count")
    call check_condition(mesh%data%n_quads == 0, &
        "L-shape mesh: pure triangular Delaunay mesh")
    call check_condition(n_boundary_vertices > 0, &
        "L-shape mesh: boundary vertices identified")

    call plot(mesh, filename="build/test_lshape_boundary_mesh.png", &
              title="L-shape Delaunay mesh")

    write(*,*) "   L-shape Delaunay plot: generated build/test_lshape_boundary_mesh.png"

    call check_summary("L-shape Delaunay Mesh Plot")

end program test_lshape_mesh_plot

