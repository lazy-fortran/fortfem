program test_basic_fem_facade
    use check, only: check_condition, check_summary
    use fortfem_core, only: function_space_t, mesh_t, unit_square_mesh
    use fortfem_feec, only: function_space
    implicit none

    type(mesh_t) :: mesh
    type(function_space_t) :: space

    mesh = unit_square_mesh(2)
    space = function_space(mesh, "Lagrange", 1)

    call check_condition(mesh%data%n_vertices == 4, &
        "core facade builds the analytical unit-square mesh")
    call check_condition(space%ndof == mesh%data%n_vertices, &
        "FEEC facade constructs the P1 space without the umbrella")
    call check_summary("basic FEM canonical facade")
end program test_basic_fem_facade
