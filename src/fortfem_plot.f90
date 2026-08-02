module fortfem_plot
    !! Canonical direct-import facade for solution and mesh visualization.
    !!
    !! Plotting implementations stay in the existing ``fortfem_api_plot*``
    !! modules so the compatibility umbrella and this migration surface share
    !! exactly one implementation.  New clients can import this module without
    !! depending on ``fortfem_api``.
    use fortfem_api_plot, only: plot, plot_function_scalar, plot_vector_function, &
        plot_mesh
    use fortfem_api_plot_interpolation, only: compute_scalar_plot_grid, &
        compute_vector_plot_grid
    use fortfem_api_plot_mesh, only: prepare_mesh_plot, save_mesh_figure, &
        build_mesh_edge_path
    implicit none
    private

    public :: plot
    public :: plot_function_scalar
    public :: plot_vector_function
    public :: plot_mesh
    public :: compute_scalar_plot_grid
    public :: compute_vector_plot_grid
    public :: prepare_mesh_plot
    public :: save_mesh_figure
    public :: build_mesh_edge_path

end module fortfem_plot
