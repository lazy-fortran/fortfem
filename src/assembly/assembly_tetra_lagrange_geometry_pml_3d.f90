module fortfem_assembly_tetra_lagrange_geometry_pml_3d
    !! Geometry-to-global-chain wrapper for scalar curvilinear tetrahedral PML.
    !!
    !! Layer geometry remains caller-owned. This module only composes the
    !! fixed-active-set normal-frame stretch builder with the existing P1/
    !! arbitrary-degree scalar Helmholtz CSC assembly and its derivatives.
    use fortfem_curvilinear_pml_geometry, only: &
        build_curvilinear_normal_pml_element_stretch, &
        build_curvilinear_normal_pml_element_stretch_jvp, &
        build_curvilinear_normal_pml_element_stretch_vjp
    use fortfem_assembly_tetra_lagrange_curvilinear_pml_3d, only: &
        assemble_tetra_lagrange_curvilinear_pml_csc, &
        assemble_tetra_lagrange_curvilinear_pml_csc_jvp, &
        assemble_tetra_lagrange_curvilinear_pml_csc_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_z_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_tetra_lagrange_geometry_pml_csc
    public :: assemble_tetra_lagrange_geometry_pml_csc_jvp
    public :: assemble_tetra_lagrange_geometry_pml_csc_vjp

contains

    subroutine assemble_tetra_lagrange_geometry_pml_csc( &
            mesh_vertices, tetrahedra, degree, layer_origins, layer_normals, &
            layer_width, wave_number, sigma_max, polynomial_degree, matrix, &
            status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, polynomial_degree
        real(dp), intent(in) :: layer_origins(:, :), layer_normals(:, :)
        real(dp), intent(in) :: layer_width(:), wave_number, sigma_max
        type(csc_z_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: stretch(:, :, :)
        integer :: local_status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Geometry-generated scalar tetrahedral PML assembly failed")
        call build_curvilinear_normal_pml_element_stretch( &
            mesh_vertices, tetrahedra, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree, stretch, local_status)
        if (local_status /= 0) return
        call assemble_tetra_lagrange_curvilinear_pml_csc( &
            mesh_vertices, tetrahedra, degree, stretch, wave_number, matrix, status)
    end subroutine assemble_tetra_lagrange_geometry_pml_csc

    subroutine assemble_tetra_lagrange_geometry_pml_csc_jvp( &
            mesh_vertices, tetrahedra, degree, layer_origins, layer_normals, &
            layer_width, wave_number, sigma_max, polynomial_degree, &
            mesh_vertices_dot, origins_dot, normals_dot, width_dot, &
            wave_number_dot, sigma_max_dot, matrix_dot, status)
        real(dp), intent(in) :: mesh_vertices(:, :), mesh_vertices_dot(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, polynomial_degree
        real(dp), intent(in) :: layer_origins(:, :), layer_normals(:, :)
        real(dp), intent(in) :: layer_width(:), wave_number, sigma_max
        real(dp), intent(in) :: origins_dot(:, :), normals_dot(:, :)
        real(dp), intent(in) :: width_dot(:), wave_number_dot, sigma_max_dot
        type(csc_z_t), intent(out) :: matrix_dot
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: stretch(:, :, :), stretch_dot(:, :, :)
        integer :: local_status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Geometry-generated scalar tetrahedral PML JVP failed")
        call build_curvilinear_normal_pml_element_stretch( &
            mesh_vertices, tetrahedra, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree, stretch, local_status)
        if (local_status /= 0) return
        call build_curvilinear_normal_pml_element_stretch_jvp( &
            mesh_vertices, tetrahedra, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree, mesh_vertices_dot, &
            origins_dot, normals_dot, width_dot, wave_number_dot, sigma_max_dot, &
            stretch_dot, local_status)
        if (local_status /= 0) return
        call assemble_tetra_lagrange_curvilinear_pml_csc_jvp( &
            mesh_vertices, tetrahedra, degree, stretch, wave_number, &
            mesh_vertices_dot, stretch_dot, wave_number_dot, matrix_dot, status)
    end subroutine assemble_tetra_lagrange_geometry_pml_csc_jvp

    subroutine assemble_tetra_lagrange_geometry_pml_csc_vjp( &
            mesh_vertices, tetrahedra, degree, layer_origins, layer_normals, &
            layer_width, wave_number, sigma_max, polynomial_degree, &
            matrix_values_bar, mesh_vertices_bar, origins_bar, normals_bar, &
            width_bar, wave_number_bar, sigma_max_bar, status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, polynomial_degree
        real(dp), intent(in) :: layer_origins(:, :), layer_normals(:,:)
        real(dp), intent(in) :: layer_width(:), wave_number, sigma_max
        complex(dp), intent(in) :: matrix_values_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :), origins_bar(:, :)
        real(dp), intent(out) :: normals_bar(:, :), width_bar(:)
        real(dp), intent(out) :: wave_number_bar, sigma_max_bar
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: stretch(:, :, :), stretch_bar(:, :, :)
        real(dp), allocatable :: mesh_vertices_from_pml(:, :)
        real(dp), allocatable :: origins_from_pml(:, :), normals_from_pml(:, :)
        real(dp), allocatable :: width_from_pml(:)
        real(dp) :: wave_number_from_pml, sigma_max_from_pml
        integer :: local_status

        mesh_vertices_bar = 0.0_dp
        origins_bar = 0.0_dp
        normals_bar = 0.0_dp
        width_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        sigma_max_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Geometry-generated scalar tetrahedral PML VJP failed")
        call build_curvilinear_normal_pml_element_stretch( &
            mesh_vertices, tetrahedra, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree, stretch, local_status)
        if (local_status /= 0) return
        allocate(stretch_bar(3, 3, size(tetrahedra, 2)))
        call assemble_tetra_lagrange_curvilinear_pml_csc_vjp( &
            mesh_vertices, tetrahedra, degree, stretch, wave_number, &
            matrix_values_bar, mesh_vertices_bar, stretch_bar, wave_number_bar, &
            status)
        if (status%code /= 0) return
        allocate( &
            mesh_vertices_from_pml(size(mesh_vertices, 1), size(mesh_vertices, 2)), &
            origins_from_pml(size(layer_origins, 1), size(layer_origins, 2)), &
            normals_from_pml(size(layer_normals, 1), size(layer_normals, 2)), &
            width_from_pml(size(layer_width)))
        call build_curvilinear_normal_pml_element_stretch_vjp( &
            mesh_vertices, tetrahedra, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree, stretch_bar, &
            mesh_vertices_from_pml, origins_from_pml, normals_from_pml, &
            width_from_pml, wave_number_from_pml, sigma_max_from_pml, local_status)
        if (local_status /= 0) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Geometry-generated scalar tetrahedral PML geometry VJP failed")
            return
        end if
        mesh_vertices_bar = mesh_vertices_bar + mesh_vertices_from_pml
        origins_bar = origins_from_pml
        normals_bar = normals_from_pml
        width_bar = width_from_pml
        wave_number_bar = wave_number_bar + wave_number_from_pml
        sigma_max_bar = sigma_max_from_pml
        call status_set(status, 0, "")
    end subroutine assemble_tetra_lagrange_geometry_pml_csc_vjp

end module fortfem_assembly_tetra_lagrange_geometry_pml_3d
