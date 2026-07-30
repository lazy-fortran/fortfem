module fortfem_helmholtz_torus_curved_fem_bem_coupling_3d
    !! Symmetric Costabel-Han Helmholtz coupling on a tetrahedral solid torus
    !! with exact-parametric curved boundary operators.
    !!
    !! Block signs follow the Helmholtz analogue of equation (28) in
    !! Aurada et al., arXiv:1211.4225.
    use fortfem_assembly_tetra_lagrange_arbitrary_order_3d, only: &
        assemble_tetra_lagrange_stiffness_element
    use fortfem_helmholtz_torus_curved_bem_3d, only: &
        assemble_helmholtz_torus_curved_calderon_3d
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    implicit none
    private

    public :: assemble_helmholtz_fem_bem_costabel_torus_curved_3d
    public :: solve_helmholtz_fem_bem_costabel_torus_curved_3d

contains

    subroutine solve_helmholtz_fem_bem_costabel_torus_curved_3d( &
            vertices, tetrahedra, parameters, boundary_triangles, &
            major_radius, minor_radius, interior_wave_number, &
            exterior_wave_number, volume_load, quadrature_degree, &
            potential, normal_flux, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        integer, intent(in) :: tetrahedra(:, :), boundary_triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: interior_wave_number, exterior_wave_number
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: quadrature_degree
        complex(dp), allocatable, intent(out) :: potential(:), normal_flux(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix(:, :), right_hand_side(:)
        complex(dp), allocatable :: solution(:)
        integer :: info, total_dof_count, vertex_count

        status = 1
        if (allocated(potential)) deallocate(potential)
        if (allocated(normal_flux)) deallocate(normal_flux)
        if (size(volume_load) /= size(vertices, 2)) return
        call assemble_helmholtz_fem_bem_costabel_torus_curved_3d( &
            vertices, tetrahedra, parameters, boundary_triangles, &
            major_radius, minor_radius, interior_wave_number, &
            exterior_wave_number, quadrature_degree, matrix, status)
        if (status /= 0) return
        total_dof_count = size(matrix, 1)
        vertex_count = size(vertices, 2)
        allocate(right_hand_side(total_dof_count), solution(total_dof_count))
        right_hand_side = cmplx(0.0_dp, 0.0_dp, dp)
        right_hand_side(:vertex_count) = volume_load
        call dense_solve(matrix, right_hand_side, solution, info)
        if (info /= 0) then
            status = 2
            return
        end if
        allocate(potential(vertex_count))
        allocate(normal_flux(total_dof_count - vertex_count))
        potential = solution(:vertex_count)
        normal_flux = solution(vertex_count + 1:)
        status = 0
    end subroutine solve_helmholtz_fem_bem_costabel_torus_curved_3d

    subroutine assemble_helmholtz_fem_bem_costabel_torus_curved_3d( &
            vertices, tetrahedra, parameters, boundary_triangles, &
            major_radius, minor_radius, interior_wave_number, &
            exterior_wave_number, quadrature_degree, matrix, status)
        real(dp), intent(in) :: vertices(:, :), parameters(:, :)
        integer, intent(in) :: tetrahedra(:, :), boundary_triangles(:, :)
        real(dp), intent(in) :: major_radius, minor_radius
        real(dp), intent(in) :: interior_wave_number, exterior_wave_number
        integer, intent(in) :: quadrature_degree
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        integer, parameter :: basis_to_vertex(4) = [4, 3, 2, 1]
        complex(dp), allocatable :: adjoint(:, :), double_layer(:, :)
        complex(dp), allocatable :: hypersingular(:, :), single_layer(:, :)
        real(dp), allocatable :: element_matrix(:, :), mass(:, :)
        real(dp) :: tetra_vertices(3, 4)
        integer :: column, local_status, row, tetrahedron
        integer :: total_dof_count, triangle_count, vertex_count

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        if (interior_wave_number < 0.0_dp) return
        if (exterior_wave_number <= 0.0_dp) return
        if (size(vertices, 1) /= 3 .or. size(vertices, 2) < 4) return
        if (size(parameters, 1) /= 2) return
        if (size(parameters, 2) /= size(vertices, 2)) return
        if (size(tetrahedra, 1) /= 4 .or. size(tetrahedra, 2) < 1) return
        if (size(boundary_triangles, 1) /= 3 .or. &
            size(boundary_triangles, 2) < 1) return
        if (any(tetrahedra < 1)) return
        if (any(tetrahedra > size(vertices, 2))) return
        if (any(boundary_triangles < 1)) return
        if (any(boundary_triangles > size(vertices, 2))) return
        call assemble_helmholtz_torus_curved_calderon_3d( &
            parameters, boundary_triangles, major_radius, minor_radius, &
            exterior_wave_number, quadrature_degree, single_layer, &
            double_layer, adjoint, hypersingular, mass, local_status)
        if (local_status /= 0) return

        vertex_count = size(vertices, 2)
        triangle_count = size(boundary_triangles, 2)
        total_dof_count = vertex_count + triangle_count
        allocate(matrix(total_dof_count, total_dof_count))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        do tetrahedron = 1, size(tetrahedra, 2)
            tetra_vertices = vertices(:, tetrahedra(:, tetrahedron))
            call assemble_tetra_lagrange_stiffness_element( &
                tetra_vertices, 1, 2, element_matrix, local_status, &
                stiffness_coefficient=1.0_dp, &
                mass_coefficient=-interior_wave_number**2)
            if (local_status /= 0) return
            do column = 1, 4
                do row = 1, 4
                    matrix( &
                        tetrahedra(basis_to_vertex(row), tetrahedron), &
                        tetrahedra(basis_to_vertex(column), tetrahedron)) = &
                        matrix( &
                        tetrahedra(basis_to_vertex(row), tetrahedron), &
                        tetrahedra(basis_to_vertex(column), tetrahedron)) + &
                        element_matrix(row, column)
                end do
            end do
        end do
        matrix(:vertex_count, :vertex_count) = &
            matrix(:vertex_count, :vertex_count) + hypersingular
        matrix(:vertex_count, vertex_count + 1:) = &
            adjoint - 0.5_dp*cmplx(transpose(mass), 0.0_dp, dp)
        matrix(vertex_count + 1:, :vertex_count) = &
            double_layer - 0.5_dp*cmplx(mass, 0.0_dp, dp)
        matrix(vertex_count + 1:, vertex_count + 1:) = -single_layer
        status = 0
    end subroutine assemble_helmholtz_fem_bem_costabel_torus_curved_3d

end module fortfem_helmholtz_torus_curved_fem_bem_coupling_3d
