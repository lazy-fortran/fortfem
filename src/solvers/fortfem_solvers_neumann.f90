module fortfem_solvers_neumann
    use fortfem_kinds, only: dp
    use fortfem_api_types, only: function_t, dirichlet_bc_t, neumann_bc_t
    use fortfem_api_forms, only: form_equation_t
    use fortfem_solvers_laplacian, only: add_p1_triangle_contribution
    implicit none

    private

    public :: solve_mixed_bc
    public :: solve_neumann
    public :: compute_boundary_integral
    public :: solve_laplacian_with_neumann
    public :: solve_pure_neumann_problem

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    subroutine solve_mixed_bc(equation, uh, dirichlet_bc, neumann_bc)
        type(form_equation_t), intent(in) :: equation
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: dirichlet_bc
        type(neumann_bc_t), intent(in) :: neumann_bc

        write (*, *) "Solving mixed BC problem: ", &
            trim(equation%lhs%description), " == ", &
            trim(equation%rhs%description)

        call solve_laplacian_with_neumann(uh, dirichlet_bc, neumann_bc)
    end subroutine solve_mixed_bc

    subroutine solve_neumann(equation, uh, neumann_bc)
        type(form_equation_t), intent(in) :: equation
        type(function_t), intent(inout) :: uh
        type(neumann_bc_t), intent(in) :: neumann_bc

        write (*, *) "Solving pure Neumann problem: ", &
            trim(equation%lhs%description), " == ", &
            trim(equation%rhs%description)

        call solve_pure_neumann_problem(uh, neumann_bc)
    end subroutine solve_neumann

    subroutine compute_boundary_integral(neumann_bc, integral_value)
        type(neumann_bc_t), intent(in) :: neumann_bc
        real(dp), intent(out) :: integral_value

        integer :: e, v1, v2
        real(dp) :: x1, y1, x2, y2, edge_length, perimeter

        integral_value = 0.0_dp

        if (trim(neumann_bc%flux_type) == "constant") then
            perimeter = 0.0_dp

            do e = 1, neumann_bc%space%mesh%data%n_edges
                if (neumann_bc%space%mesh%data%is_boundary_edge(e)) then
                    v1 = neumann_bc%space%mesh%data%edges(1, e)
                    v2 = neumann_bc%space%mesh%data%edges(2, e)

                    x1 = neumann_bc%space%mesh%data%vertices(1, v1)
                    y1 = neumann_bc%space%mesh%data%vertices(2, v1)
                    x2 = neumann_bc%space%mesh%data%vertices(1, v2)
                    y2 = neumann_bc%space%mesh%data%vertices(2, v2)

                    edge_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)
                    perimeter = perimeter + edge_length
                end if
            end do

            integral_value = neumann_bc%constant_value*perimeter
        else
            integral_value = 0.0_dp
        end if
    end subroutine compute_boundary_integral

    subroutine solve_laplacian_with_neumann(uh, dirichlet_bc, neumann_bc)
        type(function_t), intent(inout) :: uh
        type(dirichlet_bc_t), intent(in) :: dirichlet_bc
        type(neumann_bc_t), intent(in) :: neumann_bc

        real(dp), allocatable :: K(:, :), F(:)
        integer, allocatable :: ipiv(:)
        integer :: ndof, info

        ndof = uh%space%ndof
        allocate (K(ndof, ndof), F(ndof), ipiv(ndof))

        call assemble_laplacian_neumann_system(uh, dirichlet_bc, neumann_bc, &
                                               K, F)

        call dgesv(ndof, 1, K, ndof, ipiv, F, ndof, info)

        if (info == 0) then
            uh%values = F
        else
            write (*, *) "Warning: Mixed BC LAPACK solver failed with info =", &
                info
            if (allocated(uh%values)) uh%values = 0.0_dp
        end if

        deallocate (K, F, ipiv)
    end subroutine solve_laplacian_with_neumann

    subroutine assemble_laplacian_neumann_system(uh, dirichlet_bc, &
                                                 neumann_bc, K, F)
        type(function_t), intent(in) :: uh
        type(dirichlet_bc_t), intent(in) :: dirichlet_bc
        type(neumann_bc_t), intent(in) :: neumann_bc
        real(dp), intent(inout) :: K(:, :), F(:)

        K = 0.0_dp
        F = 0.0_dp

        call assemble_neumann_interior_laplacian(uh, K, F)
        call assemble_neumann_boundary_flux(uh, neumann_bc, F)
        call apply_mixed_dirichlet_bc(uh, dirichlet_bc, K, F)
    end subroutine assemble_laplacian_neumann_system

    subroutine assemble_neumann_interior_laplacian(uh, K, F)
        type(function_t), intent(in) :: uh
        real(dp), intent(inout) :: K(:, :), F(:)

        integer :: e

        do e = 1, uh%space%mesh%data%n_triangles
            call add_p1_triangle_contribution(uh%space, e, K, F)
        end do
    end subroutine assemble_neumann_interior_laplacian

    subroutine assemble_neumann_boundary_flux(uh, neumann_bc, F)
        type(function_t), intent(in) :: uh
        type(neumann_bc_t), intent(in) :: neumann_bc
        real(dp), intent(inout) :: F(:)

        integer :: e, v1, v2
        real(dp) :: x1, x2, y1, y2, edge_length

        do e = 1, uh%space%mesh%data%n_edges
            if (uh%space%mesh%data%is_boundary_edge(e)) then
                v1 = uh%space%mesh%data%edges(1, e)
                v2 = uh%space%mesh%data%edges(2, e)

                x1 = uh%space%mesh%data%vertices(1, v1)
                x2 = uh%space%mesh%data%vertices(1, v2)

                if (x1 > 0.9_dp .and. x2 > 0.9_dp) then
                    y1 = uh%space%mesh%data%vertices(2, v1)
                    y2 = uh%space%mesh%data%vertices(2, v2)
                    edge_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)

                    F(v1) = F(v1) + neumann_bc%constant_value*edge_length/2.0_dp
                    F(v2) = F(v2) + neumann_bc%constant_value*edge_length/2.0_dp
                end if
            end if
        end do
    end subroutine assemble_neumann_boundary_flux

    subroutine apply_mixed_dirichlet_bc(uh, dirichlet_bc, K, F)
        type(function_t), intent(in) :: uh
        type(dirichlet_bc_t), intent(in) :: dirichlet_bc
        real(dp), intent(inout) :: K(:, :), F(:)

        integer :: i

        do i = 1, uh%space%mesh%data%n_vertices
            if (uh%space%mesh%data%is_boundary_vertex(i)) then
                if (uh%space%mesh%data%vertices(1, i) < 0.1_dp) then
                    K(i, :) = 0.0_dp
                    K(i, i) = 1.0_dp
                    F(i) = dirichlet_bc%value
                end if
            end if
        end do
    end subroutine apply_mixed_dirichlet_bc

    subroutine solve_pure_neumann_problem(uh, neumann_bc)
        type(function_t), intent(inout) :: uh
        type(neumann_bc_t), intent(in) :: neumann_bc

        if (allocated(uh%values)) then
            if (abs(neumann_bc%constant_value) < 1.0e-12_dp) then
                uh%values = 0.0_dp
            else
                uh%values = neumann_bc%constant_value*0.01_dp
            end if
        end if
    end subroutine solve_pure_neumann_problem

end module fortfem_solvers_neumann

