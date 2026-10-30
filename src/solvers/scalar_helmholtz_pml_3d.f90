module fortfem_scalar_helmholtz_pml_3d
    !! P1 tetrahedral scalar Helmholtz with a cellwise diagonal complex stretch
    !! stretch(3, cell), with JVP and VJP; assembled through fortfem_block_csc by
    !! fortfem_scalar_helmholtz_pml_common, which states the weak form and the
    !! time convention.
    use fortfem_kinds, only: dp
    use fortfem_scalar_helmholtz_pml_common, only: pml_solve, pml_solve_jvp, &
        pml_solve_vjp
    implicit none
    private

    public :: solve_scalar_helmholtz_pml_p1_3d
    public :: solve_scalar_helmholtz_pml_p1_3d_jvp
    public :: solve_scalar_helmholtz_pml_p1_3d_vjp

contains

    subroutine solve_scalar_helmholtz_pml_p1_3d( &
            vertices, tetrahedra, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, status)
        !! status: 0 success, 1 invalid input, 2 failed solve.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        complex(dp), intent(in) :: stretch(:, :)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status

        call pml_solve(vertices, tetrahedra, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, .true., 2, solution, status)
    end subroutine solve_scalar_helmholtz_pml_p1_3d

    subroutine solve_scalar_helmholtz_pml_p1_3d_jvp( &
            vertices, tetrahedra, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, vertices_dot, stretch_dot, &
            wave_number_dot, volume_load_dot, dirichlet_values_dot, &
            solution_dot, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        complex(dp), intent(in) :: stretch(:, :)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        real(dp), intent(in) :: vertices_dot(:, :)
        complex(dp), intent(in) :: stretch_dot(:, :)
        real(dp), intent(in) :: wave_number_dot
        complex(dp), intent(in) :: volume_load_dot(:)
        complex(dp), intent(in) :: dirichlet_values_dot(:)
        complex(dp), allocatable, intent(out) :: solution_dot(:)
        integer, intent(out) :: status

        call pml_solve_jvp(vertices, tetrahedra, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, vertices_dot, stretch_dot, &
            wave_number_dot, volume_load_dot, dirichlet_values_dot, &
            solution_dot, status)
    end subroutine solve_scalar_helmholtz_pml_p1_3d_jvp

    subroutine solve_scalar_helmholtz_pml_p1_3d_vjp( &
            vertices, tetrahedra, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, solution_bar, &
            vertices_bar, stretch_bar, wave_number_bar, volume_load_bar, &
            dirichlet_values_bar, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :)
        complex(dp), intent(in) :: stretch(:, :)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(in) :: solution(:), solution_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :)
        complex(dp), intent(out) :: stretch_bar(:, :)
        real(dp), intent(out) :: wave_number_bar
        complex(dp), intent(out) :: volume_load_bar(:)
        complex(dp), intent(out) :: dirichlet_values_bar(:)
        integer, intent(out) :: status

        call pml_solve_vjp(vertices, tetrahedra, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, solution_bar, &
            vertices_bar, stretch_bar, wave_number_bar, volume_load_bar, &
            dirichlet_values_bar, status)
    end subroutine solve_scalar_helmholtz_pml_p1_3d_vjp

end module fortfem_scalar_helmholtz_pml_3d
