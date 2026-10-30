module fortfem_scalar_helmholtz_planar_dtn_2d
    !! P1 scalar Helmholtz -div grad u - k**2 u = f on a triangle mesh with a
    !! periodic planar Dirichlet-to-Neumann (DtN) boundary, with JVP and VJP.
    !!
    !! Time dependence is exp(-i omega t): the DtN multiplier of
    !! fortfem_planar_helmholtz_dtn, i k_n for propagating modes with the
    !! principal root, belongs to the outgoing wave exp(i k x). The DtN samples
    !! dtn_nodes are m equispaced points of one period of the exterior.
    !!
    !! One CSC pattern from fortfem_block_csc holds the triangles (generated P1
    !! Helmholtz cell) and one dense m x m block over dtn_nodes holding minus the
    !! DtN form. Dirichlet values go to sparse_direct_solve_constrained. Storage is
    !! O(triangles + nnz + m**2).
    !!
    !! Derivatives cover the wave number, the period, the load and the Dirichlet
    !! values. Cotangents follow sparse_direct_solve_constrained_vjp.
    use fortfem_block_csc, only: block_csc_add, block_csc_add_group, &
        block_csc_finalize, block_csc_gather, block_csc_init, &
        block_csc_new_matrix, block_csc_pattern_t
    use fortfem_dirichlet_constraints, only: dirichlet_mask, dirichlet_mask_vjp
    use fortfem_kinds, only: dp
    use fortfem_planar_helmholtz_dtn, only: &
        assemble_planar_helmholtz_dtn_form, &
        assemble_planar_helmholtz_dtn_form_jvp, &
        assemble_planar_helmholtz_dtn_form_vjp
    use fortfem_scalar_helmholtz_p1_csc, only: helmholtz_p1_assemble, &
        helmholtz_p1_assemble_jvp, helmholtz_p1_assemble_vjp
    use fortfem_sparse_direct, only: sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, sparse_direct_solve_constrained_vjp
    use fortsparse, only: csc_z_t
    implicit none
    private

    public :: solve_scalar_helmholtz_planar_dtn_p1
    public :: solve_scalar_helmholtz_planar_dtn_p1_jvp
    public :: solve_scalar_helmholtz_planar_dtn_p1_vjp

    type :: dtn_system_t
        type(block_csc_pattern_t) :: pattern
        integer :: cell_group = 0
        integer :: dtn_group = 0
        type(csc_z_t) :: matrix
        complex(dp), allocatable :: rhs(:), fixed(:)
        logical, allocatable :: constrained(:)
        complex(dp), allocatable :: gradient(:, :), mass(:)
    end type dtn_system_t

contains

    subroutine solve_scalar_helmholtz_planar_dtn_p1( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, status)
        !! status: 0 success, 1 invalid input, otherwise the failed solve.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), dtn_nodes(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        type(dtn_system_t) :: system

        solution = (0.0_dp, 0.0_dp)
        call build_system(vertices, triangles, dtn_nodes, wavenumber, period, &
            volume_load, dirichlet_nodes, dirichlet_values, size(solution), &
            system, status)
        if (status /= 0) return
        call sparse_direct_solve_constrained(system%matrix, system%rhs, &
            system%constrained, system%fixed, solution, status)
    end subroutine solve_scalar_helmholtz_planar_dtn_p1

    subroutine solve_scalar_helmholtz_planar_dtn_p1_jvp( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, wavenumber_dot, period_dot, &
            volume_load_dot, dirichlet_values_dot, solution_dot, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), dtn_nodes(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        real(dp), intent(in) :: wavenumber_dot, period_dot
        complex(dp), intent(in) :: volume_load_dot(:)
        complex(dp), intent(in) :: dirichlet_values_dot(:)
        complex(dp), intent(out) :: solution_dot(:)
        integer, intent(out) :: status

        type(dtn_system_t) :: system
        complex(dp), allocatable :: values_dot(:), form_dot(:, :), fixed_dot(:)
        logical, allocatable :: constrained(:)
        integer :: m

        solution_dot = (0.0_dp, 0.0_dp)
        status = 1
        if (size(volume_load_dot) /= size(vertices, 2)) return
        if (size(dirichlet_values_dot) /= size(dirichlet_values)) return
        call build_system(vertices, triangles, dtn_nodes, wavenumber, period, &
            volume_load, dirichlet_nodes, dirichlet_values, size(solution_dot), &
            system, status)
        if (status /= 0) return
        m = size(dtn_nodes)
        allocate(values_dot(system%pattern%nnz))
        values_dot = (0.0_dp, 0.0_dp)
        call helmholtz_p1_assemble_jvp(system%pattern, system%cell_group, &
            vertices, triangles, system%gradient, system%mass, wavenumber, &
            0.0_dp*vertices, 0.0_dp*system%gradient, 0.0_dp*system%mass, &
            wavenumber_dot, values_dot, status)
        if (status /= 0) return
        allocate(form_dot(m, m))
        call assemble_planar_helmholtz_dtn_form_jvp(m, wavenumber, period, &
            wavenumber_dot, period_dot, form_dot, status)
        if (status /= 0) return
        call block_csc_add(system%pattern, system%dtn_group, &
            reshape(-form_dot, [m, m, 1]), values_dot)
        call dirichlet_mask(size(vertices, 2), dirichlet_nodes, &
            dirichlet_values_dot, constrained, fixed_dot)
        call sparse_direct_solve_constrained_jvp(system%matrix, system%rhs, &
            system%constrained, system%fixed, values_dot, volume_load_dot, &
            fixed_dot, solution_dot, status)
    end subroutine solve_scalar_helmholtz_planar_dtn_p1_jvp

    subroutine solve_scalar_helmholtz_planar_dtn_p1_vjp( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, solution_bar, &
            volume_load_bar, dirichlet_values_bar, wavenumber_bar, period_bar, &
            status)
        !! Cotangents for the solution cotangent solution_bar at the solution
        !! `solution`; the outputs are overwritten.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), dtn_nodes(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        complex(dp), intent(in) :: solution(:), solution_bar(:)
        complex(dp), intent(out) :: volume_load_bar(:)
        complex(dp), intent(out) :: dirichlet_values_bar(:)
        real(dp), intent(out) :: wavenumber_bar, period_bar
        integer, intent(out) :: status

        type(dtn_system_t) :: system
        complex(dp), allocatable :: values_bar(:), fixed_bar(:)
        complex(dp), allocatable :: block_bar(:, :, :), gradient_bar(:, :)
        complex(dp), allocatable :: mass_bar(:)
        real(dp), allocatable :: vertices_bar(:, :)
        real(dp) :: cell_wavenumber_bar
        integer :: m

        volume_load_bar = (0.0_dp, 0.0_dp)
        dirichlet_values_bar = (0.0_dp, 0.0_dp)
        wavenumber_bar = 0.0_dp
        period_bar = 0.0_dp
        status = 1
        if (size(solution_bar) /= size(vertices, 2)) return
        if (size(volume_load_bar) /= size(vertices, 2)) return
        if (size(dirichlet_values_bar) /= size(dirichlet_values)) return
        call build_system(vertices, triangles, dtn_nodes, wavenumber, period, &
            volume_load, dirichlet_nodes, dirichlet_values, size(solution), &
            system, status)
        if (status /= 0) return
        m = size(dtn_nodes)
        allocate(values_bar(system%pattern%nnz), fixed_bar(size(solution)))
        call sparse_direct_solve_constrained_vjp(system%matrix, system%rhs, &
            system%constrained, system%fixed, solution, solution_bar, &
            values_bar, volume_load_bar, fixed_bar, status)
        if (status /= 0) return
        call dirichlet_mask_vjp(dirichlet_nodes, fixed_bar, dirichlet_values_bar)

        allocate(block_bar(m, m, 1))
        call block_csc_gather(system%pattern, system%dtn_group, values_bar, &
            block_bar)
        call assemble_planar_helmholtz_dtn_form_vjp(m, wavenumber, period, &
            -block_bar(:, :, 1), wavenumber_bar, period_bar, status)
        if (status /= 0) return
        allocate(vertices_bar(size(vertices, 1), size(vertices, 2)))
        allocate(gradient_bar(2, size(triangles, 2)), mass_bar(size(triangles, 2)))
        call helmholtz_p1_assemble_vjp(system%pattern, system%cell_group, &
            vertices, triangles, system%gradient, system%mass, wavenumber, &
            values_bar, vertices_bar, gradient_bar, mass_bar, &
            cell_wavenumber_bar, status)
        if (status /= 0) return
        wavenumber_bar = wavenumber_bar + cell_wavenumber_bar
    end subroutine solve_scalar_helmholtz_planar_dtn_p1_vjp

    subroutine build_system(vertices, triangles, dtn_nodes, wavenumber, period, &
            volume_load, dirichlet_nodes, dirichlet_values, solution_size, &
            system, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), dtn_nodes(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        integer, intent(in) :: solution_size
        type(dtn_system_t), intent(out) :: system
        integer, intent(out) :: status

        complex(dp), allocatable :: form(:, :)
        integer :: m, vertex_count

        status = 1
        if (.not. valid_solver_inputs(vertices, triangles, dtn_nodes, &
            wavenumber, period, volume_load, dirichlet_nodes, dirichlet_values, &
            solution_size)) return
        vertex_count = size(vertices, 2)
        m = size(dtn_nodes)
        call block_csc_init(system%pattern, vertex_count)
        call block_csc_add_group(system%pattern, triangles, system%cell_group)
        call block_csc_add_group(system%pattern, reshape(dtn_nodes, [m, 1]), &
            system%dtn_group)
        call block_csc_finalize(system%pattern)
        call block_csc_new_matrix(system%pattern, system%matrix)

        allocate(system%gradient(2, size(triangles, 2)), &
            system%mass(size(triangles, 2)))
        system%gradient = (1.0_dp, 0.0_dp)
        system%mass = (1.0_dp, 0.0_dp)
        call helmholtz_p1_assemble(system%pattern, system%cell_group, vertices, &
            triangles, system%gradient, system%mass, wavenumber, &
            system%matrix%val, status)
        if (status /= 0) return
        allocate(form(m, m))
        call assemble_planar_helmholtz_dtn_form(m, wavenumber, period, form, &
            status)
        if (status /= 0) return
        call block_csc_add(system%pattern, system%dtn_group, &
            reshape(-form, [m, m, 1]), system%matrix%val)
        system%rhs = volume_load
        call dirichlet_mask(vertex_count, dirichlet_nodes, dirichlet_values, &
            system%constrained, system%fixed)
    end subroutine build_system

    pure logical function valid_solver_inputs( &
            vertices, triangles, dtn_nodes, wavenumber, period, volume_load, &
            dirichlet_nodes, dirichlet_values, solution_size) result(valid)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :), dtn_nodes(:)
        real(dp), intent(in) :: wavenumber, period
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        integer, intent(in) :: solution_size

        integer :: vertex_count

        valid = .false.
        vertex_count = size(vertices, 2)
        if (size(vertices, 1) /= 2 .or. size(triangles, 1) /= 3) return
        if (vertex_count < 1 .or. size(triangles, 2) < 1) return
        if (wavenumber <= 0.0_dp .or. period <= 0.0_dp) return
        if (size(dtn_nodes) < 3) return
        if (size(volume_load) /= vertex_count) return
        if (solution_size /= vertex_count) return
        if (size(dirichlet_nodes) /= size(dirichlet_values)) return
        if (any(triangles < 1) .or. any(triangles > vertex_count)) return
        if (any(dtn_nodes < 1) .or. any(dtn_nodes > vertex_count)) return
        if (any(dirichlet_nodes < 1) .or. &
            any(dirichlet_nodes > vertex_count)) return
        valid = .true.
    end function valid_solver_inputs

end module fortfem_scalar_helmholtz_planar_dtn_2d
