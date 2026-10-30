module fortfem_scalar_helmholtz_pml_common
    !! P1 scalar Helmholtz with a cellwise diagonal complex stretch on triangles or
    !! tetrahedra, the shared engine of fortfem_scalar_helmholtz_pml_2d and _3d.
    !!
    !! With S = diag(s_1, .., s_dim) per cell the weak form is
    !! (det(S) S^-2 grad u, grad v) - k**2 (det(S) u, v) = (f, v) (see
    !! fortfem_cartesian_helmholtz_pml). The stretch is the caller's: under the
    !! exp(-i omega t) convention of this library, s = 1 + i sigma/k with
    !! sigma >= 0 damps the outgoing wave exp(i k x).
    !!
    !! The cells form one group of a fortfem_block_csc pattern (generated P1
    !! Helmholtz cell, generated stretch coefficients) and Dirichlet values go to
    !! sparse_direct_solve_constrained. Storage is O(cells + nnz).
    !! Cotangents follow sparse_direct_solve_constrained_vjp.
    use fortfem_block_csc, only: block_csc_add_group, block_csc_finalize, &
        block_csc_init, block_csc_new_matrix, block_csc_pattern_t
    use fortfem_cartesian_helmholtz_pml, only: &
        cartesian_scalar_helmholtz_pml_coefficients, &
        cartesian_scalar_helmholtz_pml_coefficients_jvp, &
        cartesian_scalar_helmholtz_pml_coefficients_vjp
    use fortfem_dirichlet_constraints, only: dirichlet_mask, dirichlet_mask_vjp
    use fortfem_kinds, only: dp
    use fortfem_scalar_helmholtz_p1_csc, only: helmholtz_p1_assemble, &
        helmholtz_p1_assemble_jvp, helmholtz_p1_assemble_vjp
    use fortfem_sparse_direct, only: sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, sparse_direct_solve_constrained_vjp
    use fortsparse, only: csc_z_t
    implicit none
    private

    public :: pml_solve
    public :: pml_solve_jvp
    public :: pml_solve_vjp

    type :: pml_system_t
        type(block_csc_pattern_t) :: pattern
        integer :: cell_group = 0
        type(csc_z_t) :: matrix
        complex(dp), allocatable :: rhs(:), fixed(:)
        logical, allocatable :: constrained(:)
        complex(dp), allocatable :: gradient(:, :), mass(:)
    end type pml_system_t

contains

    subroutine pml_solve(vertices, cells, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, reject_repeated, solve_failure, &
            solution, status)
        !! status: 0 success, 1 invalid input (a repeated Dirichlet node too when
        !! reject_repeated), solve_failure for a failed factorization. Otherwise
        !! a repeated node keeps its last value.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :)
        complex(dp), intent(in) :: stretch(:, :)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        logical, intent(in) :: reject_repeated
        integer, intent(in) :: solve_failure
        complex(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status

        type(pml_system_t) :: system

        call build_system(vertices, cells, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, reject_repeated, system, status)
        if (status /= 0) return
        allocate(solution(size(vertices, 2)))
        call sparse_direct_solve_constrained(system%matrix, system%rhs, &
            system%constrained, system%fixed, solution, status)
        if (status /= 0) status = solve_failure
    end subroutine pml_solve

    subroutine pml_solve_jvp(vertices, cells, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, vertices_dot, stretch_dot, &
            wave_number_dot, volume_load_dot, dirichlet_values_dot, &
            solution_dot, status)
        !! Directional derivative of the solution; repeated Dirichlet nodes are
        !! rejected (status 1).
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :)
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

        type(pml_system_t) :: system
        complex(dp), allocatable :: values_dot(:), fixed_dot(:)
        complex(dp), allocatable :: gradient_dot(:, :), mass_dot(:)
        logical, allocatable :: constrained(:)

        status = 1
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (any(shape(stretch_dot) /= shape(stretch))) return
        if (size(volume_load_dot) /= size(vertices, 2)) return
        if (size(dirichlet_values_dot) /= size(dirichlet_values)) return
        call build_system(vertices, cells, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, .true., system, status)
        if (status /= 0) return
        call stretch_coefficients_jvp(stretch, stretch_dot, gradient_dot, &
            mass_dot, status)
        if (status /= 0) return
        allocate(values_dot(system%pattern%nnz))
        values_dot = (0.0_dp, 0.0_dp)
        call helmholtz_p1_assemble_jvp(system%pattern, system%cell_group, &
            vertices, cells, system%gradient, system%mass, wave_number, &
            vertices_dot, gradient_dot, mass_dot, wave_number_dot, values_dot, &
            status)
        if (status /= 0) return
        call dirichlet_mask(size(vertices, 2), dirichlet_nodes, &
            dirichlet_values_dot, constrained, fixed_dot)
        allocate(solution_dot(size(vertices, 2)))
        call sparse_direct_solve_constrained_jvp(system%matrix, system%rhs, &
            system%constrained, system%fixed, values_dot, volume_load_dot, &
            fixed_dot, solution_dot, status)
    end subroutine pml_solve_jvp

    subroutine pml_solve_vjp(vertices, cells, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, solution, solution_bar, &
            vertices_bar, stretch_bar, wave_number_bar, volume_load_bar, &
            dirichlet_values_bar, status)
        !! Cotangents for solution_bar at the solution `solution`; the outputs are
        !! overwritten. Repeated Dirichlet nodes are rejected (status 1).
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :)
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

        type(pml_system_t) :: system
        complex(dp), allocatable :: values_bar(:), fixed_bar(:)
        complex(dp), allocatable :: gradient_bar(:, :), mass_bar(:)

        vertices_bar = 0.0_dp
        stretch_bar = (0.0_dp, 0.0_dp)
        wave_number_bar = 0.0_dp
        volume_load_bar = (0.0_dp, 0.0_dp)
        dirichlet_values_bar = (0.0_dp, 0.0_dp)
        status = 1
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (any(shape(stretch_bar) /= shape(stretch))) return
        if (size(solution) /= size(vertices, 2)) return
        if (size(solution_bar) /= size(vertices, 2)) return
        if (size(volume_load_bar) /= size(vertices, 2)) return
        if (size(dirichlet_values_bar) /= size(dirichlet_values)) return
        call build_system(vertices, cells, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, .true., system, status)
        if (status /= 0) return
        allocate(values_bar(system%pattern%nnz), fixed_bar(size(solution)))
        call sparse_direct_solve_constrained_vjp(system%matrix, system%rhs, &
            system%constrained, system%fixed, solution, solution_bar, &
            values_bar, volume_load_bar, fixed_bar, status)
        if (status /= 0) return
        call dirichlet_mask_vjp(dirichlet_nodes, fixed_bar, dirichlet_values_bar)
        allocate(gradient_bar, mold=system%gradient)
        allocate(mass_bar, mold=system%mass)
        call helmholtz_p1_assemble_vjp(system%pattern, system%cell_group, &
            vertices, cells, system%gradient, system%mass, wave_number, &
            values_bar, vertices_bar, gradient_bar, mass_bar, wave_number_bar, &
            status)
        if (status /= 0) return
        call stretch_coefficients_vjp(stretch, gradient_bar, mass_bar, &
            stretch_bar, status)
    end subroutine pml_solve_vjp

    subroutine build_system(vertices, cells, stretch, wave_number, volume_load, &
            dirichlet_nodes, dirichlet_values, reject_repeated, system, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :)
        complex(dp), intent(in) :: stretch(:, :)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: volume_load(:)
        integer, intent(in) :: dirichlet_nodes(:)
        complex(dp), intent(in) :: dirichlet_values(:)
        logical, intent(in) :: reject_repeated
        type(pml_system_t), intent(out) :: system
        integer, intent(out) :: status

        integer :: vertex_count

        status = 1
        if (.not. valid_inputs(vertices, cells, stretch, wave_number, &
            volume_load, dirichlet_nodes, dirichlet_values)) return
        vertex_count = size(vertices, 2)
        if (reject_repeated) then
            if (repeated(vertex_count, dirichlet_nodes)) return
        end if
        call stretch_coefficients(stretch, system%gradient, system%mass, status)
        if (status /= 0) return
        call block_csc_init(system%pattern, vertex_count)
        call block_csc_add_group(system%pattern, cells, system%cell_group)
        call block_csc_finalize(system%pattern)
        call block_csc_new_matrix(system%pattern, system%matrix)
        call helmholtz_p1_assemble(system%pattern, system%cell_group, vertices, &
            cells, system%gradient, system%mass, wave_number, system%matrix%val, &
            status)
        if (status /= 0) return
        system%rhs = volume_load
        call dirichlet_mask(vertex_count, dirichlet_nodes, dirichlet_values, &
            system%constrained, system%fixed)
    end subroutine build_system

    subroutine stretch_coefficients(stretch, gradient, mass, status)
        !! Gradient and mass coefficients of each cell's stretch (s, 1, ..).
        complex(dp), intent(in) :: stretch(:, :)
        complex(dp), allocatable, intent(out) :: gradient(:, :), mass(:)
        integer, intent(out) :: status

        complex(dp) :: full(3), gradient_full(3)
        integer :: dim, e, local_status
        logical :: failed

        dim = size(stretch, 1)
        allocate(gradient(dim, size(stretch, 2)), mass(size(stretch, 2)))
        failed = .false.
        !$omp parallel do default(none) shared(stretch, gradient, mass, dim) &
        !$omp private(e, full, gradient_full, local_status) &
        !$omp reduction(.or.: failed) schedule(static)
        do e = 1, size(stretch, 2)
            full = (1.0_dp, 0.0_dp)
            full(1:dim) = stretch(:, e)
            call cartesian_scalar_helmholtz_pml_coefficients(full, &
                gradient_full, mass(e), local_status)
            gradient(:, e) = gradient_full(1:dim)
            if (local_status /= 0) failed = .true.
        end do
        !$omp end parallel do
        status = merge(1, 0, failed)
    end subroutine stretch_coefficients

    subroutine stretch_coefficients_jvp(stretch, stretch_dot, gradient_dot, &
            mass_dot, status)
        complex(dp), intent(in) :: stretch(:, :), stretch_dot(:, :)
        complex(dp), allocatable, intent(out) :: gradient_dot(:, :), mass_dot(:)
        integer, intent(out) :: status

        complex(dp) :: full(3), full_dot(3), gradient_full(3)
        integer :: dim, e, local_status
        logical :: failed

        dim = size(stretch, 1)
        allocate(gradient_dot(dim, size(stretch, 2)), mass_dot(size(stretch, 2)))
        failed = .false.
        !$omp parallel do default(none) &
        !$omp shared(stretch, stretch_dot, gradient_dot, mass_dot, dim) &
        !$omp private(e, full, full_dot, gradient_full, local_status) &
        !$omp reduction(.or.: failed) schedule(static)
        do e = 1, size(stretch, 2)
            full = (1.0_dp, 0.0_dp)
            full(1:dim) = stretch(:, e)
            full_dot = (0.0_dp, 0.0_dp)
            full_dot(1:dim) = stretch_dot(:, e)
            call cartesian_scalar_helmholtz_pml_coefficients_jvp(full, full_dot, &
                gradient_full, mass_dot(e), local_status)
            gradient_dot(:, e) = gradient_full(1:dim)
            if (local_status /= 0) failed = .true.
        end do
        !$omp end parallel do
        status = merge(1, 0, failed)
    end subroutine stretch_coefficients_jvp

    subroutine stretch_coefficients_vjp(stretch, gradient_bar, mass_bar, &
            stretch_bar, status)
        complex(dp), intent(in) :: stretch(:, :), gradient_bar(:, :), mass_bar(:)
        complex(dp), intent(out) :: stretch_bar(:, :)
        integer, intent(out) :: status

        complex(dp) :: full(3), full_gradient_bar(3), full_bar(3)
        integer :: dim, e, local_status
        logical :: failed

        dim = size(stretch, 1)
        failed = .false.
        !$omp parallel do default(none) &
        !$omp shared(stretch, gradient_bar, mass_bar, stretch_bar, dim) &
        !$omp private(e, full, full_gradient_bar, full_bar, local_status) &
        !$omp reduction(.or.: failed) schedule(static)
        do e = 1, size(stretch, 2)
            full = (1.0_dp, 0.0_dp)
            full(1:dim) = stretch(:, e)
            full_gradient_bar = (0.0_dp, 0.0_dp)
            full_gradient_bar(1:dim) = gradient_bar(:, e)
            call cartesian_scalar_helmholtz_pml_coefficients_vjp(full, &
                full_gradient_bar, mass_bar(e), full_bar, local_status)
            stretch_bar(:, e) = full_bar(1:dim)
            if (local_status /= 0) failed = .true.
        end do
        !$omp end parallel do
        status = merge(1, 0, failed)
    end subroutine stretch_coefficients_vjp

    pure logical function repeated(vertex_count, nodes)
        integer, intent(in) :: vertex_count, nodes(:)
        logical, allocatable :: seen(:)
        integer :: k

        allocate(seen(vertex_count))
        seen = .false.
        repeated = .false.
        do k = 1, size(nodes)
            if (seen(nodes(k))) then
                repeated = .true.
                return
            end if
            seen(nodes(k)) = .true.
        end do
    end function repeated

    pure logical function valid_inputs(vertices, cells, stretch, wave_number, &
            volume_load, dirichlet_nodes, dirichlet_values) result(valid)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :), dirichlet_nodes(:)
        complex(dp), intent(in) :: stretch(:, :), volume_load(:)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: dirichlet_values(:)

        integer :: dim, vertex_count

        valid = .false.
        dim = size(vertices, 1)
        vertex_count = size(vertices, 2)
        if (dim /= 2 .and. dim /= 3) return
        if (vertex_count < dim + 1) return
        if (size(cells, 1) /= dim + 1 .or. size(cells, 2) < 1) return
        if (size(stretch, 1) /= dim .or. size(stretch, 2) /= size(cells, 2)) return
        if (size(volume_load) /= vertex_count) return
        if (size(dirichlet_nodes) /= size(dirichlet_values)) return
        if (wave_number <= 0.0_dp) return
        if (any(cells < 1) .or. any(cells > vertex_count)) return
        if (any(dirichlet_nodes < 1) .or. &
            any(dirichlet_nodes > vertex_count)) return
        valid = .true.
    end function valid_inputs

end module fortfem_scalar_helmholtz_pml_common
