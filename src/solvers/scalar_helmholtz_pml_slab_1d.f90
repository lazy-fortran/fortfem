module fortfem_scalar_helmholtz_pml_slab_1d
    !! P1 finite elements for a scalar Helmholtz slab -u'' - k**2 u = 0 on
    !! [nodes(1), nodes(n)] terminated by a polynomial complex-coordinate perfectly
    !! matched layer on (physical_end, nodes(n)), with JVP and VJP.
    !!
    !! Time dependence is exp(-i omega t). An element with midpoint x_m beyond
    !! physical_end has the stretch s = 1 + i sigma/k with
    !! sigma = sigma_max ((x_m - physical_end)/(nodes(n) - physical_end))**p, so the
    !! outgoing wave exp(i k x) decays in the layer. u(nodes(1)) = left_value and
    !! u(nodes(n)) = 0.
    !!
    !! The attenuation sigma/k and the stretch coefficients are generated
    !! (gen_cartesian_pml_products), the interval matrices too
    !! (gen_scalar_helmholtz_products); the intervals are one group of a
    !! fortfem_block_csc pattern and the end values go to
    !! sparse_direct_solve_constrained. Storage is O(n).
    !! Cotangents follow sparse_direct_solve_constrained_vjp.
    use fortfem_block_csc, only: block_csc_add_group, block_csc_finalize, &
        block_csc_init, block_csc_new_matrix, block_csc_pattern_t
    use fortfem_cartesian_helmholtz_pml, only: &
        cartesian_scalar_helmholtz_pml_coefficients, &
        cartesian_scalar_helmholtz_pml_coefficients_jvp, &
        cartesian_scalar_helmholtz_pml_coefficients_vjp
    use fortfem_generated_cartesian_pml_attenuation, only: &
        generated_cartesian_pml_attenuation
    use fortfem_generated_cartesian_pml_attenuation_jvp, only: &
        generated_cartesian_pml_attenuation_jvp
    use fortfem_generated_cartesian_pml_attenuation_vjp, only: &
        generated_cartesian_pml_attenuation_vjp
    use fortfem_kinds, only: dp
    use fortfem_scalar_helmholtz_p1_csc, only: helmholtz_p1_assemble, &
        helmholtz_p1_assemble_jvp, helmholtz_p1_assemble_vjp
    use fortfem_sparse_direct, only: sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, sparse_direct_solve_constrained_vjp
    use fortsparse, only: csc_z_t
    implicit none
    private

    public :: solve_scalar_helmholtz_pml_slab_1d
    public :: solve_scalar_helmholtz_pml_slab_1d_jvp
    public :: solve_scalar_helmholtz_pml_slab_1d_vjp

    type :: slab_system_t
        type(block_csc_pattern_t) :: pattern
        integer :: cell_group = 0
        type(csc_z_t) :: matrix
        integer, allocatable :: intervals(:, :)
        complex(dp), allocatable :: stretch(:), gradient(:, :), mass(:)
        complex(dp), allocatable :: rhs(:), fixed(:)
        logical, allocatable :: constrained(:)
    end type slab_system_t

contains

    subroutine solve_scalar_helmholtz_pml_slab_1d( &
            nodes, physical_end, wave_number, sigma_max, polynomial_degree, &
            left_value, solution, status)
        !! status: 0 success, 1 invalid input, 2 failed solve.
        real(dp), intent(in) :: nodes(:)
        real(dp), intent(in) :: physical_end, wave_number, sigma_max
        integer, intent(in) :: polynomial_degree
        complex(dp), intent(in) :: left_value
        complex(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status

        type(slab_system_t) :: system

        call build_system(nodes, physical_end, wave_number, sigma_max, &
            polynomial_degree, left_value, system, status)
        if (status /= 0) return
        allocate(solution(size(nodes)))
        call sparse_direct_solve_constrained(system%matrix, system%rhs, &
            system%constrained, system%fixed, solution, status)
        if (status /= 0) status = 2
    end subroutine solve_scalar_helmholtz_pml_slab_1d

    subroutine solve_scalar_helmholtz_pml_slab_1d_jvp( &
            nodes, physical_end, wave_number, sigma_max, polynomial_degree, &
            left_value, nodes_dot, physical_end_dot, wave_number_dot, &
            sigma_max_dot, left_value_dot, solution_dot, status)
        real(dp), intent(in) :: nodes(:)
        real(dp), intent(in) :: physical_end, wave_number, sigma_max
        integer, intent(in) :: polynomial_degree
        complex(dp), intent(in) :: left_value
        real(dp), intent(in) :: nodes_dot(:)
        real(dp), intent(in) :: physical_end_dot, wave_number_dot
        real(dp), intent(in) :: sigma_max_dot
        complex(dp), intent(in) :: left_value_dot
        complex(dp), allocatable, intent(out) :: solution_dot(:)
        integer, intent(out) :: status

        type(slab_system_t) :: system
        complex(dp), allocatable :: values_dot(:), fixed_dot(:)
        complex(dp), allocatable :: gradient_dot(:, :), mass_dot(:)
        complex(dp) :: full(3), full_dot(3), gradient_full(3)
        real(dp) :: attenuation_dot, layer_width, layer_width_dot, distance
        real(dp) :: distance_dot
        integer :: e, n, local_status

        status = 1
        n = size(nodes)
        if (size(nodes_dot) /= n) return
        call build_system(nodes, physical_end, wave_number, sigma_max, &
            polynomial_degree, left_value, system, status)
        if (status /= 0) return
        layer_width = nodes(n) - physical_end
        layer_width_dot = nodes_dot(n) - physical_end_dot
        allocate(gradient_dot(1, n - 1), mass_dot(n - 1))
        do e = 1, n - 1
            full = [system%stretch(e), (1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)]
            full_dot = (0.0_dp, 0.0_dp)
            distance = midpoint(nodes, e) - physical_end
            if (distance > 0.0_dp) then
                distance_dot = midpoint(nodes_dot, e) - physical_end_dot
                call generated_cartesian_pml_attenuation_jvp(distance, &
                    layer_width, wave_number, sigma_max, &
                    real(polynomial_degree, dp), distance_dot, layer_width_dot, &
                    wave_number_dot, sigma_max_dot, attenuation_dot)
                full_dot(1) = cmplx(0.0_dp, attenuation_dot, dp)
            end if
            call cartesian_scalar_helmholtz_pml_coefficients_jvp(full, full_dot, &
                gradient_full, mass_dot(e), local_status)
            if (local_status /= 0) return
            gradient_dot(1, e) = gradient_full(1)
        end do
        allocate(values_dot(system%pattern%nnz))
        values_dot = (0.0_dp, 0.0_dp)
        call helmholtz_p1_assemble_jvp(system%pattern, system%cell_group, &
            reshape(nodes, [1, n]), system%intervals, system%gradient, &
            system%mass, wave_number, reshape(nodes_dot, [1, n]), gradient_dot, &
            mass_dot, wave_number_dot, values_dot, status)
        if (status /= 0) return
        allocate(fixed_dot(n))
        fixed_dot = (0.0_dp, 0.0_dp)
        fixed_dot(1) = left_value_dot
        allocate(solution_dot(n))
        call sparse_direct_solve_constrained_jvp(system%matrix, system%rhs, &
            system%constrained, system%fixed, values_dot, 0.0_dp*system%rhs, &
            fixed_dot, solution_dot, status)
        if (status /= 0) status = 2
    end subroutine solve_scalar_helmholtz_pml_slab_1d_jvp

    subroutine solve_scalar_helmholtz_pml_slab_1d_vjp( &
            nodes, physical_end, wave_number, sigma_max, polynomial_degree, &
            left_value, solution, solution_bar, nodes_bar, physical_end_bar, &
            wave_number_bar, sigma_max_bar, left_value_bar, status)
        !! Cotangents for solution_bar at the solution `solution`; the outputs
        !! are overwritten.
        real(dp), intent(in) :: nodes(:)
        real(dp), intent(in) :: physical_end, wave_number, sigma_max
        integer, intent(in) :: polynomial_degree
        complex(dp), intent(in) :: left_value
        complex(dp), intent(in) :: solution(:)
        complex(dp), intent(in) :: solution_bar(:)
        real(dp), intent(out) :: nodes_bar(:)
        real(dp), intent(out) :: physical_end_bar, wave_number_bar
        real(dp), intent(out) :: sigma_max_bar
        complex(dp), intent(out) :: left_value_bar
        integer, intent(out) :: status

        type(slab_system_t) :: system
        complex(dp), allocatable :: values_bar(:), rhs_bar(:), fixed_bar(:)
        complex(dp), allocatable :: gradient_bar(:, :), mass_bar(:)
        real(dp), allocatable :: vertices_bar(:, :)
        complex(dp) :: full(3), full_gradient_bar(3), full_bar(3)
        real(dp) :: distance, layer_width, distance_bar, width_bar, k_bar
        real(dp) :: attenuation_bar, smax_bar
        integer :: e, n

        nodes_bar = 0.0_dp
        physical_end_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        sigma_max_bar = 0.0_dp
        left_value_bar = (0.0_dp, 0.0_dp)
        status = 1
        n = size(nodes)
        if (size(solution) /= n .or. size(solution_bar) /= n) return
        if (size(nodes_bar) /= n) return
        call build_system(nodes, physical_end, wave_number, sigma_max, &
            polynomial_degree, left_value, system, status)
        if (status /= 0) return
        allocate(values_bar(system%pattern%nnz), rhs_bar(n), fixed_bar(n))
        call sparse_direct_solve_constrained_vjp(system%matrix, system%rhs, &
            system%constrained, system%fixed, solution, solution_bar, &
            values_bar, rhs_bar, fixed_bar, status)
        if (status /= 0) then
            status = 2
            return
        end if
        left_value_bar = fixed_bar(1)
        allocate(vertices_bar(1, n), gradient_bar(1, n - 1), mass_bar(n - 1))
        call helmholtz_p1_assemble_vjp(system%pattern, system%cell_group, &
            reshape(nodes, [1, n]), system%intervals, system%gradient, &
            system%mass, wave_number, values_bar, vertices_bar, gradient_bar, &
            mass_bar, wave_number_bar, status)
        if (status /= 0) return
        nodes_bar = vertices_bar(1, :)
        layer_width = nodes(n) - physical_end
        ! Serial in element order: the layer-width cotangent is a reduction.
        do e = 1, n - 1
            distance = midpoint(nodes, e) - physical_end
            if (distance <= 0.0_dp) cycle
            full = [system%stretch(e), (1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)]
            full_gradient_bar = [gradient_bar(1, e), (0.0_dp, 0.0_dp), &
                (0.0_dp, 0.0_dp)]
            call cartesian_scalar_helmholtz_pml_coefficients_vjp(full, &
                full_gradient_bar, mass_bar(e), full_bar, status)
            if (status /= 0) return
            ! s = 1 + i attenuation: attenuation_bar = Re(conj(i) s_bar).
            attenuation_bar = aimag(full_bar(1))
            call generated_cartesian_pml_attenuation_vjp(distance, layer_width, &
                wave_number, sigma_max, real(polynomial_degree, dp), &
                attenuation_bar, distance_bar, width_bar, k_bar, smax_bar)
            sigma_max_bar = sigma_max_bar + smax_bar
            wave_number_bar = wave_number_bar + k_bar
            nodes_bar(e) = nodes_bar(e) + 0.5_dp*distance_bar
            nodes_bar(e + 1) = nodes_bar(e + 1) + 0.5_dp*distance_bar
            nodes_bar(n) = nodes_bar(n) + width_bar
            physical_end_bar = physical_end_bar - distance_bar - width_bar
        end do
    end subroutine solve_scalar_helmholtz_pml_slab_1d_vjp

    subroutine build_system(nodes, physical_end, wave_number, sigma_max, &
            polynomial_degree, left_value, system, status)
        real(dp), intent(in) :: nodes(:)
        real(dp), intent(in) :: physical_end, wave_number, sigma_max
        integer, intent(in) :: polynomial_degree
        complex(dp), intent(in) :: left_value
        type(slab_system_t), intent(out) :: system
        integer, intent(out) :: status

        complex(dp) :: gradient_full(3)
        real(dp) :: attenuation, distance
        integer :: e, n, local_status

        status = 1
        n = size(nodes)
        if (n < 3) return
        if (any(nodes(2:) <= nodes(:n - 1))) return
        if (physical_end <= nodes(1) .or. physical_end >= nodes(n)) return
        if (wave_number <= 0.0_dp .or. sigma_max <= 0.0_dp) return
        if (polynomial_degree < 1) return

        allocate(system%intervals(2, n - 1), system%stretch(n - 1))
        allocate(system%gradient(1, n - 1), system%mass(n - 1))
        do e = 1, n - 1
            system%intervals(:, e) = [e, e + 1]
            system%stretch(e) = (1.0_dp, 0.0_dp)
            distance = midpoint(nodes, e) - physical_end
            if (distance > 0.0_dp) then
                call generated_cartesian_pml_attenuation(distance, &
                    nodes(n) - physical_end, wave_number, sigma_max, &
                    real(polynomial_degree, dp), attenuation)
                system%stretch(e) = cmplx(1.0_dp, attenuation, dp)
            end if
            call cartesian_scalar_helmholtz_pml_coefficients( &
                [system%stretch(e), (1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)], &
                gradient_full, system%mass(e), local_status)
            if (local_status /= 0) return
            system%gradient(1, e) = gradient_full(1)
        end do
        call block_csc_init(system%pattern, n)
        call block_csc_add_group(system%pattern, system%intervals, &
            system%cell_group)
        call block_csc_finalize(system%pattern)
        call block_csc_new_matrix(system%pattern, system%matrix)
        call helmholtz_p1_assemble(system%pattern, system%cell_group, &
            reshape(nodes, [1, n]), system%intervals, system%gradient, &
            system%mass, wave_number, system%matrix%val, status)
        if (status /= 0) return
        allocate(system%rhs(n), system%fixed(n), system%constrained(n))
        system%rhs = (0.0_dp, 0.0_dp)
        system%fixed = (0.0_dp, 0.0_dp)
        system%fixed(1) = left_value
        system%constrained = .false.
        system%constrained([1, n]) = .true.
    end subroutine build_system

    pure real(dp) function midpoint(x, e)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: e

        midpoint = 0.5_dp*(x(e) + x(e + 1))
    end function midpoint

end module fortfem_scalar_helmholtz_pml_slab_1d
