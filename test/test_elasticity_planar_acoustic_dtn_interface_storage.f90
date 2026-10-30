program test_elasticity_planar_acoustic_dtn_interface_storage
    !! The planar DtN couples only the normal displacement component. For
    !! axis-aligned interface normals the solver stores the m x m acoustic
    !! form on the normal-component dofs (m*m entries), not the dense
    !! 2m x 2m block (4*m*m); tilted normals keep the dense block. The
    !! component kernel must equal the matching entry of the block kernel,
    !! with consistent JVP/VJP, and the reduced path must reproduce the
    !! manufactured longitudinal wave.
    use check, only: check_condition, check_summary
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_elasticity_planar_acoustic_dtn_2d, only: &
        planar_dtn_interface_entries, solve_elasticity_planar_acoustic_dtn_p1, &
        solve_elasticity_planar_acoustic_dtn_p1_jvp
    use fortfem_generated_interface_normal_block, only: &
        generated_interface_normal_block
    use fortfem_generated_interface_normal_block_jvp, only: &
        generated_interface_normal_block_jvp
    use fortfem_generated_interface_normal_block_vjp, only: &
        generated_interface_normal_block_vjp
    use fortfem_generated_interface_normal_component_block, only: &
        generated_interface_normal_component_block
    use fortfem_generated_interface_normal_component_block_jvp, only: &
        generated_interface_normal_component_block_jvp
    use fortfem_generated_interface_normal_component_block_vjp, only: &
        generated_interface_normal_component_block_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: m = 7
    real(dp) :: normals(2, m)
    integer(int64) :: entries
    logical :: all_passed

    all_passed = .true.
    normals = 0.0_dp
    normals(1, :) = 1.0_dp
    entries = planar_dtn_interface_entries(m, normals)
    call record_condition(entries == int(m, int64)*int(m, int64), &
        "axis-aligned normals store the m x m normal-component block")
    normals(2, :) = 0.3_dp
    entries = planar_dtn_interface_entries(m, normals)
    call record_condition(entries == 4_int64*int(m, int64)*int(m, int64), &
        "tilted normals keep the dense 2m x 2m block")

    call kernel_equivalence()
    call manufactured_longitudinal_wave()
    call tilted_jvp_matches_finite_differences()

    call check_summary("planar acoustic DtN interface-block storage")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        if (.not. condition) all_passed = .false.
        call check_condition(condition, description)
    end subroutine record_condition

    pure function pseudo_random(index) result(value)
        integer, intent(in) :: index
        real(dp) :: value

        integer(int64) :: state

        state = int(index, int64)*2654435761_int64 + 101_int64
        state = state*6364136223846793005_int64 + 1442695040888963407_int64
        value = real(iand(state, 1048575_int64), dp)/1048575.0_dp
    end function pseudo_random

    subroutine kernel_equivalence()
        !! The component kernel equals the matching entry of the block kernel
        !! and its JVP matches; the VJP pairs like the block VJP.
        real(dp) :: na(2), nb(2), value_dot(2)
        complex(dp) :: value, block(2, 2), block_dot(2, 2), zero(2, 2)
        complex(dp) :: scalar, scalar_dot, scalar_bar, value_bar, dense_bar
        complex(dp) :: left, right
        integer :: k, p, q

        left = (0.0_dp, 0.0_dp)
        right = (0.0_dp, 0.0_dp)
        do k = 1, 12
            na = [pseudo_random(k), pseudo_random(k + 31)]
            nb = [pseudo_random(k + 7), pseudo_random(k + 41)]
            value = cmplx(pseudo_random(k + 19), pseudo_random(k + 23), dp)
            value_dot = [pseudo_random(k + 29), pseudo_random(k + 37)]
            call generated_interface_normal_block(na(1), na(2), nb(1), nb(2), &
                value, block)
            call generated_interface_normal_block_jvp(na(1), na(2), nb(1), &
                nb(2), value, cmplx(value_dot(1), value_dot(2), dp), block_dot)
            do q = 1, 2
                do p = 1, 2
                    call generated_interface_normal_component_block(na(p), &
                        nb(q), value, scalar)
                    if (scalar /= block(p, q)) then
                        call record_condition(.false., &
                            "component kernel matches the block kernel entry")
                        return
                    end if
                    call generated_interface_normal_component_block_jvp(na(p), &
                        nb(q), value, cmplx(value_dot(1), value_dot(2), dp), &
                        scalar_dot)
                    if (scalar_dot /= block_dot(p, q)) then
                        call record_condition(.false., &
                            "component JVP matches the block JVP entry")
                        return
                    end if
                    ! Cotangent on one entry; both VJP kernels are the plain
                    ! transpose, so the cotangent of the form entry must agree.
                    scalar_bar = cmplx(pseudo_random(k + 53), 0.0_dp, dp)
                    call generated_interface_normal_component_block_vjp(na(p), &
                        nb(q), value, scalar_bar, value_bar)
                    left = left + value_bar*cmplx(value_dot(1), value_dot(2), dp)
                    zero = (0.0_dp, 0.0_dp)
                    zero(p, q) = scalar_bar
                    call generated_interface_normal_block_vjp(na(1), na(2), &
                        nb(1), nb(2), value, zero, dense_bar)
                    right = right + dense_bar* &
                        cmplx(value_dot(1), value_dot(2), dp)
                end do
            end do
        end do
        call record_condition(abs(left - right) <= &
            1.0e-12_dp*max(1.0_dp, abs(right)), &
            "component VJP pairs like the block VJP")
    end subroutine kernel_equivalence

    subroutine manufactured_longitudinal_wave()
        !! The reduced normal-component path reproduces the manufactured
        !! longitudinal wave of the convergence test, fixture and oracle.
        integer, parameter :: n = 9
        real(dp), parameter :: omega = 2.3_dp
        complex(dp), parameter :: young = cmplx(10.0_dp, -0.1_dp, dp)
        real(dp), parameter :: nu = 0.25_dp, rho_s = 1.0_dp, rho_f = 1.0_dp
        complex(dp), parameter :: cf = cmplx(1.0_dp, -0.01_dp, dp)
        real(dp) :: vertices(2, n*n), interface_normals(2, n), spacing
        integer :: triangles(3, 2*(n - 1)**2), interface_nodes(n)
        integer, allocatable :: dirichlet_dofs(:)
        integer :: status, node, k
        complex(dp) :: lambda, mu, modulus, k_s, incident
        complex(dp) :: load(2*n*n), incident_pressure(n), solution(2*n*n)
        complex(dp), allocatable :: dirichlet_values(:)
        real(dp) :: rms_error

        mu = young/(2.0_dp*(1.0_dp + nu))
        lambda = young*nu/((1.0_dp + nu)*(1.0_dp - 2.0_dp*nu))
        modulus = lambda + 2.0_dp*mu
        k_s = omega*sqrt(rho_s/modulus)
        call build_square_mesh(n, vertices, triangles)
        do node = 1, n
            interface_nodes(node) = node*n
        end do
        interface_normals = 0.0_dp
        interface_normals(1, :) = 1.0_dp
        incident = -modulus*k_s*cos(k_s) + cmplx(0.0_dp, 1.0_dp, dp)* &
            rho_f*omega*cf*sin(k_s)
        incident_pressure = incident
        load = (0.0_dp, 0.0_dp)
        spacing = 1.0_dp/real(n - 1, dp)
        call build_dirichlet_dofs(n, dirichlet_dofs)
        allocate(dirichlet_values(size(dirichlet_dofs)))
        do node = 1, size(dirichlet_dofs)
            k = (dirichlet_dofs(node) + 1)/2
            if (mod(dirichlet_dofs(node), 2) == 1) then
                dirichlet_values(node) = sin(k_s*vertices(1, k))
            else
                dirichlet_values(node) = (0.0_dp, 0.0_dp)
            end if
        end do
        call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, &
            interface_nodes, interface_normals, omega, cf, rho_f, &
            spacing*real(n, dp), 2, young, nu, rho_s, load, incident_pressure, &
            dirichlet_dofs, dirichlet_values, solution, status)
        call record_condition(status == 0, &
            "normal-component DtN path solves the manufactured wave")
        rms_error = 0.0_dp
        do node = 1, n*n
            rms_error = rms_error + abs(solution(2*node - 1) - &
                sin(k_s*vertices(1, node)))**2 + abs(solution(2*node))**2
        end do
        rms_error = sqrt(rms_error/real(2*n*n, dp))
        call record_condition(rms_error < 2.0e-2_dp, &
            "normal-component DtN path reproduces the manufactured wave")
    end subroutine manufactured_longitudinal_wave

    subroutine tilted_jvp_matches_finite_differences()
        !! The dense fallback path stays differentiable: the JVP of the
        !! solution with respect to omega matches a central difference when
        !! both interface components are active.
        integer, parameter :: n = 5
        real(dp) :: vertices(2, n*n), normals(2, n), vertices_dot(2, n*n)
        integer :: triangles(3, 2*(n - 1)**2), interface_nodes(n)
        integer :: dofs(2*n), status, solve_status, node, dof
        real(dp), parameter :: omega = 1.7_dp, period = 0.4_dp
        real(dp), parameter :: nu = 0.3_dp, rho_s = 2.7_dp, rho_f = 1.0_dp
        complex(dp), parameter :: young = cmplx(70.0_dp, -0.7_dp, dp)
        complex(dp), parameter :: cf = cmplx(1.48_dp, -0.01_dp, dp)
        complex(dp) :: load(2*n*n), incident_pressure(n), incident_dot(n)
        complex(dp) :: solution(2*n*n), solution_dot(2*n*n)
        complex(dp) :: forward(2*n*n), backward(2*n*n)
        complex(dp) :: dirichlet_values(2*n)
        complex(dp) :: functional_dot, slope
        real(dp) :: weight(2*n*n), h

        normals = 0.0_dp
        normals(1, :) = 0.8_dp
        normals(2, :) = 0.6_dp
        call build_square_mesh(n, vertices, triangles)
        do node = 1, n
            interface_nodes(node) = node
        end do
        incident_pressure = cmplx(1.0_dp, 0.2_dp, dp)
        incident_dot = (0.0_dp, 0.0_dp)
        load = (0.0_dp, 0.0_dp)
        vertices_dot = 0.0_dp
        do dof = 1, 2*n
            dofs(dof) = dof
            dirichlet_values(dof) = (0.0_dp, 0.0_dp)
        end do
        do node = 1, 2*n*n
            weight(node) = pseudo_random(node)
        end do
        call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, &
            interface_nodes, normals, omega, cf, rho_f, period, 2, young, nu, &
            rho_s, load, incident_pressure, dofs, dirichlet_values, solution, &
            status)
        call solve_elasticity_planar_acoustic_dtn_p1_jvp(vertices, triangles, &
            interface_nodes, normals, omega, cf, rho_f, period, 2, young, nu, &
            rho_s, load, incident_pressure, dofs, dirichlet_values, &
            vertices_dot, 1.0_dp, (0.0_dp, 0.0_dp), 0.0_dp, &
            (0.0_dp, 0.0_dp), 0.0_dp, 0.0_dp, incident_dot, solution_dot, &
            solve_status)
        call record_condition(status == 0 .and. solve_status == 0, &
            "tilted interface solves and differentiates")
        functional_dot = sum(solution_dot*weight)
        h = 1.0e-5_dp
        call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, &
            interface_nodes, normals, omega + h, cf, rho_f, period, 2, young, &
            nu, rho_s, load, incident_pressure, dofs, dirichlet_values, &
            forward, status)
        call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, &
            interface_nodes, normals, omega - h, cf, rho_f, period, 2, young, &
            nu, rho_s, load, incident_pressure, dofs, dirichlet_values, &
            backward, status)
        slope = (sum(forward*weight) - sum(backward*weight))/(2.0_dp*h)
        call record_condition(abs(functional_dot - slope) <= &
            1.0e-4_dp*max(1.0_dp, abs(slope)), &
            "tilted interface JVP matches the finite-difference slope")
    end subroutine tilted_jvp_matches_finite_differences

    subroutine build_square_mesh(n, vertices, triangles)
        integer, intent(in) :: n
        real(dp), intent(out) :: vertices(:, :)
        integer, intent(out) :: triangles(:, :)

        real(dp) :: spacing
        integer :: column, lower_left, row, triangle, vertex

        spacing = 1.0_dp/real(n - 1, dp)
        vertex = 0
        do row = 0, n - 1
            do column = 0, n - 1
                vertex = vertex + 1
                vertices(:, vertex) = spacing*real([column, row], dp)
            end do
        end do
        triangle = 0
        do row = 1, n - 1
            do column = 1, n - 1
                lower_left = column + (row - 1)*n
                triangle = triangle + 1
                triangles(:, triangle) = [lower_left, lower_left + 1, &
                    lower_left + n + 1]
                triangle = triangle + 1
                triangles(:, triangle) = [lower_left, lower_left + n + 1, &
                    lower_left + n]
            end do
        end do
    end subroutine build_square_mesh

    subroutine build_dirichlet_dofs(n, dofs)
        integer, intent(in) :: n
        integer, allocatable, intent(out) :: dofs(:)

        logical, allocatable :: constrained(:)
        integer :: column, count_dofs, node, row

        allocate(constrained(2*n*n))
        constrained = .false.
        do row = 1, n
            node = (row - 1)*n + 1
            constrained(2*node - 1:2*node) = .true.
        end do
        do column = 1, n
            node = column
            constrained(2*node - 1:2*node) = .true.
            node = (n - 1)*n + column
            constrained(2*node - 1:2*node) = .true.
        end do
        count_dofs = count(constrained)
        allocate(dofs(count_dofs))
        count_dofs = 0
        do node = 1, 2*n*n
            if (.not. constrained(node)) cycle
            count_dofs = count_dofs + 1
            dofs(count_dofs) = node
        end do
    end subroutine build_dirichlet_dofs

end program test_elasticity_planar_acoustic_dtn_interface_storage
