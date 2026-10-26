program test_equilibrated_mixed_2d_slow
    !! Equilibrated-residual bound of the RT_k-DG_k mixed Poisson solver
    !! (actual solve_mixed_poisson_rt solutions on the unit square). Oracles:
    !! (1) u = x(1-x)y(1-y) with k = 3: sigma_h = -grad u exactly, the local
    !! postprocessing recovers u, so the bound vanishes to rounding;
    !! (2) u = sin(pi x) sin(pi y), k = 0..2 on three meshes: the flux bound is
    !! >= the exact flux error ||sigma_h + grad u||, the potential bound is
    !! >= ||grad(u - s)||, the efficiency index stays below a measured cap and
    !! the bound converges at the optimal rate k + 1.
    use check, only: check_condition, check_summary
    use fortfem_core, only: mesh_t, rectangle_mesh
    use fortfem_equilibrated_nonconforming_2d, only: &
        MIXED_POTENTIAL_OSWALD, MIXED_POTENTIAL_PATCH, &
        mixed_poisson_flux_bound, nonconforming_estimate_t, NONCONFORMING_OK
    use fortfem_feec, only: build_triangle_trimmed_dof_map, &
        evaluate_triangle_lagrange_basis, evaluate_triangle_rt_interpolant, &
        initialize_triangle_lagrange_basis, initialize_triangle_raviart_thomas, &
        solve_mixed_poisson_rt, triangle_duffy_quadrature, &
        triangle_lagrange_basis_t, triangle_rt_basis_t
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: friedrichs = 1.0_dp/(pi*sqrt(2.0_dp))
    !> Measured indices: patch 1.04-1.45, Oswald 2.0-3.1.
    real(dp), parameter :: efficiency_cap(0:2, 2) = reshape( &
        [1.8_dp, 1.8_dp, 1.8_dp, 3.6_dp, 3.6_dp, 3.6_dp], [3, 2])
    type(nonconforming_estimate_t) :: est
    real(dp) :: flux_error, potential_error, bound(3, 2), index, rate
    integer :: k, level, vertex_count, method
    type(mesh_t) :: solved_mesh
    real(dp), allocatable :: flux_dofs(:), pressure_dofs(:)
    logical :: polynomial

    polynomial = .true.
    call solve_problem(5, 3, solved_mesh, flux_dofs, pressure_dofs)
    do method = MIXED_POTENTIAL_PATCH, MIXED_POTENTIAL_OSWALD
        call estimate_solution(solved_mesh%data, 3, method, flux_dofs, &
            pressure_dofs, est, flux_error, potential_error)
        call check_condition(est%bound < 1.0e-9_dp, &
            "exact polynomial solution: flux bound vanishes")
        call check_condition(est%potential_bound < 1.0e-9_dp, &
            "exact polynomial solution: potential bound vanishes")
        print '(a, i0, a, 3es10.2)', "polynomial case, method ", method, &
            " bound, flux error, potential error ", est%bound, flux_error, &
            potential_error
    end do

    polynomial = .false.
    do k = 0, 2
        do level = 1, 3
            vertex_count = 1 + 2**(level + 1)
            ! Both reconstructions use the same mixed solution. Factor once.
            call solve_problem(vertex_count, k, solved_mesh, flux_dofs, pressure_dofs)
            do method = MIXED_POTENTIAL_PATCH, MIXED_POTENTIAL_OSWALD
                call estimate_solution(solved_mesh%data, k, method, flux_dofs, &
                    pressure_dofs, est, flux_error, potential_error)
                bound(level, method) = est%bound
                index = est%bound/flux_error
                call check_condition(est%bound >= flux_error, &
                    "flux bound is above the exact flux error")
                call check_condition(est%potential_bound >= potential_error, &
                    "potential bound is above the exact potential error")
                call check_condition(index <= efficiency_cap(k, method), &
                    "flux efficiency index below the cap")
                print '(a, i1, a, i2, a, i3, 5(a, es10.2))', "method ", method, &
                    " k", k, " n", vertex_count - 1, &
                    " flux err", flux_error, " bound", est%bound, " index", index, &
                    " X", est%flux_defect, " osc+div", &
                    est%oscillation + est%divergence_defect
            end do
        end do
        do method = MIXED_POTENTIAL_PATCH, MIXED_POTENTIAL_OSWALD
            rate = log(bound(2, method)/bound(3, method))/log(2.0_dp)
            call check_condition(rate > real(k + 1, dp) - 0.2_dp, &
                "flux bound converges at the optimal rate k + 1")
            print '(a, i2, a, f6.2)', "k", k, " bound rate", rate
        end do
    end do
    call check_summary("equilibrated_mixed_2d")

contains

    subroutine solve_problem(vertex_count, degree, mesh, flux_dofs, pressure_dofs)
        integer, intent(in) :: vertex_count, degree
        type(mesh_t), intent(out) :: mesh
        real(dp), allocatable, intent(out) :: flux_dofs(:), pressure_dofs(:)
        type(fortsparse_status_t) :: status

        mesh = rectangle_mesh(vertex_count, vertex_count, &
            [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
        if (polynomial) then
            call solve_mixed_poisson_rt(mesh%data, degree, 2*degree + 6, &
                polynomial_source, flux_dofs, pressure_dofs, status)
        else
            call solve_mixed_poisson_rt(mesh%data, degree, 2*degree + 6, &
                smooth_source, flux_dofs, pressure_dofs, status)
        end if
        if (status%code /= 0) error stop "mixed solve failed"
    end subroutine solve_problem

    subroutine estimate_solution(mesh, degree, reconstruction, flux_dofs, &
            pressure_dofs, estimate, flux_err, potential_err)
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree, reconstruction
        real(dp), intent(in) :: flux_dofs(:), pressure_dofs(:)
        type(nonconforming_estimate_t), intent(out) :: estimate
        real(dp), intent(out) :: flux_err, potential_err
        integer :: info

        if (polynomial) then
            call mixed_poisson_flux_bound(mesh, degree, 2*degree + 8, &
                flux_dofs, pressure_dofs, polynomial_source, friedrichs, &
                estimate, info, reconstruction)
        else
            call mixed_poisson_flux_bound(mesh, degree, 2*degree + 8, &
                flux_dofs, pressure_dofs, smooth_source, friedrichs, &
                estimate, info, reconstruction)
        end if
        call check_condition(info == NONCONFORMING_OK, "mixed estimate succeeds")
        if (info /= NONCONFORMING_OK) then
            print '(a, i0)', "mixed estimate status ", info
            error stop "mixed estimate failed"
        end if
        call exact_errors(mesh, degree, flux_dofs, estimate%potential, &
            flux_err, potential_err)
    end subroutine estimate_solution

    subroutine exact_errors(mesh, degree, flux_dofs, potential, flux_err, &
            potential_err)
        !! ||sigma_h + grad u||_L2 and ||grad(u - s)||_L2 by a rule of degree
        !! 2k + 12 (independent of the estimator's rule and kernels).
        type(mesh_2d_t), intent(inout) :: mesh
        integer, intent(in) :: degree
        real(dp), intent(in) :: flux_dofs(:), potential(:, :)
        real(dp), intent(out) :: flux_err, potential_err
        type(triangle_lagrange_basis_t) :: pbasis
        type(triangle_rt_basis_t) :: fbasis
        integer, allocatable :: flux_map(:, :), transforms(:, :)
        real(dp), allocatable :: xi(:), eta(:), w(:), local_flux(:)
        real(dp), allocatable :: pval(:), pgrad(:, :)
        real(dp) :: v(2, 3), jac(2, 2), det, point(2), sigma(2), div, grad(2)
        real(dp) :: gref(2), gu(2), jit(2, 2)
        integer :: flux_count, st, t, q, potential_degree

        call build_triangle_trimmed_dof_map(mesh, degree + 1, flux_map, &
            transforms, flux_count, st)
        if (st /= 0) error stop "flux map"
        call initialize_triangle_raviart_thomas(degree, fbasis, st)
        potential_degree = nint(0.5_dp*(sqrt(8.0_dp*size(potential, 1) + 1.0_dp) &
            - 3.0_dp))
        call initialize_triangle_lagrange_basis(potential_degree, pbasis, st)
        call triangle_duffy_quadrature(2*degree + 12, xi, eta, w, st)
        allocate (local_flux(size(flux_map, 1)), pval(size(potential, 1)))
        allocate (pgrad(2, size(potential, 1)))
        flux_err = 0.0_dp
        potential_err = 0.0_dp
        do t = 1, mesh%n_triangles
            v = mesh%vertices(:, mesh%triangles(1:3, t))
            jac(:, 1) = v(:, 2) - v(:, 1)
            jac(:, 2) = v(:, 3) - v(:, 1)
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            jit(1, :) = [jac(2, 2), -jac(2, 1)]/det
            jit(2, :) = [-jac(1, 2), jac(1, 1)]/det
            local_flux = real(transforms(:, t), dp)*flux_dofs(flux_map(:, t))
            do q = 1, size(w)
                point = v(:, 1) + xi(q)*jac(:, 1) + eta(q)*jac(:, 2)
                call evaluate_triangle_rt_interpolant(v, fbasis, local_flux, &
                    xi(q), eta(q), sigma, div, st)
                gu = exact_gradient(point)
                flux_err = flux_err + w(q)*det*sum((sigma + gu)**2)
                call evaluate_triangle_lagrange_basis(pbasis, xi(q), eta(q), &
                    pval, pgrad, st)
                gref = matmul(pgrad, potential(:, t))
                grad = matmul(jit, gref)
                potential_err = potential_err + w(q)*det*sum((gu - grad)**2)
            end do
        end do
        flux_err = sqrt(flux_err)
        potential_err = sqrt(potential_err)
    end subroutine exact_errors

    function exact_gradient(point) result(g)
        real(dp), intent(in) :: point(2)
        real(dp) :: g(2)

        if (polynomial) then
            g = [(1.0_dp - 2.0_dp*point(1))*point(2)*(1.0_dp - point(2)), &
                point(1)*(1.0_dp - point(1))*(1.0_dp - 2.0_dp*point(2))]
        else
            g = pi*[cos(pi*point(1))*sin(pi*point(2)), &
                sin(pi*point(1))*cos(pi*point(2))]
        end if
    end function exact_gradient

    pure function polynomial_source(x, y) result(value)
        real(dp), intent(in) :: x, y
        real(dp) :: value

        value = 2.0_dp*(x*(1.0_dp - x) + y*(1.0_dp - y))
    end function polynomial_source

    pure function smooth_source(x, y) result(value)
        real(dp), intent(in) :: x, y
        real(dp) :: value

        value = 2.0_dp*pi*pi*sin(pi*x)*sin(pi*y)
    end function smooth_source

end program test_equilibrated_mixed_2d_slow
