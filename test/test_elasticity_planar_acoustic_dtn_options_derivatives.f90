program test_elasticity_planar_acoustic_dtn_options_derivatives
    !! Derivatives of the planar elastic-acoustic DtN solver with the optional
    !! inputs: complex per-edge absorbing impedances (Z_p, Z_s) and a baffled
    !! DtN period (dtn_period_factor = 2.5, so 15 DtN samples for 6 interface
    !! samples). JVP against central differences for each input alone and all
    !! together, VJP against JVP through Re<u_bar, u_dot> = <x_bar, x_dot> +
    !! Re<z_bar, z_dot>, on a jittered mesh with damping.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_elasticity_planar_acoustic_dtn_p1, &
        solve_elasticity_planar_acoustic_dtn_p1_jvp, &
        solve_elasticity_planar_acoustic_dtn_p1_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: n = 6, nv = n*n, m = n
    real(dp), parameter :: omega = 2.3_dp, rho_f = 1.1_dp, nu = 0.27_dp
    real(dp), parameter :: rho_s = 1.3_dp, period = 1.0_dp + 1.0_dp/(n - 1)
    real(dp), parameter :: lfac = 2.5_dp
    ! FortFEM uses exp(-i omega t), so loss is Im(c) < 0 and Im(E) < 0.
    complex(dp), parameter :: c = (1.0_dp, -0.02_dp), e = (10.0_dp, -0.3_dp)
    real(dp) :: vertices(2, nv), normals(2, m), anormals(2, n - 1)
    integer :: triangles(3, 2*(n - 1)**2), inodes(m), edges(2, n - 1)
    integer, allocatable :: dofs(:)
    complex(dp), allocatable :: values(:)
    complex(dp) :: load(2*nv), pinc(m)
    ! Tangent direction
    real(dp) :: v_dot(2, nv), omega_dot, rho_f_dot, nu_dot, rho_s_dot
    complex(dp) :: c_dot, e_dot, p_dot(m)
    complex(dp) :: z(2, n - 1), z_dot(2, n - 1), z_bar(2, n - 1)
    real(dp) :: weights(9)
    complex(dp) :: u_dot(2*nv), u_fd(2*nv), u_bar(2*nv)
    real(dp) :: v_bar(2, nv), omega_bar, rho_f_bar, nu_bar, rho_s_bar
    real(dp) :: lhs, rhs, error
    complex(dp) :: c_bar, e_bar, p_bar(m)
    integer :: k, status, part
    logical :: all_passed
    character(*), parameter :: part_names(10) = [character(17) :: "vertices", &
        "angular frequency", "sound speed", "fluid density", "Young modulus", &
        "Poisson ratio", "solid density", "incident", "edge impedances", &
        "all together"]

    all_passed = .true.
    call build_problem()
    do k = 1, nv
        v_dot(:, k) = 0.03_dp*[sin(1.3_dp*k), cos(2.1_dp*k)]
    end do
    omega_dot = 0.2_dp
    c_dot = (0.05_dp, -0.01_dp)
    rho_f_dot = 0.4_dp
    nu_dot = 0.05_dp
    rho_s_dot = -0.3_dp
    e_dot = (0.7_dp, -0.2_dp)
    do k = 1, m
        p_dot(k) = cmplx(cos(0.9_dp*k), sin(0.4_dp*k), dp)
    end do
    do k = 1, n - 1
        z(:, k) = [cmplx(2.0_dp + 0.1_dp*k, -0.2_dp, dp), &
            cmplx(0.3_dp*k, 0.05_dp, dp)]
        z_dot(:, k) = [cmplx(0.4_dp*sin(1.1_dp*k), 0.1_dp, dp), &
            cmplx(-0.2_dp, 0.3_dp*cos(0.6_dp*k), dp)]
    end do

    do part = 1, 10
        weights = 0.0_dp
        if (part <= 9) then
            weights(part) = 1.0_dp
        else
            weights = 1.0_dp
        end if
        call solve_elasticity_planar_acoustic_dtn_p1_jvp(vertices, triangles, &
            inodes, normals, omega, c, rho_f, period, 2, e, nu, rho_s, load, pinc, &
            dofs, values, weights(1)*v_dot, weights(2)*omega_dot, &
            weights(3)*c_dot, weights(4)*rho_f_dot, weights(5)*e_dot, &
            weights(6)*nu_dot, weights(7)*rho_s_dot, weights(8)*p_dot, u_dot, &
            status, absorbing_edges=edges, absorbing_normals=anormals, &
            absorbing_impedances=z, absorbing_impedances_dot=weights(9)*z_dot, &
            dtn_period_factor=lfac)
        call record(status == 0, "JVP succeeds: "//trim(part_names(part)))
        call central_difference(weights, u_fd)
        error = maxval(abs(u_dot - u_fd))/maxval(abs(u_fd))
        write (*, "(a,a,es10.3)") trim(part_names(part)), ": JVP vs FD ", error
        call record(error < 1.0e-6_dp, "JVP matches central differences: "// &
            trim(part_names(part)))
    end do

    do k = 1, 2*nv
        u_bar(k) = cmplx(sin(0.7_dp*k), cos(1.9_dp*k), dp)
    end do
    call solve_elasticity_planar_acoustic_dtn_p1_vjp(vertices, triangles, inodes, &
        normals, omega, c, rho_f, period, 2, e, nu, rho_s, load, pinc, dofs, &
        values, u_bar, v_bar, omega_bar, c_bar, rho_f_bar, e_bar, nu_bar, &
        rho_s_bar, p_bar, status, absorbing_edges=edges, absorbing_normals=anormals, &
        absorbing_impedances=z, absorbing_impedances_bar=z_bar, &
        dtn_period_factor=lfac)
    call record(status == 0, "VJP succeeds")
    lhs = real(sum(conjg(u_bar)*u_dot), dp)
    rhs = sum(v_bar*v_dot) + omega_bar*omega_dot + &
        real(conjg(c_bar)*c_dot, dp) + rho_f_bar*rho_f_dot + &
        real(conjg(e_bar)*e_dot, dp) + nu_bar*nu_dot + rho_s_bar*rho_s_dot + &
        real(sum(conjg(p_bar)*p_dot), dp) + real(sum(conjg(z_bar)*z_dot), dp)
    write (*, "(a,2es24.15)") "dot test ", lhs, rhs
    call record(abs(lhs - rhs) <= 1.0e-10_dp*abs(lhs), &
        "VJP is the adjoint of the JVP")

    call check_summary("Planar elastic-acoustic DtN solver derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine central_difference(weights, derivative)
        real(dp), intent(in) :: weights(9)
        complex(dp), intent(out) :: derivative(:)
        real(dp), parameter :: h = 1.0e-5_dp
        complex(dp) :: plus(2*nv), minus(2*nv)

        call perturbed(h, weights, plus)
        call perturbed(-h, weights, minus)
        derivative = (plus - minus)/(2.0_dp*h)
    end subroutine central_difference

    subroutine perturbed(step, weights, solution)
        real(dp), intent(in) :: step, weights(9)
        complex(dp), intent(out) :: solution(:)
        integer :: solve_status

        call solve_elasticity_planar_acoustic_dtn_p1( &
            vertices + step*weights(1)*v_dot, triangles, inodes, normals, &
            omega + step*weights(2)*omega_dot, c + step*weights(3)*c_dot, &
            rho_f + step*weights(4)*rho_f_dot, period, 2, &
            e + step*weights(5)*e_dot, nu + step*weights(6)*nu_dot, &
            rho_s + step*weights(7)*rho_s_dot, load, &
            pinc + step*weights(8)*p_dot, dofs, values, solution, &
            solve_status, absorbing_edges=edges, absorbing_normals=anormals, &
            absorbing_impedances=z + step*weights(9)*z_dot, &
            dtn_period_factor=lfac)
        if (solve_status /= 0) error stop "perturbed solve failed"
    end subroutine perturbed

    subroutine build_problem()
        integer :: row, col, t, v, ll

        v = 0
        do row = 0, n - 1
            do col = 0, n - 1
                v = v + 1
                vertices(:, v) = [real(col, dp), real(row, dp)]/real(n - 1, dp) + &
                    0.05_dp/real(n, dp)*[sin(3.0_dp*v), cos(5.0_dp*v)]
            end do
        end do
        t = 0
        do row = 1, n - 1
            do col = 1, n - 1
                ll = col + (row - 1)*n
                t = t + 1
                triangles(:, t) = [ll, ll + 1, ll + n + 1]
                t = t + 1
                triangles(:, t) = [ll, ll + n + 1, ll + n]
            end do
        end do
        normals = 0.0_dp
        normals(1, :) = 1.0_dp
        do k = 1, m
            inodes(k) = k*n
            pinc(k) = cmplx(1.0_dp + 0.3_dp*sin(1.7_dp*k), 0.2_dp*cos(1.0_dp*k), dp)
        end do
        do k = 1, n - 1
            edges(:, k) = [(n - 1)*n + k, (n - 1)*n + k + 1]
            anormals(:, k) = [0.0_dp, 1.0_dp]
        end do
        dofs = [(2*k - 1, k=1, nv, n), (2*k, k=1, nv, n)]
        allocate(values(size(dofs)))
        do k = 1, size(dofs)
            values(k) = cmplx(0.01_dp*k, -0.02_dp*k, dp)
        end do
        do k = 1, 2*nv
            load(k) = cmplx(0.1_dp*sin(0.3_dp*k), 0.05_dp*cos(0.7_dp*k), dp)
        end do
    end subroutine build_problem

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_elasticity_planar_acoustic_dtn_options_derivatives
