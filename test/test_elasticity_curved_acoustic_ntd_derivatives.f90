program test_elasticity_curved_acoustic_ntd_derivatives
    !! JVP of the curved elastic-acoustic NtD solver against central differences
    !! (each parameter alone and all together) and VJP against JVP, on a two-ring
    !! mesh of a perturbed ellipse. The interface is held fixed: interior vertices
    !! move, and a tangent that moves an interface vertex is rejected.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_elasticity_curved_acoustic_ntd_p1, &
        solve_elasticity_curved_acoustic_ntd_p1_jvp, &
        solve_elasticity_curved_acoustic_ntd_p1_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: np = 16, nv = 2*np + 1
    real(dp), parameter :: omega = 1.1_dp, c = 0.8_dp, rho_f = 1.7_dp
    real(dp), parameter :: nu = 0.25_dp, rho_s = 2.0_dp
    ! FortFEM uses exp(-i omega t), so loss is Im(E) < 0.
    complex(dp), parameter :: e = (10.0_dp, -0.2_dp)
    real(dp) :: vertices(2, nv)
    integer :: triangles(3, 3*np), inodes(np), pnodes(2, np)
    integer, allocatable :: dofs(:)
    complex(dp), allocatable :: values(:)
    complex(dp) :: load(2*nv), pinc(np)
    real(dp) :: v_dot(2, nv), rho_f_dot, nu_dot, rho_s_dot
    complex(dp) :: e_dot, p_dot(np)
    real(dp) :: weights(6)
    complex(dp) :: u_dot(2*nv), u_fd(2*nv), u_bar(2*nv)
    real(dp) :: v_bar(2, nv), rho_f_bar, nu_bar, rho_s_bar, lhs, rhs, error
    complex(dp) :: e_bar, p_bar(np)
    real(dp) :: moved(2, nv)
    integer :: k, status, part
    logical :: all_passed
    character(*), parameter :: part_names(7) = [character(16) :: "vertices", &
        "fluid density", "Young modulus", "Poisson ratio", "solid density", &
        "incident", "all together"]

    all_passed = .true.
    call build_problem()
    v_dot = 0.0_dp
    do k = 1, np + 1
        v_dot(:, k) = 0.03_dp*[sin(1.3_dp*k), cos(2.1_dp*k)]
    end do
    rho_f_dot = 0.4_dp
    nu_dot = 0.05_dp
    rho_s_dot = -0.3_dp
    e_dot = (0.7_dp, -0.2_dp)
    do k = 1, np
        p_dot(k) = cmplx(cos(0.9_dp*k), sin(0.4_dp*k), dp)
    end do

    do part = 1, 7
        weights = 0.0_dp
        if (part <= 6) then
            weights(part) = 1.0_dp
        else
            weights = 1.0_dp
        end if
        call solve_elasticity_curved_acoustic_ntd_p1_jvp(vertices, triangles, &
            inodes, pnodes, omega, c, rho_f, e, nu, rho_s, load, pinc, dofs, &
            values, weights(1)*v_dot, weights(2)*rho_f_dot, weights(3)*e_dot, &
            weights(4)*nu_dot, weights(5)*rho_s_dot, weights(6)*p_dot, u_dot, &
            status)
        call record(status == 0, "JVP succeeds: "//trim(part_names(part)))
        call central_difference(weights, u_fd)
        error = maxval(abs(u_dot - u_fd))/maxval(abs(u_fd))
        write (*, "(a,a,es10.3)") trim(part_names(part)), ": JVP vs FD ", error
        call record(error < 1.0e-6_dp, "JVP matches central differences: "// &
            trim(part_names(part)))
    end do

    moved = v_dot
    moved(:, np + 2) = [0.01_dp, 0.0_dp]
    call solve_elasticity_curved_acoustic_ntd_p1_jvp(vertices, triangles, inodes, &
        pnodes, omega, c, rho_f, e, nu, rho_s, load, pinc, dofs, values, moved, &
        0.0_dp, (0.0_dp, 0.0_dp), 0.0_dp, 0.0_dp, 0.0_dp*p_dot, u_fd, status)
    call record(status == 3, "JVP rejects a moving interface vertex")

    do k = 1, 2*nv
        u_bar(k) = cmplx(sin(0.7_dp*k), cos(1.9_dp*k), dp)
    end do
    call solve_elasticity_curved_acoustic_ntd_p1_vjp(vertices, triangles, inodes, &
        pnodes, omega, c, rho_f, e, nu, rho_s, load, pinc, dofs, values, u_bar, &
        v_bar, rho_f_bar, e_bar, nu_bar, rho_s_bar, p_bar, status)
    call record(status == 0, "VJP succeeds")
    call record(all(v_bar(:, np + 2:) == 0.0_dp), &
        "VJP leaves the interface vertex cotangents at zero")
    lhs = real(sum(conjg(u_bar)*u_dot), dp)
    rhs = sum(v_bar*v_dot) + rho_f_bar*rho_f_dot + &
        real(conjg(e_bar)*e_dot, dp) + nu_bar*nu_dot + rho_s_bar*rho_s_dot + &
        real(sum(conjg(p_bar)*p_dot), dp)
    write (*, "(a,2es24.15)") "dot test ", lhs, rhs
    call record(abs(lhs - rhs) <= 1.0e-10_dp*abs(lhs), &
        "VJP is the adjoint of the JVP")

    call check_summary("Curved elastic-acoustic NtD solver derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine central_difference(weights, derivative)
        real(dp), intent(in) :: weights(6)
        complex(dp), intent(out) :: derivative(:)
        real(dp), parameter :: h = 1.0e-5_dp
        complex(dp) :: plus(2*nv), minus(2*nv)

        call perturbed(h, weights, plus)
        call perturbed(-h, weights, minus)
        derivative = (plus - minus)/(2.0_dp*h)
    end subroutine central_difference

    subroutine perturbed(step, weights, solution)
        real(dp), intent(in) :: step, weights(6)
        complex(dp), intent(out) :: solution(:)
        integer :: solve_status

        call solve_elasticity_curved_acoustic_ntd_p1( &
            vertices + step*weights(1)*v_dot, triangles, inodes, pnodes, &
            omega, c, rho_f + step*weights(2)*rho_f_dot, &
            e + step*weights(3)*e_dot, nu + step*weights(4)*nu_dot, &
            rho_s + step*weights(5)*rho_s_dot, load, &
            pinc + step*weights(6)*p_dot, dofs, values, solution, solve_status)
        if (solve_status /= 0) error stop "perturbed solve failed"
    end subroutine perturbed

    subroutine build_problem()
        real(dp) :: angle, pi
        integer :: next

        pi = acos(-1.0_dp)
        vertices(:, 1) = [0.02_dp, -0.01_dp]
        do k = 1, np
            angle = 2.0_dp*pi*real(k - 1, dp)/real(np, dp)
            vertices(:, 1 + k) = 0.4_dp*[cos(angle + 0.1_dp), &
                1.1_dp*sin(angle + 0.1_dp)]
            vertices(:, 1 + np + k) = (0.8_dp + 0.04_dp*cos(3.0_dp*angle))* &
                [cos(angle), 1.1_dp*sin(angle)]
            inodes(k) = 1 + np + k
            pnodes(:, k) = [k, modulo(k, np) + 1]
        end do
        do k = 1, np
            next = modulo(k, np) + 1
            triangles(:, 3*k - 2) = [1, 1 + k, 1 + next]
            triangles(:, 3*k - 1) = [1 + k, 1 + np + k, 1 + np + next]
            triangles(:, 3*k) = [1 + k, 1 + np + next, 1 + next]
            pinc(k) = cmplx(cos(0.5_dp*k), 0.1_dp*sin(1.0_dp*k), dp)
        end do
        do k = 1, 2*nv
            load(k) = cmplx(0.1_dp*sin(0.3_dp*k), 0.05_dp, dp)
        end do
        dofs = [2]
        values = [(0.01_dp, 0.0_dp)]
    end subroutine build_problem

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_elasticity_curved_acoustic_ntd_derivatives
