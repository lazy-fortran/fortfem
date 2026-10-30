program test_elasticity_planar_acoustic_plate_impedance
    !! Steel plate between two water half-spaces at normal incidence: water on
    !! the wet face x = d through the DtN map, water behind the back face x = 0
    !! through per-edge absorbing impedances (Z_p, Z_s) = (rho_w c_w, 0). At
    !! normal incidence that edge is the exact water half-space, so the plate
    !! must reflect and transmit as the analytic three-layer solution.
    !!
    !! Oracle (time convention exp(-i omega t), incident pressure exp(-i k_w x)
    !! of unit amplitude at x = d, Z_1 = rho_w c_w, Z_2 = rho_s c_p,
    !! phi = k_p d): continuity of displacement and normal stress at both faces
    !! gives
    !!   T = 1/(cos(phi) - i/2 (Z_2/Z_1 + Z_1/Z_2) sin(phi)),
    !!   R = -i/2 (Z_2/Z_1 - Z_1/Z_2) sin(phi) T,
    !! with R the reflected pressure at x = d and T the transmitted pressure at
    !! x = 0. From the solution, with the blocked pressure P_b = 2 p_inc,
    !! 1 + R = (P_b - i omega Z_1 u_x(d))/p_inc and T = i omega Z_1 u_x(0)/p_inc.
    !! The lateral faces are rollers, sampled as the even extension of the wet
    !! face (period 2W, fluid density and load halved) as in
    !! test_elasticity_planar_acoustic_absorbing_halfspace. With 120 cells
    !! (330 per P wavelength) |R - R_exact| = 1.2e-6 and |T/T_exact - 1| = 2.4e-6.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_elasticity_planar_acoustic_dtn_p1
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 120, ny = 4
    integer, parameter :: nv = (nx + 1)*(ny + 1)
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: young_modulus = 200.0e9_dp, poisson_ratio = 0.3_dp
    real(dp), parameter :: solid_density = 7850.0_dp
    real(dp), parameter :: sound_speed = 1480.0_dp, fluid_density = 1000.0_dp
    real(dp), parameter :: thickness = 6.0e-3_dp, h = thickness/nx
    real(dp), parameter :: blocked_pressure = 1.0_dp
    real(dp) :: vertices(2, nv), normals(2, 2*ny), anormals(2, ny)
    integer :: triangles(3, 2*nx*ny), inodes(2*ny), aedges(2, ny)
    integer :: dofs(2*(nx + 1))
    complex(dp) :: impedances(2, ny), values(2*(nx + 1)), load(2*nv)
    complex(dp) :: pinc(2*ny), solution(2*nv)
    real(dp) :: modulus, p_speed, omega, z_1, z_2, phi, p_inc
    complex(dp) :: r_exact, t_exact, r_fem, t_fem, u_wet, u_back
    real(dp) :: r_error, t_error, energy_error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    modulus = young_modulus*(1.0_dp - poisson_ratio)/ &
        ((1.0_dp + poisson_ratio)*(1.0_dp - 2.0_dp*poisson_ratio))
    p_speed = sqrt(modulus/solid_density)
    phi = 2.0_dp*pi*0.3_dp
    omega = phi*p_speed/thickness
    z_1 = fluid_density*sound_speed
    z_2 = solid_density*p_speed
    t_exact = 1.0_dp/(cos(phi) - (0.0_dp, 0.5_dp)*(z_2/z_1 + z_1/z_2)*sin(phi))
    r_exact = -(0.0_dp, 0.5_dp)*(z_2/z_1 - z_1/z_2)*sin(phi)*t_exact
    p_inc = blocked_pressure/2.0_dp

    call build_problem()
    impedances(1, :) = cmplx(z_1, 0.0_dp, dp)
    impedances(2, :) = (0.0_dp, 0.0_dp)
    call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, inodes, &
        normals, omega, sound_speed, fluid_density/2.0_dp, 2.0_dp*ny*h, 3, &
        young_modulus, poisson_ratio, solid_density, load, pinc, dofs, values, &
        solution, status, absorbing_edges=aedges, absorbing_normals=anormals, &
        absorbing_impedances=impedances)
    call record(status == 0, "plate with water behind the back face solves")

    u_wet = solution(2*node(nx, ny/2) - 1)
    u_back = solution(2*node(0, ny/2) - 1)
    r_fem = (blocked_pressure - (0.0_dp, 1.0_dp)*omega*z_1*u_wet)/p_inc - 1.0_dp
    t_fem = (0.0_dp, 1.0_dp)*omega*z_1*u_back/p_inc
    r_error = abs(r_fem - r_exact)
    t_error = abs(t_fem - t_exact)/abs(t_exact)
    energy_error = abs(abs(r_fem)**2 + abs(t_fem)**2 - 1.0_dp)
    write (*, "(a,2es14.6,a,2es14.6)") "R fem ", r_fem, "  exact ", r_exact
    write (*, "(a,2es14.6,a,2es14.6)") "T fem ", t_fem, "  exact ", t_exact
    write (*, "(a,3es10.3)") "|R error|, |T error|/|T|, energy ", r_error, &
        t_error, energy_error
    call record(r_error < 1.0e-4_dp, &
        "plate reflection matches the analytic three-layer coefficient")
    call record(t_error < 1.0e-4_dp, &
        "plate transmission matches the analytic three-layer coefficient")
    call record(energy_error < 1.0e-6_dp, &
        "lossless plate conserves energy, |R|**2 + |T|**2 = 1")

    call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, inodes, &
        normals, omega, sound_speed, fluid_density/2.0_dp, 2.0_dp*ny*h, 3, &
        young_modulus, poisson_ratio, solid_density, load, pinc, dofs, values, &
        solution, status, absorbing_edges=aedges, absorbing_normals=anormals)
    t_fem = (0.0_dp, 1.0_dp)*omega*z_1*solution(2*node(0, ny/2) - 1)/p_inc
    call record(status == 0 .and. abs(t_fem - t_exact) > 0.1_dp*abs(t_exact), &
        "default steel impedance on the back face is a different problem")

    impedances(1, 1) = (1.0_dp, 0.0_dp)
    call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, inodes, &
        normals, omega, sound_speed, fluid_density/2.0_dp, 2.0_dp*ny*h, 3, &
        young_modulus, poisson_ratio, solid_density, load, pinc, dofs, values, &
        solution, status, absorbing_edges=aedges, absorbing_normals=anormals, &
        absorbing_impedances=impedances(:, 1:ny - 1))
    call record(status == 1, "impedances must match the absorbing edges")

    call check_summary("Steel plate with water behind its back face")
    if (.not. all_passed) error stop 1

contains

    pure integer function node(col, row)
        integer, intent(in) :: col, row

        node = 1 + row + (ny + 1)*col
    end function node

    subroutine build_problem()
        integer :: col, row, t, ll, k

        do col = 0, nx
            do row = 0, ny
                vertices(:, node(col, row)) = h*real([col, row], dp)
            end do
        end do
        t = 0
        do col = 0, nx - 1
            do row = 0, ny - 1
                ll = node(col, row)
                t = t + 1
                triangles(:, t) = [ll, ll + ny + 1, ll + ny + 2]
                t = t + 1
                triangles(:, t) = [ll, ll + ny + 2, ll + 1]
            end do
        end do
        do row = 0, ny
            inodes(row + 1) = node(nx, row)
        end do
        do row = ny - 1, 1, -1
            inodes(2*ny - row + 1) = node(nx, row)
        end do
        normals = 0.0_dp
        normals(1, :) = 1.0_dp
        pinc = blocked_pressure/2.0_dp
        do row = 0, ny - 1
            aedges(:, row + 1) = [node(0, row), node(0, row + 1)]
        end do
        anormals = 0.0_dp
        anormals(1, :) = -1.0_dp
        k = 0
        do col = 0, nx
            k = k + 1
            dofs(k) = 2*node(col, 0)
            k = k + 1
            dofs(k) = 2*node(col, ny)
        end do
        values = (0.0_dp, 0.0_dp)
        load = (0.0_dp, 0.0_dp)
    end subroutine build_problem

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_elasticity_planar_acoustic_plate_impedance
