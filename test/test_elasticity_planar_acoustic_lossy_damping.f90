program test_elasticity_planar_acoustic_lossy_damping
    !! Lossy water on a lossy steel strip at normal incidence, exp(-i omega t).
    !! Loss means Im(E) < 0 and Im(c) < 0, so k = omega sqrt(rho_s/M) and
    !! k_f = omega/c have Im(k) > 0: a wave travelling away from the wet face
    !! decays with distance in both media.
    !!
    !! Solid: strip 0 <= x <= L, wet face x = L, absorbing face x = 0, rollers on
    !! the lateral faces. The transmitted wave is u_x = A exp(-i k x), so
    !! |u_x(0)|/|u_x(L)| = exp(-Im(k) L) < 1. Reversing the time convention of
    !! the absorbing edge makes x = 0 a source (the wave grows with depth);
    !! reversing it everywhere turns the solution into the conjugate of a gain
    !! problem (also growing). Reversing the DtN map alone changes the wet-face
    !! displacement U = P_b/(i omega (Z_s + Z_f)) to P_b/(i omega (Z_s - Z_f))
    !! and the fluid wavenumber recovered from the map to -k_f.
    !!
    !! The lateral faces are mirror planes, so the interface is the even
    !! extension of the wet face over period 2W: rows 0..ny, then ny-1..1. Each
    !! interior row appears twice, so fluid density and pressure are halved.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: apply_planar_acoustic_displacement_dtn, &
        solve_elasticity_planar_acoustic_dtn_p1
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 340, ny = 20, m = 2*ny
    integer, parameter :: nv = (nx + 1)*(ny + 1)
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: loss = 0.02_dp
    complex(dp), parameter :: young_modulus = &
        200.0e9_dp*cmplx(1.0_dp, -loss, dp)
    complex(dp), parameter :: sound_speed = 1480.0_dp*cmplx(1.0_dp, -loss, dp)
    real(dp), parameter :: poisson_ratio = 0.3_dp, solid_density = 7850.0_dp
    real(dp), parameter :: fluid_density = 1000.0_dp
    real(dp), parameter :: h = 0.5e-4_dp, thickness = nx*h, width = ny*h
    real(dp), parameter :: blocked_pressure = 1.0_dp
    real(dp), allocatable :: vertices(:, :), normals(:, :), anormals(:, :)
    integer, allocatable :: triangles(:, :), inodes(:), dofs(:), aedges(:, :)
    complex(dp), allocatable :: values(:), load(:), pinc(:), solution(:)
    complex(dp) :: modulus, wavenumber, z_s, z_f, exact, wet, back
    complex(dp) :: unit_displacement(m), pressure(m), fluid_wavenumber
    real(dp) :: omega, decay, predicted_decay, wet_error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    modulus = young_modulus*(1.0_dp - poisson_ratio)/ &
        ((1.0_dp + poisson_ratio)*(1.0_dp - 2.0_dp*poisson_ratio))
    omega = 2.0_dp*pi*2.125_dp/thickness*real(sqrt(modulus/solid_density), dp)
    wavenumber = omega*sqrt(solid_density/modulus)
    z_s = sqrt(solid_density*modulus)
    z_f = fluid_density*sound_speed
    call build_problem()

    call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, inodes, &
        normals, omega, sound_speed, fluid_density/2.0_dp, 2.0_dp*width, 3, &
        young_modulus, poisson_ratio, solid_density, load, pinc, dofs, values, &
        solution, status, absorbing_edges=aedges, absorbing_normals=anormals)
    call record(status == 0, "lossy steel strip solve succeeds")

    wet = solution(2*node(nx, ny/2) - 1)
    back = solution(2*node(0, ny/2) - 1)
    decay = abs(back)/abs(wet)
    predicted_decay = exp(-aimag(wavenumber)*thickness)
    exact = blocked_pressure/(cmplx(0.0_dp, omega, dp)*(z_s + z_f))
    wet_error = abs(wet - exact)/abs(exact)
    write (*, "(a,f8.5,a,f8.5)") "solid |u(0)|/|u(L)| ", decay, &
        ", exp(-Im(k) L) ", predicted_decay
    write (*, "(a,es10.3)") "wet-face displacement relative error ", wet_error
    call record(decay < 1.0_dp, "lossy solid amplitude decays with depth")
    call record(abs(decay/predicted_decay - 1.0_dp) < 1.0e-2_dp, &
        "solid decay over the strip is exp(-Im(k) L)")
    call record(wet_error < 1.0e-2_dp, &
        "wet-face displacement is P_b/(i omega (Z_s + Z_f)) with lossy Z")

    ! The map sends u_n = exp(i k_f z) at z = 0 to p = -i rho omega**2 u_n/k_f.
    unit_displacement = (1.0_dp, 0.0_dp)
    call apply_planar_acoustic_displacement_dtn(unit_displacement, omega, &
        sound_speed, fluid_density, 2.0_dp*width, 3, pressure, status)
    call record(status == 0, "uniform DtN application succeeds")
    fluid_wavenumber = (0.0_dp, -1.0_dp)*fluid_density*omega**2/pressure(1)
    write (*, "(a,2es12.4)") "fluid wavenumber from the DtN map ", &
        fluid_wavenumber
    call record(real(fluid_wavenumber, dp) > 0.0_dp .and. &
        aimag(fluid_wavenumber) > 0.0_dp, &
        "lossy fluid wave is outgoing and decays away from the interface")
    call record(abs(fluid_wavenumber - omega/sound_speed) < &
        1.0e-12_dp*abs(omega/sound_speed), "fluid wavenumber is omega/c")

    call check_summary("Lossy steel strip under lossy water")
    if (.not. all_passed) error stop 1

contains

    pure integer function node(col, row)
        integer, intent(in) :: col, row

        node = 1 + row + (ny + 1)*col
    end function node

    subroutine build_problem()
        integer :: col, row, t, ll, k

        allocate(vertices(2, nv), triangles(3, 2*nx*ny))
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
        allocate(inodes(m), normals(2, m), pinc(m))
        do row = 0, ny
            inodes(row + 1) = node(nx, row)
        end do
        do row = ny - 1, 1, -1
            inodes(2*ny - row + 1) = node(nx, row)
        end do
        normals = 0.0_dp
        normals(1, :) = 1.0_dp
        pinc = blocked_pressure/2.0_dp
        allocate(aedges(2, ny), anormals(2, ny))
        do row = 0, ny - 1
            aedges(:, row + 1) = [node(0, row), node(0, row + 1)]
        end do
        anormals = 0.0_dp
        anormals(1, :) = -1.0_dp
        allocate(dofs(2*(nx + 1)))
        k = 0
        do col = 0, nx
            k = k + 1
            dofs(k) = 2*node(col, 0)
            k = k + 1
            dofs(k) = 2*node(col, ny)
        end do
        allocate(values(size(dofs)), load(2*nv), solution(2*nv))
        values = (0.0_dp, 0.0_dp)
        load = (0.0_dp, 0.0_dp)
    end subroutine build_problem

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_elasticity_planar_acoustic_lossy_damping
