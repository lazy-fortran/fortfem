program test_elasticity_planar_acoustic_absorbing_halfspace
    !! Water plane wave at normal incidence on a steel half-space, truncated by
    !! an absorbing face. The DtN map fixes the time convention exp(-i omega t)
    !! (outgoing exp(i k_n z), multiplier -i rho omega**2/k_n), under which the
    !! Lysmer-Kuhlemeyer condition sigma n = -Z du/dt reads sigma n = i omega Z u
    !! and enters the weak form as -i omega Z (u, v)_edge.
    !!
    !! Strip 0 <= x <= L, wet face x = L, absorbing face x = 0, rollers on the
    !! lateral faces. In the solid u_x = A exp(-i k x) + B exp(i k x); the exact
    !! absorbing face gives B = 0, and with the blocked pressure P_b the wet-face
    !! displacement is U = P_b/(i omega (Z_s + Z_f)), Z = density*speed (P wave
    !! in steel). The fit of A, B uses the 1D P1 wavenumber
    !! cos(k_h h) = (1 - (kh)**2/3)/(1 + (kh)**2/6), so |B/A| measures only the
    !! reflection of the discrete absorbing face. With 80, 160 and 320 cells per
    !! P wavelength, |B/A| was 1.8e-4, 4.3e-5, 9.8e-6 and the wet-face error
    !! 3.0e-2, 7.2e-3, 1.7e-3 (second order); this test uses 160. The former
    !! +i omega Z edge gave |B/A| = 5.5e3 and a wet-face error of 2.07.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_elasticity_planar_acoustic_dtn_p1
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 340, ny = 20
    integer, parameter :: nv = (nx + 1)*(ny + 1)
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: young_modulus = 200.0e9_dp, poisson_ratio = 0.3_dp
    real(dp), parameter :: solid_density = 7850.0_dp
    real(dp), parameter :: sound_speed = 1480.0_dp, fluid_density = 1000.0_dp
    real(dp), parameter :: h = 0.5e-4_dp, thickness = nx*h
    real(dp), parameter :: blocked_pressure = 1.0_dp
    real(dp), allocatable :: vertices(:, :), normals(:, :), anormals(:, :)
    integer, allocatable :: triangles(:, :), inodes(:), dofs(:), aedges(:, :)
    complex(dp), allocatable :: values(:), load(:), pinc(:), solution(:)
    real(dp) :: modulus, p_speed, wavenumber, z_s, z_f, omega
    real(dp) :: displacement_error, reflection, interface_reflection_error
    real(dp) :: amplitude
    complex(dp) :: exact, centre, interface_reflection
    integer :: status
    logical :: all_passed

    all_passed = .true.
    modulus = young_modulus*(1.0_dp - poisson_ratio)/ &
        ((1.0_dp + poisson_ratio)*(1.0_dp - 2.0_dp*poisson_ratio))
    p_speed = sqrt(modulus/solid_density)
    wavenumber = 2.0_dp*pi*2.125_dp/thickness
    omega = wavenumber*p_speed
    z_s = solid_density*p_speed
    z_f = fluid_density*sound_speed
    call build_problem()

    call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, inodes, &
        normals, omega, sound_speed, fluid_density, ny*h, 3, young_modulus, &
        poisson_ratio, solid_density, load, pinc, dofs, values, solution, &
        status, absorbing_edges=aedges, absorbing_normals=anormals)
    call record(status == 0, "steel half-space solve succeeds")

    call fit_waves(reflection, amplitude)
    exact = blocked_pressure/(cmplx(0.0_dp, omega*(z_s + z_f), dp))
    centre = solution(2*inodes(ny/2) - 1)
    displacement_error = abs(centre - exact)/abs(exact)
    interface_reflection = (blocked_pressure - &
        (0.0_dp, 1.0_dp)*omega*z_f*centre)/(0.5_dp*blocked_pressure) - 1.0_dp
    interface_reflection_error = abs(interface_reflection - &
        (z_s - z_f)/(z_s + z_f))
    write (*, "(a,es10.3)") "absorbing-face reflection |B/A| ", reflection
    write (*, "(a,es10.3)") "transmitted amplitude |A|/|U| - 1 ", amplitude - 1
    write (*, "(a,es10.3)") "wet-face displacement relative error ", &
        displacement_error
    write (*, "(a,es10.3)") "water-steel reflection coefficient error ", &
        interface_reflection_error
    call record(reflection < 1.0e-3_dp, &
        "absorbing face reflects less than 0.1 % at normal incidence")
    call record(abs(amplitude - 1.0_dp) < 1.0e-3_dp, &
        "transmitted wave carries the half-space amplitude")
    call record(displacement_error < 1.0e-2_dp, &
        "wet-face displacement matches the half-space solution")
    call record(interface_reflection_error < 1.0e-3_dp, &
        "water-steel reflection is (Z_s - Z_f)/(Z_s + Z_f)")

    call check_summary("Steel half-space with absorbing face")
    if (.not. all_passed) error stop 1

contains

    subroutine fit_waves(ratio, amplitude)
        !! Least-squares fit of u_x along the middle row to A e^{-ikx} + B e^{ikx};
        !! ratio = |B/A|, amplitude = |A|/|U|.
        real(dp), intent(out) :: ratio, amplitude
        complex(dp) :: basis(0:nx, 2), gram(2, 2), rhs(2), samples(0:nx)
        complex(dp) :: a, b, det
        real(dp) :: q, discrete_k, x
        integer :: col

        q = (wavenumber*h)**2
        discrete_k = acos((1.0_dp - q/3.0_dp)/(1.0_dp + q/6.0_dp))/h
        do col = 0, nx
            x = h*real(col, dp)
            basis(col, 1) = exp(cmplx(0.0_dp, -discrete_k*x, dp))
            basis(col, 2) = exp(cmplx(0.0_dp, discrete_k*x, dp))
            samples(col) = solution(2*(1 + ny/2 + (ny + 1)*col) - 1)
        end do
        gram = matmul(conjg(transpose(basis)), basis)
        rhs = matmul(conjg(transpose(basis)), samples)
        det = gram(1, 1)*gram(2, 2) - gram(1, 2)*gram(2, 1)
        a = (gram(2, 2)*rhs(1) - gram(1, 2)*rhs(2))/det
        b = (gram(1, 1)*rhs(2) - gram(2, 1)*rhs(1))/det
        ratio = abs(b)/abs(a)
        amplitude = abs(a)*omega*(z_s + z_f)/blocked_pressure
    end subroutine fit_waves

    subroutine build_problem()
        !! Vertex (column, row) is 1 + row + (ny + 1) column.
        integer :: col, row, t, ll, k

        allocate(vertices(2, nv), triangles(3, 2*nx*ny))
        do col = 0, nx
            do row = 0, ny
                vertices(:, 1 + row + (ny + 1)*col) = h*real([col, row], dp)
            end do
        end do
        t = 0
        do col = 0, nx - 1
            do row = 0, ny - 1
                ll = 1 + row + (ny + 1)*col
                t = t + 1
                triangles(:, t) = [ll, ll + ny + 1, ll + ny + 2]
                t = t + 1
                triangles(:, t) = [ll, ll + ny + 2, ll + 1]
            end do
        end do
        allocate(inodes(ny), normals(2, ny), pinc(ny))
        do row = 0, ny - 1
            inodes(row + 1) = 1 + row + (ny + 1)*nx
        end do
        normals = 0.0_dp
        normals(1, :) = 1.0_dp
        pinc = blocked_pressure
        allocate(aedges(2, ny), anormals(2, ny))
        do row = 0, ny - 1
            aedges(:, row + 1) = [1 + row, 2 + row]
        end do
        anormals = 0.0_dp
        anormals(1, :) = -1.0_dp
        allocate(dofs(2*(nx + 1)))
        k = 0
        do col = 0, nx
            k = k + 1
            dofs(k) = 2*(1 + (ny + 1)*col)
            k = k + 1
            dofs(k) = 2*(1 + ny + (ny + 1)*col)
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

end program test_elasticity_planar_acoustic_absorbing_halfspace
