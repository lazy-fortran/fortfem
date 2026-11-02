program test_elasticity_planar_acoustic_dtn_steel_slow
    !! Normal incidence of a water plane wave on a steel layer of 2.125 P-wave
    !! lengths with a traction-free back face, solved with 202202 dofs (a dense
    !! matrix would need 650 GB). With u_x = A cos(k x), x = 0 the back face and
    !! x = L the wet face, and the blocked pressure P_b (twice the incident
    !! amplitude) the interface displacement is
    !!     U = P_b cos(kL)/(omega (Z_s sin(kL) + i Z_f cos(kL))),
    !! Z = density*speed (P wave in steel). The layer is lossless, so the
    !! reflection coefficient R = (P_b - i omega Z_f U)/(P_b/2) - 1 has |R| = 1.
    !! Rollers (u_y = 0) on the lateral faces keep the exact solution
    !! one-dimensional; the interface samples rows 0..ny-1 of the periodic cell.
    !! The error budget is the FE phase error, about 2 kL (kh)**2/24 = 2e-4, plus
    !! the lateral-edge mismatch of the periodic DtN sampling (the top row has no
    !! fluid load), measured as 0.7e-3 to 1.3e-3 for ny = 50, 100, 200.
    !! The vertices are numbered column by column, so the natural ordering is a
    !! band of 2(ny + 1) dofs (the fpm build of fortsparse has no fill-reducing
    !! ordering) and the interface dofs are contiguous.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_elasticity_planar_acoustic_dtn_p1
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 1000, ny = 100
    integer, parameter :: nv = (nx + 1)*(ny + 1)
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: young_modulus = 200.0e9_dp, poisson_ratio = 0.3_dp
    real(dp), parameter :: solid_density = 7850.0_dp
    real(dp), parameter :: sound_speed = 1480.0_dp, fluid_density = 1000.0_dp
    real(dp), parameter :: h = 1.0e-4_dp, thickness = nx*h
    real(dp), parameter :: blocked_pressure = 1.0_dp
    real(dp), allocatable :: vertices(:, :), normals(:, :)
    integer, allocatable :: triangles(:, :), inodes(:), dofs(:)
    complex(dp), allocatable :: values(:), load(:), pinc(:), solution(:)
    real(dp) :: modulus, p_speed, wavenumber, z_s, z_f, omega, seconds
    real(dp) :: peak_bytes, displacement_error, spread
    complex(dp) :: exact, centre, reflection
    integer(8) :: start, finish, rate
    integer :: status
    logical :: all_passed, measured

    all_passed = .true.
    modulus = young_modulus*(1.0_dp - poisson_ratio)/ &
        ((1.0_dp + poisson_ratio)*(1.0_dp - 2.0_dp*poisson_ratio))
    p_speed = sqrt(modulus/solid_density)
    wavenumber = 2.0_dp*pi*2.125_dp/thickness
    omega = wavenumber*p_speed
    z_s = solid_density*p_speed
    z_f = fluid_density*sound_speed
    call build_problem()

    call system_clock(start, rate)
    ! All modes |m| >= 1 are evanescent (period < fluid wavelength).
    call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, inodes, &
        normals, omega, sound_speed, fluid_density, ny*h, 3, young_modulus, &
        poisson_ratio, solid_density, load, pinc, dofs, values, solution, status)
    call system_clock(finish)
    seconds = real(finish - start, dp)/real(rate, dp)
    call record(status == 0, "steel-water solve with 202202 dofs succeeds")

    exact = blocked_pressure*cos(wavenumber*thickness)/(omega*cmplx( &
        z_s*sin(wavenumber*thickness), z_f*cos(wavenumber*thickness), dp))
    centre = solution(2*inodes(ny/2) - 1)
    displacement_error = abs(centre - exact)/abs(exact)
    spread = maxval(abs(solution(2*inodes - 1) - centre))/abs(exact)
    reflection = (blocked_pressure - (0.0_dp, 1.0_dp)*omega*z_f*centre)/ &
        (0.5_dp*blocked_pressure) - 1.0_dp
    call peak_resident_bytes(peak_bytes, measured)
    write (*, "(a,i0,a,f8.2,a)") "dofs ", 2*nv, ", solve ", seconds, " s"
    if (measured) write (*, "(a,f8.1,a)") "peak resident memory ", &
        peak_bytes/1024.0_dp**2, " MiB"
    write (*, "(a,2es14.6)") "interface displacement ", centre
    write (*, "(a,2es14.6)") "analytic displacement  ", exact
    write (*, "(a,es10.3,a,es10.3,a,es10.3)") "relative error ", &
        displacement_error, ", spread along interface ", spread, &
        ", |R| - 1 ", abs(reflection) - 1.0_dp
    call record(displacement_error < 2.5e-3_dp, &
        "interface displacement matches the layer solution")
    call record(abs(abs(reflection) - 1.0_dp) < 1.0e-4_dp, &
        "lossless layer reflects all energy (|R| = 1)")
    if (measured) call record(peak_bytes < 4.0_dp*1024.0_dp**3, &
        "peak resident memory below 4 GiB (no N^2 storage)")

    call check_summary("Steel layer under normal plane-wave incidence")
    if (.not. all_passed) error stop 1

contains

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

    subroutine peak_resident_bytes(bytes, measured)
        !! VmHWM from /proc/self/status; unavailable outside Linux.
        real(dp), intent(out) :: bytes
        logical, intent(out) :: measured
        character(256) :: line
        integer :: unit, ios, kib

        bytes = 0.0_dp
        measured = .false.
        open (newunit=unit, file="/proc/self/status", action="read", &
            status="old", iostat=ios)
        if (ios /= 0) return
        do
            read (unit, "(a)", iostat=ios) line
            if (ios /= 0) exit
            if (line(1:6) /= "VmHWM:") cycle
            read (line(7:), *, iostat=ios) kib
            if (ios == 0) then
                bytes = 1024.0_dp*real(kib, dp)
                measured = .true.
            end if
            exit
        end do
        close (unit)
    end subroutine peak_resident_bytes

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_elasticity_planar_acoustic_dtn_steel_slow
