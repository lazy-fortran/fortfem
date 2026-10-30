program test_scalar_helmholtz_2d_scale
    !! The plane-wave cases of test_scalar_helmholtz_pml_2d (stretched plane wave,
    !! exact Dirichlet values on the boundary) and of
    !! test_scalar_helmholtz_planar_dtn_convergence (DtN boundary at x = 1) on
    !! about 8e4 vertices, where a dense complex matrix would need over 100 GiB.
    !! Vertices are numbered column by column (y fastest), so the natural ordering
    !! is a band of about ny entries (the fpm build of fortsparse has no
    !! fill-reducing ordering).
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_scalar_helmholtz_planar_dtn_p1, &
        solve_scalar_helmholtz_pml_p1_2d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: triangles(:, :)
    complex(dp), allocatable :: solution(:)
    real(dp) :: seconds, peak_bytes
    integer(8) :: start, finish, rate
    logical :: measured

    call pml_case()
    call dtn_case()
    call peak_resident_bytes(peak_bytes, measured)
    if (measured) then
        write (*, "(a,f8.1,a)") "peak resident memory ", &
            peak_bytes/1024.0_dp**2, " MiB"
        call check_condition(peak_bytes < 2.0_dp*1024.0_dp**3, &
            "peak resident memory below 2 GiB (no N^2 storage)")
    end if
    call check_summary("scalar_helmholtz_2d_scale")

contains

    subroutine pml_case()
        integer, parameter :: nx = 2001, ny = 41
        real(dp), parameter :: k = 2.0_dp, sigma_max = 3.0_dp
        complex(dp), allocatable :: stretch(:, :), exact(:), load(:)
        integer, allocatable :: boundary(:)
        real(dp) :: centroid, error
        integer :: e, status, v

        call strip(nx, ny, 2.0_dp, vertices, triangles)
        allocate(stretch(2, size(triangles, 2)), exact(nx*ny), load(nx*ny))
        stretch = (1.0_dp, 0.0_dp)
        do e = 1, size(triangles, 2)
            centroid = sum(vertices(1, triangles(:, e)))/3.0_dp
            if (centroid > 1.0_dp) stretch(1, e) = cmplx(1.0_dp, &
                sigma_max*(centroid - 1.0_dp)**2/k, dp)
        end do
        do v = 1, nx*ny
            exact(v) = exp((0.0_dp, 1.0_dp)*k*cmplx(vertices(1, v), &
                sigma_max*max(vertices(1, v) - 1.0_dp, 0.0_dp)**3/(3.0_dp*k), dp))
        end do
        boundary = pack([(v, v=1, nx*ny)], on_boundary(nx, ny))
        load = (0.0_dp, 0.0_dp)
        call system_clock(start, rate)
        call solve_scalar_helmholtz_pml_p1_2d(vertices, triangles, stretch, k, &
            load, boundary, exact(boundary), solution, status)
        call system_clock(finish)
        seconds = real(finish - start, dp)/real(rate, dp)
        error = sqrt(sum(abs(solution - exact)**2)/sum(abs(exact)**2))
        call report("PML", nx*ny, seconds, error)
        call check_condition(status == 0 .and. error < 1.0e-3_dp, &
            "2D PML on 82041 vertices matches the stretched plane wave")
    end subroutine pml_case

    subroutine dtn_case()
        integer, parameter :: nx = 4001, ny = 21
        real(dp), parameter :: k = 2.0_dp
        complex(dp), allocatable :: exact(:), load(:)
        integer, allocatable :: dirichlet(:), dtn(:)
        logical, allocatable :: fixed(:)
        real(dp) :: error
        integer :: status, v

        call strip(nx, ny, 1.0_dp, vertices, triangles)
        ! Height (ny - 1) h with h = 1/(nx - 1): the DtN period is ny h.
        vertices(2, :) = vertices(2, :)*real(ny - 1, dp)/real(nx - 1, dp)
        exact = exp((0.0_dp, 1.0_dp)*k*vertices(1, :))
        dtn = [((nx - 1)*ny + v, v=1, ny)]
        fixed = on_boundary(nx, ny)
        fixed(dtn(2:ny - 1)) = .false.
        dirichlet = pack([(v, v=1, nx*ny)], fixed)
        if (allocated(solution)) deallocate(solution)
        allocate(load(nx*ny), solution(nx*ny))
        load = (0.0_dp, 0.0_dp)
        call system_clock(start, rate)
        call solve_scalar_helmholtz_planar_dtn_p1(vertices, triangles, dtn, k, &
            real(ny, dp)/real(nx - 1, dp), load, dirichlet, exact(dirichlet), &
            solution, status)
        call system_clock(finish)
        seconds = real(finish - start, dp)/real(rate, dp)
        error = sqrt(sum(abs(solution - exact)**2)/sum(abs(exact)**2))
        call report("DtN", nx*ny, seconds, error)
        call check_condition(status == 0 .and. error < 1.0e-3_dp, &
            "2D DtN on 84021 vertices matches the plane wave")
    end subroutine dtn_case

    subroutine report(name, n, time, error)
        character(*), intent(in) :: name
        integer, intent(in) :: n
        real(dp), intent(in) :: time, error

        write (*, "(a,a,i0,a,f6.1,a,f7.2,a,es9.2)") name, ": vertices ", n, &
            ", dense would need ", 16.0_dp*real(n, dp)**2/1024.0_dp**3, &
            " GiB, solve ", time, " s, relative error ", error
    end subroutine report

    function on_boundary(nx, ny) result(mask)
        integer, intent(in) :: nx, ny
        logical :: mask(nx*ny)
        integer :: i, j

        do i = 0, nx - 1
            do j = 0, ny - 1
                mask(1 + j + ny*i) = i == 0 .or. i == nx - 1 .or. j == 0 .or. &
                    j == ny - 1
            end do
        end do
    end function on_boundary

    subroutine strip(nx, ny, width, points, cells)
        !! [0, width] x [0, 1]; vertex 1 + j + ny*i at grid point (i, j).
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: width
        real(dp), allocatable, intent(out) :: points(:, :)
        integer, allocatable, intent(out) :: cells(:, :)
        integer :: c, i, j, t

        allocate(points(2, nx*ny), cells(3, 2*(nx - 1)*(ny - 1)))
        do i = 0, nx - 1
            do j = 0, ny - 1
                points(:, 1 + j + ny*i) = [width*real(i, dp)/real(nx - 1, dp), &
                    real(j, dp)/real(ny - 1, dp)]
            end do
        end do
        t = 0
        do i = 0, nx - 2
            do j = 1, ny - 1
                c = j + ny*i
                cells(:, t + 1) = [c, c + ny, c + ny + 1]
                cells(:, t + 2) = [c, c + ny + 1, c + 1]
                t = t + 2
            end do
        end do
    end subroutine strip

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

end program test_scalar_helmholtz_2d_scale
