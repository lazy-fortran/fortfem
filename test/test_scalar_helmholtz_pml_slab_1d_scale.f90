program test_scalar_helmholtz_pml_slab_1d_scale
    !! The PML slab of test_scalar_helmholtz_pml_slab_1d on 400001 nodes, where a
    !! dense complex matrix would need 2.6 TB: the solution matches the exact
    !! stretched-coordinate solution and the peak memory stays O(n).
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_scalar_helmholtz_pml_slab_1d
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: n = 400001
    real(dp), parameter :: k = 3.0_dp, physical_end = 1.0_dp
    real(dp), parameter :: outer_end = 2.0_dp, sigma_max = 4.0_dp
    real(dp), allocatable :: nodes(:)
    complex(dp), allocatable :: solution(:), exact(:)
    complex(dp) :: stretched_outer, stretched
    real(dp) :: error, seconds, peak_bytes, dense_bytes
    integer(8) :: start, finish, rate
    integer :: i, status
    logical :: measured

    nodes = [(outer_end*real(i - 1, dp)/real(n - 1, dp), i=1, n)]
    call system_clock(start, rate)
    call solve_scalar_helmholtz_pml_slab_1d(nodes, physical_end, k, sigma_max, &
        2, (1.0_dp, 0.0_dp), solution, status)
    call system_clock(finish)
    seconds = real(finish - start, dp)/real(rate, dp)
    call check_condition(status == 0, "slab solve on 400001 nodes succeeds")

    stretched_outer = cmplx(outer_end, sigma_max*(outer_end - physical_end)/ &
        (3.0_dp*k), dp)
    allocate(exact(n))
    do i = 1, n
        stretched = cmplx(nodes(i), 0.0_dp, dp)
        if (nodes(i) > physical_end) stretched = cmplx(nodes(i), sigma_max* &
            (nodes(i) - physical_end)**3/(3.0_dp*k*(outer_end - physical_end)**2), &
            dp)
        exact(i) = sin(k*(stretched_outer - stretched))/sin(k*stretched_outer)
    end do
    error = sqrt(sum(abs(solution - exact)**2)/sum(abs(exact)**2))
    dense_bytes = 16.0_dp*real(n, dp)**2
    call peak_resident_bytes(peak_bytes, measured)
    write (*, "(a,i0,a,es9.2,a,f7.3,a,es9.2)") "nodes ", n, &
        ", dense would need ", dense_bytes/1024.0_dp**3, " GiB, solve ", &
        seconds, " s, relative error ", error
    if (measured) write (*, "(a,f8.1,a)") "peak resident memory ", &
        peak_bytes/1024.0_dp**2, " MiB"
    ! Discretization error ~ (k h)**2 ~ 1e-10; roundoff ~ cond*eps ~ 1e-6.
    call check_condition(error < 1.0e-5_dp, "solution matches the exact slab")
    if (measured) call check_condition(peak_bytes < 1024.0_dp**3, &
        "peak resident memory below 1 GiB (no N^2 storage)")
    call check_summary("scalar_helmholtz_pml_slab_1d_scale")

contains

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

end program test_scalar_helmholtz_pml_slab_1d_scale
