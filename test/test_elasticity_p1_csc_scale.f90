program test_elasticity_p1_csc_scale
    !! Vector P1 elasticity on 2.0e5 dofs, where a dense complex matrix would need
    !! about 650 GB: the sparse assembly runs in O(N) memory and nnz grows linearly
    !! with the number of dofs.
    use, intrinsic :: iso_fortran_env, only: int64
    use check, only: check_condition, check_summary
    use fortfem_block_csc, only: block_csc_pattern_t
    use fortfem_elasticity_p1_csc_2d, only: elasticity_p1_assemble, &
        elasticity_p1_pattern
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx_small = 160, nx_large = 318
    real(dp) :: dense_bytes, ratio_small, ratio_large, peak_bytes
    integer(int64) :: start, finish, rate
    integer :: n_large, nnz_large, nnz_small, status(2)
    logical :: measured

    call system_clock(start, rate)
    call assemble(nx_small, nnz_small, status(1))
    call assemble(nx_large, nnz_large, status(2))
    call system_clock(finish)
    n_large = 2*nx_large**2
    dense_bytes = 16.0_dp*real(n_large, dp)**2
    ratio_small = real(nnz_small, dp)/real(2*nx_small**2, dp)
    ratio_large = real(nnz_large, dp)/real(n_large, dp)
    write (*, "(a,i0,a,i0,a,es9.2,a)") "dofs ", n_large, ", nnz ", nnz_large, &
        ", dense would need ", dense_bytes/1024.0_dp**3, " GiB"
    write (*, "(a,f8.3,a,f8.3,a,f6.2,a)") "nnz per dof ", ratio_small, " -> ", &
        ratio_large, ", both assemblies ", real(finish - start, dp)/real(rate, dp), &
        " s"
    call check_condition(all(status == 0), "both assemblies succeed")
    call check_condition(n_large >= 200000 .and. dense_bytes > 100.0_dp*1024.0_dp**3, &
        "problem size: >= 2e5 dofs, dense matrix above 100 GiB")
    ! Interior rows of this mesh hold 14 entries; boundary rows fewer.
    call check_condition(ratio_large <= 14.0_dp .and. &
        abs(ratio_large - ratio_small) <= 0.02_dp*ratio_large, &
        "nnz grows linearly with the number of dofs")
    call peak_resident_bytes(peak_bytes, measured)
    if (measured) then
        write (*, "(a,f8.1,a)") "peak resident memory ", peak_bytes/1024.0_dp**2, &
            " MiB"
        call check_condition(peak_bytes < 2.0_dp*1024.0_dp**3, &
            "peak resident memory below 2 GiB (no N^2 storage)")
    end if
    call check_summary("elasticity_p1_csc_scale")

contains

    subroutine assemble(nx, nnz, status)
        integer, intent(in) :: nx
        integer, intent(out) :: nnz, status
        type(block_csc_pattern_t) :: pattern
        real(dp), allocatable :: vertices(:, :)
        integer, allocatable :: triangles(:, :)
        complex(dp), allocatable :: values(:)
        integer :: cell_group, i, j, k, nt

        allocate(vertices(2, nx*nx), triangles(3, 2*(nx - 1)**2))
        do j = 1, nx
            do i = 1, nx
                vertices(:, i + (j - 1)*nx) = [real(i - 1, dp), real(j - 1, dp)]/ &
                    real(nx - 1, dp)
            end do
        end do
        k = 0
        do j = 1, nx - 1
            do i = 1, nx - 1
                triangles(:, k + 1) = [i + (j - 1)*nx, i + 1 + (j - 1)*nx, i + 1 + j*nx]
                triangles(:, k + 2) = [i + (j - 1)*nx, i + 1 + j*nx, i + j*nx]
                k = k + 2
            end do
        end do
        nt = size(triangles, 2)
        call elasticity_p1_pattern(nx*nx, triangles, pattern, cell_group)
        allocate(values(pattern%nnz))
        values = (0.0_dp, 0.0_dp)
        call elasticity_p1_assemble(pattern, cell_group, vertices, triangles, &
            spread((1.0_dp, 0.1_dp), 1, nt), spread((0.5_dp, 0.05_dp), 1, nt), &
            spread((3.0_dp, 0.0_dp), 1, nt), values, status)
        nnz = pattern%nnz
    end subroutine assemble

    subroutine peak_resident_bytes(bytes, measured)
        !! VmHWM from /proc/self/status; unavailable outside Linux.
        real(dp), intent(out) :: bytes
        logical, intent(out) :: measured
        character(256) :: line
        integer :: ios, kib, unit

        bytes = 0.0_dp
        measured = .false.
        open (newunit=unit, file="/proc/self/status", status="old", &
            action="read", iostat=ios)
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

end program test_elasticity_p1_csc_scale
