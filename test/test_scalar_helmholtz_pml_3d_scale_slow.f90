program test_scalar_helmholtz_pml_3d_scale_slow
    !! The stretched plane wave of test_scalar_helmholtz_pml_3d (constant stretch
    !! 1 + 0.4 i in x, exact Dirichlet values on the boundary) on a
    !! [0, 2] x [0, 0.25]**2 box with 98196 vertices, where a dense complex
    !! matrix would need 144 GiB. Vertices are numbered with x slowest, so the
    !! natural ordering is a band of about (ny + 1)**2 entries (the fpm build of
    !! fortsparse has no fill-reducing ordering).
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_scalar_helmholtz_pml_p1_3d
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 501, ny = 13, m = ny + 1
    real(dp), parameter :: k = 1.6_dp
    complex(dp), parameter :: x_stretch = (1.0_dp, 0.4_dp)
    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: tetrahedra(:, :), boundary(:)
    complex(dp), allocatable :: stretch(:, :), exact(:), load(:), solution(:)
    real(dp) :: error, seconds, peak_bytes, dense_bytes
    integer(8) :: start, finish, rate
    integer :: n, status, v
    logical :: measured

    call box(vertices, tetrahedra)
    n = size(vertices, 2)
    exact = exp((0.0_dp, 1.0_dp)*k*x_stretch*vertices(1, :))
    boundary = pack([(v, v=1, n)], vertices(1, :) < 1.0e-12_dp .or. &
        vertices(1, :) > 2.0_dp - 1.0e-12_dp .or. &
        any(vertices(2:3, :) < 1.0e-12_dp, dim=1) .or. &
        any(vertices(2:3, :) > 0.25_dp - 1.0e-12_dp, dim=1))
    allocate(stretch(3, size(tetrahedra, 2)), load(n))
    stretch(1, :) = x_stretch
    stretch(2:3, :) = (1.0_dp, 0.0_dp)
    load = (0.0_dp, 0.0_dp)
    call system_clock(start, rate)
    call solve_scalar_helmholtz_pml_p1_3d(vertices, tetrahedra, stretch, k, &
        load, boundary, exact(boundary), solution, status)
    call system_clock(finish)
    seconds = real(finish - start, dp)/real(rate, dp)
    call check_condition(status == 0, "3D PML solve on 98196 vertices succeeds")
    error = sqrt(sum(abs(solution - exact)**2)/sum(abs(exact)**2))
    dense_bytes = 16.0_dp*real(n, dp)**2
    call peak_resident_bytes(peak_bytes, measured)
    write (*, "(a,i0,a,i0,a,f6.1,a,f7.2,a,es9.2)") "vertices ", n, &
        ", tetrahedra ", size(tetrahedra, 2), ", dense would need ", &
        dense_bytes/1024.0_dp**3, " GiB, solve ", seconds, &
        " s, relative error ", error
    if (measured) write (*, "(a,f8.1,a)") "peak resident memory ", &
        peak_bytes/1024.0_dp**2, " MiB"
    call check_condition(error < 1.0e-3_dp, &
        "solution matches the stretched plane wave")
    if (measured) call check_condition(peak_bytes < 8.0_dp*1024.0_dp**3, &
        "peak resident memory below 8 GiB (no N^2 storage)")
    call check_summary("scalar_helmholtz_pml_3d_scale_slow")

contains

    subroutine box(points, cells)
        !! Vertex 1 + l + m*(j + m*i) at grid point (i, j, l); six positive
        !! tetrahedra per hexahedron.
        real(dp), allocatable, intent(out) :: points(:, :)
        integer, allocatable, intent(out) :: cells(:, :)
        integer, parameter :: kuhn(4, 6) = reshape([1, 2, 4, 8, 1, 4, 3, 8, &
            1, 3, 7, 8, 1, 7, 5, 8, 1, 5, 6, 8, 1, 6, 2, 8], [4, 6])
        integer :: corner(8), i, j, l, c, t, q
        real(dp) :: e(3, 3)

        allocate(points(3, nx*m*m), cells(4, 6*(nx - 1)*ny*ny))
        do i = 0, nx - 1
            do j = 0, ny
                do l = 0, ny
                    points(:, 1 + l + m*(j + m*i)) = [2.0_dp*real(i, dp)/ &
                        real(nx - 1, dp), 0.25_dp*real(j, dp)/real(ny, dp), &
                        0.25_dp*real(l, dp)/real(ny, dp)]
                end do
            end do
        end do
        t = 0
        do i = 0, nx - 2
            do j = 0, ny - 1
                do l = 0, ny - 1
                    c = 1 + l + m*(j + m*i)
                    corner = [c, c + m*m, c + m, c + m*m + m, c + 1, c + m*m + 1, &
                        c + m + 1, c + m*m + m + 1]
                    do q = 1, 6
                        t = t + 1
                        cells(:, t) = corner(kuhn(:, q))
                        e(:, 1) = points(:, cells(2, t)) - points(:, cells(1, t))
                        e(:, 2) = points(:, cells(3, t)) - points(:, cells(1, t))
                        e(:, 3) = points(:, cells(4, t)) - points(:, cells(1, t))
                        if (dot_product(e(:, 1), [e(2, 2)*e(3, 3) - e(3, 2)* &
                            e(2, 3), e(3, 2)*e(1, 3) - e(1, 2)*e(3, 3), &
                            e(1, 2)*e(2, 3) - e(2, 2)*e(1, 3)]) < 0.0_dp) &
                            cells(3:4, t) = cells([4, 3], t)
                    end do
                end do
            end do
        end do
    end subroutine box

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

end program test_scalar_helmholtz_pml_3d_scale_slow
