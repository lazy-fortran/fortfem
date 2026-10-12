program fortfem_csr_spmv_scaling
    !! Strong scaling of the row-parallel CSR product and of the CSR assembly on a
    !! synthetic P3 triangle lattice (exchange form, random zero-row-sum blocks).
    !! Usage: fortfem_csr_spmv_scaling [cells_per_side] [max_threads]
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_element_to_csr, only: csr_assemble_blocks, csr_build_pattern
    use fortfem_kinds, only: dp
    use fortfem_sparse_matrix, only: csr_matrix_t, spmv
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    implicit none

    integer, parameter :: degree = 3, nloc = 10
    integer, allocatable :: conn(:, :)
    real(dp), allocatable :: blocks(:, :, :), x(:), y(:), y_ref(:)
    type(csr_matrix_t) :: matrix
    integer(int64) :: t0, t1, rate
    real(dp) :: assemble, bandwidth, pattern, seconds, serial
    character(len=32) :: text
    integer :: i, n_dofs, nx, repeat, repeats, threads, top

    nx = 300
    top = 1
    !$ top = omp_get_max_threads()
    if (command_argument_count() >= 1) then
        call get_command_argument(1, text)
        read (text, *) nx
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, text)
        read (text, *) top
    end if
    call build_lattice(nx, conn, blocks)
    n_dofs = maxval(conn)
    allocate(x(n_dofs), y(n_dofs), y_ref(n_dofs))
    do i = 1, n_dofs
        x(i) = sin(0.001_dp*real(i, dp)) + 0.3_dp*cos(0.37_dp*real(i, dp))
    end do

    write (*, "(a,i0,a,i0)") "P3 lattice nx=", nx, " dofs=", n_dofs
    threads = 1
    do
        !$ call omp_set_num_threads(threads)
        call system_clock(t0, rate)
        call csr_build_pattern(n_dofs, conn, .true., matrix)
        call system_clock(t1)
        pattern = real(t1 - t0, dp)/real(rate, dp)
        call system_clock(t0)
        call csr_assemble_blocks(matrix, conn, blocks)
        call system_clock(t1)
        assemble = real(t1 - t0, dp)/real(rate, dp)
        if (threads == 1) then
            write (*, "(a,i0,a,f8.1,a)") "nnz=", size(matrix%val), "  matrix ", &
                real(size(matrix%val), dp)*12.0e-6_dp, " MB"
            write (*, "(a)") "threads  pattern[s]  assemble[s]  spmv[ms]  speedup  GB/s" &
                //"  bitwise_equal_to_1_thread"
            repeats = max(5, min(500, int(4.0e9_dp/real(size(matrix%val), dp))))
        end if
        call spmv(matrix, x, y)
        seconds = huge(1.0_dp)
        do repeat = 1, repeats
            call system_clock(t0)
            call spmv(matrix, x, y)
            call system_clock(t1)
            seconds = min(seconds, real(t1 - t0, dp)/real(rate, dp))
        end do
        if (threads == 1) then
            serial = seconds
            y_ref = y
        end if
        bandwidth = 1.0e-9_dp*(12.0_dp*real(size(matrix%val), dp) &
            + 16.0_dp*real(n_dofs, dp))/seconds
        write (*, "(i7,f11.3,f13.3,f10.3,f9.2,f7.1,l8)") threads, pattern, assemble, &
            1.0e3_dp*seconds, serial/seconds, bandwidth, all(y == y_ref)
        if (threads >= top) exit
        threads = min(2*threads, top)
        if (threads == 6) threads = 8
    end do

contains

    subroutine build_lattice(nx, conn, blocks)
        integer, intent(in) :: nx
        integer, allocatable, intent(out) :: conn(:, :)
        real(dp), allocatable, intent(out) :: blocks(:, :, :)

        integer :: a, b, e, half, i, ia, ib, j, lattice, m, p2, p3, q2, q3, side
        real(dp) :: weight

        lattice = degree*nx + 1
        allocate(conn(nloc, 2*nx*nx), blocks(nloc, nloc, 2*nx*nx))
        e = 0
        do j = 0, nx - 1
            do i = 0, nx - 1
                do half = 0, 1
                    e = e + 1
                    side = 1 - 2*half
                    m = 0
                    do p2 = 0, degree
                        do p3 = 0, degree - p2
                            m = m + 1
                            ia = i*degree + half*degree + side*p2
                            ib = j*degree + half*degree + side*p3
                            conn(m, e) = ia + lattice*ib + 1
                        end do
                    end do
                    ! Symmetric block with zero row sums (a graph Laplacian).
                    blocks(:, :, e) = 0.0_dp
                    do a = 1, nloc
                        do b = a + 1, nloc
                            q2 = mod(7*a + 13*b + e, 11)
                            q3 = mod(3*a + 5*b + e, 7)
                            weight = 0.1_dp + 0.05_dp*real(q2, dp) + 0.01_dp*real(q3, dp)
                            blocks(a, b, e) = -weight
                            blocks(b, a, e) = -weight
                            blocks(a, a, e) = blocks(a, a, e) + weight
                            blocks(b, b, e) = blocks(b, b, e) + weight
                        end do
                    end do
                end do
            end do
        end do
    end subroutine build_lattice

end program fortfem_csr_spmv_scaling
