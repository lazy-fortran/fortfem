program fortfem_element_precision_scaling
    use, intrinsic :: iso_fortran_env, only: real32, real64, int64
    use fortfem_element_batch, only: p1_element_batch
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, &
                                      csr_assemble_batched
!$  use omp_lib, only: omp_get_max_threads
    implicit none
    integer, parameter :: dp = real64, sp = real32, nel = 524288, repeats = 30
    real(dp), allocatable :: jac(:, :), diffusion(:), reaction(:), blocks(:, :, :)
    real(sp), allocatable :: jac32(:, :), d32(:), r32(:), blocks32(:, :, :)
    integer, allocatable :: status(:), conn(:, :)
    type(csr_matrix_t) :: matrix
    integer(int64) :: first, last, rate
    real(dp) :: elapsed, checksum, expected, seconds(5), error
    integer :: e, trial, pass, variant, threads, i, j, iterations
    real(dp) :: temporary

    allocate (jac(nel, 4), diffusion(nel), reaction(nel), blocks(nel, 3, 3))
    allocate (jac32(nel, 4), d32(nel), r32(nel), blocks32(nel, 3, 3), status(nel))
    !$omp parallel do
    do e = 1, nel
        jac(e, :) = [1.0_dp, 0.125_dp, -0.25_dp, 0.75_dp]
        diffusion(e) = 1; reaction(e) = 0.25_dp
        jac32(e, :) = real(jac(e, :), sp); d32(e) = 1; r32(e) = 0.25_sp
    end do
    !$omp end parallel do
    allocate (conn(3, nel))
    do e = 1, nel
        conn(:, e) = [e, e + 1, e + 2]
    end do
    call csr_build_pattern(nel + 2, conn, .false., matrix)
    threads = 1
!$  threads = omp_get_max_threads()
    ! Sum of all entries is the exact integral of the reaction coefficient.
    expected = real(nel, dp)*0.25_dp*0.78125_dp/2
  write (*, '(a)') 'variant,threads,elements,median_seconds,Mcells_per_second,rel_error'
    do variant = 1, 5
        iterations = repeats
        if (variant >= 4) iterations = 3
        do trial = 1, 5
            call system_clock(first, rate)
            do pass = 1, iterations
                select case (variant)
                case (1)
                    call p1_element_batch(jac, diffusion, reaction, blocks, status)
                case (2)
                    call p1_element_batch(jac32, d32, r32, blocks32, status)
                case (3)
                    call p1_element_batch(jac32, d32, r32, blocks, status)
                case (4)
                    call csr_assemble_batched(matrix, conn, blocks)
                case (5)
                    call csr_assemble_batched(matrix, conn, blocks32)
                end select
            end do
            call system_clock(last)
            seconds(trial) = real(last - first, dp)/real(rate, dp)/iterations
            if (any(status /= 0)) error stop 'benchmark invalid cells'
        end do
        do i = 2, 5
            temporary = seconds(i); j = i - 1
            do while (j >= 1)
                if (seconds(j) <= temporary) exit
                seconds(j + 1) = seconds(j); j = j - 1
            end do
            seconds(j + 1) = temporary
        end do
        elapsed = seconds(3)
        if (variant >= 4) then
            checksum = sum(matrix%val)
        else if (variant == 2) then
            checksum = sum(real(blocks32, dp))
        else
            checksum = sum(blocks)
        end if
        error = abs(checksum - expected)/expected
        if (error > 1.0e-4_dp) error stop 'benchmark checksum failed'
        write (*, '(i0,a,i0,a,i0,3(a,es14.6))') variant, ',', threads, ',', nel, &
            ',', elapsed, ',', real(nel, dp)/elapsed/1.0e6_dp, ',', error
    end do
end program fortfem_element_precision_scaling
