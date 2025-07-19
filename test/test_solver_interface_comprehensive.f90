program test_solver_interface_comprehensive
    ! Comprehensive unit tests for solver interface
    use fortfem_kinds
    use fortfem_solver_interface
    use fortfem_sparse_matrix
    implicit none
    
    logical :: all_passed = .true.
    
    print *, "=== Testing Solver Interface (Comprehensive) ==="
    print *, ""
    
    call test_lapack_dense_solver()
    call test_sparse_solver_conversion()
    call test_solver_accuracy()
    call test_solver_edge_cases()
    call test_multiple_rhs()
    call test_solver_performance()
    
    print *, ""
    if (all_passed) then
        print *, "✅ ALL SOLVER INTERFACE TESTS PASSED!"
    else
        print *, "❌ SOME SOLVER INTERFACE TESTS FAILED!"
        stop 1
    end if
    
contains

    subroutine test_lapack_dense_solver()
        type(lapack_dense_solver_t) :: solver
        type(csr_matrix_t) :: A_csr
        type(triplet_matrix_t) :: A_triplet
        real(dp), allocatable :: b(:), x(:), x_exact(:)
        integer :: n, i, info
        real(dp) :: error
        real(dp), parameter :: tol = 1e-12_dp
        logical :: passed = .true.
        
        print *, "Testing LAPACK dense solver..."
        
        n = 5
        allocate(b(n), x(n), x_exact(n))
        
        ! Create a simple SPD system: tridiagonal with diagonal dominance
        call A_triplet%init(n, 3*n)
        do i = 1, n
            call A_triplet%add(i, i, 4.0_dp)
            if (i > 1) call A_triplet%add(i, i-1, -1.0_dp)
            if (i < n) call A_triplet%add(i, i+1, -1.0_dp)
        end do
        
        ! Convert to CSR
        call A_triplet%to_csr(A_csr)
        
        ! Set exact solution and compute RHS
        x_exact = [(real(i, dp), i=1,n)]
        call A_csr%matvec(x_exact, b)
        
        ! Solve
        call solver%init()
        call solver%solve(A_csr, b, x, info)
        
        if (info /= 0) then
            print *, "  ❌ LAPACK solver failed with info =", info
            passed = .false.
        else
            ! Check solution
            error = sqrt(sum((x - x_exact)**2))
            if (error > tol) then
                print *, "  ❌ Solution error too large:", error
                passed = .false.
            end if
        end if
        
        if (passed) then
            print *, "  ✅ LAPACK dense solver test passed"
        else
            all_passed = .false.
        end if
        
        ! Clean up
        deallocate(b, x, x_exact)
        call solver%destroy()
        call A_csr%destroy()
        call A_triplet%destroy()
    end subroutine test_lapack_dense_solver
    
    subroutine test_sparse_solver_conversion()
        type(csr_matrix_t) :: A_csr
        type(triplet_matrix_t) :: A_triplet
        real(dp), allocatable :: values(:)
        integer, allocatable :: row_ptr(:), col_idx(:)
        integer :: n, nnz, i
        logical :: passed = .true.
        
        print *, "Testing sparse matrix conversion..."
        
        n = 4
        call A_triplet%init(n, 10)
        
        ! Add entries in random order
        call A_triplet%add(2, 1, 2.0_dp)
        call A_triplet%add(1, 1, 1.0_dp)
        call A_triplet%add(3, 3, 3.0_dp)
        call A_triplet%add(1, 2, 1.5_dp)
        call A_triplet%add(4, 4, 4.0_dp)
        call A_triplet%add(2, 2, 2.5_dp)
        call A_triplet%add(3, 2, 3.5_dp)
        
        ! Convert to CSR
        call A_triplet%to_csr(A_csr)
        
        ! Check CSR structure
        if (A_csr%n /= n) then
            print *, "  ❌ Wrong matrix size:", A_csr%n
            passed = .false.
        end if
        
        if (A_csr%nnz /= 7) then
            print *, "  ❌ Wrong number of non-zeros:", A_csr%nnz
            passed = .false.
        end if
        
        ! Check row pointers are increasing
        do i = 1, n
            if (A_csr%row_ptr(i) > A_csr%row_ptr(i+1)) then
                print *, "  ❌ Row pointers not increasing!"
                passed = .false.
                exit
            end if
        end do
        
        if (passed) then
            print *, "  ✅ Sparse matrix conversion test passed"
        else
            all_passed = .false.
        end if
        
        call A_csr%destroy()
        call A_triplet%destroy()
    end subroutine test_sparse_solver_conversion
    
    subroutine test_solver_accuracy()
        type(lapack_dense_solver_t) :: solver
        type(csr_matrix_t) :: A_csr
        type(triplet_matrix_t) :: A_triplet
        real(dp), allocatable :: b(:), x(:), r(:)
        integer :: n, i, j, info
        real(dp) :: residual_norm, b_norm
        real(dp), parameter :: tol = 1e-10_dp
        logical :: passed = .true.
        
        print *, "Testing solver accuracy..."
        
        ! Test with Hilbert matrix (ill-conditioned)
        n = 6
        allocate(b(n), x(n), r(n))
        
        call A_triplet%init(n, n*n)
        do i = 1, n
            do j = 1, n
                call A_triplet%add(i, j, 1.0_dp / real(i + j - 1, dp))
            end do
        end do
        
        call A_triplet%to_csr(A_csr)
        
        ! RHS with known solution x = [1,1,...,1]
        b = 0.0_dp
        do i = 1, n
            do j = 1, n
                b(i) = b(i) + 1.0_dp / real(i + j - 1, dp)
            end do
        end do
        b_norm = sqrt(sum(b**2))
        
        ! Solve
        call solver%init()
        call solver%solve(A_csr, b, x, info)
        
        if (info /= 0) then
            print *, "  ❌ Solver failed on Hilbert matrix!"
            passed = .false.
        else
            ! Compute residual: r = b - Ax
            call A_csr%matvec(x, r)
            r = b - r
            residual_norm = sqrt(sum(r**2))
            
            if (residual_norm / b_norm > tol) then
                print *, "  ❌ Relative residual too large:", residual_norm / b_norm
                passed = .false.
            end if
        end if
        
        if (passed) then
            print *, "  ✅ Solver accuracy test passed"
        else
            all_passed = .false.
        end if
        
        deallocate(b, x, r)
        call solver%destroy()
        call A_csr%destroy()
        call A_triplet%destroy()
    end subroutine test_solver_accuracy
    
    subroutine test_solver_edge_cases()
        type(lapack_dense_solver_t) :: solver
        type(csr_matrix_t) :: A_csr
        type(triplet_matrix_t) :: A_triplet
        real(dp), allocatable :: b(:), x(:)
        integer :: info
        logical :: passed = .true.
        
        print *, "Testing solver edge cases..."
        
        ! Test 1: 1x1 system
        allocate(b(1), x(1))
        call A_triplet%init(1, 1)
        call A_triplet%add(1, 1, 2.0_dp)
        call A_triplet%to_csr(A_csr)
        
        b(1) = 4.0_dp
        call solver%init()
        call solver%solve(A_csr, b, x, info)
        
        if (info /= 0 .or. abs(x(1) - 2.0_dp) > 1e-14_dp) then
            print *, "  ❌ 1x1 system test failed!"
            passed = .false.
        end if
        
        call A_csr%destroy()
        call A_triplet%destroy()
        call solver%destroy()
        deallocate(b, x)
        
        ! Test 2: Diagonal system
        allocate(b(5), x(5))
        call A_triplet%init(5, 5)
        do info = 1, 5
            call A_triplet%add(info, info, real(info, dp))
            b(info) = real(info**2, dp)
        end do
        call A_triplet%to_csr(A_csr)
        
        call solver%init()
        call solver%solve(A_csr, b, x, info)
        
        if (info /= 0) then
            print *, "  ❌ Diagonal system test failed!"
            passed = .false.
        else
            ! Check x(i) = i
            if (any(abs(x - [(real(info, dp), info=1,5)]) > 1e-14_dp)) then
                print *, "  ❌ Diagonal system solution incorrect!"
                passed = .false.
            end if
        end if
        
        if (passed) then
            print *, "  ✅ Solver edge cases test passed"
        else
            all_passed = .false.
        end if
        
        deallocate(b, x)
        call solver%destroy()
        call A_csr%destroy()
        call A_triplet%destroy()
    end subroutine test_solver_edge_cases
    
    subroutine test_multiple_rhs()
        type(lapack_dense_solver_t) :: solver
        type(csr_matrix_t) :: A_csr
        type(triplet_matrix_t) :: A_triplet
        real(dp), allocatable :: b1(:), b2(:), x1(:), x2(:)
        integer :: n, i, info1, info2
        logical :: passed = .true.
        
        print *, "Testing multiple right-hand sides..."
        
        n = 4
        allocate(b1(n), b2(n), x1(n), x2(n))
        
        ! Create SPD matrix
        call A_triplet%init(n, 3*n)
        do i = 1, n
            call A_triplet%add(i, i, 3.0_dp)
            if (i > 1) call A_triplet%add(i, i-1, -0.5_dp)
            if (i < n) call A_triplet%add(i, i+1, -0.5_dp)
        end do
        call A_triplet%to_csr(A_csr)
        
        ! Two different RHS vectors
        b1 = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        b2 = [4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
        
        ! Solve both systems
        call solver%init()
        call solver%solve(A_csr, b1, x1, info1)
        call solver%solve(A_csr, b2, x2, info2)
        
        if (info1 /= 0 .or. info2 /= 0) then
            print *, "  ❌ Multiple RHS solving failed!"
            passed = .false.
        end if
        
        ! Solutions should be different
        if (sqrt(sum((x1 - x2)**2)) < 1e-10_dp) then
            print *, "  ❌ Solutions incorrectly identical!"
            passed = .false.
        end if
        
        if (passed) then
            print *, "  ✅ Multiple RHS test passed"
        else
            all_passed = .false.
        end if
        
        deallocate(b1, b2, x1, x2)
        call solver%destroy()
        call A_csr%destroy()
        call A_triplet%destroy()
    end subroutine test_multiple_rhs
    
    subroutine test_solver_performance()
        type(lapack_dense_solver_t) :: solver
        type(csr_matrix_t) :: A_csr
        type(triplet_matrix_t) :: A_triplet
        real(dp), allocatable :: b(:), x(:)
        integer :: n, i, info
        real(dp) :: t1, t2, assembly_time, solve_time
        logical :: passed = .true.
        
        print *, "Testing solver performance..."
        
        ! Test with medium-sized system
        n = 100
        allocate(b(n), x(n))
        
        ! Time assembly
        call cpu_time(t1)
        call A_triplet%init(n, 3*n)
        do i = 1, n
            call A_triplet%add(i, i, 4.0_dp)
            if (i > 1) call A_triplet%add(i, i-1, -1.0_dp)
            if (i < n) call A_triplet%add(i, i+1, -1.0_dp)
        end do
        call A_triplet%to_csr(A_csr)
        call cpu_time(t2)
        assembly_time = t2 - t1
        
        ! Set RHS
        b = 1.0_dp
        
        ! Time solve
        call cpu_time(t1)
        call solver%init()
        call solver%solve(A_csr, b, x, info)
        call cpu_time(t2)
        solve_time = t2 - t1
        
        if (info /= 0) then
            print *, "  ❌ Performance test solver failed!"
            passed = .false.
        else
            print *, "    Assembly time:", assembly_time, "seconds"
            print *, "    Solve time:", solve_time, "seconds"
            
            ! Basic sanity check - solve shouldn't take too long
            if (solve_time > 1.0_dp) then
                print *, "  ⚠️  Solve time seems excessive for n =", n
            end if
        end if
        
        if (passed) then
            print *, "  ✅ Solver performance test passed"
        else
            all_passed = .false.
        end if
        
        deallocate(b, x)
        call solver%destroy()
        call A_csr%destroy()
        call A_triplet%destroy()
    end subroutine test_solver_performance

end program test_solver_interface_comprehensive