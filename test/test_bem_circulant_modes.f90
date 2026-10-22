program test_bem_circulant_modes
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_bem_circulant, only: circulant_matrix_t
    use fortnum_linalg, only: dense_solve
    use check, only: check_condition, check_summary
    implicit none
    type(circulant_matrix_t) :: op
    real(dp), allocatable :: dense(:, :), rhs(:), x(:), reference(:), y(:)
    integer :: n, i, s, status
    logical :: passed
    passed = .true.
    do n = 1, 13
        allocate (op%blocks(2, 2, n), op%row_dofs(2, n), op%column_dofs(2, n))
        op%blocks = 0._dp
        do s = 1, n
            ! Nonsymmetric blocks make the Fourier sign independently visible.
            op%blocks(1, 1, s) = .2_dp*s
            op%blocks(2, 1, s) = -.1_dp*s
            op%blocks(1, 2, s) = .03_dp*s
            op%blocks(2, 2, s) = .1_dp*s
            op%row_dofs(1, s) = s
            op%row_dofs(2, s) = n + s
        end do
        op%blocks(1, 1, 1) = op%blocks(1, 1, 1) + 20._dp
        op%blocks(2, 2, 1) = op%blocks(2, 2, 1) + 20._dp
        op%column_dofs = op%row_dofs
        call op%expand(dense, status)
        call record(status == 0, 'synthetic expansion succeeds')
        allocate (rhs(2*n), x(2*n), reference(2*n), y(2*n))
        do i = 1, 2*n
            rhs(i) = sin(real(i*i, dp))
        end do
        call dense_solve(dense, rhs, reference, status)
        call record(status == 0, 'independent dense solve succeeds')
        call op%solve(rhs, x, status)
        call record(status == 0, 'arbitrary-length harmonic solve succeeds')
        call record(norm2(x - reference) < 1.e-12_dp, &
                    'all harmonics match the nonsymmetric dense oracle')
        call op%apply(x, y, status)
        call record(status == 0, 'circulant apply succeeds')
        call record(norm2(y - rhs) < 1.e-12_dp, 'solve residual matches original RHS')
        deallocate (op%blocks, op%row_dofs, op%column_dofs, dense, rhs, x, reference, y)
    end do
    allocate (rhs(1), x(1))
    rhs = 1._dp
    call op%solve(rhs, x, status)
    call record(status == 1, 'uninitialized operators are rejected')
    allocate (op%blocks(1, 1, 1), op%row_dofs(1, 1), op%column_dofs(1, 1))
    op%blocks = 0._dp; op%row_dofs = 1; op%column_dofs = 1
    call op%solve(rhs, x, status)
    call record(status /= 0, 'singular harmonic systems return failure')
    op%row_dofs = 0
    call op%apply(rhs, x, status)
    call record(status == 1, 'invalid DOF maps are rejected before indexing')
    call check_summary('Circulant harmonic systems')
    if (.not. passed) error stop 1
contains
    subroutine record(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok, label)
    end subroutine
end program
