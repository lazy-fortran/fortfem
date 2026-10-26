program tensor_setup_scaling
    ! Strong scaling of build_tensor_operator and assemble_tensor_operator.
    ! usage: fortfem_lab_setup1 dim p n nphi nq repetitions
    ! The thread count comes from OMP_NUM_THREADS. Times are the minimum over the
    ! repetitions; the checksums must not change with the thread count.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_tensor_mesh, only: tensor_mesh_t, ogrid_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                                       assemble_tensor_operator
    implicit none
    type(tensor_mesh_t) :: mesh
    type(tensor_operator_t) :: op
    type(csr_matrix_t) :: matrix
    character(32) :: argument
    integer :: arg(6), i, rep, status, mode
    real(dp) :: t0, best(2, 2)
    arg = [2, 6, 4, 16, 15, 3]
    do i = 1, 6
        call get_command_argument(i, argument)
        if (len_trim(argument) > 0) read (argument, *) arg(i)
    end do
    call ogrid_tensor_mesh(arg(1), arg(3), arg(2), mesh, status, nphi=arg(4), &
                           q=arg(2) + 1)
    if (status /= 0) error stop 'mesh'
    best = huge(1.0_dp)
    ! Untimed warm-up: page-faults the large arrays before either mode is timed.
    call build_tensor_operator(mesh, coefficient, op, status, nq=arg(5))
    if (status /= 0) error stop 'operator'
    call assemble_tensor_operator(op, matrix)
    do mode = 1, 2 ! 1: serial path (parallel=.false.), 2: threaded
        do rep = 1, arg(6)
            t0 = omp_get_wtime()
            call build_tensor_operator(mesh, coefficient, op, status, nq=arg(5), &
                                       parallel=mode == 2)
            best(1, mode) = min(best(1, mode), omp_get_wtime() - t0)
            if (status /= 0) error stop 'operator'
            t0 = omp_get_wtime()
            call assemble_tensor_operator(op, matrix, parallel=mode == 2)
            best(2, mode) = min(best(2, mode), omp_get_wtime() - t0)
        end do
    end do
    print '(a,6(i0,1x))', 'dim,p,n,nphi,nq,threads=', arg(1:5), omp_get_max_threads()
    print '(a,3(i0,1x))', 'nelem,points,dofs=', mesh%nelem, &
        mesh%nelem*product(op%no), mesh%ndof
    print '(a,2es13.5)', 'build_s serial,threaded=', best(1, :)
    print '(a,2es13.5)', 'assemble_s serial,threaded=', best(2, :)
    print '(a,4es24.16)', 'checksums points,jac,metric,matrix=', sum(op%points), &
        sum(op%jac), sum(op%metric), sum(matrix%val)
contains
    subroutine coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(x)
            k(i, i) = 2.0_dp + sin(3.0_dp*x(1) + x(i))
        end do
        k(1, 2) = 0.3_dp*cos(x(1)*x(2)); k(2, 1) = k(1, 2)
        c = 1.0_dp + x(1)**2
    end subroutine
end program
