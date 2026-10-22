program fortfem_lab_mf1
    use, intrinsic :: iso_fortran_env, only: dp => real64, compiler_version, &
                                                                        compiler_options
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    use fortfem_tensor_mesh, only: tensor_mesh_t, structured_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                                       apply_tensor_operator, assemble_tensor_operator
    use fortfem_tensor_solver, only: solve_tensor_diffusion, tensor_solve_info_t
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    type(tensor_mesh_t) :: mesh
    type(tensor_operator_t) :: op, stream
    type(csr_matrix_t) :: matrix
    type(tensor_solve_info_t) :: info
    integer :: p, h, repeats, status, i, j, mode, iteration
    real(dp), allocatable :: x(:), y(:), z(:), rhs(:)
    real(dp) :: start, elapsed, flops, error
    character(40) :: arg
    p = 4; h = 4; repeats = 20; mode = 0
    call get_command_argument(1, arg); if (len_trim(arg) > 0) read (arg, *) p
    call get_command_argument(2, arg); if (len_trim(arg) > 0) read (arg, *) h
    call get_command_argument(3, arg); if (len_trim(arg) > 0) read (arg, *) repeats
    call get_command_argument(4, arg); if (len_trim(arg) > 0) read (arg, *) mode
    print '(a)', compiler_version()
    print '(a)', compiler_options()
    call structured_tensor_mesh(3, [h, h, h], p, mesh, status, map=curved)
    if (status /= 0) error stop 'benchmark mesh'
    call build_tensor_operator(mesh, coefficient, op, status)
    if (status /= 0) error stop 'benchmark operator'
    allocate (x(mesh%ndof), y(mesh%ndof), z(mesh%ndof), rhs(mesh%ndof))
    do i = 1, mesh%ndof
        x(i) = product(sin(acos(-1.0_dp)*mesh%coordinates(:, i)))
    end do
    call apply_tensor_operator(op, x, rhs)
    if (mode == 1) then
        x = 0
        start = omp_get_wtime()
        call solve_tensor_diffusion(mesh, coefficient, rhs, 1.e-10_dp, x, info, &
                                    maxiter=500)
        elapsed = omp_get_wtime() - start
        call apply_tensor_operator(op, x, y)
        error = norm2(y - rhs)/norm2(rhs)
        call report('solve', elapsed, 0.0_dp, error, info%iterations)
        if (info%status /= 0) error stop 'benchmark solve'
    else
       call build_tensor_operator(mesh, coefficient, stream, status, precompute=.false.)
        if (status /= 0) error stop 'benchmark streaming operator'
        call apply_tensor_operator(op, x, y)
        call apply_tensor_operator(stream, x, z)
        error = norm2(y - z)/norm2(y)
        if (error > 1.e-12_dp) error stop 'benchmark recompute result'
        flops = real(mesh%nelem, dp)*(48.0_dp*(p + 1)**4 + 20.0_dp*(p + 1)**3)
        start = omp_get_wtime()
        do iteration = 1, repeats
            call apply_tensor_operator(op, x, y)
        end do
        elapsed = (omp_get_wtime() - start)/repeats
        call report('precomputed', elapsed, flops, error, 0)
        start = omp_get_wtime()
        do iteration = 1, repeats
            call apply_tensor_operator(stream, x, z)
        end do
        elapsed = (omp_get_wtime() - start)/repeats
        call report('recompute', elapsed, flops, error, 0)
        call assemble_tensor_operator(op, matrix)
        z = 0
        do i = 1, matrix%n
            do j = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                z(i) = z(i) + matrix%val(j)*x(matrix%colidx(j))
            end do
        end do
        error = norm2(y - z)/norm2(y)
        if (error > 1.e-12_dp) error stop 'benchmark CSR result'
        start = omp_get_wtime()
        do iteration = 1, repeats
            !$omp parallel do default(shared) private(i,j) schedule(static)
            do i = 1, matrix%n
                z(i) = 0
                do j = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                    z(i) = z(i) + matrix%val(j)*x(matrix%colidx(j))
                end do
            end do
            !$omp end parallel do
        end do
        elapsed = (omp_get_wtime() - start)/repeats
        call report('CSR', elapsed, 2.0_dp*size(matrix%val), error, 0)
    end if
contains
    subroutine report(name, seconds, operations, err, iterations)
        character(*), intent(in) :: name
        real(dp), intent(in) :: seconds, operations, err
        integer, intent(in) :: iterations
        print '(a,1x,5i9,4es16.7)', trim(name), p, h, mesh%ndof, &
            omp_get_max_threads(), iterations, seconds, mesh%ndof/seconds, &
            operations/seconds/1.e9_dp, err
    end subroutine
    subroutine curved(xi, x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        x = xi
        x(2) = xi(2) + 0.03_dp*sin(2*acos(-1.0_dp)*xi(1))*xi(2)*(1 - xi(2))
    end subroutine
    subroutine coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(x)
            k(i, i) = 1
        end do
        c = 0
    end subroutine
end program
