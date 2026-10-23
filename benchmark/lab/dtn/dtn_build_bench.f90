program dtn_build_bench
    !! Time and verify build_exterior_dtn on synthetic dense operators.
    !! Usage: dtn_build_bench N [legacy]
    !!   N       trace size (panels = vertices = N)
    !!   legacy  1 also times the pre-fix stages (dense_solve, intrinsic matmul)
    !!           and compares both Steklov matrices (reference oracle)
    !! Operators are Fibonacci-sphere kernels: V ~ 1/(4 pi r), K ~ double-layer,
    !! M diagonal, W symmetric positive. They exercise the dense algebra only.
    !! The thread count comes from OPENBLAS_NUM_THREADS (BLAS/LAPACK) and
    !! OMP_NUM_THREADS (OpenMP loops); both are printed in the row.
    use fortfem_kinds, only: dp
    use fortfem_parametric_laplace_bem_3d, only: build_exterior_dtn, &
        bem_operators_t, bem_exterior_dtn_t
    use fortnum_linalg, only: dense_solve
    implicit none
    type(bem_operators_t) :: op
    type(bem_exterior_dtn_t) :: exterior
    real(dp), allocatable :: x(:, :), jump(:, :), flux(:, :), ref(:, :), s(:, :)
    real(dp) :: t0, t_build, t_solve, t_mat, rel, sym
    character(len=32) :: arg, nthreads_blas, nthreads_omp
    character(len=256) :: row
    integer :: n, status, legacy, length

    n = 1536
    legacy = 0
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read (arg, *) n
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read (arg, *) legacy
    end if
    call build_operators(n, op)

    call tic(t0)
    call build_exterior_dtn(op, exterior, status)
    t_build = toc(t0)
    if (status /= 0) error stop 'build_exterior_dtn failed'
    call exterior%matrix(s, status)
    if (status /= 0) error stop 'matrix failed'
    sym = maxval(abs(s - transpose(s)))/maxval(abs(s))

    t_solve = -1.0_dp
    t_mat = -1.0_dp
    rel = -1.0_dp
    if (legacy /= 0) then
        jump = 0.5_dp*op%mass - op%k
        allocate (flux(n, n))
        call tic(t0)
        call dense_solve(op%v, -jump, flux, status)
        t_solve = toc(t0)
        if (status /= 0) error stop 'legacy dense_solve failed'
        call tic(t0)
        ref = op%w - matmul(transpose(jump), flux)
        t_mat = toc(t0)
        ref = 0.5_dp*(ref + transpose(ref))
        rel = maxval(abs(s - ref))/maxval(abs(ref))
    end if

    call get_environment_variable('OPENBLAS_NUM_THREADS', nthreads_blas, length)
    if (length == 0) nthreads_blas = 'unset'
    call get_environment_variable('OMP_NUM_THREADS', nthreads_omp, length)
    if (length == 0) nthreads_omp = 'unset'
    write (*, '(a)') 'DTN_BENCH n,openblas_threads,omp_threads,t_build,' // &
        't_legacy_solve,t_legacy_matmul,rel_diff_vs_legacy,symmetry_defect'
    write (row, '(i0,5a,4(es12.4,a),es12.4)') n, ',', trim(nthreads_blas), ',', &
        trim(nthreads_omp), ',', t_build, ',', t_solve, ',', t_mat, ',', rel, &
        ',', sym
    write (*, '(2a)') 'DTN_BENCH ', trim(row)

contains

    subroutine build_operators(n, op)
        integer, intent(in) :: n
        type(bem_operators_t), intent(out) :: op
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp), allocatable :: nrm(:, :)
        real(dp) :: z, r, h, d(3), dist
        integer :: i, j

        allocate (x(3, n), nrm(3, n))
        do i = 1, n
            z = 1.0_dp - 2.0_dp*(real(i, dp) - 0.5_dp)/real(n, dp)
            r = sqrt(1.0_dp - z*z)
            nrm(:, i) = [r*cos(pi*(3.0_dp - sqrt(5.0_dp))*real(i, dp)), &
                r*sin(pi*(3.0_dp - sqrt(5.0_dp))*real(i, dp)), z]
            x(:, i) = nrm(:, i)
        end do
        h = sqrt(4.0_dp*pi/real(n, dp))
        allocate (op%v(n, n), op%k(n, n), op%mass(n, n), op%w(n, n))
        allocate (op%adjoint(n, n), op%area(n))
        op%mass = 0.0_dp
        op%area = h*h
        do j = 1, n
            do i = 1, n
                if (i == j) then
                    op%v(i, j) = h/(4.0_dp*pi)*3.0_dp
                    op%k(i, j) = 0.0_dp
                    op%w(i, j) = 1.0_dp
                    op%mass(i, j) = h*h
                else
                    d = x(:, i) - x(:, j)
                    dist = norm2(d)
                    op%v(i, j) = h*h/(4.0_dp*pi*dist)
                    op%k(i, j) = h*h*dot_product(nrm(:, j), d)/ &
                        (4.0_dp*pi*dist**3)
                    op%w(i, j) = -h*h/(4.0_dp*pi*dist**3)*0.01_dp
                end if
            end do
        end do
        op%adjoint = transpose(op%k)
        op%genus = 0
    end subroutine

    subroutine tic(t)
        real(dp), intent(out) :: t
        integer :: c, rate
        call system_clock(c, rate)
        t = real(c, dp)/real(rate, dp)
    end subroutine

    function toc(t) result(dt)
        real(dp), intent(in) :: t
        real(dp) :: dt
        integer :: c, rate
        call system_clock(c, rate)
        dt = real(c, dp)/real(rate, dp) - t
    end function

end program
