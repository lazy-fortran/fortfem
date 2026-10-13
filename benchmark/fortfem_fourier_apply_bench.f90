module bench_tensor_module
    use fortfem_kinds, only: dp
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t
    implicit none
    private
    public :: bench_tensor_t

    type, extends(meridian_fourier_tensor_t) :: bench_tensor_t
    contains
        procedure :: tensor => bench_eval
    end type

contains

    subroutine bench_eval(self, t, x, k, status)
        class(bench_tensor_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: status
        integer :: p
        do p = 1, size(x, 2)
            k(:, :, p) = 0.0_dp
            k(1, 1, p) = 2.0_dp + 0.2_dp*x(1, p) + 0.5_dp*cos(x(3, p))
            k(2, 2, p) = 1.5_dp + 0.1_dp*x(2, p) + 0.4_dp*cos(2*x(3, p))
            k(3, 3, p) = 1.0_dp + 0.1_dp*x(1, p) + 0.3_dp*sin(x(3, p))
            k(1, 2, p) = 0.3_dp*sin(x(3, p)); k(2, 1, p) = k(1, 2, p)
            k(1, 3, p) = 0.1_dp + 0.05_dp*x(1, p); k(3, 1, p) = k(1, 3, p)
            k(2, 3, p) = 0.2_dp*cos(x(3, p)); k(3, 2, p) = k(2, 3, p)
        end do
        status = 0
    end subroutine bench_eval

end module bench_tensor_module

program fortfem_fourier_apply_bench
    !! Memory and time of the matrix-free FFT apply of the mode-coupled V2 mass
    !! against the assembled matrix, N = 4..16, nphi = 4N + 4 (the KIN6D mapped
    !! NN1 choice).  Usage: fortfem_fourier_apply_bench [cells_per_side] [N...]
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_from_map
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build
    use fortfem_meridian_fourier_apply, only: meridian_fourier_op_t, &
        meridian_fourier_op_init, meridian_fourier_mass_apply, &
        meridian_fourier_op_bytes
    use fortfem_meridian_fourier_curl, only: meridian_fourier_mass
    use bench_tensor_module, only: bench_tensor_t
    use fortsparse, only: csc_t, csc_matvec
    !$ use omp_lib, only: omp_get_max_threads
    implicit none

    type(mesh_2d_t) :: pmesh, topo
    type(meridian_mesh_t) :: mm
    type(meridian_feec_t) :: fs
    type(bench_tensor_t) :: coef
    type(csc_t) :: m
    type(meridian_fourier_op_t) :: op
    real(dp), allocatable :: pv(:, :, :), u(:), ya(:), ym(:)
    integer(int64) :: c0, c1, rate
    real(dp) :: t_asm, t_setup, t_app, t_mat
    character(len=32) :: text
    integer :: ncell, modes, nphi, st, t, rep, nreps, ig, nth
    integer :: list(16), nlist

    ncell = 12
    if (command_argument_count() >= 1) then
        call get_command_argument(1, text); read (text, *) ncell
    end if
    nlist = 0
    do ig = 2, command_argument_count()
        call get_command_argument(ig, text); nlist = nlist + 1
        read (text, *) list(nlist)
    end do
    if (nlist == 0) then
        list(1:4) = [4, 8, 12, 16]; nlist = 4
    end if
    nth = 1
    !$ nth = omp_get_max_threads()
    call pmesh%create_rectangular(ncell, ncell, 2.0_dp, 3.0_dp, -0.5_dp, 0.5_dp)
    topo = pmesh
    allocate (pv(2, 3, topo%n_triangles))
    do t = 1, topo%n_triangles
        pv(:, :, t) = pmesh%vertices(:, pmesh%triangles(:, t))
    end do
    call meridian_mesh_from_map(topo, pv, 1, mm, st, ident_map)
    call meridian_feec_build(mm, 2, fs, st)
    write (*, "(a,i0,a,i0,a,i0,a)") "triangles=", topo%n_triangles, " V2 ndof=", &
        fs%sp(2)%ndof, " nloc=", fs%sp(2)%nloc, " (P2, V2 mass)"
    write (*, "(a)") "  N nphi   asm_MB  apply_MB    t_asm[s] t_setup[s]  t_matvec[s]" &
        //" t_apply[s]  relerr threads"
    do ig = 1, nlist
        modes = list(ig); nphi = 4*modes + 4
        call system_clock(c0, rate)
        call meridian_fourier_mass(fs, 2, modes, nphi, coef, m, st)
        call system_clock(c1); t_asm = real(c1 - c0, dp)/rate
        call system_clock(c0)
        call meridian_fourier_op_init(fs, 2, modes, nphi, coef, op, st)
        call system_clock(c1); t_setup = real(c1 - c0, dp)/rate
        allocate (u(m%ncol), ya(m%ncol), ym(m%ncol))
        call random_number(u)
        ym = csc_matvec(m, u)
        call meridian_fourier_mass_apply(op, u, ya)
        nreps = 5
        call system_clock(c0)
        do rep = 1, nreps
            ym = csc_matvec(m, u)
        end do
        call system_clock(c1); t_mat = real(c1 - c0, dp)/rate/nreps
        call system_clock(c0)
        do rep = 1, nreps
            call meridian_fourier_mass_apply(op, u, ya)
        end do
        call system_clock(c1); t_app = real(c1 - c0, dp)/rate/nreps
        write (*, "(2i4,2f10.2,4es12.3,es9.1,i4)") modes, nphi, &
            (12.0_dp*m%nnz + 4.0_dp*(m%ncol + 1))/1.0e6_dp, &
            real(meridian_fourier_op_bytes(op), dp)/1.0e6_dp, t_asm, t_setup, &
            t_mat, t_app, maxval(abs(ya - ym))/maxval(abs(ym)), nth
        deallocate (u, ya, ym)
    end do

contains

    subroutine ident_map(p, x)
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: x(2)
        x = p
    end subroutine ident_map


end program fortfem_fourier_apply_bench
