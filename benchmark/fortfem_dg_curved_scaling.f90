module fortfem_dg_curved_bench_coef
    !! Smooth complex anisotropic coefficients of the DG scaling benchmark.
    use fortfem_dg_curved_2d, only: dg_coefficients_t
    use fortfem_kinds, only: dp
    implicit none
    private

    type, extends(dg_coefficients_t), public :: bench_coef_t
    contains
        procedure :: evaluate => bench_coefficients
    end type bench_coef_t

contains

    subroutine bench_coefficients(self, t, x, kk, b, c, status)
        class(bench_coef_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        complex(dp), intent(out) :: kk(:, :, :)
        real(dp), intent(out) :: b(:, :)
        complex(dp), intent(out) :: c(:)
        integer, intent(out) :: status

        integer :: q

        associate (unused => self)
        end associate
        do q = 1, size(x, 2)
            kk(:, 1, q) = [(1.0_dp, 0.2_dp), (0.3_dp, 0.0_dp)]
            kk(:, 2, q) = [(0.3_dp, 0.0_dp), (1.0e-3_dp, 0.0_dp)]
            b(:, q) = [1.0_dp + x(2, q), 0.5_dp - x(1, q)]
            c(q) = cmplx(0.1_dp, real(t, dp)*1.0e-6_dp, dp)
        end do
        status = 0
    end subroutine bench_coefficients

end module fortfem_dg_curved_bench_coef

program fortfem_dg_curved_scaling
    !! Strong scaling of the curved DG advection-diffusion assembly (P3 on a
    !! warped square, geometry order 3, complex anisotropic tensor).
    !! Usage: fortfem_dg_curved_scaling [cells_per_side] [max_threads]
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_dg_curved_2d, only: DG_FORM_CONSERVATIVE, dg_curved_assemble, &
        dg_curved_space_init, dg_curved_space_t
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_from_map, &
        meridian_mesh_t
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortsparse, only: csc_z_t
    use fortfem_dg_curved_bench_coef, only: bench_coef_t
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    implicit none

    type(mesh_2d_t) :: pmesh, topo
    type(meridian_mesh_t) :: mm
    type(dg_curved_space_t) :: space
    type(bench_coef_t) :: coef
    type(csc_z_t) :: a, a1
    real(dp), allocatable :: pv(:, :, :)
    integer(int64) :: t0, t1, rate
    real(dp) :: seconds, serial, best
    character(len=32) :: text
    integer :: nx, top, threads, st, t, v, rep

    nx = 160
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
    call pmesh%create_rectangular(nx + 1, nx + 1, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    topo = pmesh
    allocate (pv(2, 3, topo%n_triangles))
    do t = 1, topo%n_triangles
        pv(:, :, t) = pmesh%vertices(:, pmesh%triangles(:, t))
    end do
    do v = 1, topo%n_vertices
        call warp(pmesh%vertices(:, v), topo%vertices(:, v))
    end do
    call meridian_mesh_from_map(topo, pv, 3, mm, st, warp)
    if (st /= 0) error stop "mesh failed"
    call dg_curved_space_init(mm, 3, space, st)
    if (st /= 0) error stop "space failed"
    write (*, "(a,i0,a,i0,a,i0)") "P3 DG warped square nx=", nx, " triangles=", &
        space%nt, " dofs=", space%ndof
    write (*, "(a)") "threads  assemble[s]  speedup  max|A-A_1|/max|A_1|"
    threads = 1
    serial = 1.0_dp
    do
        !$ call omp_set_num_threads(threads)
        best = huge(1.0_dp)
        do rep = 1, 3
            call system_clock(t0, rate)
            call dg_curved_assemble(space, coef, DG_FORM_CONSERVATIVE, a, st)
            call system_clock(t1)
            if (st /= 0) error stop "assembly failed"
            best = min(best, real(t1 - t0, dp)/real(rate, dp))
        end do
        seconds = best
        if (threads == 1) then
            serial = seconds
            a1 = a
            write (*, "(a,i0)") "nnz=", a%nnz
        end if
        write (*, "(i7, f13.4, f9.2, es12.3)") threads, seconds, serial/seconds, &
            maxval(abs(a%val - a1%val))/maxval(abs(a1%val))
        if (threads >= top) exit
        threads = min(2*threads, top)
    end do

contains

    subroutine warp(p, x)
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: x(2)

        x(1) = p(1) + 0.05_dp*sin(6.0_dp*p(2))
        x(2) = p(2) + 0.04_dp*sin(5.0_dp*p(1))
    end subroutine warp

end program fortfem_dg_curved_scaling
