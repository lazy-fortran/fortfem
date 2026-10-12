program fortfem_toroidal_mode_bem_scaling
    !! Strong scaling of the panel-parity coloured assembly of the per-n
    !! boundary operators on a D-shaped meridional curve.
    !! Usage: fortfem_toroidal_mode_bem_scaling [panels] [n] [max_threads]
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_kinds, only: dp
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_curve_from_nodes
    use fortfem_toroidal_mode_boundary_operators_2d, only: tm_bem_assemble
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    implicit none

    type(gs_bem_curve_t) :: curve
    type(gs_bem_operators_t) :: ops, ref
    real(dp), allocatable :: nodes(:, :)
    real(dp) :: theta, seconds, serial, diff
    integer(int64) :: t0, t1, rate
    character(len=32) :: text
    integer :: np, n, top, threads, i, status

    np = 96
    n = 3
    top = 1
    !$ top = omp_get_max_threads()
    if (command_argument_count() >= 1) then
        call get_command_argument(1, text)
        read (text, *) np
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, text)
        read (text, *) n
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, text)
        read (text, *) top
    end if
    allocate(nodes(2, 2*np))
    do i = 1, 2*np
        theta = 2.0_dp*acos(-1.0_dp)*real(i - 1, dp)/real(2*np, dp)
        nodes(:, i) = [3.0_dp + cos(theta + 0.4_dp*sin(theta)), &
            0.5_dp + 1.4_dp*sin(theta)]
    end do
    call gs_bem_curve_from_nodes(nodes, curve, status)
    write (*, "(a,i0,a,i0)") "toroidal-mode BEM assembly, panels=", np, " n=", n
    write (*, "(a)") "threads  assemble[s]  speedup  max|M - M_1thread|/max|M|"
    threads = 1
    do
        !$ call omp_set_num_threads(threads)
        call system_clock(t0, rate)
        call tm_bem_assemble(curve, n, 1, ops, status)
        call system_clock(t1)
        seconds = real(t1 - t0, dp)/real(rate, dp)
        if (threads == 1) then
            serial = seconds
            ref = ops
        end if
        diff = max(maxval(abs(ops%w - ref%w))/maxval(abs(ref%w)), &
            maxval(abs(ops%v - ref%v))/maxval(abs(ref%v)))
        write (*, "(i7,f13.3,f9.2,es14.2)") threads, seconds, serial/seconds, diff
        if (threads >= top) exit
        threads = min(2*threads, top)
    end do
end program fortfem_toroidal_mode_bem_scaling
