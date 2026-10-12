program test_toroidal_mode_boundary_operators_threads_2d
    !! The panel-parity coloured assembly of V_n, K_n, K'_n, W_n and the mass
    !! matrices, and the coil flux table, do not depend on the number of
    !! OpenMP threads to roundoff (vacuous in a build without OpenMP).
    !$ use omp_lib, only: omp_set_num_threads
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_toroidal_mode_kernel, only: tm_rect_coil_t, &
        tm_rect_coil_flux_table
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: gs_bem_curve_t, &
        gs_bem_operators_t, gs_bem_curve_from_nodes
    use fortfem_toroidal_mode_boundary_operators_2d, only: tm_bem_assemble
    implicit none

    integer, parameter :: thread_counts(4) = [1, 2, 8, 48]
    real(dp), parameter :: tolerance = 1.0e-12_dp
    logical :: all_passed

    all_passed = .true.
    call run_case(24, 1, 0)
    call run_case(25, 2, 1)
    call run_case(40, 0, 3)
    call run_case(31, 1, 6)
    call coil_table_case()
    call check_summary("Toroidal-mode boundary operators, thread invariance")
    if (.not. all_passed) error stop 1

contains

    subroutine make_curve(np, curve)
        integer, intent(in) :: np
        type(gs_bem_curve_t), intent(out) :: curve
        real(dp) :: nodes(2, 2*np), theta
        integer :: i, status

        do i = 1, 2*np
            theta = 2.0_dp*acos(-1.0_dp)*real(i - 1, dp)/real(2*np, dp)
            nodes(:, i) = [3.0_dp + cos(theta + 0.4_dp*sin(theta)), &
                0.5_dp + 1.4_dp*sin(theta)]
        end do
        call gs_bem_curve_from_nodes(nodes, curve, status)
        if (status /= 0) error stop "curve construction failed"
    end subroutine make_curve

    function close_to(a, b) result(same)
        real(dp), intent(in) :: a(:, :), b(:, :)
        logical :: same

        same = all(shape(a) == shape(b))
        if (same) same = maxval(abs(a - b)) <= tolerance*max(maxval(abs(b)), &
            tiny(1.0_dp))
    end function close_to

    subroutine run_case(np, deg, n)
        integer, intent(in) :: np, deg, n
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ref, ops
        character(len=64) :: tag
        integer :: it, status

        call make_curve(np, curve)
        write (tag, '(a,i0,a,i0,a,i0)') "n=", n, " np=", np, " flux degree=", deg
        !$ call omp_set_num_threads(thread_counts(1))
        call tm_bem_assemble(curve, n, deg, ref, status)
        do it = 2, size(thread_counts)
            !$ call omp_set_num_threads(thread_counts(it))
            call tm_bem_assemble(curve, n, deg, ops, status)
            call record(close_to(ops%v, ref%v) .and. close_to(ops%k, ref%k) &
                .and. close_to(ops%kt, ref%kt) .and. close_to(ops%w, ref%w) &
                .and. close_to(ops%m, ref%m) .and. close_to(ops%mxx, ref%mxx) &
                .and. close_to(ops%kxx, ref%kxx), &
                "operators independent of thread count, "//trim(tag))
        end do
    end subroutine run_case

    subroutine coil_table_case()
        type(tm_rect_coil_t) :: coils(5)
        real(dp) :: r(20), z(20), ref(20, 5), table(20, 5)
        integer :: i, it

        do i = 1, 5
            coils(i) = tm_rect_coil_t(r=1.5_dp + 0.4_dp*i, z=0.3_dp*i - 0.6_dp, &
                dr=0.1_dp, dz=0.15_dp, turns=2.0_dp)
        end do
        do i = 1, 20
            r(i) = 0.8_dp + 0.05_dp*i
            z(i) = 0.05_dp*i - 0.5_dp
        end do
        !$ call omp_set_num_threads(thread_counts(1))
        call tm_rect_coil_flux_table(3, coils, r, z, 6, 6, ref)
        do it = 2, size(thread_counts)
            !$ call omp_set_num_threads(thread_counts(it))
            call tm_rect_coil_flux_table(3, coils, r, z, 6, 6, table)
            call record(close_to(table, ref), &
                "coil flux table independent of thread count")
        end do
    end subroutine coil_table_case

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_toroidal_mode_boundary_operators_threads_2d
