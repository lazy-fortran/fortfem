program test_axisymmetric_gs_boundary_operators_threads_2d
    !! The panel-parity coloured assembly of V, K, K', W and the mass
    !! matrices is bitwise independent of the number of OpenMP threads (this
    !! is vacuous in a build without OpenMP).
    !$ use omp_lib, only: omp_set_num_threads
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_axisymmetric_gs_kernel, only: gs_rect_coil_t, &
        gs_rect_coil_flux, gs_rect_coil_flux_table
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: &
        gs_bem_curve_t, gs_bem_operators_t, gs_bem_curve_from_nodes, &
        gs_bem_assemble
    implicit none

    integer, parameter :: thread_counts(4) = [1, 2, 8, 48]
    logical :: all_passed

    all_passed = .true.
    call run_case(24, 1)
    call run_case(25, 2)
    call run_case(40, 0)
    call coil_table_case()
    call check_summary("Axisymmetric GS boundary operators, thread invariance")
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

    function identical(a, b) result(same)
        real(dp), intent(in) :: a(:, :), b(:, :)
        logical :: same

        same = all(shape(a) == shape(b))
        if (same) same = all(a == b)
    end function identical

    subroutine run_case(np, deg)
        integer, intent(in) :: np, deg
        type(gs_bem_curve_t) :: curve
        type(gs_bem_operators_t) :: ref, ops
        character(len=64) :: tag
        integer :: it, status

        call make_curve(np, curve)
        write (tag, '(a,i0,a,i0)') "np=", np, " flux degree=", deg
        !$ call omp_set_num_threads(thread_counts(1))
        call gs_bem_assemble(curve, deg, ref, status)
        all_passed = all_passed .and. status == 0
        call check_condition(status == 0, "assembly status, "//trim(tag))
        do it = 2, size(thread_counts)
            !$ call omp_set_num_threads(thread_counts(it))
            call gs_bem_assemble(curve, deg, ops, status)
            call record(identical(ref%v, ops%v) .and. &
                identical(ref%k, ops%k) .and. identical(ref%kt, ops%kt) .and. &
                identical(ref%w, ops%w) .and. identical(ref%m, ops%m) .and. &
                identical(ref%mxx, ops%mxx) .and. &
                identical(ref%kxx, ops%kxx), &
                "bitwise equal operators, "//trim(tag))
        end do
    end subroutine run_case

    subroutine coil_table_case()
        type(gs_rect_coil_t) :: coils(3)
        real(dp) :: r(37), z(37), ref(37, 3), psi(37, 3)
        integer :: i, it

        do i = 1, 37
            r(i) = 2.0_dp + 0.05_dp*i
            z(i) = -1.0_dp + 0.06_dp*i
        end do
        do i = 1, 3
            coils(i) = gs_rect_coil_t(r=6.0_dp + i, z=2.0_dp*i, dr=0.3_dp, &
                dz=0.4_dp, turns=10.0_dp*i)
        end do
        do i = 1, 3
            ref(:, i) = [(gs_rect_coil_flux(coils(i), r(it), z(it), 6, 6), &
                it=1, 37)]
        end do
        do it = 1, size(thread_counts)
            !$ call omp_set_num_threads(thread_counts(it))
            call gs_rect_coil_flux_table(coils, r, z, 6, 6, psi)
            call record(identical(ref, psi), "coil flux table equals scalar")
        end do
    end subroutine coil_table_case

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_axisymmetric_gs_boundary_operators_threads_2d
