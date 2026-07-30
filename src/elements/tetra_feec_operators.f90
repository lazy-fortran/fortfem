module fortfem_tetra_feec_operators
    use fortfem_kinds, only: dp
    use fortfem_tetra_lagrange_arbitrary_order, only: &
        evaluate_tetra_lagrange, initialize_tetra_lagrange, &
        tetra_lagrange_dof_count, tetra_lagrange_t
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, tetra_nedelec_dof_count, &
        tetra_nedelec_first_kind_t
    use fortfem_tetra_nedelec_interpolation, only: &
        interpolate_reference_tetra_nedelec
    use fortfem_tetra_rt_arbitrary_order, only: &
        initialize_tetra_rt, tetra_rt_dof_count, tetra_rt_t
    use fortfem_tetra_rt_interpolation, only: interpolate_reference_tetra_rt
    implicit none

    private

    public :: build_tetra_discrete_gradient
    public :: build_tetra_discrete_curl

contains

    subroutine build_tetra_discrete_curl(order, matrix, status)
        integer, intent(in) :: order
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        type(tetra_nedelec_first_kind_t) :: nedelec_basis
        type(tetra_rt_t) :: rt_basis
        real(dp), allocatable :: dofs(:)
        integer :: nedelec_dof_count, selected_dof

        status = 1
        if (order < 1 .or. order > 4) return
        call initialize_tetra_nedelec_first_kind( &
            order, nedelec_basis, status)
        if (status /= 0) return
        call initialize_tetra_rt(order - 1, rt_basis, status)
        if (status /= 0) return
        nedelec_dof_count = tetra_nedelec_dof_count(nedelec_basis)
        allocate(matrix(tetra_rt_dof_count(rt_basis), nedelec_dof_count))
        do selected_dof = 1, nedelec_dof_count
            call interpolate_reference_tetra_rt( &
                rt_basis, selected_curl, dofs, status)
            if (status /= 0) return
            matrix(:, selected_dof) = dofs
        end do
        status = 0

    contains

        subroutine selected_curl(point, value)
            real(dp), intent(in) :: point(3)
            real(dp), intent(out) :: value(3)

            real(dp) :: curls(3, nedelec_dof_count)
            real(dp) :: values(3, nedelec_dof_count)
            integer :: local_status

            call evaluate_tetra_nedelec_first_kind( &
                nedelec_basis, point, values, curls, local_status)
            if (local_status /= 0) then
                value = 0.0_dp
            else
                value = curls(:, selected_dof)
            end if
        end subroutine selected_curl

    end subroutine build_tetra_discrete_curl

    subroutine build_tetra_discrete_gradient(order, matrix, status)
        integer, intent(in) :: order
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        type(tetra_lagrange_t) :: lagrange_basis
        type(tetra_nedelec_first_kind_t) :: nedelec_basis
        real(dp), allocatable :: dofs(:)
        integer :: lagrange_dof_count, selected_dof

        status = 1
        if (order < 1 .or. order > 4) return
        call initialize_tetra_lagrange(order, lagrange_basis, status)
        if (status /= 0) return
        call initialize_tetra_nedelec_first_kind( &
            order, nedelec_basis, status)
        if (status /= 0) return
        lagrange_dof_count = tetra_lagrange_dof_count(lagrange_basis)
        allocate(matrix(tetra_nedelec_dof_count(nedelec_basis), &
            lagrange_dof_count))
        do selected_dof = 1, lagrange_dof_count
            call interpolate_reference_tetra_nedelec( &
                nedelec_basis, selected_gradient, dofs, status)
            if (status /= 0) return
            matrix(:, selected_dof) = dofs
        end do
        status = 0

    contains

        pure subroutine selected_gradient(point, value)
            real(dp), intent(in) :: point(3)
            real(dp), intent(out) :: value(3)

            real(dp) :: gradients(3, lagrange_dof_count)
            real(dp) :: values(lagrange_dof_count)
            integer :: local_status

            call evaluate_tetra_lagrange( &
                lagrange_basis, point, values, gradients, local_status)
            if (local_status /= 0) then
                value = 0.0_dp
            else
                value = gradients(:, selected_dof)
            end if
        end subroutine selected_gradient

    end subroutine build_tetra_discrete_gradient

end module fortfem_tetra_feec_operators
