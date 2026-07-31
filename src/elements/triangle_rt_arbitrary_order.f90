module fortfem_triangle_rt_arbitrary_order
    use fortfem_kinds, only: dp
    use fortfem_triangle_nedelec_arbitrary_order, only: &
        assignment(=), evaluate_triangle_nedelec_first_kind, &
        evaluate_triangle_nedelec_first_kind_jvp, &
        initialize_triangle_nedelec_first_kind, triangle_nedelec_first_kind_t
    implicit none

    private

    public :: assignment(=)
    public :: evaluate_triangle_raviart_thomas
    public :: evaluate_triangle_raviart_thomas_jvp
    public :: evaluate_triangle_raviart_thomas_vjp
    public :: initialize_triangle_raviart_thomas
    public :: triangle_rt_basis_t
    public :: triangle_rt_dof_count

    type :: triangle_rt_basis_t
        private
        integer :: degree = -1
        integer :: dof_count = 0
        type(triangle_nedelec_first_kind_t) :: rotated_basis
    end type triangle_rt_basis_t

    interface assignment(=)
        module procedure assign_triangle_rt_basis
    end interface

contains

    subroutine initialize_triangle_raviart_thomas(degree, basis, status)
        integer, intent(in) :: degree
        type(triangle_rt_basis_t), intent(out) :: basis
        integer, intent(out) :: status

        status = 1
        if (degree < 0) return
        basis%degree = degree
        basis%dof_count = (degree + 1) * (degree + 3)
        call initialize_triangle_nedelec_first_kind( &
            degree + 1, basis%rotated_basis, status)
    end subroutine initialize_triangle_raviart_thomas

    subroutine evaluate_triangle_raviart_thomas( &
            basis, xi, eta, values, divergences, status)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(:, :), divergences(:)
        integer, intent(out) :: status

        real(dp), allocatable :: curls(:), nedelec_values(:, :)
        integer :: edge_dof_count, interior_dof_count
        integer :: local_dof, source_dof, target_dof

        values = 0.0_dp
        divergences = 0.0_dp
        status = 1
        if (basis%degree < 0 .or. basis%dof_count < 1) return
        if (size(values, 1) /= 2 .or. &
            size(values, 2) /= basis%dof_count) return
        if (size(divergences) /= basis%dof_count) return

        allocate(curls(basis%dof_count))
        allocate(nedelec_values(2, basis%dof_count))
        call evaluate_triangle_nedelec_first_kind( &
            basis%rotated_basis, xi, eta, nedelec_values, curls, status)
        if (status /= 0) return

        edge_dof_count = 3 * (basis%degree + 1)
        interior_dof_count = basis%degree * (basis%degree + 1) / 2
        do local_dof = 1, edge_dof_count
            values(1, local_dof) = nedelec_values(2, local_dof)
            values(2, local_dof) = -nedelec_values(1, local_dof)
            divergences(local_dof) = curls(local_dof)
        end do
        do local_dof = 1, interior_dof_count
            target_dof = edge_dof_count + local_dof
            source_dof = edge_dof_count + interior_dof_count + local_dof
            values(1, target_dof) = nedelec_values(2, source_dof)
            values(2, target_dof) = -nedelec_values(1, source_dof)
            divergences(target_dof) = curls(source_dof)

            target_dof = edge_dof_count + interior_dof_count + local_dof
            source_dof = edge_dof_count + local_dof
            values(1, target_dof) = -nedelec_values(2, source_dof)
            values(2, target_dof) = nedelec_values(1, source_dof)
            divergences(target_dof) = -curls(source_dof)
        end do
        status = 0
    end subroutine evaluate_triangle_raviart_thomas

    subroutine evaluate_triangle_raviart_thomas_jvp( &
            basis, xi, eta, xi_dot, eta_dot, values_dot, divergences_dot, &
            status)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta, xi_dot, eta_dot
        real(dp), intent(out) :: values_dot(:, :), divergences_dot(:)
        integer, intent(out) :: status

        real(dp), allocatable :: curls_dot(:), nedelec_values_dot(:, :)
        integer :: edge_dof_count, interior_dof_count
        integer :: local_dof, source_dof, target_dof

        values_dot = 0.0_dp
        divergences_dot = 0.0_dp
        status = 1
        if (.not. valid_evaluation_shapes( &
            basis, values_dot, divergences_dot)) return
        allocate(curls_dot(basis%dof_count))
        allocate(nedelec_values_dot(2, basis%dof_count))
        call evaluate_triangle_nedelec_first_kind_jvp( &
            basis%rotated_basis, xi, eta, xi_dot, eta_dot, &
            nedelec_values_dot, curls_dot, status)
        if (status /= 0) return
        edge_dof_count = 3*(basis%degree + 1)
        interior_dof_count = basis%degree*(basis%degree + 1)/2
        do local_dof = 1, edge_dof_count
            values_dot(1, local_dof) = nedelec_values_dot(2, local_dof)
            values_dot(2, local_dof) = -nedelec_values_dot(1, local_dof)
            divergences_dot(local_dof) = curls_dot(local_dof)
        end do
        do local_dof = 1, interior_dof_count
            target_dof = edge_dof_count + local_dof
            source_dof = edge_dof_count + interior_dof_count + local_dof
            values_dot(1, target_dof) = nedelec_values_dot(2, source_dof)
            values_dot(2, target_dof) = -nedelec_values_dot(1, source_dof)
            divergences_dot(target_dof) = curls_dot(source_dof)
            target_dof = edge_dof_count + interior_dof_count + local_dof
            source_dof = edge_dof_count + local_dof
            values_dot(1, target_dof) = -nedelec_values_dot(2, source_dof)
            values_dot(2, target_dof) = nedelec_values_dot(1, source_dof)
            divergences_dot(target_dof) = -curls_dot(source_dof)
        end do
        status = 0
    end subroutine evaluate_triangle_raviart_thomas_jvp

    subroutine evaluate_triangle_raviart_thomas_vjp( &
            basis, xi, eta, values_bar, divergences_bar, xi_bar, eta_bar, &
            status)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(in) :: values_bar(:, :), divergences_bar(:)
        real(dp), intent(out) :: xi_bar, eta_bar
        integer, intent(out) :: status

        real(dp), allocatable :: divergences_eta(:), divergences_xi(:)
        real(dp), allocatable :: values_eta(:, :), values_xi(:, :)

        xi_bar = 0.0_dp
        eta_bar = 0.0_dp
        status = 1
        if (.not. valid_evaluation_shapes( &
            basis, values_bar, divergences_bar)) return
        allocate(values_xi, mold=values_bar)
        allocate(values_eta, mold=values_bar)
        allocate(divergences_xi, mold=divergences_bar)
        allocate(divergences_eta, mold=divergences_bar)
        call evaluate_triangle_raviart_thomas_jvp( &
            basis, xi, eta, 1.0_dp, 0.0_dp, values_xi, divergences_xi, status)
        if (status /= 0) return
        call evaluate_triangle_raviart_thomas_jvp( &
            basis, xi, eta, 0.0_dp, 1.0_dp, values_eta, divergences_eta, status)
        if (status /= 0) return
        xi_bar = sum(values_bar*values_xi) + &
            dot_product(divergences_bar, divergences_xi)
        eta_bar = sum(values_bar*values_eta) + &
            dot_product(divergences_bar, divergences_eta)
    end subroutine evaluate_triangle_raviart_thomas_vjp

    pure logical function valid_evaluation_shapes( &
            basis, values, divergences) result(valid)
        type(triangle_rt_basis_t), intent(in) :: basis
        real(dp), intent(in) :: values(:, :), divergences(:)

        valid = basis%degree >= 0
        if (.not. valid) return
        valid = basis%dof_count >= 1
        if (.not. valid) return
        valid = size(values, 1) == 2
        if (.not. valid) return
        valid = size(values, 2) == basis%dof_count
        if (.not. valid) return
        valid = size(divergences) == basis%dof_count
    end function valid_evaluation_shapes

    pure function triangle_rt_dof_count(basis) result(dof_count)
        type(triangle_rt_basis_t), intent(in) :: basis
        integer :: dof_count

        dof_count = basis%dof_count
    end function triangle_rt_dof_count

    subroutine assign_triangle_rt_basis(left, right)
        type(triangle_rt_basis_t), intent(out) :: left
        type(triangle_rt_basis_t), intent(in) :: right

        left%degree = right%degree
        left%dof_count = right%dof_count
        left%rotated_basis = right%rotated_basis
    end subroutine assign_triangle_rt_basis

end module fortfem_triangle_rt_arbitrary_order
