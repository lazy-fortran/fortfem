module fortfem_triangle_nedelec_arbitrary_order
    !! First-kind Nedelec basis of order k on the reference triangle. The basis
    !! is dual to edge moments against shifted Legendre polynomials and to
    !! interior moments against L2-orthogonal polynomials q_m of P_{k-2}
    !! (evaluate_triangle_nedelec_interior_tests). FortSym solves these moment
    !! systems in exact
    !! rational arithmetic (tools/codegen/app/gen_triangle_feec_tables.f90),
    !! and the basis is stored in degree-k Bernstein form with every
    !! coefficient rounded once. Bernstein polynomials are nonnegative, so the
    !! evaluation is well conditioned for all supported orders.
    use fortfem_kinds, only: dp
    use fortfem_generated_triangle_feec_tables, only: &
        triangle_feec_table_max_order, triangle_nedelec_bernstein_table, &
        triangle_nedelec_interior_test_table
    use fortfem_triangle_bernstein_jets, only: triangle_bernstein_jet
    implicit none

    private

    public :: assignment(=)
    public :: evaluate_triangle_nedelec_first_kind
    public :: evaluate_triangle_nedelec_first_kind_jvp
    public :: evaluate_triangle_nedelec_first_kind_vjp
    public :: evaluate_triangle_nedelec_interior_tests
    public :: initialize_triangle_nedelec_first_kind
    public :: triangle_nedelec_dof_count
    public :: triangle_nedelec_first_kind_t
    public :: triangle_nedelec_max_order

    integer, parameter :: triangle_nedelec_max_order = &
        triangle_feec_table_max_order

    type :: triangle_nedelec_first_kind_t
        private
        integer :: order = 0
        integer :: dof_count = 0
        integer :: bernstein_count = 0
        real(dp), allocatable :: coefficients(:, :, :) !! (2, beta, dof)
        real(dp), allocatable :: interior_tests(:, :) !! (beta, m), degree k-2
    end type triangle_nedelec_first_kind_t

    interface assignment(=)
        module procedure assign_triangle_nedelec_first_kind
    end interface

contains

    subroutine initialize_triangle_nedelec_first_kind(order, basis, status)
        !! status = 1 for order < 1 or order > triangle_nedelec_max_order.
        integer, intent(in) :: order
        type(triangle_nedelec_first_kind_t), intent(out) :: basis
        integer, intent(out) :: status

        status = 1
        if (order < 1 .or. order > triangle_nedelec_max_order) return
        call triangle_nedelec_bernstein_table(order, basis%coefficients, &
            status)
        if (status /= 0) return
        call triangle_nedelec_interior_test_table(order, basis%interior_tests, &
            status)
        if (status /= 0) return
        basis%order = order
        basis%dof_count = order*(order + 2)
        basis%bernstein_count = (order + 1)*(order + 2)/2
        if (size(basis%coefficients, 3) /= basis%dof_count) status = 2
    end subroutine initialize_triangle_nedelec_first_kind

    subroutine evaluate_triangle_nedelec_first_kind( &
            basis, xi, eta, values, curls, status)
        type(triangle_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(:, :), curls(:)
        integer, intent(out) :: status

        real(dp) :: b(basis%bernstein_count)
        real(dp) :: db(2, basis%bernstein_count), d2b(3, basis%bernstein_count)
        integer :: dof, beta

        values = 0.0_dp
        curls = 0.0_dp
        status = 1
        if (.not. valid_evaluation_shapes(basis, values, curls)) return
        if (.not. valid_reference_point(xi, eta)) return
        call triangle_bernstein_jet(basis%order, xi, eta, b, db, d2b)
        do dof = 1, basis%dof_count
            do beta = 1, basis%bernstein_count
                values(:, dof) = values(:, dof) + &
                    basis%coefficients(:, beta, dof)*b(beta)
                curls(dof) = curls(dof) + &
                    basis%coefficients(2, beta, dof)*db(1, beta) - &
                    basis%coefficients(1, beta, dof)*db(2, beta)
            end do
        end do
        status = 0
    end subroutine evaluate_triangle_nedelec_first_kind

    subroutine evaluate_triangle_nedelec_interior_tests(basis, xi, eta, &
            values, status)
        !! Values of the interior test polynomials q_m; the interior dofs are
        !! int_T q_m v_1 and int_T q_m v_2 (component 1 first), m in graded order.
        type(triangle_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: status

        real(dp), allocatable :: b(:), db(:, :), d2b(:, :)
        integer :: nb

        values = 0.0_dp
        status = 1
        if (basis%order < 1) return
        if (size(values) /= size(basis%interior_tests, 2)) return
        if (size(values) == 0) then
            status = 0
            return
        end if
        nb = size(basis%interior_tests, 1)
        allocate (b(nb), db(2, nb), d2b(3, nb))
        call triangle_bernstein_jet(basis%order - 2, xi, eta, b, db, d2b)
        values = matmul(b, basis%interior_tests)
        status = 0
    end subroutine evaluate_triangle_nedelec_interior_tests

    subroutine evaluate_triangle_nedelec_first_kind_jvp( &
            basis, xi, eta, xi_dot, eta_dot, values_dot, curls_dot, status)
        type(triangle_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta, xi_dot, eta_dot
        real(dp), intent(out) :: values_dot(:, :), curls_dot(:)
        integer, intent(out) :: status

        real(dp) :: b(basis%bernstein_count)
        real(dp) :: db(2, basis%bernstein_count), d2b(3, basis%bernstein_count)
        real(dp) :: b_dot, bx_dot, by_dot
        integer :: dof, beta

        values_dot = 0.0_dp
        curls_dot = 0.0_dp
        status = 1
        if (.not. valid_evaluation_shapes(basis, values_dot, curls_dot)) return
        if (.not. valid_reference_point(xi, eta)) return
        call triangle_bernstein_jet(basis%order, xi, eta, b, db, d2b)
        do dof = 1, basis%dof_count
            do beta = 1, basis%bernstein_count
                b_dot = db(1, beta)*xi_dot + db(2, beta)*eta_dot
                bx_dot = d2b(1, beta)*xi_dot + d2b(2, beta)*eta_dot
                by_dot = d2b(2, beta)*xi_dot + d2b(3, beta)*eta_dot
                values_dot(:, dof) = values_dot(:, dof) + &
                    basis%coefficients(:, beta, dof)*b_dot
                curls_dot(dof) = curls_dot(dof) + &
                    basis%coefficients(2, beta, dof)*bx_dot - &
                    basis%coefficients(1, beta, dof)*by_dot
            end do
        end do
        status = 0
    end subroutine evaluate_triangle_nedelec_first_kind_jvp

    subroutine evaluate_triangle_nedelec_first_kind_vjp( &
            basis, xi, eta, values_bar, curls_bar, xi_bar, eta_bar, status)
        type(triangle_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: xi, eta
        real(dp), intent(in) :: values_bar(:, :), curls_bar(:)
        real(dp), intent(out) :: xi_bar, eta_bar
        integer, intent(out) :: status

        real(dp), allocatable :: curls_eta(:), curls_xi(:)
        real(dp), allocatable :: values_eta(:, :), values_xi(:, :)

        xi_bar = 0.0_dp
        eta_bar = 0.0_dp
        status = 1
        if (.not. valid_evaluation_shapes(basis, values_bar, curls_bar)) return
        allocate(values_xi, mold=values_bar)
        allocate(values_eta, mold=values_bar)
        allocate(curls_xi, mold=curls_bar)
        allocate(curls_eta, mold=curls_bar)
        call evaluate_triangle_nedelec_first_kind_jvp( &
            basis, xi, eta, 1.0_dp, 0.0_dp, values_xi, curls_xi, status)
        if (status /= 0) return
        call evaluate_triangle_nedelec_first_kind_jvp( &
            basis, xi, eta, 0.0_dp, 1.0_dp, values_eta, curls_eta, status)
        if (status /= 0) return
        xi_bar = sum(values_bar*values_xi) + dot_product(curls_bar, curls_xi)
        eta_bar = sum(values_bar*values_eta) + &
            dot_product(curls_bar, curls_eta)
    end subroutine evaluate_triangle_nedelec_first_kind_vjp

    pure function triangle_nedelec_dof_count(basis) result(dof_count)
        type(triangle_nedelec_first_kind_t), intent(in) :: basis
        integer :: dof_count

        dof_count = basis%dof_count
    end function triangle_nedelec_dof_count

    pure logical function valid_evaluation_shapes( &
            basis, values, curls) result(valid)
        type(triangle_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: values(:, :), curls(:)

        valid = basis%order >= 1
        if (.not. valid) return
        valid = basis%dof_count >= 1
        if (.not. valid) return
        valid = size(values, 1) == 2
        if (.not. valid) return
        valid = size(values, 2) == basis%dof_count
        if (.not. valid) return
        valid = size(curls) == basis%dof_count
    end function valid_evaluation_shapes

    pure logical function valid_reference_point(xi, eta) result(valid)
        real(dp), intent(in) :: xi, eta

        valid = xi >= -64.0_dp*epsilon(1.0_dp)
        if (.not. valid) return
        valid = eta >= -64.0_dp*epsilon(1.0_dp)
        if (.not. valid) return
        valid = xi + eta <= 1.0_dp + 64.0_dp*epsilon(1.0_dp)
    end function valid_reference_point

    subroutine assign_triangle_nedelec_first_kind(left, right)
        type(triangle_nedelec_first_kind_t), intent(out) :: left
        type(triangle_nedelec_first_kind_t), intent(in) :: right

        left%order = right%order
        left%dof_count = right%dof_count
        left%bernstein_count = right%bernstein_count
        if (allocated(right%coefficients)) left%coefficients = right%coefficients
        if (allocated(right%interior_tests)) &
            left%interior_tests = right%interior_tests
    end subroutine assign_triangle_nedelec_first_kind

end module fortfem_triangle_nedelec_arbitrary_order
