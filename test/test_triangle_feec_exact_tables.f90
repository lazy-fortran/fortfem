program test_triangle_feec_exact_tables
    !! Independent oracles for the FortSym reference tables of the triangle
    !! complex P_k -> N_k -> RT_{k-1} -> P_{k-1}, k = 1..8:
    !! - table-level exactness curl(grad) = 0 and div(rot grad) = 0, relative
    !!   to the absolute sum of the products (the rounding scale);
    !! - the tables agree with the pointwise evaluators (Bernstein jets for
    !!   N_k and RT_{k-1}, the nodal product form for the Lagrange gradients);
    !! - Nedelec interpolation reproduces every polynomial of N_k.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_generated_triangle_feec_tables, only: &
        triangle_feec_operator_tables, triangle_feec_table_max_order
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        evaluate_triangle_lagrange_basis, initialize_triangle_lagrange_basis, &
        triangle_lagrange_basis_t, triangle_lagrange_nodes
    use fortfem_triangle_nedelec_arbitrary_order, only: &
        evaluate_triangle_nedelec_first_kind, &
        initialize_triangle_nedelec_first_kind, &
        triangle_nedelec_first_kind_t
    use fortfem_triangle_rt_arbitrary_order, only: &
        evaluate_triangle_raviart_thomas, initialize_triangle_raviart_thomas, &
        triangle_rt_basis_t
    use fortfem_triangle_vector_interpolation, only: &
        evaluate_triangle_nedelec_interpolant, interpolate_triangle_nedelec
    use fortfem_triangle_feec_operators, only: build_triangle_discrete_gradient
    implicit none

    real(dp), parameter :: eps = epsilon(1.0_dp)
    real(dp), parameter :: reference_vertices(2, 3) = reshape( &
        [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
    real(dp), parameter :: points(2, 6) = reshape([ &
        0.1_dp, 0.2_dp, 0.7_dp, 0.05_dp, 0.3_dp, 0.6_dp, &
        1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp, 0.45_dp, 0.9_dp, 0.1_dp], [2, 6])

    integer :: field_component, field_x, field_y
    logical :: field_rotational
    logical :: all_passed
    integer :: k, status
    real(dp), allocatable :: g(:, :), c(:, :), d(:, :)
    type(triangle_nedelec_first_kind_t) :: unused

    all_passed = .true.
    do k = 1, triangle_feec_table_max_order
        call triangle_feec_operator_tables(k, g, c, d, status)
        call record(status == 0 .and. size(g, 1) == k*(k + 2) .and. &
            size(g, 2) == (k + 1)*(k + 2)/2 .and. size(c, 1) == k*(k + 1)/2 &
            .and. size(c, 2) == k*(k + 2) .and. all(shape(d) == shape(c)), &
            "reference tables have the FEEC dimensions", k)
        call check_table_exactness(k, g, c, d)
        call check_tables_against_evaluators(k, g, c, d)
        call check_polynomial_reproduction(k)
    end do

    call triangle_feec_operator_tables(0, g, c, d, status)
    call record(status /= 0, "order 0 is rejected", 0)
    call triangle_feec_operator_tables(triangle_feec_table_max_order + 1, &
        g, c, d, status)
    call record(status /= 0, "orders above the table range are rejected", &
        triangle_feec_table_max_order + 1)
    call initialize_triangle_nedelec_first_kind( &
        triangle_feec_table_max_order + 1, unused, status)
    call record(status /= 0, "Nedelec orders above the table range are rejected", &
        triangle_feec_table_max_order + 1)
    call build_triangle_discrete_gradient(4, g, status)
    call record(status == 0 .and. size(g, 1) == 24, &
        "discrete gradient is the generated table", 4)

    call check_summary("Exact triangle FEEC reference tables")
    if (.not. all_passed) error stop 1

contains

    subroutine check_table_exactness(k, g, c, d)
        integer, intent(in) :: k
        real(dp), intent(in) :: g(:, :), c(:, :), d(:, :)

        real(dp), allocatable :: pg(:, :), residual(:, :), scale(:, :)
        integer :: j, pj, sj, half, e3

        residual = matmul(c, g)
        scale = matmul(abs(c), abs(g))
        call report(k, "curl(grad)", maxval(abs(residual))/maxval(scale))
        call record(maxval(abs(residual)) <= 4.0_dp*eps*maxval(scale), &
            "curl(grad) = 0 at the table level to rounding", k)

        ! rot90 N_j = sj RT_pj, so div(rot grad) uses the signed permutation.
        e3 = 3*k
        half = k*(k - 1)/2
        allocate (pg(size(g, 1), size(g, 2)))
        do j = 1, size(g, 1)
            if (j <= e3) then
                pj = j
                sj = 1
            else if (j <= e3 + half) then
                pj = j + half
                sj = -1
            else
                pj = j - half
                sj = 1
            end if
            pg(pj, :) = real(sj, dp)*g(j, :)
        end do
        residual = matmul(d, pg)
        scale = matmul(abs(d), abs(pg))
        call report(k, "div(rot grad)", maxval(abs(residual))/maxval(scale))
        call record(maxval(abs(residual)) <= 4.0_dp*eps*maxval(scale), &
            "div(rot grad) = 0 at the table level to rounding", k)
        call record(maxval(abs(matmul(g, [(1.0_dp, j=1, size(g, 2))]))) <= &
            4.0_dp*eps*maxval(abs(g)), "gradient annihilates constants", k)
    end subroutine check_table_exactness

    subroutine check_tables_against_evaluators(k, g, c, d)
        integer, intent(in) :: k
        real(dp), intent(in) :: g(:, :), c(:, :), d(:, :)

        type(triangle_nedelec_first_kind_t) :: nb
        type(triangle_rt_basis_t) :: rb
        type(triangle_lagrange_basis_t) :: lb, db
        real(dp), allocatable :: nodes(:, :), values(:, :), curls(:)
        real(dp), allocatable :: lv(:), lg(:, :), rebuilt(:, :)
        real(dp) :: error_c, error_d, error_g, scale_g
        integer :: l, p, status, ndof

        ndof = k*(k + 2)
        call initialize_triangle_nedelec_first_kind(k, nb, status)
        call record(status == 0, "Nedelec basis initialises", k)
        call initialize_triangle_raviart_thomas(k - 1, rb, status)
        call record(status == 0, "RT basis initialises", k)
        call initialize_triangle_lagrange_basis(k - 1, db, status)
        call triangle_lagrange_nodes(db, nodes)
        allocate (values(2, ndof), curls(ndof))
        error_c = 0.0_dp
        error_d = 0.0_dp
        do l = 1, size(nodes, 2)
            call evaluate_triangle_nedelec_first_kind(nb, nodes(1, l), &
                nodes(2, l), values, curls, status)
            error_c = max(error_c, maxval(abs(curls - c(l, :))))
            call evaluate_triangle_raviart_thomas(rb, nodes(1, l), &
                nodes(2, l), values, curls, status)
            error_d = max(error_d, maxval(abs(curls - d(l, :))))
        end do
        call report(k, "curl table vs Bernstein evaluator", &
            error_c/maxval(abs(c)))
        call report(k, "div table vs Bernstein evaluator", &
            error_d/maxval(abs(d)))
        call record(error_c <= 64.0_dp*eps*maxval(abs(c)), &
            "curl table matches the Nedelec evaluator", k)
        call record(error_d <= 64.0_dp*eps*maxval(abs(d)), &
            "divergence table matches the RT evaluator", k)

        call initialize_triangle_lagrange_basis(k, lb, status)
        allocate (lv((k + 1)*(k + 2)/2), lg(2, (k + 1)*(k + 2)/2))
        error_g = 0.0_dp
        scale_g = 0.0_dp
        do p = 1, size(points, 2)
            call evaluate_triangle_lagrange_basis(lb, points(1, p), &
                points(2, p), lv, lg, status)
            call evaluate_triangle_nedelec_first_kind(nb, points(1, p), &
                points(2, p), values, curls, status)
            rebuilt = matmul(values, g)
            error_g = max(error_g, maxval(abs(rebuilt - lg)))
            scale_g = max(scale_g, maxval(matmul(abs(values), abs(g))))
        end do
        call report(k, "sum gref N - grad phi", error_g/scale_g)
        call record(error_g <= 64.0_dp*eps*scale_g, &
            "gradient table reproduces the Lagrange gradients", k)
    end subroutine check_tables_against_evaluators

    subroutine check_polynomial_reproduction(k)
        !! Interpolate every monomial field of P_{k-1}^2 and the rotational
        !! fields (-y, x) x^a y^(k-1-a); N_k contains all of them. The dofs come
        !! from floating quadrature, each with an error of order eps max|u|
        !! (max|u| <= 1 here), so the reproduction error is measured against
        !! eps sum_j |N_j(x)| and eps sum_j |curl N_j(x)|: the rounding scale of
        !! the interpolant, which grows with k because the interior dofs are
        !! monomial moments.
        integer, intent(in) :: k

        type(triangle_nedelec_first_kind_t) :: nb
        real(dp), allocatable :: dofs(:), basis_values(:, :), basis_curls(:)
        real(dp) :: value(2), exact(2), curl, exact_curl, worst, worst_absolute
        real(dp) :: value_scale, curl_scale
        integer :: total, a, p, status

        call initialize_triangle_nedelec_first_kind(k, nb, status)
        allocate (basis_values(2, k*(k + 2)), basis_curls(k*(k + 2)))
        worst = 0.0_dp
        worst_absolute = 0.0_dp
        do field_component = 1, 3
            do total = 0, k - 1
                if (field_component == 3 .and. total /= k - 1) cycle
                do a = 0, total
                    field_x = a
                    field_y = total - a
                    field_rotational = field_component == 3
                    call interpolate_triangle_nedelec(reference_vertices, k, &
                        2*k + 2, polynomial_field, dofs, status)
                    if (status /= 0) then
                        call record(.false., "interpolation succeeds", k)
                        return
                    end if
                    do p = 1, size(points, 2)
                        call evaluate_triangle_nedelec_interpolant( &
                            reference_vertices, nb, dofs, points(1, p), &
                            points(2, p), value, curl, status)
                        call evaluate_triangle_nedelec_first_kind(nb, &
                            points(1, p), points(2, p), basis_values, &
                            basis_curls, status)
                        value_scale = eps*sum(abs(basis_values))
                        curl_scale = eps*sum(abs(basis_curls))
                        call polynomial_field(points(1, p), points(2, p), exact)
                        exact_curl = polynomial_curl(points(1, p), points(2, p))
                        worst = max(worst, &
                            maxval(abs(value - exact))/value_scale, &
                            abs(curl - exact_curl)/curl_scale)
                        worst_absolute = max(worst_absolute, &
                            maxval(abs(value - exact)))
                    end do
                end do
            end do
        end do
        call report(k, "interpolation of N_k polynomials, max abs", &
            worst_absolute)
        call report(k, "interpolation error / (eps sum|N_j|)", worst)
        call record(worst <= 16.0_dp, &
            "Nedelec interpolation reproduces N_k polynomials to rounding", k)
    end subroutine check_polynomial_reproduction

    subroutine polynomial_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        real(dp) :: monomial

        monomial = x**field_x*y**field_y
        value = 0.0_dp
        if (field_rotational) then
            value = [-y*monomial, x*monomial]
        else
            value(field_component) = monomial
        end if
    end subroutine polynomial_field

    real(dp) function polynomial_curl(x, y) result(curl)
        real(dp), intent(in) :: x, y

        real(dp) :: dmx, dmy

        dmx = 0.0_dp
        dmy = 0.0_dp
        if (field_x > 0) dmx = real(field_x, dp)*x**(field_x - 1)*y**field_y
        if (field_y > 0) dmy = real(field_y, dp)*x**field_x*y**(field_y - 1)
        if (field_rotational) then
            curl = real(field_x + field_y + 2, dp)*x**field_x*y**field_y
        else if (field_component == 1) then
            curl = -dmy
        else
            curl = dmx
        end if
    end function polynomial_curl

    subroutine report(k, label, value)
        integer, intent(in) :: k
        character(*), intent(in) :: label
        real(dp), intent(in) :: value

        write (*, "(a,i0,a,a,a,es10.3)") "  k=", k, " ", label, ": ", value
    end subroutine report

    subroutine record(condition, label, k)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        integer, intent(in) :: k

        character(16) :: suffix

        write (suffix, "(a,i0,a)") " (k=", k, ")"
        call check_condition(condition, label//trim(suffix))
        if (.not. condition) all_passed = .false.
    end subroutine record

end program test_triangle_feec_exact_tables
