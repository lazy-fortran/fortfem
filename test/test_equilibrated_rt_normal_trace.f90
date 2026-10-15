program test_equilibrated_rt_normal_trace
    !! The edge degrees of freedom of the reference RT_k basis, on which the
    !! equilibrated-flux estimator builds its Neumann data and edge orientation
    !! rule: basis function (e - 1)(k + 1) + j + 1 has the shifted Legendre
    !! moment j of the normal trace sigma . (t_y, -t_x) on reference edge e
    !! (counter-clockwise from local vertex e to e + 1, parameter s in [0, 1])
    !! equal to one and every other edge moment zero; interior functions have no
    !! normal trace. The Piola map leaves this flux per unit edge parameter
    !! unchanged on a curved isoparametric cell.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortfem_triangle_rt_arbitrary_order, only: triangle_rt_basis_t, &
        initialize_triangle_raviart_thomas, evaluate_triangle_raviart_thomas
    use equilibrated_test_mesh, only: test_mesh_t, build_test_mesh, cell_point
    implicit none
    integer :: k

    do k = 0, 4
        call reference_moments(k)
    end do
    call curved_cell_trace(3)
    call check_summary("equilibrated_rt_normal_trace")

contains

    pure real(dp) function legendre(j, s) result(p)
        !! Shifted Legendre polynomial of degree j on [0, 1], P_j(1) = 1.
        integer, intent(in) :: j
        real(dp), intent(in) :: s
        real(dp) :: x, previous, current
        integer :: l

        x = 2*s - 1
        previous = 1.0_dp
        current = x
        if (j == 0) then
            p = 1.0_dp
            return
        end if
        do l = 1, j - 1
            p = ((2*l + 1)*x*current - l*previous)/(l + 1)
            previous = current
            current = p
        end do
        p = current
    end function legendre

    pure subroutine edge_point(e, s, point, normal)
        integer, intent(in) :: e
        real(dp), intent(in) :: s
        real(dp), intent(out) :: point(2), normal(2)

        select case (e)
        case (1)
            point = [s, 0.0_dp]
            normal = [0.0_dp, -1.0_dp]
        case (2)
            point = [1.0_dp - s, s]
            normal = [1.0_dp, 1.0_dp]
        case default
            point = [0.0_dp, 1.0_dp - s]
            normal = [-1.0_dp, 0.0_dp]
        end select
    end subroutine edge_point

    subroutine reference_moments(k)
        integer, intent(in) :: k
        type(triangle_rt_basis_t) :: basis
        real(dp), allocatable :: nodes(:), weights(:), values(:, :), div(:)
        real(dp) :: point(2), normal(2), moment, expected, worst
        integer :: stat, e, j, i, q, nrt, edge_dof

        call initialize_triangle_raviart_thomas(k, basis, stat)
        if (stat /= 0) error stop "rt"
        nrt = (k + 1)*(k + 3)
        allocate (nodes(k + 4), weights(k + 4), values(2, nrt), div(nrt))
        call gauss_legendre_ab(k + 4, 0.0_dp, 1.0_dp, nodes, weights)
        worst = 0.0_dp
        do i = 1, nrt
            do e = 1, 3
                do j = 0, k
                    moment = 0.0_dp
                    do q = 1, size(nodes)
                        call edge_point(e, nodes(q), point, normal)
                        call evaluate_triangle_raviart_thomas(basis, point(1), &
                            point(2), values, div, stat)
                        moment = moment + weights(q)*legendre(j, nodes(q)) &
                            *dot_product(values(:, i), normal)
                    end do
                    edge_dof = (e - 1)*(k + 1) + j + 1
                    expected = merge(1.0_dp, 0.0_dp, i == edge_dof)
                    worst = max(worst, abs(moment - expected))
                end do
            end do
        end do
        call check_condition(worst < 1.0e-8_dp, &
            "RT edge degrees of freedom are the Legendre normal-trace moments")
        print '(a, i2, a, es10.2)', "RT degree", k, " worst moment deviation", worst
    end subroutine reference_moments

    subroutine curved_cell_trace(k)
        !! Physical flux per unit edge parameter, sigma . (t_y, -t_x) with
        !! t = J t^ and sigma = J sigma^/det J, equals sigma^ . (t^_y, -t^_x).
        integer, intent(in) :: k
        type(test_mesh_t) :: mesh
        type(triangle_rt_basis_t) :: basis
        real(dp), allocatable :: values(:, :), div(:), phi(:), gref(:, :)
        real(dp) :: point(2), normal(2), tref(2), x(2), jac(2, 2), det, sigma(2)
        real(dp) :: tphys(2), reference, physical, worst
        integer :: stat, t, e, q, i, nrt

        call build_test_mesh(2, 2, 0.08_dp, mesh, stat)
        if (stat /= 0) error stop "mesh"
        call initialize_triangle_raviart_thomas(k, basis, stat)
        nrt = (k + 1)*(k + 3)
        allocate (values(2, nrt), div(nrt), phi(mesh%nb), gref(2, mesh%nb))
        worst = 0.0_dp
        do t = 1, mesh%nt
            do e = 1, 3
                do q = 1, 3
                    call edge_point(e, 0.2_dp*q, point, normal)
                    select case (e)
                    case (1)
                        tref = [1.0_dp, 0.0_dp]
                    case (2)
                        tref = [-1.0_dp, 1.0_dp]
                    case default
                        tref = [0.0_dp, -1.0_dp]
                    end select
                    call cell_point(mesh, t, point(1), point(2), x, jac, det, phi, &
                        gref)
                    call evaluate_triangle_raviart_thomas(basis, point(1), point(2), &
                        values, div, stat)
                    tphys = matmul(jac, tref)
                    do i = 1, nrt
                        sigma = matmul(jac, values(:, i))/det
                        physical = sigma(1)*tphys(2) - sigma(2)*tphys(1)
                        reference = values(1, i)*tref(2) - values(2, i)*tref(1)
                        worst = max(worst, abs(physical - reference))
                    end do
                end do
            end do
        end do
        call check_condition(det > 0.0_dp .and. worst < 1.0e-12_dp, &
            "Piola map preserves the flux per unit edge parameter on curved cells")
        print '(a, es10.2)', "curved cell Piola trace deviation", worst
    end subroutine curved_cell_trace

end program test_equilibrated_rt_normal_trace
