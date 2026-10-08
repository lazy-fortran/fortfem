program test_p2_hessian_reproduction
    use fortfem_kinds, only: dp
    use basis_p2_2d_module, only: basis_p2_2d_t
    use check, only: check_condition, check_summary
    implicit none
    type(basis_p2_2d_t) :: basis
    real(dp) :: coefficients(6), actual(2, 2), exact(2, 2)
    real(dp) :: samples(2, 3), x, y
    integer :: polynomial, i, point

    samples(:, 1) = [0.2_dp, 0.3_dp]
    samples(:, 2) = [0.0_dp, 0.0_dp]
    samples(:, 3) = [0.5_dp, 0.5_dp]
    do polynomial = 1, 6
        exact = 0.0_dp
        select case (polynomial)
        case (4)
            exact(1, 1) = 2.0_dp
        case (5)
            exact(2, 2) = 2.0_dp
        case (6)
            exact(1, 2) = 1.0_dp
            exact(2, 1) = 1.0_dp
        end select
        do i = 1, 6
            x = basis%nodes(1, i)
            y = basis%nodes(2, i)
            select case (polynomial)
            case (1)
                coefficients(i) = 1.0_dp
            case (2)
                coefficients(i) = x
            case (3)
                coefficients(i) = y
            case (4)
                coefficients(i) = x*x
            case (5)
                coefficients(i) = y*y
            case (6)
                coefficients(i) = x*y
            end select
        end do
        do point = 1, size(samples, 2)
            actual = 0.0_dp
            do i = 1, 6
                actual = actual + coefficients(i)*basis%hessian( &
                    i, samples(1, point), samples(2, point))
            end do
            call check_condition(maxval(abs(actual - exact)) < 1.0e-14_dp, &
                'P2 Hessian reproduces constant, affine and all quadratics')
        end do
    end do
    call check_summary('P2 Hessian polynomial reproduction')
end program test_p2_hessian_reproduction
