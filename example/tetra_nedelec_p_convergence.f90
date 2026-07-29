program tetra_nedelec_p_convergence
    use fortfem_api, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, &
        interpolate_reference_tetra_nedelec, tetra_duffy_quadrature, &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    implicit none

    type(tetra_nedelec_first_kind_t) :: basis
    real(dp) :: errors(4)
    integer :: order, status

    print "(a)", "order  L2 error for grad(x^4)"
    do order = 1, 4
        call initialize_tetra_nedelec_first_kind(order, basis, status)
        if (status /= 0) error stop "basis initialization failed"
        call interpolation_error(basis, errors(order), status)
        if (status /= 0) error stop "interpolation failed"
        print "(i5,2x,es14.6)", order, errors(order)
    end do
    if (any(errors(2:4) >= errors(1:3))) then
        error stop "interpolation error did not decrease"
    end if
    if (errors(4) >= 2.0e-11_dp) then
        error stop "order four did not reproduce the cubic field"
    end if

contains

    subroutine interpolation_error(basis, error, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(out) :: error
        integer, intent(out) :: status

        real(dp), allocatable :: curls(:, :), dofs(:), values(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: difference(3), point(3)
        integer :: node

        call interpolate_reference_tetra_nedelec( &
            basis, cubic_gradient, dofs, status)
        if (status /= 0) return
        allocate( &
            values(3, tetra_nedelec_dof_count(basis)), &
            curls(3, tetra_nedelec_dof_count(basis)))
        call tetra_duffy_quadrature(12, x, y, z, weights, status)
        if (status /= 0) return
        error = 0.0_dp
        do node = 1, size(weights)
            point = [x(node), y(node), z(node)]
            call evaluate_tetra_nedelec_first_kind( &
                basis, point, values, curls, status)
            if (status /= 0) return
            difference = matmul(values, dofs) - cubic_gradient_value(point)
            error = error + weights(node) * dot_product(difference, difference)
        end do
        error = sqrt(error)
    end subroutine interpolation_error

    pure subroutine cubic_gradient(point, value)
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: value(3)

        value = cubic_gradient_value(point)
    end subroutine cubic_gradient

    pure function cubic_gradient_value(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)

        value = [4.0_dp * point(1)**3, 0.0_dp, 0.0_dp]
    end function cubic_gradient_value

end program tetra_nedelec_p_convergence
