module fortfem_nested_geometry_differential_jet
    !! Closure-neutral differential jet for a sampled three-dimensional map.
    !!
    !! The caller owns the sampled map value and derivatives.  This module
    !! only applies tensor calculus to those arrays: it does not identify an
    !! axis, select a coordinate convention, or impose an equilibrium model.
    !! An ``axis_limit`` entry says that the caller has supplied a finite
    !! limiting jet at that sample.  At such a sample the inverse Jacobian is
    !! not inferred when the Jacobian is singular; its returned value and
    !! derivatives are zero.  A caller needing a contravariant limit must
    !! supply that limit separately.  This explicit branch keeps a hidden
    !! division by zero out of generated/client code.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
                          FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_nested_geometry_differential_jet
    public :: evaluate_nested_geometry_differential_jet_jvp
    public :: evaluate_nested_geometry_differential_jet_vjp
    public :: evaluate_nested_geometry_polynomial_jet
    public :: validate_nested_geometry_axis_flags

contains

    subroutine evaluate_nested_geometry_differential_jet( &
        value, jacobian, hessian, third_derivative, axis_limit, metric, &
        metric_gradient, metric_hessian, determinant, determinant_gradient, &
        determinant_hessian, inverse_jacobian, inverse_jacobian_gradient, &
        status)
        !! Evaluate metric, volume, and inverse-map diagnostics from a jet.
        real(dp), intent(in) :: value(:, :), jacobian(:, :, :)
        real(dp), intent(in) :: hessian(:, :, :, :), third_derivative(:, :, :, :, :)
        logical, intent(in) :: axis_limit(:)
        real(dp), intent(out) :: metric(:, :, :), metric_gradient(:, :, :, :)
        real(dp), intent(out) :: metric_hessian(:, :, :, :, :)
        real(dp), intent(out) :: determinant(:), determinant_gradient(:, :)
        real(dp), intent(out) :: determinant_hessian(:, :, :)
        real(dp), intent(out) :: inverse_jacobian(:, :, :)
        real(dp), intent(out) :: inverse_jacobian_gradient(:, :, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample, first, second, direction, other
        real(dp) :: inv(3, 3), inv_gradient(3, 3, 3)
        real(dp) :: cofactor(3, 3), h_first(3, 3), h_second(3, 3)

        metric = 0.0_dp
        metric_gradient = 0.0_dp
        metric_hessian = 0.0_dp
        determinant = 0.0_dp
        determinant_gradient = 0.0_dp
        determinant_hessian = 0.0_dp
        inverse_jacobian = 0.0_dp
        inverse_jacobian_gradient = 0.0_dp
        call validate_jet_inputs( &
            value, jacobian, hessian, third_derivative, axis_limit, status)
        if (status%code /= FORTSPARSE_OK) return

        do sample = 1, size(value, 2)
            call evaluate_metric_part( &
                jacobian(:, :, sample), hessian(:, :, :, sample), &
                third_derivative(:, :, :, :, sample), metric(:, :, sample), &
                metric_gradient(:, :, :, sample), metric_hessian(:, :, :, :, sample))
            determinant(sample) = determinant3(jacobian(:, :, sample))
            call determinant3_gradient(jacobian(:, :, sample), cofactor)
            do direction = 1, 3
                h_first = hessian(:, :, direction, sample)
                determinant_gradient(direction, sample) = sum(cofactor*h_first)
                do other = 1, 3
                    h_second = hessian(:, :, other, sample)
                    determinant_hessian(direction, other, sample) = &
                        determinant_second_directional( &
                        jacobian(:, :, sample), h_first, h_second) + &
                        determinant_directional( &
                        jacobian(:, :, sample), &
                        third_derivative(:, :, direction, other, sample))
                end do
            end do
            if (abs(determinant(sample)) <= epsilon(1.0_dp)) then
                if (.not. axis_limit(sample)) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                                 "nested geometry jet has a singular non-axis Jacobian")
                    return
                end if
                inverse_jacobian(:, :, sample) = 0.0_dp
                inverse_jacobian_gradient(:, :, :, sample) = 0.0_dp
            else
                call invert3(jacobian(:, :, sample), inv)
                inverse_jacobian(:, :, sample) = inv
                do direction = 1, 3
                    inv_gradient(:, :, direction) = -matmul( &
                                     matmul(inv, hessian(:, :, direction, sample)), inv)
                    inverse_jacobian_gradient(:, :, direction, sample) = &
                        inv_gradient(:, :, direction)
                end do
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_nested_geometry_differential_jet

    subroutine evaluate_nested_geometry_differential_jet_jvp( &
        value, jacobian, hessian, third_derivative, axis_limit, value_dot, &
        jacobian_dot, hessian_dot, third_derivative_dot, metric_dot, &
        metric_gradient_dot, metric_hessian_dot, determinant_dot, &
        determinant_gradient_dot, determinant_hessian_dot, inverse_dot, &
        inverse_gradient_dot, status)
        !! Apply the exact fixed-topology tangent of the jet diagnostics.
        real(dp), intent(in) :: value(:, :), jacobian(:, :, :)
        real(dp), intent(in) :: hessian(:, :, :, :), third_derivative(:, :, :, :, :)
        logical, intent(in) :: axis_limit(:)
        real(dp), intent(in) :: value_dot(:, :), jacobian_dot(:, :, :)
    real(dp), intent(in) :: hessian_dot(:, :, :, :), third_derivative_dot(:, :, :, :, :)
        real(dp), intent(out) :: metric_dot(:, :, :), metric_gradient_dot(:, :, :, :)
        real(dp), intent(out) :: metric_hessian_dot(:, :, :, :, :)
        real(dp), intent(out) :: determinant_dot(:), determinant_gradient_dot(:, :)
        real(dp), intent(out) :: determinant_hessian_dot(:, :, :)
        real(dp), intent(out) :: inverse_dot(:, :, :)
        real(dp), intent(out) :: inverse_gradient_dot(:, :, :, :)
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: metric(3, 3), metric_gradient(3, 3, 3)
        real(dp) :: metric_hessian(3, 3, 3, 3)
        real(dp) :: determinant, determinant_gradient(3), determinant_hessian(3, 3)
        real(dp) :: inverse(3, 3), inverse_gradient(3, 3, 3)
        integer :: sample, first, second, direction, other
        real(dp) :: h_first(3, 3), h_second(3, 3)
        real(dp) :: h_first_dot(3, 3), h_second_dot(3, 3)

        metric_dot = 0.0_dp
        metric_gradient_dot = 0.0_dp
        metric_hessian_dot = 0.0_dp
        determinant_dot = 0.0_dp
        determinant_gradient_dot = 0.0_dp
        determinant_hessian_dot = 0.0_dp
        inverse_dot = 0.0_dp
        inverse_gradient_dot = 0.0_dp
        call validate_jet_inputs( &
            value, jacobian, hessian, third_derivative, axis_limit, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_jet_shapes( &
            value_dot, jacobian_dot, hessian_dot, third_derivative_dot, &
            size(value, 2))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "nested geometry jet JVP has incompatible tangent shapes")
            return
        end if
        do sample = 1, size(value, 2)
            call local_jet_quantities( &
                jacobian(:, :, sample), hessian(:, :, :, sample), &
                third_derivative(:, :, :, :, sample), metric, metric_gradient, &
                metric_hessian, determinant, determinant_gradient, &
                determinant_hessian, inverse, inverse_gradient, axis_limit(sample), &
                status)
            if (status%code /= FORTSPARSE_OK) return
            call local_jet_tangent( &
                jacobian(:, :, sample), hessian(:, :, :, sample), &
                third_derivative(:, :, :, :, sample), jacobian_dot(:, :, sample), &
               hessian_dot(:, :, :, sample), third_derivative_dot(:, :, :, :, sample), &
                metric_dot(:, :, sample), metric_gradient_dot(:, :, :, sample), &
                metric_hessian_dot(:, :, :, :, sample), determinant_dot(sample), &
           determinant_gradient_dot(:, sample), determinant_hessian_dot(:, :, sample), &
                inverse_dot(:, :, sample), inverse_gradient_dot(:, :, :, sample), &
                inverse, inverse_gradient, axis_limit(sample))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_nested_geometry_differential_jet_jvp

    subroutine evaluate_nested_geometry_differential_jet_vjp( &
        value, jacobian, hessian, third_derivative, axis_limit, metric_bar, &
        metric_gradient_bar, metric_hessian_bar, determinant_bar, &
        determinant_gradient_bar, determinant_hessian_bar, inverse_bar, &
        inverse_gradient_bar, value_bar, jacobian_bar, hessian_bar, &
        third_derivative_bar, status)
        !! Apply the real transpose of the fixed-topology diagnostic map.
        real(dp), intent(in) :: value(:, :), jacobian(:, :, :)
        real(dp), intent(in) :: hessian(:, :, :, :), third_derivative(:, :, :, :, :)
        logical, intent(in) :: axis_limit(:)
        real(dp), intent(in) :: metric_bar(:, :, :), metric_gradient_bar(:, :, :, :)
        real(dp), intent(in) :: metric_hessian_bar(:, :, :, :, :)
        real(dp), intent(in) :: determinant_bar(:), determinant_gradient_bar(:, :)
        real(dp), intent(in) :: determinant_hessian_bar(:, :, :)
        real(dp), intent(in) :: inverse_bar(:, :, :), inverse_gradient_bar(:, :, :, :)
        real(dp), intent(out) :: value_bar(:, :), jacobian_bar(:, :, :)
   real(dp), intent(out) :: hessian_bar(:, :, :, :), third_derivative_bar(:, :, :, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample, first, second, direction, other
        real(dp) :: metric(3, 3), metric_gradient(3, 3, 3)
        real(dp) :: metric_hessian(3, 3, 3, 3), determinant, det_gradient(3)
        real(dp) :: det_hessian(3, 3), inverse(3, 3), inverse_gradient(3, 3, 3)
        real(dp) :: cofactor(3, 3), jbar(3, 3), hbar_local(3, 3, 3)
        real(dp) :: tbar_local(3, 3, 3, 3), qbar
        real(dp) :: first_h(3, 3), second_h(3, 3), first_t(3, 3)

        value_bar = 0.0_dp
        jacobian_bar = 0.0_dp
        hessian_bar = 0.0_dp
        third_derivative_bar = 0.0_dp
        call validate_jet_inputs( &
            value, jacobian, hessian, third_derivative, axis_limit, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_vjp_shapes( &
            metric_bar, metric_gradient_bar, metric_hessian_bar, determinant_bar, &
            determinant_gradient_bar, determinant_hessian_bar, inverse_bar, &
            inverse_gradient_bar, size(value, 2))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "nested geometry jet VJP has incompatible cotangent shapes")
            return
        end if
      if (.not. all_finite(metric_bar) .or. .not. all_finite(metric_gradient_bar) .or. &
      .not. all_finite(metric_hessian_bar) .or. .not. all_finite(determinant_bar) .or. &
            .not. all_finite(determinant_gradient_bar) .or. &
     .not. all_finite(determinant_hessian_bar) .or. .not. all_finite(inverse_bar) .or. &
            .not. all_finite(inverse_gradient_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "nested geometry jet VJP has non-finite cotangents")
            return
        end if

        do sample = 1, size(value, 2)
            call local_jet_quantities( &
                jacobian(:, :, sample), hessian(:, :, :, sample), &
                third_derivative(:, :, :, :, sample), metric, metric_gradient, &
                metric_hessian, determinant, det_gradient, det_hessian, inverse, &
                inverse_gradient, axis_limit(sample), status)
            if (status%code /= FORTSPARSE_OK) return
            jbar = 0.0_dp
            hbar_local = 0.0_dp
            tbar_local = 0.0_dp
            if (.not. axis_limit(sample)) then
                call reverse_inverse( &
                    inverse, inverse_bar(:, :, sample), inverse_gradient, &
                    inverse_gradient_bar(:, :, :, sample), jacobian(:, :, sample), &
                    hessian(:, :, :, sample), jbar, hbar_local)
            end if
            do first = 1, 3
                do second = 1, 3
                    qbar = metric_bar(first, second, sample)
                    jbar(:, first) = jbar(:, first) + &
                                     qbar*jacobian(:, second, sample)
                    jbar(:, second) = jbar(:, second) + &
                                      qbar*jacobian(:, first, sample)
                    do direction = 1, 3
                        qbar = metric_gradient_bar(first, second, direction, sample)
                   hbar_local(:, first, direction) = hbar_local(:, first, direction) + &
                                                        qbar*jacobian(:, second, sample)
                        jbar(:, second) = jbar(:, second) + &
                                          qbar*hessian(:, first, direction, sample)
                        jbar(:, first) = jbar(:, first) + &
                                         qbar*hessian(:, second, direction, sample)
                        hbar_local(:, second, direction) = &
                      hbar_local(:, second, direction) + qbar*jacobian(:, first, sample)
                        do other = 1, 3
                      qbar = metric_hessian_bar(first, second, direction, other, sample)
                            tbar_local(:, first, direction, other) = &
                                tbar_local(:, first, direction, other) + &
                                qbar*jacobian(:, second, sample)
                            jbar(:, second) = jbar(:, second) + &
                               qbar*third_derivative(:, first, direction, other, sample)
                            hbar_local(:, first, direction) = &
                                hbar_local(:, first, direction) + &
                                qbar*hessian(:, second, other, sample)
                            hbar_local(:, second, other) = &
                                hbar_local(:, second, other) + &
                                qbar*hessian(:, first, direction, sample)
                            hbar_local(:, first, other) = &
                                hbar_local(:, first, other) + &
                                qbar*hessian(:, second, direction, sample)
                            hbar_local(:, second, direction) = &
                                hbar_local(:, second, direction) + &
                                qbar*hessian(:, first, other, sample)
                            jbar(:, first) = jbar(:, first) + &
                              qbar*third_derivative(:, second, direction, other, sample)
                            tbar_local(:, second, direction, other) = &
                                tbar_local(:, second, direction, other) + &
                                qbar*jacobian(:, first, sample)
                        end do
                    end do
                end do
            end do
            call determinant3_gradient(jacobian(:, :, sample), cofactor)
            jbar = jbar + determinant_bar(sample)*cofactor
            do direction = 1, 3
                hbar_local(:, :, direction) = hbar_local(:, :, direction) + &
                                    determinant_gradient_bar(direction, sample)*cofactor
                call add_det_second_directional_gradient( &
                    jacobian(:, :, sample), hessian(:, :, direction, sample), &
                    determinant_gradient_bar(direction, sample), jbar)
                do other = 1, 3
                    qbar = determinant_hessian_bar(direction, other, sample)
                    call add_det_second_base_gradient( &
                        hessian(:, :, direction, sample), &
                        hessian(:, :, other, sample), qbar, jbar)
                    call add_det_second_first_gradient( &
                        jacobian(:, :, sample), hessian(:, :, other, sample), qbar, &
                        hbar_local(:, :, direction))
                    call add_det_second_first_gradient( &
                       jacobian(:, :, sample), hessian(:, :, direction, sample), qbar, &
                        hbar_local(:, :, other))
                    call add_det_second_directional_gradient( &
             jacobian(:, :, sample), third_derivative(:, :, direction, other, sample), &
                        qbar, jbar)
                    tbar_local(:, :, direction, other) = &
                        tbar_local(:, :, direction, other) + qbar*cofactor
                end do
            end do
            jacobian_bar(:, :, sample) = jbar
            hessian_bar(:, :, :, sample) = hbar_local
            third_derivative_bar(:, :, :, :, sample) = tbar_local
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_nested_geometry_differential_jet_vjp

    subroutine evaluate_nested_geometry_polynomial_jet( &
        constant, linear, quadratic, cubic, coordinates, axis_limit, value, &
        jacobian, hessian, third_derivative, metric, metric_gradient, &
        metric_hessian, determinant, determinant_gradient, determinant_hessian, &
        inverse_jacobian, inverse_jacobian_gradient, status)
        !! Evaluate a cubic polynomial map and its complete sampled geometry jet.
        real(dp), intent(in) :: constant(3), linear(3, 3)
        real(dp), intent(in) :: quadratic(3, 3, 3), cubic(3, 3, 3, 3)
        real(dp), intent(in) :: coordinates(:, :)
        logical, intent(in) :: axis_limit(:)
        real(dp), intent(out) :: value(:, :), jacobian(:, :, :)
        real(dp), intent(out) :: hessian(:, :, :, :), third_derivative(:, :, :, :, :)
        real(dp), intent(out) :: metric(:, :, :), metric_gradient(:, :, :, :)
        real(dp), intent(out) :: metric_hessian(:, :, :, :, :)
        real(dp), intent(out) :: determinant(:), determinant_gradient(:, :)
        real(dp), intent(out) :: determinant_hessian(:, :, :)
        real(dp), intent(out) :: inverse_jacobian(:, :, :)
        real(dp), intent(out) :: inverse_jacobian_gradient(:, :, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample, component, first, second, third
        real(dp) :: q(3)
        call clear_outputs( &
            value, jacobian, hessian, third_derivative, metric, metric_gradient, &
            metric_hessian, determinant, determinant_gradient, determinant_hessian, &
            inverse_jacobian, inverse_jacobian_gradient)
        if (size(coordinates, 1) /= 3 .or. size(coordinates, 2) < 1 .or. &
            size(axis_limit) /= size(coordinates, 2) .or. &
            .not. all_finite(constant) .or. .not. all_finite(linear) .or. &
            .not. all_finite(quadratic) .or. .not. all_finite(cubic) .or. &
            .not. all_finite(coordinates)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "nested geometry polynomial jet has invalid inputs")
            return
        end if
        do sample = 1, size(coordinates, 2)
            q = coordinates(:, sample)
            do component = 1, 3
                value(component, sample) = constant(component)
                do first = 1, 3
                    value(component, sample) = value(component, sample) + &
                                               linear(component, first)*q(first)
                    do second = 1, 3
                        value(component, sample) = value(component, sample) + 0.5_dp* &
                                  quadratic(component, first, second)*q(first)*q(second)
                        do third = 1, 3
                value(component, sample) = value(component, sample) + (1.0_dp/6.0_dp)* &
                      cubic(component, first, second, third)*q(first)*q(second)*q(third)
                        end do
                    end do
                end do
                do first = 1, 3
                    jacobian(component, first, sample) = linear(component, first)
                    do second = 1, 3
             jacobian(component, first, sample) = jacobian(component, first, sample) + &
                                           quadratic(component, first, second)*q(second)
                        do third = 1, 3
     jacobian(component, first, sample) = jacobian(component, first, sample) + 0.5_dp* &
                               cubic(component, first, second, third)*q(second)*q(third)
                        end do
                    end do
                    do second = 1, 3
         hessian(component, first, second, sample) = quadratic(component, first, second)
                        do third = 1, 3
                            hessian(component, first, second, sample) = &
                                hessian(component, first, second, sample) + &
                                cubic(component, first, second, third)*q(third)
                        end do
                        do third = 1, 3
                           third_derivative(component, first, second, third, sample) = &
                                cubic(component, first, second, third)
                        end do
                    end do
                end do
            end do
        end do
        call evaluate_nested_geometry_differential_jet( &
            value, jacobian, hessian, third_derivative, axis_limit, metric, &
            metric_gradient, metric_hessian, determinant, determinant_gradient, &
            determinant_hessian, inverse_jacobian, inverse_jacobian_gradient, status)
    end subroutine evaluate_nested_geometry_polynomial_jet

    subroutine validate_nested_geometry_axis_flags(axis_limit, sample_count, status)
        !! Validate only fixed-topology axis metadata; no physics is inferred.
        logical, intent(in) :: axis_limit(:)
        integer, intent(in) :: sample_count
        type(fortsparse_status_t), intent(out) :: status
        if (sample_count < 1 .or. size(axis_limit) /= sample_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                         "nested geometry axis-limit flags have the wrong sample count")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_nested_geometry_axis_flags

    subroutine validate_jet_inputs(value, jacobian, hessian, third_derivative, &
                                   axis_limit, status)
        real(dp), intent(in) :: value(:, :), jacobian(:, :, :)
        real(dp), intent(in) :: hessian(:, :, :, :), third_derivative(:, :, :, :, :)
        logical, intent(in) :: axis_limit(:)
        type(fortsparse_status_t), intent(out) :: status
        if (size(value, 1) /= 3 .or. size(value, 2) < 1 .or. &
            size(jacobian, 1) /= 3 .or. size(jacobian, 2) /= 3 .or. &
            size(jacobian, 3) /= size(value, 2) .or. &
            size(hessian, 1) /= 3 .or. size(hessian, 2) /= 3 .or. &
            size(hessian, 3) /= 3 .or. size(hessian, 4) /= size(value, 2) .or. &
            size(third_derivative, 1) /= 3 .or. size(third_derivative, 2) /= 3 .or. &
            size(third_derivative, 3) /= 3 .or. size(third_derivative, 4) /= 3 .or. &
            size(third_derivative, 5) /= size(value, 2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "nested geometry jet arrays have incompatible dimensions")
            return
        end if
        call validate_nested_geometry_axis_flags(axis_limit, size(value, 2), status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. all_finite(value) .or. .not. all_finite(jacobian) .or. &
            .not. all_finite(hessian) .or. .not. all_finite(third_derivative)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "nested geometry jet arrays contain non-finite values")
            return
        end if
    end subroutine validate_jet_inputs

    subroutine local_jet_quantities(jacobian, hessian, third_derivative, metric, &
                   metric_gradient, metric_hessian, determinant, determinant_gradient, &
                     determinant_hessian, inverse, inverse_gradient, axis_limit, status)
  real(dp), intent(in) :: jacobian(3, 3), hessian(3, 3, 3), third_derivative(3, 3, 3, 3)
        real(dp), intent(out) :: metric(3, 3), metric_gradient(3, 3, 3)
        real(dp), intent(out) :: metric_hessian(3, 3, 3, 3), determinant
        real(dp), intent(out) :: determinant_gradient(3), determinant_hessian(3, 3)
        real(dp), intent(out) :: inverse(3, 3), inverse_gradient(3, 3, 3)
        logical, intent(in) :: axis_limit
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: cofactor(3, 3), h_first(3, 3), h_second(3, 3)
        integer :: first, second, direction, other
        call evaluate_metric_part(jacobian, hessian, third_derivative, metric, &
                                  metric_gradient, metric_hessian)
        determinant = determinant3(jacobian)
        call determinant3_gradient(jacobian, cofactor)
        do direction = 1, 3
            h_first = hessian(:, :, direction)
            determinant_gradient(direction) = sum(cofactor*h_first)
            do other = 1, 3
                h_second = hessian(:, :, other)
               determinant_hessian(direction, other) = determinant_second_directional( &
                               jacobian, h_first, h_second) + determinant_directional( &
                                     jacobian, third_derivative(:, :, direction, other))
            end do
        end do
        if (abs(determinant) <= epsilon(1.0_dp)) then
            if (.not. axis_limit) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                                "nested geometry jet has a singular non-axis Jacobian")
                return
            end if
            inverse = 0.0_dp
            inverse_gradient = 0.0_dp
        else
            call invert3(jacobian, inverse)
            do direction = 1, 3
                inverse_gradient(:, :, direction) = -matmul( &
                                     matmul(inverse, hessian(:, :, direction)), inverse)
            end do
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine local_jet_quantities

    subroutine local_jet_tangent(jacobian, hessian, third_derivative, jacobian_dot, &
                   hessian_dot, third_derivative_dot, metric_dot, metric_gradient_dot, &
                        metric_hessian_dot, determinant_dot, determinant_gradient_dot, &
                  determinant_hessian_dot, inverse_dot, inverse_gradient_dot, inverse, &
                                 inverse_gradient, axis_limit)
  real(dp), intent(in) :: jacobian(3, 3), hessian(3, 3, 3), third_derivative(3, 3, 3, 3)
        real(dp), intent(in) :: jacobian_dot(3, 3), hessian_dot(3, 3, 3)
        real(dp), intent(in) :: third_derivative_dot(3, 3, 3, 3)
        real(dp), intent(out) :: metric_dot(3, 3), metric_gradient_dot(3, 3, 3)
        real(dp), intent(out) :: metric_hessian_dot(3, 3, 3, 3)
        real(dp), intent(out) :: determinant_dot, determinant_gradient_dot(3)
        real(dp), intent(out) :: determinant_hessian_dot(3, 3)
        real(dp), intent(out) :: inverse_dot(3, 3), inverse_gradient_dot(3, 3, 3)
        real(dp), intent(in) :: inverse(3, 3), inverse_gradient(3, 3, 3)
        logical, intent(in) :: axis_limit
        integer :: first, second, direction, other
        real(dp) :: cofactor(3, 3), h_first(3, 3), h_second(3, 3)
        real(dp) :: h_first_dot(3, 3), h_second_dot(3, 3)
        real(dp) :: third_slice(3, 3), third_slice_dot(3, 3)
        call evaluate_metric_tangent(jacobian, hessian, third_derivative, &
                          jacobian_dot, hessian_dot, third_derivative_dot, metric_dot, &
                                     metric_gradient_dot, metric_hessian_dot)
        call determinant3_gradient(jacobian, cofactor)
        determinant_dot = sum(cofactor*jacobian_dot)
        do direction = 1, 3
            h_first = hessian(:, :, direction)
            h_first_dot = hessian_dot(:, :, direction)
            determinant_gradient_dot(direction) = &
                determinant_second_directional(jacobian, jacobian_dot, h_first) + &
                sum(cofactor*h_first_dot)
            do other = 1, 3
                h_second = hessian(:, :, other)
                h_second_dot = hessian_dot(:, :, other)
                third_slice = third_derivative(:, :, direction, other)
                third_slice_dot = third_derivative_dot(:, :, direction, other)
                determinant_hessian_dot(direction, other) = &
                    determinant_second_directional(jacobian_dot, h_first, h_second) + &
                    determinant_second_directional(jacobian, h_first_dot, h_second) + &
                    determinant_second_directional(jacobian, h_first, h_second_dot) + &
                 determinant_second_directional(jacobian, jacobian_dot, third_slice) + &
                    determinant_directional(jacobian, third_slice_dot)
            end do
        end do
        if (axis_limit .and. maxval(abs(inverse)) == 0.0_dp) then
            inverse_dot = 0.0_dp
            inverse_gradient_dot = 0.0_dp
        else
            inverse_dot = -matmul(matmul(inverse, jacobian_dot), inverse)
            do direction = 1, 3
                inverse_gradient_dot(:, :, direction) = -matmul( &
                             matmul(inverse_dot, hessian(:, :, direction)), inverse) - &
                      matmul(matmul(inverse, hessian_dot(:, :, direction)), inverse) - &
                          matmul(matmul(inverse, hessian(:, :, direction)), inverse_dot)
            end do
        end if
    end subroutine local_jet_tangent

    subroutine evaluate_metric_part(jacobian, hessian, third_derivative, metric, &
                                    metric_gradient, metric_hessian)
  real(dp), intent(in) :: jacobian(3, 3), hessian(3, 3, 3), third_derivative(3, 3, 3, 3)
        real(dp), intent(out) :: metric(3, 3), metric_gradient(3, 3, 3)
        real(dp), intent(out) :: metric_hessian(3, 3, 3, 3)
        integer :: first, second, direction, other
        metric = 0.0_dp
        metric_gradient = 0.0_dp
        metric_hessian = 0.0_dp
        do first = 1, 3
            do second = 1, 3
                metric(first, second) = dot_product( &
                                        jacobian(:, first), jacobian(:, second))
                do direction = 1, 3
                    metric_gradient(first, second, direction) = &
                      dot_product(hessian(:, first, direction), jacobian(:, second)) + &
                        dot_product(jacobian(:, first), hessian(:, second, direction))
                    do other = 1, 3
                        metric_hessian(first, second, direction, other) = &
      dot_product(third_derivative(:, first, direction, other), jacobian(:, second)) + &
                dot_product(hessian(:, first, direction), hessian(:, second, other)) + &
                dot_product(hessian(:, first, other), hessian(:, second, direction)) + &
          dot_product(jacobian(:, first), third_derivative(:, second, direction, other))
                    end do
                end do
            end do
        end do
    end subroutine evaluate_metric_part

 subroutine evaluate_metric_tangent(jacobian, hessian, third_derivative, jacobian_dot, &
                   hessian_dot, third_derivative_dot, metric_dot, metric_gradient_dot, &
                                       metric_hessian_dot)
  real(dp), intent(in) :: jacobian(3, 3), hessian(3, 3, 3), third_derivative(3, 3, 3, 3)
        real(dp), intent(in) :: jacobian_dot(3, 3), hessian_dot(3, 3, 3)
        real(dp), intent(in) :: third_derivative_dot(3, 3, 3, 3)
        real(dp), intent(out) :: metric_dot(3, 3), metric_gradient_dot(3, 3, 3)
        real(dp), intent(out) :: metric_hessian_dot(3, 3, 3, 3)
        integer :: first, second, direction, other
        metric_dot = 0.0_dp
        metric_gradient_dot = 0.0_dp
        metric_hessian_dot = 0.0_dp
        do first = 1, 3
            do second = 1, 3
metric_dot(first, second) = dot_product(jacobian_dot(:, first), jacobian(:, second)) + &
                                dot_product(jacobian(:, first), jacobian_dot(:, second))
                do direction = 1, 3
                    metric_gradient_dot(first, second, direction) = &
                  dot_product(hessian_dot(:, first, direction), jacobian(:, second)) + &
                  dot_product(hessian(:, first, direction), jacobian_dot(:, second)) + &
                  dot_product(jacobian_dot(:, first), hessian(:, second, direction)) + &
                      dot_product(jacobian(:, first), hessian_dot(:, second, direction))
                    do other = 1, 3
                        metric_hessian_dot(first, second, direction, other) = &
  dot_product(third_derivative_dot(:, first, direction, other), jacobian(:, second)) + &
  dot_product(third_derivative(:, first, direction, other), jacobian_dot(:, second)) + &
            dot_product(hessian_dot(:, first, direction), hessian(:, second, other)) + &
            dot_product(hessian(:, first, direction), hessian_dot(:, second, other)) + &
            dot_product(hessian_dot(:, first, other), hessian(:, second, direction)) + &
            dot_product(hessian(:, first, other), hessian_dot(:, second, direction)) + &
  dot_product(jacobian_dot(:, first), third_derivative(:, second, direction, other)) + &
      dot_product(jacobian(:, first), third_derivative_dot(:, second, direction, other))
                    end do
                end do
            end do
        end do
    end subroutine evaluate_metric_tangent

    subroutine reverse_inverse(inverse, inverse_bar, inverse_gradient, &
                     inverse_gradient_bar, jacobian, hessian, jacobian_bar, hessian_bar)
        real(dp), intent(in) :: inverse(3, 3), inverse_bar(3, 3)
        real(dp), intent(in) :: inverse_gradient(3, 3, 3), inverse_gradient_bar(3, 3, 3)
        real(dp), intent(in) :: jacobian(3, 3), hessian(3, 3, 3)
        real(dp), intent(inout) :: jacobian_bar(3, 3), hessian_bar(3, 3, 3)
        real(dp) :: inverse_bar_total(3, 3), hbar(3, 3), k(3, 3)
        integer :: direction
        inverse_bar_total = inverse_bar
        do direction = 1, 3
            hbar = -matmul(matmul(transpose(inverse), &
                             inverse_gradient_bar(:, :, direction)), transpose(inverse))
            hessian_bar(:, :, direction) = hessian_bar(:, :, direction) + hbar
            k = -matmul(inverse_gradient_bar(:, :, direction), &
                        transpose(matmul(hessian(:, :, direction), inverse))) - &
                matmul(transpose(hessian(:, :, direction)), matmul(transpose(inverse), &
                                                 inverse_gradient_bar(:, :, direction)))
            inverse_bar_total = inverse_bar_total + k
        end do
        jacobian_bar = jacobian_bar - matmul(transpose(inverse), &
                                          matmul(inverse_bar_total, transpose(inverse)))
    end subroutine reverse_inverse

    pure real(dp) function determinant3(matrix) result(value)
        real(dp), intent(in) :: matrix(3, 3)
        value = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) - &
                matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
                matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
    end function determinant3

    pure subroutine determinant3_gradient(matrix, gradient)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp), intent(out) :: gradient(3, 3)
        gradient(1, 1) = matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)
        gradient(1, 2) = matrix(2, 3)*matrix(3, 1) - matrix(2, 1)*matrix(3, 3)
        gradient(1, 3) = matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1)
        gradient(2, 1) = matrix(1, 3)*matrix(3, 2) - matrix(1, 2)*matrix(3, 3)
        gradient(2, 2) = matrix(1, 1)*matrix(3, 3) - matrix(1, 3)*matrix(3, 1)
        gradient(2, 3) = matrix(1, 2)*matrix(3, 1) - matrix(1, 1)*matrix(3, 2)
        gradient(3, 1) = matrix(1, 2)*matrix(2, 3) - matrix(1, 3)*matrix(2, 2)
        gradient(3, 2) = matrix(1, 3)*matrix(2, 1) - matrix(1, 1)*matrix(2, 3)
        gradient(3, 3) = matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1)
    end subroutine determinant3_gradient

    pure real(dp) function determinant_directional(matrix, direction) result(value)
        real(dp), intent(in) :: matrix(3, 3), direction(3, 3)
        real(dp) :: cofactor(3, 3)
        call determinant3_gradient(matrix, cofactor)
        value = sum(cofactor*direction)
    end function determinant_directional

    pure real(dp) function determinant_second_directional(matrix, first, second) result(value)
        real(dp), intent(in) :: matrix(3, 3), first(3, 3), second(3, 3)
        integer, parameter :: permutations(3, 6) = reshape([ &
                          1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
        real(dp), parameter :: signs(6) = [1.0_dp, -1.0_dp, -1.0_dp, &
                                           1.0_dp, 1.0_dp, -1.0_dp]
        integer :: permutation, i, j, k, pi, pj, pk
        value = 0.0_dp
        do permutation = 1, 6
            pi = permutations(1, permutation)
            pj = permutations(2, permutation)
            pk = permutations(3, permutation)
            value = value + signs(permutation)*( &
                    first(1, pi)*second(2, pj)*matrix(3, pk) + &
                    second(1, pi)*first(2, pj)*matrix(3, pk) + &
                    first(1, pi)*matrix(2, pj)*second(3, pk) + &
                    second(1, pi)*matrix(2, pj)*first(3, pk) + &
                    matrix(1, pi)*first(2, pj)*second(3, pk) + &
                    matrix(1, pi)*second(2, pj)*first(3, pk))
        end do
    end function determinant_second_directional

    subroutine add_det_second_directional_gradient(matrix, direction, scale, gradient)
        real(dp), intent(in) :: matrix(3, 3), direction(3, 3), scale
        real(dp), intent(inout) :: gradient(3, 3)
        integer, parameter :: permutations(3, 6) = reshape([ &
                          1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
        real(dp), parameter :: signs(6) = [1.0_dp, -1.0_dp, -1.0_dp, &
                                           1.0_dp, 1.0_dp, -1.0_dp]
        integer :: permutation, pi, pj, pk
        do permutation = 1, 6
            pi = permutations(1, permutation)
            pj = permutations(2, permutation)
            pk = permutations(3, permutation)
            gradient(1, pi) = gradient(1, pi) + scale*signs(permutation)*( &
                        direction(2, pj)*matrix(3, pk) + matrix(2, pj)*direction(3, pk))
            gradient(2, pj) = gradient(2, pj) + scale*signs(permutation)*( &
                        direction(1, pi)*matrix(3, pk) + matrix(1, pi)*direction(3, pk))
            gradient(3, pk) = gradient(3, pk) + scale*signs(permutation)*( &
                        direction(1, pi)*matrix(2, pj) + matrix(1, pi)*direction(2, pj))
        end do
    end subroutine add_det_second_directional_gradient

    subroutine add_det_second_base_gradient(first, second, scale, gradient)
        real(dp), intent(in) :: first(3, 3), second(3, 3), scale
        real(dp), intent(inout) :: gradient(3, 3)
        integer, parameter :: permutations(3, 6) = reshape([ &
                          1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
        real(dp), parameter :: signs(6) = [1.0_dp, -1.0_dp, -1.0_dp, &
                                           1.0_dp, 1.0_dp, -1.0_dp]
        integer :: permutation, pi, pj, pk
        do permutation = 1, 6
            pi = permutations(1, permutation)
            pj = permutations(2, permutation)
            pk = permutations(3, permutation)
            gradient(1, pi) = gradient(1, pi) + scale*signs(permutation)*( &
                              first(2, pj)*second(3, pk) + second(2, pj)*first(3, pk))
            gradient(2, pj) = gradient(2, pj) + scale*signs(permutation)*( &
                              first(1, pi)*second(3, pk) + second(1, pi)*first(3, pk))
            gradient(3, pk) = gradient(3, pk) + scale*signs(permutation)*( &
                              first(1, pi)*second(2, pj) + second(1, pi)*first(2, pj))
        end do
    end subroutine add_det_second_base_gradient

    subroutine add_det_second_first_gradient(matrix, second, scale, gradient)
        real(dp), intent(in) :: matrix(3, 3), second(3, 3), scale
        real(dp), intent(inout) :: gradient(3, 3)
        integer, parameter :: permutations(3, 6) = reshape([ &
                          1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
        real(dp), parameter :: signs(6) = [1.0_dp, -1.0_dp, -1.0_dp, &
                                           1.0_dp, 1.0_dp, -1.0_dp]
        integer :: permutation, pi, pj, pk
        do permutation = 1, 6
            pi = permutations(1, permutation)
            pj = permutations(2, permutation)
            pk = permutations(3, permutation)
            gradient(1, pi) = gradient(1, pi) + scale*signs(permutation)*( &
                              second(2, pj)*matrix(3, pk) + matrix(2, pj)*second(3, pk))
            gradient(2, pj) = gradient(2, pj) + scale*signs(permutation)*( &
                              second(1, pi)*matrix(3, pk) + matrix(1, pi)*second(3, pk))
            gradient(3, pk) = gradient(3, pk) + scale*signs(permutation)*( &
                              second(1, pi)*matrix(2, pj) + matrix(1, pi)*second(2, pj))
        end do
    end subroutine add_det_second_first_gradient

    pure subroutine invert3(matrix, inverse)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp), intent(out) :: inverse(3, 3)
        real(dp) :: cofactor(3, 3), determinant
        call determinant3_gradient(matrix, cofactor)
        determinant = determinant3(matrix)
        inverse = transpose(cofactor)/determinant
    end subroutine invert3

    logical function valid_jet_shapes(value_dot, jacobian_dot, hessian_dot, &
                                      third_derivative_dot, sample_count) result(valid)
        real(dp), intent(in) :: value_dot(:, :), jacobian_dot(:, :, :)
    real(dp), intent(in) :: hessian_dot(:, :, :, :), third_derivative_dot(:, :, :, :, :)
        integer, intent(in) :: sample_count
        valid = size(value_dot, 1) == 3 .and. size(value_dot, 2) == sample_count .and. &
                size(jacobian_dot, 1) == 3 .and. size(jacobian_dot, 2) == 3 .and. &
           size(jacobian_dot, 3) == sample_count .and. size(hessian_dot, 1) == 3 .and. &
                size(hessian_dot, 2) == 3 .and. size(hessian_dot, 3) == 3 .and. &
   size(hessian_dot, 4) == sample_count .and. size(third_derivative_dot, 1) == 3 .and. &
     size(third_derivative_dot, 2) == 3 .and. size(third_derivative_dot, 3) == 3 .and. &
  size(third_derivative_dot, 4) == 3 .and. size(third_derivative_dot, 5) == sample_count
        if (valid) valid = all_finite(value_dot) .and. all_finite(jacobian_dot) .and. &
                          all_finite(hessian_dot) .and. all_finite(third_derivative_dot)
    end function valid_jet_shapes

logical function valid_vjp_shapes(metric_bar, metric_gradient_bar, metric_hessian_bar, &
      determinant_bar, determinant_gradient_bar, determinant_hessian_bar, inverse_bar, &
                                      inverse_gradient_bar, sample_count) result(valid)
        real(dp), intent(in) :: metric_bar(:, :, :), metric_gradient_bar(:, :, :, :)
        real(dp), intent(in) :: metric_hessian_bar(:, :, :, :, :), determinant_bar(:)
real(dp), intent(in) :: determinant_gradient_bar(:, :), determinant_hessian_bar(:, :, :)
        real(dp), intent(in) :: inverse_bar(:, :, :), inverse_gradient_bar(:, :, :, :)
        integer, intent(in) :: sample_count
        valid = size(metric_bar, 1) == 3 .and. size(metric_bar, 2) == 3 .and. &
     size(metric_bar, 3) == sample_count .and. size(metric_gradient_bar, 1) == 3 .and. &
       size(metric_gradient_bar, 2) == 3 .and. size(metric_gradient_bar, 3) == 3 .and. &
            size(metric_gradient_bar, 4) == sample_count .and. size(metric_hessian_bar, 1) == 3 .and. &
         size(metric_hessian_bar, 2) == 3 .and. size(metric_hessian_bar, 3) == 3 .and. &
            size(metric_hessian_bar, 4) == 3 .and. size(metric_hessian_bar, 5) == sample_count .and. &
            size(determinant_bar) == sample_count .and. size(determinant_gradient_bar, 1) == 3 .and. &
            size(determinant_gradient_bar, 2) == sample_count .and. size(determinant_hessian_bar, 1) == 3 .and. &
            size(determinant_hessian_bar, 2) == 3 .and. size(determinant_hessian_bar, 3) == sample_count .and. &
                size(inverse_bar, 1) == 3 .and. size(inverse_bar, 2) == 3 .and. &
   size(inverse_bar, 3) == sample_count .and. size(inverse_gradient_bar, 1) == 3 .and. &
     size(inverse_gradient_bar, 2) == 3 .and. size(inverse_gradient_bar, 3) == 3 .and. &
                size(inverse_gradient_bar, 4) == sample_count
    end function valid_vjp_shapes

    subroutine clear_outputs(value, jacobian, hessian, third_derivative, metric, &
                   metric_gradient, metric_hessian, determinant, determinant_gradient, &
                       determinant_hessian, inverse_jacobian, inverse_jacobian_gradient)
        real(dp), intent(out) :: value(:, :), jacobian(:, :, :), hessian(:, :, :, :)
        real(dp), intent(out) :: third_derivative(:, :, :, :, :), metric(:, :, :)
     real(dp), intent(out) :: metric_gradient(:, :, :, :), metric_hessian(:, :, :, :, :)
        real(dp), intent(out) :: determinant(:), determinant_gradient(:, :)
        real(dp), intent(out) :: determinant_hessian(:, :, :), inverse_jacobian(:, :, :)
        real(dp), intent(out) :: inverse_jacobian_gradient(:, :, :, :)
        value = 0.0_dp
        jacobian = 0.0_dp
        hessian = 0.0_dp
        third_derivative = 0.0_dp
        metric = 0.0_dp
        metric_gradient = 0.0_dp
        metric_hessian = 0.0_dp
        determinant = 0.0_dp
        determinant_gradient = 0.0_dp
        determinant_hessian = 0.0_dp
        inverse_jacobian = 0.0_dp
        inverse_jacobian_gradient = 0.0_dp
    end subroutine clear_outputs

    logical function all_finite(values) result(valid)
        real(dp), intent(in) :: values(..)
        select rank (values)
        rank (1)
            valid = all(ieee_is_finite(values))
        rank (2)
            valid = all(ieee_is_finite(values))
        rank (3)
            valid = all(ieee_is_finite(values))
        rank (4)
            valid = all(ieee_is_finite(values))
        rank (5)
            valid = all(ieee_is_finite(values))
        rank default
            valid = .false.
        end select
    end function all_finite

end module fortfem_nested_geometry_differential_jet
