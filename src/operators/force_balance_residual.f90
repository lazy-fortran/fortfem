module fortfem_force_balance_residual
    !! Closure-neutral weak force-balance composition.
    !!
    !! For each test function i and component c, the residual is
    !!
    !!   R(i,c) = sum_q wv(q) V(q,i,c) Fv(q,c)
    !!          + sum_b wb(b) B(b,i,c) Fb(b,c)
    !!          + sum_s ws(s) S(s,i,c) Fs(s,c).
    !!
    !! Volume force samples may be supplied by magnetic, tensor-pressure,
    !! inertial, or body-force blocks.  Boundary and sheet terms are explicit
    !! so an Ampere current sheet is never silently smeared into a volume
    !! source.  This module owns only the weak composition and its derivatives;
    !! constitutive laws, geometry, and physical units remain caller-owned.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_force_balance_residual
    public :: assemble_force_balance_residual_jvp
    public :: assemble_force_balance_residual_vjp

contains

    subroutine assemble_force_balance_residual( &
            volume_test, volume_force, volume_weights, boundary_test, &
            boundary_force, boundary_weights, sheet_test, sheet_force, &
            sheet_weights, residual, status)
        real(dp), intent(in) :: volume_test(:, :, :), volume_force(:, :)
        real(dp), intent(in) :: volume_weights(:)
        real(dp), intent(in) :: boundary_test(:, :, :), boundary_force(:, :)
        real(dp), intent(in) :: boundary_weights(:)
        real(dp), intent(in) :: sheet_test(:, :, :), sheet_force(:, :)
        real(dp), intent(in) :: sheet_weights(:)
        real(dp), intent(out) :: residual(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample, test_function, component

        residual = 0.0_dp
        if (.not. validate_value_inputs( &
            volume_test, volume_force, volume_weights, boundary_test, &
            boundary_force, boundary_weights, sheet_test, sheet_force, &
            sheet_weights, residual, status)) return
        do sample = 1, size(volume_weights)
            do test_function = 1, size(residual, 1)
                do component = 1, size(residual, 2)
                    residual(test_function, component) = &
                        residual(test_function, component) + &
                        volume_weights(sample)*volume_test( &
                        sample, test_function, component)*volume_force( &
                        sample, component)
                end do
            end do
        end do
        call add_force_term( &
            boundary_test, boundary_force, boundary_weights, residual)
        call add_force_term(sheet_test, sheet_force, sheet_weights, residual)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_force_balance_residual

    subroutine assemble_force_balance_residual_jvp( &
            volume_test, volume_force, volume_weights, boundary_test, &
            boundary_force, boundary_weights, sheet_test, sheet_force, &
            sheet_weights, volume_test_dot, volume_force_dot, volume_weights_dot, &
            boundary_test_dot, boundary_force_dot, boundary_weights_dot, &
            sheet_test_dot, sheet_force_dot, sheet_weights_dot, residual_dot, status)
        real(dp), intent(in) :: volume_test(:, :, :), volume_force(:, :)
        real(dp), intent(in) :: volume_weights(:)
        real(dp), intent(in) :: boundary_test(:, :, :), boundary_force(:, :)
        real(dp), intent(in) :: boundary_weights(:)
        real(dp), intent(in) :: sheet_test(:, :, :), sheet_force(:, :)
        real(dp), intent(in) :: sheet_weights(:)
        real(dp), intent(in) :: volume_test_dot(:, :, :), volume_force_dot(:, :)
        real(dp), intent(in) :: volume_weights_dot(:)
        real(dp), intent(in) :: boundary_test_dot(:, :, :), boundary_force_dot(:, :)
        real(dp), intent(in) :: boundary_weights_dot(:)
        real(dp), intent(in) :: sheet_test_dot(:, :, :), sheet_force_dot(:, :)
        real(dp), intent(in) :: sheet_weights_dot(:)
        real(dp), intent(out) :: residual_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample, test_function, component

        residual_dot = 0.0_dp
        if (.not. validate_value_inputs( &
            volume_test, volume_force, volume_weights, boundary_test, &
            boundary_force, boundary_weights, sheet_test, sheet_force, &
            sheet_weights, residual_dot, status)) return
        if (.not. validate_direction_inputs( &
            volume_test_dot, volume_force_dot, volume_weights_dot, &
            boundary_test_dot, boundary_force_dot, boundary_weights_dot, &
            sheet_test_dot, sheet_force_dot, sheet_weights_dot, &
            size(residual_dot, 1), size(residual_dot, 2), size(volume_weights), &
            size(boundary_weights), size(sheet_weights), status)) return

        do sample = 1, size(volume_weights)
            do test_function = 1, size(residual_dot, 1)
                do component = 1, size(residual_dot, 2)
                    residual_dot(test_function, component) = &
                        residual_dot(test_function, component) + &
                        volume_weights_dot(sample)*volume_test( &
                        sample, test_function, component)*volume_force(sample, component) + &
                        volume_weights(sample)*volume_test_dot( &
                        sample, test_function, component)*volume_force(sample, component) + &
                        volume_weights(sample)*volume_test( &
                        sample, test_function, component)*volume_force_dot(sample, component)
                end do
            end do
        end do
        call add_force_term_jvp( &
            boundary_test, boundary_force, boundary_weights, boundary_test_dot, &
            boundary_force_dot, boundary_weights_dot, residual_dot)
        call add_force_term_jvp( &
            sheet_test, sheet_force, sheet_weights, sheet_test_dot, &
            sheet_force_dot, sheet_weights_dot, residual_dot)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_force_balance_residual_jvp

    subroutine assemble_force_balance_residual_vjp( &
            volume_test, volume_force, volume_weights, boundary_test, &
            boundary_force, boundary_weights, sheet_test, sheet_force, &
            sheet_weights, residual_bar, volume_test_bar, volume_force_bar, &
            volume_weights_bar, boundary_test_bar, boundary_force_bar, &
            boundary_weights_bar, sheet_test_bar, sheet_force_bar, &
            sheet_weights_bar, status)
        real(dp), intent(in) :: volume_test(:, :, :), volume_force(:, :)
        real(dp), intent(in) :: volume_weights(:)
        real(dp), intent(in) :: boundary_test(:, :, :), boundary_force(:, :)
        real(dp), intent(in) :: boundary_weights(:)
        real(dp), intent(in) :: sheet_test(:, :, :), sheet_force(:, :)
        real(dp), intent(in) :: sheet_weights(:), residual_bar(:, :)
        real(dp), intent(out) :: volume_test_bar(:, :, :), volume_force_bar(:, :)
        real(dp), intent(out) :: volume_weights_bar(:)
        real(dp), intent(out) :: boundary_test_bar(:, :, :), boundary_force_bar(:, :)
        real(dp), intent(out) :: boundary_weights_bar(:)
        real(dp), intent(out) :: sheet_test_bar(:, :, :), sheet_force_bar(:, :)
        real(dp), intent(out) :: sheet_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample, test_function, component
        real(dp) :: cotangent

        volume_test_bar = 0.0_dp
        volume_force_bar = 0.0_dp
        volume_weights_bar = 0.0_dp
        boundary_test_bar = 0.0_dp
        boundary_force_bar = 0.0_dp
        boundary_weights_bar = 0.0_dp
        sheet_test_bar = 0.0_dp
        sheet_force_bar = 0.0_dp
        sheet_weights_bar = 0.0_dp
        if (.not. validate_value_inputs( &
            volume_test, volume_force, volume_weights, boundary_test, &
            boundary_force, boundary_weights, sheet_test, sheet_force, &
            sheet_weights, residual_bar, status)) return
        if (.not. validate_adjoint_inputs( &
            volume_test, boundary_test, sheet_test, &
            volume_test_bar, volume_force_bar, volume_weights_bar, &
            boundary_test_bar, boundary_force_bar, boundary_weights_bar, &
            sheet_test_bar, sheet_force_bar, sheet_weights_bar, residual_bar, &
            status)) return

        do sample = 1, size(volume_weights)
            do test_function = 1, size(residual_bar, 1)
                do component = 1, size(residual_bar, 2)
                    cotangent = residual_bar(test_function, component)
                    volume_test_bar(sample, test_function, component) = &
                        volume_test_bar(sample, test_function, component) + &
                        volume_weights(sample)*cotangent*volume_force(sample, component)
                    volume_force_bar(sample, component) = &
                        volume_force_bar(sample, component) + &
                        volume_weights(sample)*cotangent*volume_test( &
                        sample, test_function, component)
                    volume_weights_bar(sample) = volume_weights_bar(sample) + &
                        cotangent*volume_test(sample, test_function, component)* &
                        volume_force(sample, component)
                end do
            end do
        end do
        call add_force_term_vjp( &
            boundary_test, boundary_force, boundary_weights, residual_bar, &
            boundary_test_bar, boundary_force_bar, boundary_weights_bar)
        call add_force_term_vjp( &
            sheet_test, sheet_force, sheet_weights, residual_bar, sheet_test_bar, &
            sheet_force_bar, sheet_weights_bar)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_force_balance_residual_vjp

    subroutine add_force_term(test, force, weights, residual)
        real(dp), intent(in) :: test(:, :, :), force(:, :), weights(:)
        real(dp), intent(inout) :: residual(:, :)
        integer :: sample, test_function, component

        do sample = 1, size(weights)
            do test_function = 1, size(residual, 1)
                do component = 1, size(residual, 2)
                    residual(test_function, component) = &
                        residual(test_function, component) + weights(sample)* &
                        test(sample, test_function, component)*force(sample, component)
                end do
            end do
        end do
    end subroutine add_force_term

    subroutine add_force_term_jvp( &
            test, force, weights, test_dot, force_dot, weights_dot, residual_dot)
        real(dp), intent(in) :: test(:, :, :), force(:, :), weights(:)
        real(dp), intent(in) :: test_dot(:, :, :), force_dot(:, :), weights_dot(:)
        real(dp), intent(inout) :: residual_dot(:, :)
        integer :: sample, test_function, component

        do sample = 1, size(weights)
            do test_function = 1, size(residual_dot, 1)
                do component = 1, size(residual_dot, 2)
                    residual_dot(test_function, component) = &
                        residual_dot(test_function, component) + &
                        weights_dot(sample)*test(sample, test_function, component)* &
                        force(sample, component) + weights(sample)*test_dot( &
                        sample, test_function, component)*force(sample, component) + &
                        weights(sample)*test(sample, test_function, component)* &
                        force_dot(sample, component)
                end do
            end do
        end do
    end subroutine add_force_term_jvp

    subroutine add_force_term_vjp( &
            test, force, weights, residual_bar, test_bar, force_bar, weights_bar)
        real(dp), intent(in) :: test(:, :, :), force(:, :), weights(:)
        real(dp), intent(in) :: residual_bar(:, :)
        real(dp), intent(out) :: test_bar(:, :, :), force_bar(:, :), weights_bar(:)
        integer :: sample, test_function, component
        real(dp) :: cotangent

        test_bar = 0.0_dp
        force_bar = 0.0_dp
        weights_bar = 0.0_dp
        do sample = 1, size(weights)
            do test_function = 1, size(residual_bar, 1)
                do component = 1, size(residual_bar, 2)
                    cotangent = residual_bar(test_function, component)
                    test_bar(sample, test_function, component) = &
                        test_bar(sample, test_function, component) + &
                        weights(sample)*cotangent*force(sample, component)
                    force_bar(sample, component) = force_bar(sample, component) + &
                        weights(sample)*cotangent*test(sample, test_function, component)
                    weights_bar(sample) = weights_bar(sample) + cotangent* &
                        test(sample, test_function, component)*force(sample, component)
                end do
            end do
        end do
    end subroutine add_force_term_vjp

    logical function validate_value_inputs( &
            volume_test, volume_force, volume_weights, boundary_test, &
            boundary_force, boundary_weights, sheet_test, sheet_force, &
            sheet_weights, residual, status) result(valid)
        real(dp), intent(in) :: volume_test(:, :, :), volume_force(:, :)
        real(dp), intent(in) :: volume_weights(:)
        real(dp), intent(in) :: boundary_test(:, :, :), boundary_force(:, :)
        real(dp), intent(in) :: boundary_weights(:)
        real(dp), intent(in) :: sheet_test(:, :, :), sheet_force(:, :)
        real(dp), intent(in) :: sheet_weights(:), residual(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: test_count, component_count

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "force-balance composition received incompatible arrays")
        if (size(volume_test, 1) < 1 .or. size(volume_test, 2) < 1 .or. &
            size(volume_test, 3) < 1) return
        test_count = size(volume_test, 2)
        component_count = size(volume_test, 3)
        if (size(volume_force, 1) /= size(volume_test, 1) .or. &
            size(volume_force, 2) /= component_count .or. &
            size(volume_weights) /= size(volume_test, 1)) return
        if (.not. valid_optional_term( &
            boundary_test, boundary_force, boundary_weights, test_count, &
            component_count)) return
        if (.not. valid_optional_term( &
            sheet_test, sheet_force, sheet_weights, test_count, component_count)) return
        if (size(residual, 1) /= test_count .or. &
            size(residual, 2) /= component_count) return
        if (.not. finite_term(volume_test, volume_force, volume_weights)) return
        if (any(volume_weights <= 0.0_dp)) return
        if (.not. finite_term(boundary_test, boundary_force, boundary_weights)) return
        if (.not. finite_term(sheet_test, sheet_force, sheet_weights)) return
        if (any(.not. ieee_is_finite(residual))) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_value_inputs

    logical function valid_optional_term(test, force, weights, test_count, &
            component_count) result(valid)
        real(dp), intent(in) :: test(:, :, :), force(:, :), weights(:)
        integer, intent(in) :: test_count, component_count

        valid = size(test, 2) == test_count .and. size(test, 3) == component_count
        if (.not. valid) return
        valid = size(force, 1) == size(test, 1) .and. &
            size(force, 2) == component_count .and. size(weights) == size(test, 1)
        if (.not. valid) return
        valid = all(weights > 0.0_dp)
    end function valid_optional_term

    logical function finite_term(test, force, weights) result(valid)
        real(dp), intent(in) :: test(:, :, :), force(:, :), weights(:)

        valid = all(ieee_is_finite(test)) .and. all(ieee_is_finite(force)) .and. &
            all(ieee_is_finite(weights))
    end function finite_term

    logical function validate_direction_inputs( &
            volume_test_dot, volume_force_dot, volume_weights_dot, &
            boundary_test_dot, boundary_force_dot, boundary_weights_dot, &
            sheet_test_dot, sheet_force_dot, sheet_weights_dot, test_count, &
            component_count, volume_count, boundary_count, sheet_count, status) &
            result(valid)
        real(dp), intent(in) :: volume_test_dot(:, :, :), volume_force_dot(:, :)
        real(dp), intent(in) :: volume_weights_dot(:)
        real(dp), intent(in) :: boundary_test_dot(:, :, :), boundary_force_dot(:, :)
        real(dp), intent(in) :: boundary_weights_dot(:)
        real(dp), intent(in) :: sheet_test_dot(:, :, :), sheet_force_dot(:, :)
        real(dp), intent(in) :: sheet_weights_dot(:)
        integer, intent(in) :: test_count, component_count
        integer, intent(in) :: volume_count, boundary_count, sheet_count
        type(fortsparse_status_t), intent(out) :: status

        valid = size(volume_test_dot, 1) == volume_count .and. &
            size(volume_test_dot, 2) == test_count .and. &
            size(volume_test_dot, 3) == component_count .and. &
            size(volume_force_dot, 1) == volume_count .and. &
            size(volume_force_dot, 2) == component_count .and. &
            size(volume_weights_dot) == volume_count
        if (valid) valid = valid_direction_term( &
            boundary_test_dot, boundary_force_dot, boundary_weights_dot, &
            test_count, component_count, boundary_count)
        if (valid) valid = valid_direction_term( &
            sheet_test_dot, sheet_force_dot, sheet_weights_dot, test_count, &
            component_count, sheet_count)
        if (valid) valid = all(ieee_is_finite(volume_test_dot)) .and. &
            all(ieee_is_finite(volume_force_dot)) .and. &
            all(ieee_is_finite(volume_weights_dot)) .and. &
            all(ieee_is_finite(boundary_test_dot)) .and. &
            all(ieee_is_finite(boundary_force_dot)) .and. &
            all(ieee_is_finite(boundary_weights_dot)) .and. &
            all(ieee_is_finite(sheet_test_dot)) .and. &
            all(ieee_is_finite(sheet_force_dot)) .and. &
            all(ieee_is_finite(sheet_weights_dot))
        if (.not. valid) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "force-balance JVP received incompatible increments")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_direction_inputs

    logical function valid_direction_term( &
            test, force, weights, test_count, component_count, sample_count) &
            result(valid)
        real(dp), intent(in) :: test(:, :, :), force(:, :), weights(:)
        integer, intent(in) :: test_count, component_count, sample_count

        valid = size(test, 1) == sample_count .and. &
            size(test, 2) == test_count .and. size(test, 3) == component_count
        if (.not. valid) return
        valid = size(force, 1) == sample_count .and. &
            size(force, 2) == component_count .and. size(weights) == size(test, 1)
    end function valid_direction_term

    logical function validate_adjoint_inputs( &
            volume_test, boundary_test, sheet_test, volume_test_bar, &
            volume_force_bar, volume_weights_bar, &
            boundary_test_bar, boundary_force_bar, boundary_weights_bar, &
            sheet_test_bar, sheet_force_bar, sheet_weights_bar, residual_bar, &
            status) result(valid)
        real(dp), intent(in) :: volume_test(:, :, :), boundary_test(:, :, :)
        real(dp), intent(in) :: sheet_test(:, :, :)
        real(dp), intent(in) :: volume_test_bar(:, :, :), volume_force_bar(:, :)
        real(dp), intent(in) :: volume_weights_bar(:)
        real(dp), intent(in) :: boundary_test_bar(:, :, :), boundary_force_bar(:, :)
        real(dp), intent(in) :: boundary_weights_bar(:)
        real(dp), intent(in) :: sheet_test_bar(:, :, :), sheet_force_bar(:, :)
        real(dp), intent(in) :: sheet_weights_bar(:), residual_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: test_count, component_count

        valid = size(volume_test_bar, 1) == size(volume_test, 1) .and. &
            size(volume_test_bar, 2) == size(volume_test, 2) .and. &
            size(volume_test_bar, 3) == size(volume_test, 3)
        if (.not. valid) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "force-balance VJP received empty cotangents")
            return
        end if
        test_count = size(residual_bar, 1)
        component_count = size(residual_bar, 2)
        valid = test_count == size(volume_test, 2) .and. &
            component_count == size(volume_test, 3) .and. &
            size(volume_force_bar, 1) == size(volume_test, 1) .and. &
            size(volume_force_bar, 2) == component_count .and. &
            size(volume_weights_bar) == size(volume_test, 1)
        if (valid) valid = valid_adjoint_term( &
            boundary_test_bar, boundary_force_bar, boundary_weights_bar, &
            test_count, component_count, size(boundary_test, 1))
        if (valid) valid = valid_adjoint_term( &
            sheet_test_bar, sheet_force_bar, sheet_weights_bar, test_count, &
            component_count, size(sheet_test, 1))
        if (valid) valid = all(ieee_is_finite(residual_bar))
        if (.not. valid) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "force-balance VJP received incompatible cotangents")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_adjoint_inputs

    logical function valid_adjoint_term( &
            test, force, weights, test_count, component_count, sample_count) &
            result(valid)
        real(dp), intent(in) :: test(:, :, :), force(:, :), weights(:)
        integer, intent(in) :: test_count, component_count, sample_count

        valid = size(test, 1) == sample_count .and. &
            size(test, 2) == test_count .and. size(test, 3) == component_count
        if (.not. valid) return
        valid = size(force, 1) == sample_count .and. &
            size(force, 2) == component_count .and. size(weights) == size(test, 1)
    end function valid_adjoint_term

end module fortfem_force_balance_residual
