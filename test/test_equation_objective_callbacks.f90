program test_equation_objective_callbacks
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        equation_objective_block_t, equation_objective_registry_t, &
        initialize_equation_objective_registry, &
        evaluate_equation_objective_callbacks, &
        evaluate_equation_objective_callbacks_jvp, &
        evaluate_equation_objective_callbacks_vjp, &
        REGISTRY_KIND_EQUATION, REGISTRY_KIND_OBJECTIVE, &
        REGISTRY_KIND_CONSTRAINT
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    type(equation_objective_registry_t) :: registry
    type(fortsparse_status_t) :: status
    real(dp), parameter :: state(3) = [1.0_dp, -2.0_dp, 0.5_dp]
    real(dp), parameter :: state_dot(3) = [-0.3_dp, 0.7_dp, 0.2_dp]
    real(dp), parameter :: eps_fd = 1.0e-7_dp
    real(dp), allocatable :: packed(:), packed_plus(:), packed_minus(:)
    real(dp), allocatable :: packed_dot(:), state_bar(:)
    real(dp) :: packed_bar(4), lhs, rhs
    logical :: all_passed

    all_passed = .true.
    call initialize_equation_objective_registry( &
        registry, &
        [character(len=16) :: "force", "quasisymmetry", "flux"], &
        [REGISTRY_KIND_EQUATION, REGISTRY_KIND_OBJECTIVE, &
        REGISTRY_KIND_CONSTRAINT], [2, 1, 1], &
        [character(len=8) :: "N", "1", "Wb"], &
        [character(len=8) :: "L2", "mean", "target"], &
        [character(len=12) :: "manufactured", "external", "external"], &
        [.true., .true., .true.], [.false., .false., .true.], status)
    call record(status%code == 0, "callback registry metadata is accepted")

    call evaluate_equation_objective_callbacks( &
        registry, state, packed, value_callback, status)
    call record(status%code == 0 .and. size(packed) == 4 .and. &
        maxval(abs(packed - [-3.0_dp, -2.5_dp, 1.5_dp, -0.5_dp])) < 1.0e-14_dp, &
        "named callbacks dispatch into deterministic packed rows")

    call evaluate_equation_objective_callbacks_jvp( &
        registry, state, state_dot, packed_dot, jvp_callback, status)
    call evaluate_equation_objective_callbacks( &
        registry, state + eps_fd*state_dot, packed_plus, value_callback, status)
    call evaluate_equation_objective_callbacks( &
        registry, state - eps_fd*state_dot, packed_minus, value_callback, status)
    call record(status%code == 0 .and. maxval(abs(packed_dot - &
        (packed_plus - packed_minus)/(2.0_dp*eps_fd))) < 1.0e-6_dp, &
        "callback registry JVP matches an independent finite difference")

    packed_bar = [1.0_dp, -2.0_dp, 0.4_dp, 3.0_dp]
    call evaluate_equation_objective_callbacks_vjp( &
        registry, state, packed_bar, state_bar, vjp_callback, status)
    lhs = dot_product(packed_bar, packed_dot)
    rhs = dot_product(state_bar, state_dot)
    call record(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "callback registry VJP satisfies the independent transpose oracle")

    call evaluate_equation_objective_callbacks( &
        registry, state, packed, failing_callback, status)
    call record(status%code /= 0 .and. size(packed) == 0, &
        "callback failures are surfaced without partial packed output")

    call check_summary("Equation/objective callback registry")
    if (.not. all_passed) error stop 1

contains

    subroutine value_callback(block, x, values, callback_status)
        type(equation_objective_block_t), intent(in) :: block
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: callback_status

        values = 0.0_dp
        callback_status = 1
        select case (trim(block%name))
        case ("force")
            if (size(values) /= 2) return
            values = [x(1) + 2.0_dp*x(2), x(2) - x(3)]
        case ("quasisymmetry")
            if (size(values) /= 1) return
            values(1) = x(1)**2 + x(3)
        case ("flux")
            if (size(values) /= 1) return
            values(1) = sum(x)
        case default
            return
        end select
        callback_status = 0
    end subroutine value_callback

    subroutine jvp_callback(block, x, x_dot, values_dot, callback_status)
        type(equation_objective_block_t), intent(in) :: block
        real(dp), intent(in) :: x(:), x_dot(:)
        real(dp), intent(out) :: values_dot(:)
        integer, intent(out) :: callback_status

        values_dot = 0.0_dp
        callback_status = 1
        select case (trim(block%name))
        case ("force")
            if (size(values_dot) /= 2) return
            values_dot = [x_dot(1) + 2.0_dp*x_dot(2), x_dot(2) - x_dot(3)]
        case ("quasisymmetry")
            if (size(values_dot) /= 1) return
            values_dot(1) = 2.0_dp*x(1)*x_dot(1) + x_dot(3)
        case ("flux")
            if (size(values_dot) /= 1) return
            values_dot(1) = sum(x_dot)
        case default
            return
        end select
        callback_status = 0
    end subroutine jvp_callback

    subroutine vjp_callback(block, x, values_bar, x_bar, callback_status)
        type(equation_objective_block_t), intent(in) :: block
        real(dp), intent(in) :: x(:), values_bar(:)
        real(dp), intent(out) :: x_bar(:)
        integer, intent(out) :: callback_status

        x_bar = 0.0_dp
        callback_status = 1
        if (size(x_bar) /= 3) return
        select case (trim(block%name))
        case ("force")
            if (size(values_bar) /= 2) return
            x_bar = [values_bar(1), 2.0_dp*values_bar(1) + values_bar(2), &
                -values_bar(2)]
        case ("quasisymmetry")
            if (size(values_bar) /= 1) return
            x_bar = [2.0_dp*x(1)*values_bar(1), 0.0_dp, values_bar(1)]
        case ("flux")
            if (size(values_bar) /= 1) return
            x_bar = values_bar(1)
        case default
            return
        end select
        callback_status = 0
    end subroutine vjp_callback

    subroutine failing_callback(block, x, values, callback_status)
        type(equation_objective_block_t), intent(in) :: block
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: callback_status

        associate (unused_block => block, unused_x => x)
            values = 0.0_dp
            callback_status = 11
        end associate
    end subroutine failing_callback

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_equation_objective_callbacks
