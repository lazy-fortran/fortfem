program test_linear_perturbation_composition
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        LINEAR_PASSIVITY_HERMITIAN_NONNEGATIVE, &
        LINEAR_RECIPROCITY_TRANSPOSE, &
        assemble_linear_perturbation_operator, &
        assemble_linear_perturbation_operator_jvp, &
        assemble_linear_perturbation_operator_vjp, &
        initialize_linear_perturbation_metadata, &
        linear_perturbation_metadata_t, &
        validate_linear_perturbation_metadata
    implicit none

    integer, parameter :: n = 2
    complex(dp) :: lorentz(n, n), pressure_stress(n, n), inertia(n, n)
    complex(dp) :: vacuum(n, n), wall(n, n), resistive(n, n), singular(n, n)
    complex(dp) :: lorentz_dot(n, n), pressure_stress_dot(n, n)
    complex(dp) :: inertia_dot(n, n), vacuum_dot(n, n), wall_dot(n, n)
    complex(dp) :: resistive_dot(n, n), singular_dot(n, n)
    complex(dp) :: operator(n, n), expected(n, n), operator_dot(n, n)
    complex(dp) :: operator_plus(n, n), operator_minus(n, n)
    complex(dp) :: operator_bar(n, n), lorentz_bar(n, n)
    complex(dp) :: pressure_stress_bar(n, n), inertia_bar(n, n)
    complex(dp) :: vacuum_bar(n, n), wall_bar(n, n), resistive_bar(n, n)
    complex(dp) :: singular_bar(n, n), frequency, frequency_dot, frequency_bar
    type(linear_perturbation_metadata_t) :: metadata, invalid_metadata
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp) :: lhs, rhs
    integer :: row, column, status
    logical :: all_passed

    all_passed = .true.
    call build_inputs( &
        lorentz, pressure_stress, inertia, vacuum, wall, resistive, singular, &
        frequency)
    call build_tangents( &
        lorentz_dot, pressure_stress_dot, inertia_dot, vacuum_dot, wall_dot, &
        resistive_dot, singular_dot, frequency_dot)

    call initialize_linear_perturbation_metadata( &
        metadata, "SI magnetic-energy norm", 2.5_dp, &
        "manufactured MHD-08 composition", LINEAR_RECIPROCITY_TRANSPOSE, &
        LINEAR_PASSIVITY_HERMITIAN_NONNEGATIVE, status)
    call record_condition(status == 0, "perturbation metadata initializes")
    call record_condition( &
        validate_linear_perturbation_metadata(metadata, status), &
        "perturbation metadata validates")
    call record_condition( &
        metadata%normalization_scale == 2.5_dp .and. &
        trim(metadata%normalization_label) == "SI magnetic-energy norm", &
        "normalization metadata is retained")

    call assemble_linear_perturbation_operator( &
        lorentz, pressure_stress, inertia, vacuum, wall, resistive, singular, &
        frequency, operator, status)
    call record_condition(status == 0, "seven perturbation blocks compose")
    do column = 1, n
        do row = 1, n
            expected(row, column) = &
                lorentz(row, column) + pressure_stress(row, column) + &
                vacuum(row, column) + wall(row, column) + &
                singular(row, column) - frequency*frequency*inertia(row, column) + &
                cmplx(0.0_dp, 1.0_dp, dp)*frequency*resistive(row, column)
        end do
    end do
    call record_condition( &
        maxval(abs(operator - expected)) < 1.0e-14_dp, &
        "composition matches an independent elementwise oracle")

    call assemble_linear_perturbation_operator_jvp( &
        lorentz, pressure_stress, inertia, vacuum, wall, resistive, singular, &
        frequency, lorentz_dot, pressure_stress_dot, inertia_dot, vacuum_dot, &
        wall_dot, resistive_dot, singular_dot, frequency_dot, operator_dot, &
        status)
    call record_condition(status == 0, "perturbation composition JVP assembles")
    call assemble_linear_perturbation_operator( &
        lorentz + epsilon*lorentz_dot, &
        pressure_stress + epsilon*pressure_stress_dot, &
        inertia + epsilon*inertia_dot, vacuum + epsilon*vacuum_dot, &
        wall + epsilon*wall_dot, resistive + epsilon*resistive_dot, &
        singular + epsilon*singular_dot, frequency + epsilon*frequency_dot, &
        operator_plus, status)
    call assemble_linear_perturbation_operator( &
        lorentz - epsilon*lorentz_dot, &
        pressure_stress - epsilon*pressure_stress_dot, &
        inertia - epsilon*inertia_dot, vacuum - epsilon*vacuum_dot, &
        wall - epsilon*wall_dot, resistive - epsilon*resistive_dot, &
        singular - epsilon*singular_dot, frequency - epsilon*frequency_dot, &
        operator_minus, status)
    call record_condition( &
        maxval(abs(operator_dot - &
            (operator_plus - operator_minus)/(2.0_dp*epsilon))) < 2.0e-8_dp, &
        "composition JVP matches a central difference")

    operator_bar = reshape([ &
        cmplx(0.31_dp, -0.17_dp, dp), cmplx(-0.22_dp, 0.09_dp, dp), &
        cmplx(0.14_dp, 0.25_dp, dp), cmplx(-0.06_dp, -0.19_dp, dp)], &
        shape(operator_bar))
    call assemble_linear_perturbation_operator_vjp( &
        lorentz, pressure_stress, inertia, vacuum, wall, resistive, singular, &
        frequency, operator_bar, lorentz_bar, pressure_stress_bar, inertia_bar, &
        vacuum_bar, wall_bar, resistive_bar, singular_bar, frequency_bar, status)
    call record_condition(status == 0, "perturbation composition VJP assembles")
    lhs = real(sum(conjg(operator_bar)*operator_dot), dp)
    rhs = real( &
        sum(conjg(lorentz_bar)*lorentz_dot) + &
        sum(conjg(pressure_stress_bar)*pressure_stress_dot) + &
        sum(conjg(inertia_bar)*inertia_dot) + &
        sum(conjg(vacuum_bar)*vacuum_dot) + &
        sum(conjg(wall_bar)*wall_dot) + &
        sum(conjg(resistive_bar)*resistive_dot) + &
        sum(conjg(singular_bar)*singular_dot) + &
        conjg(frequency_bar)*frequency_dot, dp)
    call record_condition( &
        abs(lhs - rhs) < 2.0e-13_dp, &
        "VJP satisfies the real-complex adjoint identity")

    call assemble_linear_perturbation_operator( &
        lorentz(1:1, :), pressure_stress, inertia, vacuum, wall, resistive, &
        singular, frequency, operator, status)
    call record_condition( &
        status /= 0 .and. maxval(abs(operator)) == 0.0_dp, &
        "incompatible blocks are rejected without partial output")
    invalid_metadata = metadata
    invalid_metadata%normalization_scale = 0.0_dp
    call record_condition( &
        .not. validate_linear_perturbation_metadata(invalid_metadata, status), &
        "nonpositive normalization is rejected")

    call check_summary("linear perturbation composition")
    if (.not. all_passed) error stop 1

contains

    subroutine build_inputs( &
            lorentz, pressure_stress, inertia, vacuum, wall, resistive, singular, &
            frequency)
        complex(dp), intent(out) :: lorentz(:, :), pressure_stress(:, :)
        complex(dp), intent(out) :: inertia(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(out) :: resistive(:, :), singular(:, :), frequency

        lorentz = reshape([ &
            cmplx(1.4_dp, 0.2_dp, dp), cmplx(-0.1_dp, 0.3_dp, dp), &
            cmplx(0.2_dp, -0.4_dp, dp), cmplx(1.1_dp, -0.1_dp, dp)], &
            shape(lorentz))
        pressure_stress = reshape([ &
            cmplx(0.4_dp, -0.1_dp, dp), cmplx(0.05_dp, 0.02_dp, dp), &
            cmplx(-0.08_dp, 0.03_dp, dp), cmplx(0.3_dp, 0.07_dp, dp)], &
            shape(pressure_stress))
        inertia = reshape([ &
            cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.03_dp, -0.01_dp, dp), &
            cmplx(0.02_dp, 0.04_dp, dp), cmplx(0.8_dp, 0.0_dp, dp)], &
            shape(inertia))
        vacuum = reshape([ &
            cmplx(0.2_dp, 0.0_dp, dp), cmplx(0.01_dp, 0.05_dp, dp), &
            cmplx(0.01_dp, -0.05_dp, dp), cmplx(0.15_dp, 0.0_dp, dp)], &
            shape(vacuum))
        wall = reshape([ &
            cmplx(0.06_dp, 0.02_dp, dp), cmplx(0.0_dp, 0.01_dp, dp), &
            cmplx(-0.02_dp, 0.0_dp, dp), cmplx(0.04_dp, 0.03_dp, dp)], &
            shape(wall))
        resistive = reshape([ &
            cmplx(0.3_dp, 0.0_dp, dp), cmplx(-0.02_dp, 0.01_dp, dp), &
            cmplx(0.04_dp, -0.03_dp, dp), cmplx(0.2_dp, 0.0_dp, dp)], &
            shape(resistive))
        singular = reshape([ &
            cmplx(0.07_dp, -0.04_dp, dp), cmplx(0.02_dp, 0.01_dp, dp), &
            cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.09_dp, 0.01_dp, dp)], &
            shape(singular))
        frequency = cmplx(0.8_dp, -0.15_dp, dp)
    end subroutine build_inputs

    subroutine build_tangents( &
            lorentz_dot, pressure_stress_dot, inertia_dot, vacuum_dot, wall_dot, &
            resistive_dot, singular_dot, frequency_dot)
        complex(dp), intent(out) :: lorentz_dot(:, :), pressure_stress_dot(:, :)
        complex(dp), intent(out) :: inertia_dot(:, :), vacuum_dot(:, :)
        complex(dp), intent(out) :: wall_dot(:, :), resistive_dot(:, :)
        complex(dp), intent(out) :: singular_dot(:, :), frequency_dot
        complex(dp) :: seed(size(lorentz_dot, 1), size(lorentz_dot, 2))

        seed = reshape([ &
            cmplx(0.03_dp, -0.02_dp, dp), cmplx(-0.01_dp, 0.04_dp, dp), &
            cmplx(0.02_dp, 0.01_dp, dp), cmplx(-0.05_dp, 0.03_dp, dp)], &
            shape(seed))
        lorentz_dot = seed
        pressure_stress_dot = cmplx(0.7_dp, -0.1_dp, dp)*seed
        inertia_dot = cmplx(-0.3_dp, 0.2_dp, dp)*seed
        vacuum_dot = cmplx(0.2_dp, 0.4_dp, dp)*seed
        wall_dot = cmplx(-0.6_dp, -0.1_dp, dp)*seed
        resistive_dot = cmplx(0.5_dp, 0.3_dp, dp)*seed
        singular_dot = cmplx(-0.2_dp, 0.6_dp, dp)*seed
        frequency_dot = cmplx(0.04_dp, -0.03_dp, dp)
    end subroutine build_tangents

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_linear_perturbation_composition
