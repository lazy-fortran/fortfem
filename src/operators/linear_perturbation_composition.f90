module fortfem_linear_perturbation_composition
    !! Closure-neutral linear ideal/resistive perturbation block composition.
    !!
    !! All matrices and conventions are supplied by the caller.  This module
    !! separates the reusable algebra from plasma models and file formats.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    complex(dp), parameter :: imaginary_unit = cmplx(0.0_dp, 1.0_dp, dp)

    integer, parameter, public :: LINEAR_RECIPROCITY_UNDECLARED = 0
    integer, parameter, public :: LINEAR_RECIPROCITY_TRANSPOSE = 1
    integer, parameter, public :: LINEAR_PASSIVITY_UNDECLARED = 0
    integer, parameter, public :: LINEAR_PASSIVITY_HERMITIAN_NONNEGATIVE = 1

    type, public :: linear_perturbation_metadata_t
        character(len=32) :: schema_version = "fortfem-linear-perturbation-1"
        character(len=64) :: normalization_label = ""
        character(len=128) :: provenance = ""
        real(dp) :: normalization_scale = 1.0_dp
        integer :: reciprocity_convention = LINEAR_RECIPROCITY_UNDECLARED
        integer :: passivity_convention = LINEAR_PASSIVITY_UNDECLARED
    end type linear_perturbation_metadata_t

    public :: initialize_linear_perturbation_metadata
    public :: validate_linear_perturbation_metadata
    public :: assemble_linear_perturbation_operator
    public :: assemble_linear_perturbation_operator_jvp
    public :: assemble_linear_perturbation_operator_vjp

    interface finite_complex
        module procedure finite_complex_scalar
        module procedure finite_complex_matrix
    end interface finite_complex

contains

    subroutine initialize_linear_perturbation_metadata( &
            metadata, normalization_label, normalization_scale, provenance, &
            reciprocity_convention, passivity_convention, status)
        type(linear_perturbation_metadata_t), intent(out) :: metadata
        character(len=*), intent(in) :: normalization_label, provenance
        real(dp), intent(in) :: normalization_scale
        integer, intent(in) :: reciprocity_convention, passivity_convention
        integer, intent(out) :: status

        metadata = linear_perturbation_metadata_t()
        status = 1
        if (len_trim(normalization_label) < 1) return
        if (len_trim(normalization_label) > len(metadata%normalization_label)) return
        if (len_trim(provenance) < 1) return
        if (len_trim(provenance) > len(metadata%provenance)) return
        if (.not. ieee_is_finite(normalization_scale)) return
        if (normalization_scale <= 0.0_dp) return
        if (.not. valid_reciprocity_convention(reciprocity_convention)) return
        if (.not. valid_passivity_convention(passivity_convention)) return

        metadata%normalization_label = trim(normalization_label)
        metadata%normalization_scale = normalization_scale
        metadata%provenance = trim(provenance)
        metadata%reciprocity_convention = reciprocity_convention
        metadata%passivity_convention = passivity_convention
        if (.not. validate_linear_perturbation_metadata(metadata, status)) return
    end subroutine initialize_linear_perturbation_metadata

    logical function validate_linear_perturbation_metadata(metadata, status) &
            result(valid)
        type(linear_perturbation_metadata_t), intent(in) :: metadata
        integer, intent(out) :: status

        valid = .false.
        status = 1
        if (metadata%schema_version /= "fortfem-linear-perturbation-1") return
        if (len_trim(metadata%normalization_label) < 1) return
        if (len_trim(metadata%provenance) < 1) return
        if (.not. ieee_is_finite(metadata%normalization_scale)) return
        if (metadata%normalization_scale <= 0.0_dp) return
        if (.not. valid_reciprocity_convention( &
            metadata%reciprocity_convention)) return
        if (.not. valid_passivity_convention(metadata%passivity_convention)) return
        valid = .true.
        status = 0
    end function validate_linear_perturbation_metadata

    subroutine assemble_linear_perturbation_operator( &
            lorentz, pressure_stress, inertia, vacuum, wall, resistive, singular, &
            frequency, operator, status)
        complex(dp), intent(in) :: lorentz(:, :), pressure_stress(:, :)
        complex(dp), intent(in) :: inertia(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(in) :: resistive(:, :), singular(:, :), frequency
        complex(dp), intent(out) :: operator(:, :)
        integer, intent(out) :: status

        operator = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_blocks( &
            lorentz, pressure_stress, inertia, vacuum, wall, resistive, &
            singular)) return
        if (.not. finite_complex(frequency)) return
        if (.not. same_shape(operator, lorentz)) return

        operator = lorentz + pressure_stress + vacuum + wall + singular - &
            frequency**2*inertia + imaginary_unit*frequency*resistive
        status = 0
    end subroutine assemble_linear_perturbation_operator

    subroutine assemble_linear_perturbation_operator_jvp( &
            lorentz, pressure_stress, inertia, vacuum, wall, resistive, singular, &
            frequency, lorentz_dot, pressure_stress_dot, inertia_dot, vacuum_dot, &
            wall_dot, resistive_dot, singular_dot, frequency_dot, operator_dot, &
            status)
        complex(dp), intent(in) :: lorentz(:, :), pressure_stress(:, :)
        complex(dp), intent(in) :: inertia(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(in) :: resistive(:, :), singular(:, :), frequency
        complex(dp), intent(in) :: lorentz_dot(:, :), pressure_stress_dot(:, :)
        complex(dp), intent(in) :: inertia_dot(:, :), vacuum_dot(:, :)
        complex(dp), intent(in) :: wall_dot(:, :), resistive_dot(:, :)
        complex(dp), intent(in) :: singular_dot(:, :), frequency_dot
        complex(dp), intent(out) :: operator_dot(:, :)
        integer, intent(out) :: status

        operator_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_blocks( &
            lorentz, pressure_stress, inertia, vacuum, wall, resistive, &
            singular)) return
        if (.not. finite_complex(frequency)) return
        if (.not. valid_block_tangents( &
            lorentz, pressure_stress, inertia, vacuum, wall, resistive, singular, &
            lorentz_dot, pressure_stress_dot, inertia_dot, vacuum_dot, wall_dot, &
            resistive_dot, singular_dot)) return
        if (.not. finite_complex(frequency_dot)) return
        if (.not. same_shape(operator_dot, lorentz)) return

        operator_dot = lorentz_dot + pressure_stress_dot + vacuum_dot + &
            wall_dot + singular_dot - frequency**2*inertia_dot - &
            2.0_dp*frequency*frequency_dot*inertia + &
            imaginary_unit*(frequency*resistive_dot + frequency_dot*resistive)
        status = 0
    end subroutine assemble_linear_perturbation_operator_jvp

    subroutine assemble_linear_perturbation_operator_vjp( &
            lorentz, pressure_stress, inertia, vacuum, wall, resistive, singular, &
            frequency, operator_bar, lorentz_bar, pressure_stress_bar, &
            inertia_bar, vacuum_bar, wall_bar, resistive_bar, singular_bar, &
            frequency_bar, status)
        complex(dp), intent(in) :: lorentz(:, :), pressure_stress(:, :)
        complex(dp), intent(in) :: inertia(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(in) :: resistive(:, :), singular(:, :), frequency
        complex(dp), intent(in) :: operator_bar(:, :)
        complex(dp), intent(out) :: lorentz_bar(:, :), pressure_stress_bar(:, :)
        complex(dp), intent(out) :: inertia_bar(:, :), vacuum_bar(:, :)
        complex(dp), intent(out) :: wall_bar(:, :), resistive_bar(:, :)
        complex(dp), intent(out) :: singular_bar(:, :), frequency_bar
        integer, intent(out) :: status
        complex(dp) :: frequency_coefficient(size(inertia, 1), size(inertia, 2))

        lorentz_bar = cmplx(0.0_dp, 0.0_dp, dp)
        pressure_stress_bar = cmplx(0.0_dp, 0.0_dp, dp)
        inertia_bar = cmplx(0.0_dp, 0.0_dp, dp)
        vacuum_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wall_bar = cmplx(0.0_dp, 0.0_dp, dp)
        resistive_bar = cmplx(0.0_dp, 0.0_dp, dp)
        singular_bar = cmplx(0.0_dp, 0.0_dp, dp)
        frequency_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_blocks( &
            lorentz, pressure_stress, inertia, vacuum, wall, resistive, &
            singular)) return
        if (.not. finite_complex(frequency)) return
        if (.not. same_shape(operator_bar, lorentz)) return
        if (.not. finite_complex(operator_bar)) return
        if (.not. valid_vjp_outputs( &
            lorentz, lorentz_bar, pressure_stress_bar, inertia_bar, vacuum_bar, &
            wall_bar, resistive_bar, singular_bar)) return

        lorentz_bar = operator_bar
        pressure_stress_bar = operator_bar
        inertia_bar = conjg(-frequency**2)*operator_bar
        vacuum_bar = operator_bar
        wall_bar = operator_bar
        resistive_bar = conjg(imaginary_unit*frequency)*operator_bar
        singular_bar = operator_bar
        frequency_coefficient = -2.0_dp*frequency*inertia + &
            imaginary_unit*resistive
        frequency_bar = sum(operator_bar*conjg(frequency_coefficient))
        status = 0
    end subroutine assemble_linear_perturbation_operator_vjp

    logical function valid_blocks( &
            lorentz, pressure_stress, inertia, vacuum, wall, resistive, singular) &
            result(valid)
        complex(dp), intent(in) :: lorentz(:, :), pressure_stress(:, :)
        complex(dp), intent(in) :: inertia(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(in) :: resistive(:, :), singular(:, :)

        valid = .false.
        if (size(lorentz, 1) < 1) return
        if (size(lorentz, 1) /= size(lorentz, 2)) return
        if (.not. same_shape(pressure_stress, lorentz)) return
        if (.not. same_shape(inertia, lorentz)) return
        if (.not. same_shape(vacuum, lorentz)) return
        if (.not. same_shape(wall, lorentz)) return
        if (.not. same_shape(resistive, lorentz)) return
        if (.not. same_shape(singular, lorentz)) return
        if (.not. finite_complex(lorentz)) return
        if (.not. finite_complex(pressure_stress)) return
        if (.not. finite_complex(inertia)) return
        if (.not. finite_complex(vacuum)) return
        if (.not. finite_complex(wall)) return
        if (.not. finite_complex(resistive)) return
        if (.not. finite_complex(singular)) return
        valid = .true.
    end function valid_blocks

    logical function valid_block_tangents( &
            lorentz, pressure_stress, inertia, vacuum, wall, resistive, singular, &
            lorentz_dot, pressure_stress_dot, inertia_dot, vacuum_dot, wall_dot, &
            resistive_dot, singular_dot) result(valid)
        complex(dp), intent(in) :: lorentz(:, :), pressure_stress(:, :)
        complex(dp), intent(in) :: inertia(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(in) :: resistive(:, :), singular(:, :)
        complex(dp), intent(in) :: lorentz_dot(:, :), pressure_stress_dot(:, :)
        complex(dp), intent(in) :: inertia_dot(:, :), vacuum_dot(:, :)
        complex(dp), intent(in) :: wall_dot(:, :), resistive_dot(:, :)
        complex(dp), intent(in) :: singular_dot(:, :)

        valid = .false.
        if (.not. same_shape(lorentz_dot, lorentz)) return
        if (.not. same_shape(pressure_stress_dot, pressure_stress)) return
        if (.not. same_shape(inertia_dot, inertia)) return
        if (.not. same_shape(vacuum_dot, vacuum)) return
        if (.not. same_shape(wall_dot, wall)) return
        if (.not. same_shape(resistive_dot, resistive)) return
        if (.not. same_shape(singular_dot, singular)) return
        if (.not. finite_complex(lorentz_dot)) return
        if (.not. finite_complex(pressure_stress_dot)) return
        if (.not. finite_complex(inertia_dot)) return
        if (.not. finite_complex(vacuum_dot)) return
        if (.not. finite_complex(wall_dot)) return
        if (.not. finite_complex(resistive_dot)) return
        if (.not. finite_complex(singular_dot)) return
        valid = .true.
    end function valid_block_tangents

    logical function valid_vjp_outputs( &
            reference, lorentz_bar, pressure_stress_bar, inertia_bar, vacuum_bar, &
            wall_bar, resistive_bar, singular_bar) result(valid)
        complex(dp), intent(in) :: reference(:, :)
        complex(dp), intent(in) :: lorentz_bar(:, :), pressure_stress_bar(:, :)
        complex(dp), intent(in) :: inertia_bar(:, :), vacuum_bar(:, :)
        complex(dp), intent(in) :: wall_bar(:, :), resistive_bar(:, :)
        complex(dp), intent(in) :: singular_bar(:, :)

        valid = same_shape(lorentz_bar, reference)
        if (.not. valid) return
        valid = same_shape(pressure_stress_bar, reference)
        if (.not. valid) return
        valid = same_shape(inertia_bar, reference)
        if (.not. valid) return
        valid = same_shape(vacuum_bar, reference)
        if (.not. valid) return
        valid = same_shape(wall_bar, reference)
        if (.not. valid) return
        valid = same_shape(resistive_bar, reference)
        if (.not. valid) return
        valid = same_shape(singular_bar, reference)
    end function valid_vjp_outputs

    pure logical function same_shape(first, second) result(same)
        complex(dp), intent(in) :: first(:, :), second(:, :)

        same = all(shape(first) == shape(second))
    end function same_shape

    pure logical function valid_reciprocity_convention(convention) result(valid)
        integer, intent(in) :: convention

        valid = convention == LINEAR_RECIPROCITY_UNDECLARED .or. &
            convention == LINEAR_RECIPROCITY_TRANSPOSE
    end function valid_reciprocity_convention

    pure logical function valid_passivity_convention(convention) result(valid)
        integer, intent(in) :: convention

        valid = convention == LINEAR_PASSIVITY_UNDECLARED .or. &
            convention == LINEAR_PASSIVITY_HERMITIAN_NONNEGATIVE
    end function valid_passivity_convention

    pure logical function finite_complex_scalar(value) result(valid)
        complex(dp), intent(in) :: value

        valid = ieee_is_finite(real(value, dp)) .and. &
            ieee_is_finite(aimag(value))
    end function finite_complex_scalar

    logical function finite_complex_matrix(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_matrix

end module fortfem_linear_perturbation_composition
