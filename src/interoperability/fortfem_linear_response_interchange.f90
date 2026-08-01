module fortfem_linear_response_interchange
    !! Neutral interchange and block-composition contracts for linear response.
    !!
    !! The record carries only externally supplied modal metadata and complex
    !! linear blocks.  It deliberately does not select an equilibrium,
    !! constitutive law, singular-layer model, or plasma normalization.
    !! External GPEC, MARS-F, GLISS, or STARWALL adapters can populate the
    !! record without importing their file formats or application physics.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    complex(dp), parameter :: imaginary_unit = cmplx(0.0_dp, 1.0_dp, dp)

    type, public :: linear_response_interchange_t
        character(len=32) :: schema_version = "fortfem-linear-response-1"
        character(len=64) :: producer = ""
        character(len=128) :: provenance = ""
        integer :: state_count = 0
        integer :: mode_count = 0
        integer :: response_count = 0
        complex(dp) :: frequency = cmplx(0.0_dp, 0.0_dp, dp)
        real(dp) :: normalization_scale = 1.0_dp
        integer, allocatable :: mode_numbers(:, :)
        character(len=64), allocatable :: mode_labels(:)
        complex(dp), allocatable :: equilibrium_block(:, :)
        complex(dp), allocatable :: inertia_block(:, :)
        complex(dp), allocatable :: resistive_block(:, :)
        complex(dp), allocatable :: vacuum_block(:, :)
        complex(dp), allocatable :: wall_block(:, :)
        complex(dp), allocatable :: response_matrix(:, :)
        character(len=64), allocatable :: response_labels(:)
    contains
        procedure, private :: assign_linear_response_interchange
        generic :: assignment(=) => assign_linear_response_interchange
    end type linear_response_interchange_t

    public :: initialize_linear_response_interchange
    public :: validate_linear_response_interchange
    public :: evaluate_linear_response_diagnostics
    public :: assemble_linear_response_operator
    public :: assemble_linear_response_operator_jvp
    public :: assemble_linear_response_operator_vjp

    interface finite_complex
        module procedure finite_complex_scalar
        module procedure finite_complex_matrix
    end interface finite_complex

contains

    subroutine initialize_linear_response_interchange( &
            interchange, frequency, mode_numbers, mode_labels, equilibrium, &
            inertia, resistive, vacuum, wall, response_matrix, response_labels, &
            normalization_scale, status)
        type(linear_response_interchange_t), intent(inout) :: interchange
        complex(dp), intent(in) :: frequency
        integer, intent(in) :: mode_numbers(:, :)
        character(len=*), intent(in) :: mode_labels(:)
        complex(dp), intent(in) :: equilibrium(:, :), inertia(:, :)
        complex(dp), intent(in) :: resistive(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(in) :: response_matrix(:, :)
        character(len=*), intent(in) :: response_labels(:)
        real(dp), intent(in) :: normalization_scale
        integer, intent(out) :: status

        call clear_linear_response_interchange(interchange)
        status = 1
        if (.not. compatible_inputs( &
            frequency, mode_numbers, mode_labels, equilibrium, inertia, &
            resistive, vacuum, wall, response_matrix, response_labels, &
            normalization_scale)) return

        interchange%frequency = frequency
        interchange%normalization_scale = normalization_scale
        interchange%state_count = size(equilibrium, 1)
        interchange%mode_count = size(mode_numbers, 2)
        interchange%response_count = size(response_matrix, 1)
        allocate(interchange%mode_numbers, source=mode_numbers)
        allocate(interchange%mode_labels(interchange%mode_count))
        allocate(interchange%equilibrium_block, source=equilibrium)
        allocate(interchange%inertia_block, source=inertia)
        allocate(interchange%resistive_block, source=resistive)
        allocate(interchange%vacuum_block, source=vacuum)
        allocate(interchange%wall_block, source=wall)
        allocate(interchange%response_matrix, source=response_matrix)
        allocate(interchange%response_labels(interchange%response_count))
        if (interchange%mode_count > 0) interchange%mode_labels = mode_labels
        if (interchange%response_count > 0) then
            interchange%response_labels = response_labels
        end if

        if (.not. validate_linear_response_interchange(interchange, status)) then
            call clear_linear_response_interchange(interchange)
            return
        end if
    end subroutine initialize_linear_response_interchange

    logical function validate_linear_response_interchange(interchange, status) &
            result(valid)
        type(linear_response_interchange_t), intent(in) :: interchange
        integer, intent(out) :: status

        valid = .false.
        status = 1
        if (interchange%schema_version /= "fortfem-linear-response-1") return
        if (interchange%state_count < 1 .or. interchange%mode_count < 1 .or. &
            interchange%response_count < 0) return
        if (.not. finite_complex(interchange%frequency) .or. &
            .not. ieee_is_finite(interchange%normalization_scale) .or. &
            interchange%normalization_scale <= 0.0_dp) return
        if (.not. allocated(interchange%mode_numbers) .or. &
            .not. allocated(interchange%mode_labels) .or. &
            .not. allocated(interchange%equilibrium_block) .or. &
            .not. allocated(interchange%inertia_block) .or. &
            .not. allocated(interchange%resistive_block) .or. &
            .not. allocated(interchange%vacuum_block) .or. &
            .not. allocated(interchange%wall_block) .or. &
            .not. allocated(interchange%response_matrix) .or. &
            .not. allocated(interchange%response_labels)) return
        if (any(shape(interchange%mode_numbers) /= [2, interchange%mode_count]) .or. &
            size(interchange%mode_labels) /= interchange%mode_count) return
        if (.not. all_square_blocks(interchange%state_count, &
            interchange%equilibrium_block, interchange%inertia_block, &
            interchange%resistive_block, interchange%vacuum_block, &
            interchange%wall_block)) return
        if (any(shape(interchange%response_matrix) /= &
            [interchange%response_count, interchange%response_count]) .or. &
            size(interchange%response_labels) /= interchange%response_count) return
        if (.not. finite_complex(interchange%equilibrium_block) .or. &
            .not. finite_complex(interchange%inertia_block) .or. &
            .not. finite_complex(interchange%resistive_block) .or. &
            .not. finite_complex(interchange%vacuum_block) .or. &
            .not. finite_complex(interchange%wall_block) .or. &
            .not. finite_complex(interchange%response_matrix)) return
        if (.not. unique_integer_pairs(interchange%mode_numbers) .or. &
            .not. unique_labels(interchange%mode_labels) .or. &
            .not. unique_labels(interchange%response_labels)) return
        valid = .true.
        status = 0
    end function validate_linear_response_interchange

    subroutine evaluate_linear_response_diagnostics( &
            response_matrix, reciprocity_error, passivity_lower_bound, status)
        complex(dp), intent(in) :: response_matrix(:, :)
        real(dp), intent(out) :: reciprocity_error, passivity_lower_bound
        integer, intent(out) :: status

        complex(dp), allocatable :: hermitian_part(:, :)
        real(dp) :: scale, radius
        integer :: row

        reciprocity_error = huge(1.0_dp)
        passivity_lower_bound = -huge(1.0_dp)
        status = 1
        if (size(response_matrix, 1) < 1 .or. &
            size(response_matrix, 1) /= size(response_matrix, 2) .or. &
            .not. finite_complex(response_matrix)) return
        scale = max(1.0_dp, maxval(abs(response_matrix)))
        reciprocity_error = maxval(abs( &
            response_matrix - transpose(response_matrix)))/scale
        allocate(hermitian_part(size(response_matrix, 1), size(response_matrix, 2)))
        hermitian_part = 0.5_dp*(response_matrix + &
            conjg(transpose(response_matrix)))
        passivity_lower_bound = huge(1.0_dp)
        do row = 1, size(response_matrix, 1)
            radius = sum(abs(hermitian_part(row, :))) - &
                abs(hermitian_part(row, row))
            passivity_lower_bound = min(passivity_lower_bound, &
                real(hermitian_part(row, row), dp) - radius)
        end do
        status = 0
    end subroutine evaluate_linear_response_diagnostics

    subroutine assemble_linear_response_operator( &
            equilibrium, inertia, resistive, vacuum, wall, frequency, &
            operator, status)
        complex(dp), intent(in) :: equilibrium(:, :), inertia(:, :)
        complex(dp), intent(in) :: resistive(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(in) :: frequency
        complex(dp), intent(out) :: operator(:, :)
        integer, intent(out) :: status

        status = 1
        if (.not. compatible_blocks(equilibrium, inertia, resistive, vacuum, wall) .or. &
            .not. finite_complex(frequency)) return
        if (any(shape(operator) /= shape(equilibrium))) return
        operator = equilibrium - frequency**2*inertia + &
            imaginary_unit*frequency*resistive + vacuum + wall
        status = 0
    end subroutine assemble_linear_response_operator

    subroutine assemble_linear_response_operator_jvp( &
            equilibrium, inertia, resistive, vacuum, wall, frequency, &
            equilibrium_dot, inertia_dot, resistive_dot, vacuum_dot, wall_dot, &
            frequency_dot, operator_dot, status)
        complex(dp), intent(in) :: equilibrium(:, :), inertia(:, :)
        complex(dp), intent(in) :: resistive(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(in) :: frequency
        complex(dp), intent(in) :: equilibrium_dot(:, :), inertia_dot(:, :)
        complex(dp), intent(in) :: resistive_dot(:, :), vacuum_dot(:, :)
        complex(dp), intent(in) :: wall_dot(:, :), frequency_dot
        complex(dp), intent(out) :: operator_dot(:, :)
        integer, intent(out) :: status

        status = 1
        if (.not. compatible_blocks(equilibrium, inertia, resistive, vacuum, wall) .or. &
            .not. finite_complex(frequency)) return
        if (.not. same_shapes(equilibrium_dot, equilibrium, inertia_dot, inertia, &
            resistive_dot, resistive, vacuum_dot, vacuum, wall_dot, wall) .or. &
            .not. finite_complex(frequency_dot)) return
        if (any(shape(operator_dot) /= shape(equilibrium))) return
        operator_dot = equilibrium_dot - frequency**2*inertia_dot - &
            2.0_dp*frequency*frequency_dot*inertia + &
            imaginary_unit*(frequency*resistive_dot + frequency_dot*resistive) + &
            vacuum_dot + wall_dot
        status = 0
    end subroutine assemble_linear_response_operator_jvp

    subroutine assemble_linear_response_operator_vjp( &
            inertia, resistive, frequency, operator_bar, equilibrium_bar, &
            inertia_bar, resistive_bar, vacuum_bar, wall_bar, frequency_bar, &
            status)
        complex(dp), intent(in) :: inertia(:, :), resistive(:, :), frequency
        complex(dp), intent(in) :: operator_bar(:, :)
        complex(dp), allocatable, intent(out) :: equilibrium_bar(:, :)
        complex(dp), allocatable, intent(out) :: inertia_bar(:, :)
        complex(dp), allocatable, intent(out) :: resistive_bar(:, :)
        complex(dp), allocatable, intent(out) :: vacuum_bar(:, :)
        complex(dp), allocatable, intent(out) :: wall_bar(:, :)
        complex(dp), intent(out) :: frequency_bar
        integer, intent(out) :: status

        complex(dp) :: frequency_coefficient(size(inertia, 1), size(inertia, 2))

        allocate(equilibrium_bar(size(operator_bar, 1), size(operator_bar, 2)))
        allocate(inertia_bar(size(operator_bar, 1), size(operator_bar, 2)))
        allocate(resistive_bar(size(operator_bar, 1), size(operator_bar, 2)))
        allocate(vacuum_bar(size(operator_bar, 1), size(operator_bar, 2)))
        allocate(wall_bar(size(operator_bar, 1), size(operator_bar, 2)))
        equilibrium_bar = cmplx(0.0_dp, 0.0_dp, dp)
        inertia_bar = cmplx(0.0_dp, 0.0_dp, dp)
        resistive_bar = cmplx(0.0_dp, 0.0_dp, dp)
        vacuum_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wall_bar = cmplx(0.0_dp, 0.0_dp, dp)
        frequency_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. compatible_blocks(inertia, inertia, resistive, resistive, &
            inertia) .or. .not. finite_complex(frequency) .or. &
            any(shape(operator_bar) /= shape(inertia))) return
        equilibrium_bar = operator_bar
        inertia_bar = conjg(-frequency**2)*operator_bar
        resistive_bar = conjg(imaginary_unit*frequency)*operator_bar
        vacuum_bar = operator_bar
        wall_bar = operator_bar
        frequency_coefficient = -2.0_dp*frequency*inertia + &
            imaginary_unit*resistive
        frequency_bar = sum(operator_bar*conjg(frequency_coefficient))
        status = 0
    end subroutine assemble_linear_response_operator_vjp

    subroutine assign_linear_response_interchange(lhs, rhs)
        class(linear_response_interchange_t), intent(out) :: lhs
        type(linear_response_interchange_t), intent(in) :: rhs

        call clear_linear_response_interchange(lhs)
        lhs%schema_version = rhs%schema_version
        lhs%producer = rhs%producer
        lhs%provenance = rhs%provenance
        lhs%state_count = rhs%state_count
        lhs%mode_count = rhs%mode_count
        lhs%response_count = rhs%response_count
        lhs%frequency = rhs%frequency
        lhs%normalization_scale = rhs%normalization_scale
        if (allocated(rhs%mode_numbers)) allocate(lhs%mode_numbers, source=rhs%mode_numbers)
        if (allocated(rhs%mode_labels)) allocate(lhs%mode_labels, source=rhs%mode_labels)
        if (allocated(rhs%equilibrium_block)) allocate(lhs%equilibrium_block, source=rhs%equilibrium_block)
        if (allocated(rhs%inertia_block)) allocate(lhs%inertia_block, source=rhs%inertia_block)
        if (allocated(rhs%resistive_block)) allocate(lhs%resistive_block, source=rhs%resistive_block)
        if (allocated(rhs%vacuum_block)) allocate(lhs%vacuum_block, source=rhs%vacuum_block)
        if (allocated(rhs%wall_block)) allocate(lhs%wall_block, source=rhs%wall_block)
        if (allocated(rhs%response_matrix)) allocate(lhs%response_matrix, source=rhs%response_matrix)
        if (allocated(rhs%response_labels)) allocate(lhs%response_labels, source=rhs%response_labels)
    end subroutine assign_linear_response_interchange

    subroutine clear_linear_response_interchange(interchange)
        type(linear_response_interchange_t), intent(inout) :: interchange

        if (allocated(interchange%mode_numbers)) deallocate(interchange%mode_numbers)
        if (allocated(interchange%mode_labels)) deallocate(interchange%mode_labels)
        if (allocated(interchange%equilibrium_block)) deallocate(interchange%equilibrium_block)
        if (allocated(interchange%inertia_block)) deallocate(interchange%inertia_block)
        if (allocated(interchange%resistive_block)) deallocate(interchange%resistive_block)
        if (allocated(interchange%vacuum_block)) deallocate(interchange%vacuum_block)
        if (allocated(interchange%wall_block)) deallocate(interchange%wall_block)
        if (allocated(interchange%response_matrix)) deallocate(interchange%response_matrix)
        if (allocated(interchange%response_labels)) deallocate(interchange%response_labels)
        interchange%state_count = 0
        interchange%mode_count = 0
        interchange%response_count = 0
    end subroutine clear_linear_response_interchange

    logical function compatible_inputs( &
            frequency, mode_numbers, mode_labels, equilibrium, inertia, &
            resistive, vacuum, wall, response_matrix, response_labels, scale) &
            result(valid)
        complex(dp), intent(in) :: frequency
        integer, intent(in) :: mode_numbers(:, :)
        character(len=*), intent(in) :: mode_labels(:)
        complex(dp), intent(in) :: equilibrium(:, :), inertia(:, :)
        complex(dp), intent(in) :: resistive(:, :), vacuum(:, :), wall(:, :)
        complex(dp), intent(in) :: response_matrix(:, :)
        character(len=*), intent(in) :: response_labels(:)
        real(dp), intent(in) :: scale

        valid = finite_complex(frequency) .and. ieee_is_finite(scale) .and. &
            scale > 0.0_dp .and. size(mode_numbers, 1) == 2 .and. &
            size(mode_numbers, 2) > 0 .and. size(mode_labels) == size(mode_numbers, 2) .and. &
            compatible_blocks(equilibrium, inertia, resistive, vacuum, wall) .and. &
            size(response_matrix, 1) == size(response_matrix, 2) .and. &
            size(response_labels) == size(response_matrix, 1) .and. &
            unique_integer_pairs(mode_numbers) .and. unique_labels(mode_labels) .and. &
            unique_labels(response_labels) .and. finite_complex(equilibrium) .and. &
            finite_complex(inertia) .and. finite_complex(resistive) .and. &
            finite_complex(vacuum) .and. finite_complex(wall) .and. &
            finite_complex(response_matrix)
    end function compatible_inputs

    logical function compatible_blocks(first, second, third, fourth, fifth) result(valid)
        complex(dp), intent(in) :: first(:, :), second(:, :), third(:, :)
        complex(dp), intent(in) :: fourth(:, :), fifth(:, :)

        valid = size(first, 1) > 0 .and. size(first, 1) == size(first, 2) .and. &
            all(shape(second) == shape(first)) .and. &
            all(shape(third) == shape(first)) .and. &
            all(shape(fourth) == shape(first)) .and. &
            all(shape(fifth) == shape(first)) .and. finite_complex(first) .and. &
            finite_complex(second) .and. finite_complex(third) .and. &
            finite_complex(fourth) .and. finite_complex(fifth)
    end function compatible_blocks

    logical function same_shapes(first, first_reference, second, second_reference, &
            third, third_reference, fourth, fourth_reference, fifth, fifth_reference) &
            result(valid)
        complex(dp), intent(in) :: first(:, :), first_reference(:, :)
        complex(dp), intent(in) :: second(:, :), second_reference(:, :)
        complex(dp), intent(in) :: third(:, :), third_reference(:, :)
        complex(dp), intent(in) :: fourth(:, :), fourth_reference(:, :)
        complex(dp), intent(in) :: fifth(:, :), fifth_reference(:, :)

        valid = all(shape(first) == shape(first_reference)) .and. &
            all(shape(second) == shape(second_reference)) .and. &
            all(shape(third) == shape(third_reference)) .and. &
            all(shape(fourth) == shape(fourth_reference)) .and. &
            all(shape(fifth) == shape(fifth_reference)) .and. &
            finite_complex(first) .and. finite_complex(second) .and. &
            finite_complex(third) .and. finite_complex(fourth) .and. &
            finite_complex(fifth)
    end function same_shapes

    pure logical function finite_complex_scalar(value) result(valid)
        complex(dp), intent(in) :: value

        valid = ieee_is_finite(real(value, dp)) .and. ieee_is_finite(aimag(value))
    end function finite_complex_scalar

    logical function finite_complex_matrix(value) result(valid)
        complex(dp), intent(in) :: value(:, :)

        valid = all(ieee_is_finite(real(value, dp))) .and. &
            all(ieee_is_finite(aimag(value)))
    end function finite_complex_matrix

    logical function all_square_blocks(state_count, equilibrium, inertia, &
            resistive, vacuum, wall) result(valid)
        integer, intent(in) :: state_count
        complex(dp), intent(in) :: equilibrium(:, :), inertia(:, :)
        complex(dp), intent(in) :: resistive(:, :), vacuum(:, :), wall(:, :)

        valid = state_count > 0 .and. size(equilibrium, 1) == state_count .and. &
            compatible_blocks(equilibrium, inertia, resistive, vacuum, wall)
    end function all_square_blocks

    logical function unique_integer_pairs(values) result(unique)
        integer, intent(in) :: values(:, :)
        integer :: first, second

        unique = size(values, 1) == 2
        if (.not. unique) return
        do first = 1, size(values, 2)
            do second = first + 1, size(values, 2)
                if (all(values(:, first) == values(:, second))) then
                    unique = .false.
                    return
                end if
            end do
        end do
    end function unique_integer_pairs

    logical function unique_labels(labels) result(unique)
        character(len=*), intent(in) :: labels(:)
        integer :: first, second

        unique = .true.
        do first = 1, size(labels)
            if (len_trim(labels(first)) == 0) then
                unique = .false.
                return
            end if
            do second = first + 1, size(labels)
                if (trim(labels(first)) == trim(labels(second))) then
                    unique = .false.
                    return
                end if
            end do
        end do
    end function unique_labels

end module fortfem_linear_response_interchange
