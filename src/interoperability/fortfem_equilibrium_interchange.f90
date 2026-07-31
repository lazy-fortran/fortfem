module fortfem_equilibrium_interchange
    !! Neutral interchange metadata for externally produced equilibria.
    !!
    !! The record deliberately contains samples and provenance only.  It does
    !! not parse an equilibrium file, select a coordinate convention, or
    !! implement a plasma closure.  External adapters may populate it from a
    !! license-compatible source and FortFEM clients can then consume the same
    !! mapped/physical samples, coefficients, profiles, and boundaries.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    type, public :: equilibrium_normalization_t
        character(len=32) :: length_unit = "1"
        character(len=32) :: magnetic_field_unit = "1"
        character(len=32) :: pressure_unit = "1"
        character(len=32) :: current_unit = "1"
        real(dp) :: length_scale = 1.0_dp
        real(dp) :: magnetic_field_scale = 1.0_dp
        real(dp) :: pressure_scale = 1.0_dp
        real(dp) :: current_scale = 1.0_dp
    end type equilibrium_normalization_t

    type, public :: equilibrium_interchange_t
        character(len=32) :: schema_version = "fortfem-equilibrium-1"
        character(len=64) :: producer = ""
        character(len=128) :: provenance = ""
        integer :: spatial_dimension = 0
        integer :: sample_count = 0
        integer :: coefficient_count = 0
        integer :: coefficient_component_count = 0
        integer :: profile_count = 0
        integer :: profile_sample_count = 0
        integer :: boundary_count = 0
        type(equilibrium_normalization_t) :: normalization
        real(dp), allocatable :: mapped_coordinates(:, :)
        real(dp), allocatable :: physical_coordinates(:, :)
        character(len=64), allocatable :: coefficient_names(:)
        integer, allocatable :: coefficient_ranks(:)
        integer, allocatable :: coefficient_offsets(:)
        character(len=64), allocatable :: coefficient_component_names(:)
        real(dp), allocatable :: coefficient_values(:, :)
        real(dp), allocatable :: profile_coordinates(:)
        character(len=64), allocatable :: profile_names(:)
        real(dp), allocatable :: profile_values(:, :)
        integer, allocatable :: boundary_offsets(:)
        character(len=64), allocatable :: boundary_names(:)
        real(dp), allocatable :: boundary_coordinates(:, :)
    contains
        procedure, private :: assign_equilibrium_interchange
        generic :: assignment(=) => assign_equilibrium_interchange
    end type equilibrium_interchange_t

    public :: initialize_equilibrium_interchange
    public :: validate_equilibrium_interchange

contains

    subroutine initialize_equilibrium_interchange( &
            interchange, mapped_coordinates, physical_coordinates, &
            coefficient_names, coefficient_values, profile_coordinates, &
            profile_names, profile_values, boundary_offsets, boundary_names, &
            boundary_coordinates, normalization, status, coefficient_ranks, &
            coefficient_offsets, coefficient_component_names)
        type(equilibrium_interchange_t), intent(inout) :: interchange
        real(dp), intent(in) :: mapped_coordinates(:, :)
        real(dp), intent(in) :: physical_coordinates(:, :)
        character(len=*), intent(in) :: coefficient_names(:)
        real(dp), intent(in) :: coefficient_values(:, :)
        real(dp), intent(in) :: profile_coordinates(:)
        character(len=*), intent(in) :: profile_names(:)
        real(dp), intent(in) :: profile_values(:, :)
        integer, intent(in) :: boundary_offsets(:)
        character(len=*), intent(in) :: boundary_names(:)
        real(dp), intent(in) :: boundary_coordinates(:, :)
        type(equilibrium_normalization_t), intent(in) :: normalization
        integer, intent(out) :: status
        integer, intent(in), optional :: coefficient_ranks(:)
        integer, intent(in), optional :: coefficient_offsets(:)
        character(len=*), intent(in), optional :: coefficient_component_names(:)

        integer :: component_count, field

        call clear_equilibrium_interchange(interchange)
        status = 1
        if (.not. input_shapes_are_compatible( &
            mapped_coordinates, physical_coordinates, coefficient_names, &
            coefficient_values, profile_coordinates, profile_names, &
            profile_values, boundary_offsets, boundary_names, &
            boundary_coordinates)) return
        component_count = size(coefficient_values, 1)
        if (present(coefficient_ranks)) then
            if (size(coefficient_ranks) /= size(coefficient_names)) return
        end if
        if (present(coefficient_offsets)) then
            if (size(coefficient_offsets) /= size(coefficient_names) + 1) return
        end if
        if (present(coefficient_component_names)) then
            if (size(coefficient_component_names) /= component_count) return
        end if

        interchange%normalization = normalization
        interchange%spatial_dimension = size(mapped_coordinates, 1)
        interchange%sample_count = size(mapped_coordinates, 2)
        interchange%coefficient_count = size(coefficient_names)
        interchange%coefficient_component_count = component_count
        interchange%profile_count = size(profile_names)
        interchange%profile_sample_count = size(profile_coordinates)
        interchange%boundary_count = size(boundary_names)
        allocate(interchange%mapped_coordinates, source=mapped_coordinates)
        allocate(interchange%physical_coordinates, source=physical_coordinates)
        allocate(interchange%coefficient_names(interchange%coefficient_count))
        allocate(interchange%coefficient_ranks(interchange%coefficient_count))
        allocate(interchange%coefficient_offsets(interchange%coefficient_count + 1))
        allocate(interchange%coefficient_component_names(component_count))
        allocate(interchange%coefficient_values, source=coefficient_values)
        allocate(interchange%profile_coordinates, source=profile_coordinates)
        allocate(interchange%profile_names(interchange%profile_count))
        allocate(interchange%profile_values, source=profile_values)
        allocate(interchange%boundary_offsets, source=boundary_offsets)
        allocate(interchange%boundary_names(interchange%boundary_count))
        allocate(interchange%boundary_coordinates, source=boundary_coordinates)
        if (interchange%coefficient_count > 0) then
            interchange%coefficient_names = coefficient_names
            interchange%coefficient_ranks = 0
            interchange%coefficient_offsets = [ &
                (field, field=1, interchange%coefficient_count + 1)]
        else
            interchange%coefficient_ranks = 0
            interchange%coefficient_offsets = 1
        end if
        if (present(coefficient_ranks)) then
            interchange%coefficient_ranks = coefficient_ranks
        end if
        if (present(coefficient_offsets)) then
            interchange%coefficient_offsets = coefficient_offsets
        end if
        if (present(coefficient_component_names)) then
            interchange%coefficient_component_names = coefficient_component_names
        else if (component_count == interchange%coefficient_count) then
            interchange%coefficient_component_names = interchange%coefficient_names
        end if
        if (interchange%profile_count > 0) then
            interchange%profile_names = profile_names
        end if
        if (interchange%boundary_count > 0) then
            interchange%boundary_names = boundary_names
        end if

        if (.not. validate_equilibrium_interchange(interchange, status)) then
            call clear_equilibrium_interchange(interchange)
            return
        end if
    end subroutine initialize_equilibrium_interchange

    logical function validate_equilibrium_interchange(interchange, status) &
            result(valid)
        type(equilibrium_interchange_t), intent(in) :: interchange
        integer, intent(out) :: status
        integer :: field

        valid = .false.
        status = 1
        if (interchange%schema_version /= "fortfem-equilibrium-1") return
        if (interchange%spatial_dimension < 1 .or. &
            interchange%spatial_dimension > 3) return
        if (interchange%sample_count < 1 .or. &
            interchange%coefficient_count < 0 .or. &
            interchange%coefficient_component_count < 0 .or. &
            interchange%profile_count < 0 .or. &
            interchange%profile_sample_count < 1 .or. &
            interchange%boundary_count < 1) return
        if (.not. allocated(interchange%mapped_coordinates) .or. &
            .not. allocated(interchange%physical_coordinates) .or. &
            .not. allocated(interchange%coefficient_names) .or. &
            .not. allocated(interchange%coefficient_ranks) .or. &
            .not. allocated(interchange%coefficient_offsets) .or. &
            .not. allocated(interchange%coefficient_component_names) .or. &
            .not. allocated(interchange%coefficient_values) .or. &
            .not. allocated(interchange%profile_coordinates) .or. &
            .not. allocated(interchange%profile_names) .or. &
            .not. allocated(interchange%profile_values) .or. &
            .not. allocated(interchange%boundary_offsets) .or. &
            .not. allocated(interchange%boundary_names) .or. &
            .not. allocated(interchange%boundary_coordinates)) return
        if (size(interchange%mapped_coordinates, 1) /= &
            interchange%spatial_dimension .or. &
            size(interchange%mapped_coordinates, 2) /= interchange%sample_count .or. &
            any(shape(interchange%physical_coordinates) /= &
            shape(interchange%mapped_coordinates))) return
        if (any(shape(interchange%coefficient_values) /= &
            [interchange%coefficient_component_count, interchange%sample_count]) .or. &
            size(interchange%coefficient_names) /= interchange%coefficient_count .or. &
            size(interchange%coefficient_ranks) /= interchange%coefficient_count .or. &
            size(interchange%coefficient_offsets) /= interchange%coefficient_count + 1 .or. &
            size(interchange%coefficient_component_names) /= &
            interchange%coefficient_component_count) return
        if (interchange%coefficient_offsets(1) /= 1 .or. &
            interchange%coefficient_offsets(interchange%coefficient_count + 1) /= &
            interchange%coefficient_component_count + 1) return
        if (interchange%coefficient_count > 0) then
            if (any(interchange%coefficient_offsets(2:) < &
                interchange%coefficient_offsets(:interchange%coefficient_count))) return
        end if
        do field = 1, interchange%coefficient_count
            if (interchange%coefficient_ranks(field) < 0 .or. &
                interchange%coefficient_ranks(field) > 2) return
            if (interchange%coefficient_offsets(field + 1) - &
                interchange%coefficient_offsets(field) /= &
                expected_component_count(interchange%coefficient_ranks(field), &
                interchange%spatial_dimension)) return
        end do
        if (size(interchange%profile_coordinates) /= &
            interchange%profile_sample_count .or. &
            any(shape(interchange%profile_values) /= &
            [interchange%profile_count, interchange%profile_sample_count]) .or. &
            size(interchange%profile_names) /= interchange%profile_count) return
        if (size(interchange%boundary_offsets) /= interchange%boundary_count + 1 .or. &
            size(interchange%boundary_names) /= interchange%boundary_count .or. &
            size(interchange%boundary_coordinates, 1) /= interchange%spatial_dimension) return
        if (interchange%boundary_offsets(1) /= 1 .or. &
            interchange%boundary_offsets(interchange%boundary_count + 1) /= &
            size(interchange%boundary_coordinates, 2) + 1) return
        if (any(interchange%boundary_offsets < 1) .or. &
            any(interchange%boundary_offsets > size(interchange%boundary_coordinates, 2) + 1) .or. &
            any(interchange%boundary_offsets(2:) < interchange%boundary_offsets(: &
            size(interchange%boundary_offsets) - 1))) return
        if (.not. finite_normalization(interchange%normalization)) return
        if (.not. all(ieee_is_finite(interchange%mapped_coordinates)) .or. &
            .not. all(ieee_is_finite(interchange%physical_coordinates)) .or. &
            .not. all(ieee_is_finite(interchange%coefficient_values)) .or. &
            .not. all(ieee_is_finite(interchange%profile_coordinates)) .or. &
            .not. all(ieee_is_finite(interchange%profile_values)) .or. &
            .not. all(ieee_is_finite(interchange%boundary_coordinates))) return
        if (.not. labels_are_unique_and_nonempty(interchange%coefficient_names) .or. &
            .not. labels_are_unique_and_nonempty(interchange%profile_names) .or. &
            .not. labels_are_unique_and_nonempty(interchange%boundary_names) .or. &
            .not. labels_are_nonempty(interchange%coefficient_component_names)) return
        valid = .true.
        status = 0
    end function validate_equilibrium_interchange

    subroutine assign_equilibrium_interchange(lhs, rhs)
        class(equilibrium_interchange_t), intent(out) :: lhs
        type(equilibrium_interchange_t), intent(in) :: rhs

        call clear_equilibrium_interchange(lhs)
        lhs%schema_version = rhs%schema_version
        lhs%producer = rhs%producer
        lhs%provenance = rhs%provenance
        lhs%spatial_dimension = rhs%spatial_dimension
        lhs%sample_count = rhs%sample_count
        lhs%coefficient_count = rhs%coefficient_count
        lhs%coefficient_component_count = rhs%coefficient_component_count
        lhs%profile_count = rhs%profile_count
        lhs%profile_sample_count = rhs%profile_sample_count
        lhs%boundary_count = rhs%boundary_count
        lhs%normalization = rhs%normalization
        if (allocated(rhs%mapped_coordinates)) allocate( &
            lhs%mapped_coordinates, source=rhs%mapped_coordinates)
        if (allocated(rhs%physical_coordinates)) allocate( &
            lhs%physical_coordinates, source=rhs%physical_coordinates)
        if (allocated(rhs%coefficient_names)) allocate( &
            lhs%coefficient_names, source=rhs%coefficient_names)
        if (allocated(rhs%coefficient_ranks)) allocate( &
            lhs%coefficient_ranks, source=rhs%coefficient_ranks)
        if (allocated(rhs%coefficient_offsets)) allocate( &
            lhs%coefficient_offsets, source=rhs%coefficient_offsets)
        if (allocated(rhs%coefficient_component_names)) allocate( &
            lhs%coefficient_component_names, source= &
            rhs%coefficient_component_names)
        if (allocated(rhs%coefficient_values)) allocate( &
            lhs%coefficient_values, source=rhs%coefficient_values)
        if (allocated(rhs%profile_coordinates)) allocate( &
            lhs%profile_coordinates, source=rhs%profile_coordinates)
        if (allocated(rhs%profile_names)) allocate( &
            lhs%profile_names, source=rhs%profile_names)
        if (allocated(rhs%profile_values)) allocate( &
            lhs%profile_values, source=rhs%profile_values)
        if (allocated(rhs%boundary_offsets)) allocate( &
            lhs%boundary_offsets, source=rhs%boundary_offsets)
        if (allocated(rhs%boundary_names)) allocate( &
            lhs%boundary_names, source=rhs%boundary_names)
        if (allocated(rhs%boundary_coordinates)) allocate( &
            lhs%boundary_coordinates, source=rhs%boundary_coordinates)
    end subroutine assign_equilibrium_interchange

    logical function input_shapes_are_compatible( &
            mapped, physical, coefficient_names, coefficients, profile_x, &
            profile_names, profiles, offsets, boundary_names, boundary) &
            result(compatible)
        real(dp), intent(in) :: mapped(:, :), physical(:, :)
        character(len=*), intent(in) :: coefficient_names(:), profile_names(:)
        real(dp), intent(in) :: coefficients(:, :), profile_x(:), profiles(:, :)
        integer, intent(in) :: offsets(:)
        character(len=*), intent(in) :: boundary_names(:)
        real(dp), intent(in) :: boundary(:, :)

        compatible = .false.
        if (size(mapped, 1) < 1 .or. size(mapped, 1) > 3 .or. &
            any(shape(physical) /= shape(mapped)) .or. size(mapped, 2) < 1) return
        if (size(coefficients, 2) /= size(mapped, 2)) return
        if (size(profile_x) < 1 .or. size(profiles, 1) /= size(profile_names) .or. &
            size(profiles, 2) /= size(profile_x)) return
        if (size(offsets) < 2 .or. size(boundary_names) /= size(offsets) - 1 .or. &
            size(boundary, 1) /= size(mapped, 1) .or. size(boundary, 2) < 1) return
        compatible = .true.
    end function input_shapes_are_compatible

    logical function finite_normalization(normalization) result(valid)
        type(equilibrium_normalization_t), intent(in) :: normalization
        real(dp) :: scales(4)

        scales = [normalization%length_scale, &
            normalization%magnetic_field_scale, normalization%pressure_scale, &
            normalization%current_scale]
        valid = len_trim(normalization%length_unit) > 0 .and. &
            len_trim(normalization%magnetic_field_unit) > 0 .and. &
            len_trim(normalization%pressure_unit) > 0 .and. &
            len_trim(normalization%current_unit) > 0 .and. &
            all(ieee_is_finite(scales)) .and. all(scales > 0.0_dp)
    end function finite_normalization

    logical function labels_are_unique_and_nonempty(labels) result(valid)
        character(len=*), intent(in) :: labels(:)
        integer :: first, second

        valid = .true.
        do first = 1, size(labels)
            if (len_trim(labels(first)) == 0) then
                valid = .false.
                return
            end if
            do second = first + 1, size(labels)
                if (trim(labels(first)) == trim(labels(second))) then
                    valid = .false.
                    return
                end if
            end do
        end do
    end function labels_are_unique_and_nonempty

    logical function labels_are_nonempty(labels) result(valid)
        character(len=*), intent(in) :: labels(:)
        integer :: label

        valid = .true.
        do label = 1, size(labels)
            if (len_trim(labels(label)) == 0) then
                valid = .false.
                return
            end if
        end do
    end function labels_are_nonempty

    integer function expected_component_count(rank, spatial_dimension) result(count)
        integer, intent(in) :: rank, spatial_dimension

        select case (rank)
        case (0)
            count = 1
        case (1)
            count = spatial_dimension
        case (2)
            count = spatial_dimension*spatial_dimension
        case default
            count = 0
        end select
    end function expected_component_count

    subroutine clear_equilibrium_interchange(interchange)
        type(equilibrium_interchange_t), intent(inout) :: interchange

        if (allocated(interchange%mapped_coordinates)) deallocate( &
            interchange%mapped_coordinates)
        if (allocated(interchange%physical_coordinates)) deallocate( &
            interchange%physical_coordinates)
        if (allocated(interchange%coefficient_names)) deallocate( &
            interchange%coefficient_names)
        if (allocated(interchange%coefficient_ranks)) deallocate( &
            interchange%coefficient_ranks)
        if (allocated(interchange%coefficient_offsets)) deallocate( &
            interchange%coefficient_offsets)
        if (allocated(interchange%coefficient_component_names)) deallocate( &
            interchange%coefficient_component_names)
        if (allocated(interchange%coefficient_values)) deallocate( &
            interchange%coefficient_values)
        if (allocated(interchange%profile_coordinates)) deallocate( &
            interchange%profile_coordinates)
        if (allocated(interchange%profile_names)) deallocate( &
            interchange%profile_names)
        if (allocated(interchange%profile_values)) deallocate( &
            interchange%profile_values)
        if (allocated(interchange%boundary_offsets)) deallocate( &
            interchange%boundary_offsets)
        if (allocated(interchange%boundary_names)) deallocate( &
            interchange%boundary_names)
        if (allocated(interchange%boundary_coordinates)) deallocate( &
            interchange%boundary_coordinates)
        interchange%schema_version = "fortfem-equilibrium-1"
        interchange%producer = ""
        interchange%provenance = ""
        interchange%spatial_dimension = 0
        interchange%sample_count = 0
        interchange%coefficient_count = 0
        interchange%coefficient_component_count = 0
        interchange%profile_count = 0
        interchange%profile_sample_count = 0
        interchange%boundary_count = 0
        interchange%normalization = equilibrium_normalization_t()
    end subroutine clear_equilibrium_interchange

end module fortfem_equilibrium_interchange
