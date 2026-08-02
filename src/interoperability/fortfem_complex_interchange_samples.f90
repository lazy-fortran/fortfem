module fortfem_complex_interchange_samples
    !! Neutral complex physical sampling contract for frequency-domain codes.
    !!
    !! Coordinates and positive quadrature weights are real, while sampled
    !! scalar, vector, or tensor components may be complex.  The record does
    !! not define a PDE, coordinate convention, equilibrium reader, or solver;
    !! it is only a common physical-grid boundary for FEM/BEM/DtN/PML and
    !! external ideal/resistive-response comparisons.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    type, public :: complex_interchange_sample_set_t
        character(len=32) :: schema_version = "fortfem-complex-samples-1"
        character(len=64) :: producer = ""
        character(len=128) :: provenance = ""
        integer :: spatial_dimension = 0
        integer :: component_count = 0
        integer :: sample_count = 0
        real(dp), allocatable :: coordinates(:, :)
        complex(dp), allocatable :: values(:, :)
        real(dp), allocatable :: weights(:)
    contains
        procedure, private :: assign_complex_interchange_samples
        generic :: assignment(=) => assign_complex_interchange_samples
    end type complex_interchange_sample_set_t

    public :: compare_complex_interchange_samples
    public :: initialize_complex_interchange_samples
    public :: validate_complex_interchange_samples

contains

    subroutine initialize_complex_interchange_samples( &
            samples, coordinates, values, weights, producer, provenance, status)
        type(complex_interchange_sample_set_t), intent(inout) :: samples
        real(dp), intent(in) :: coordinates(:, :), weights(:)
        complex(dp), intent(in) :: values(:, :)
        character(len=*), intent(in) :: producer, provenance
        integer, intent(out) :: status

        call clear_complex_interchange_samples(samples)
        status = 1
        if (size(coordinates, 1) < 1 .or. size(coordinates, 1) > 3 .or. &
            size(coordinates, 2) < 1 .or. size(values, 2) /= size(coordinates, 2) .or. &
            size(values, 1) < 1 .or. size(weights) /= size(coordinates, 2) .or. &
            len_trim(producer) == 0 .or. len_trim(provenance) == 0) return
        if (.not. all(ieee_is_finite(coordinates)) .or. &
            .not. finite_complex(values) .or. .not. all(ieee_is_finite(weights)) .or. &
            any(weights <= 0.0_dp)) return

        samples%producer = producer
        samples%provenance = provenance
        samples%spatial_dimension = size(coordinates, 1)
        samples%component_count = size(values, 1)
        samples%sample_count = size(coordinates, 2)
        allocate(samples%coordinates, source=coordinates)
        allocate(samples%values, source=values)
        allocate(samples%weights, source=weights)
        if (.not. validate_complex_interchange_samples(samples, status)) then
            call clear_complex_interchange_samples(samples)
            return
        end if
    end subroutine initialize_complex_interchange_samples

    logical function validate_complex_interchange_samples(samples, status) &
            result(valid)
        type(complex_interchange_sample_set_t), intent(in) :: samples
        integer, intent(out) :: status

        valid = .false.
        status = 1
        if (samples%schema_version /= "fortfem-complex-samples-1" .or. &
            samples%spatial_dimension < 1 .or. samples%spatial_dimension > 3 .or. &
            samples%component_count < 1 .or. samples%sample_count < 1 .or. &
            len_trim(samples%producer) == 0 .or. len_trim(samples%provenance) == 0) return
        if (.not. allocated(samples%coordinates) .or. &
            .not. allocated(samples%values) .or. .not. allocated(samples%weights)) return
        if (any(shape(samples%coordinates) /= [ &
            samples%spatial_dimension, samples%sample_count]) .or. &
            any(shape(samples%values) /= [ &
            samples%component_count, samples%sample_count]) .or. &
            size(samples%weights) /= samples%sample_count) return
        if (.not. all(ieee_is_finite(samples%coordinates)) .or. &
            .not. finite_complex(samples%values) .or. &
            .not. all(ieee_is_finite(samples%weights)) .or. &
            any(samples%weights <= 0.0_dp)) return
        valid = .true.
        status = 0
    end function validate_complex_interchange_samples

    subroutine compare_complex_interchange_samples( &
            reference, candidate, coordinate_tolerance, absolute_error, &
            relative_error, maximum_error, status)
        type(complex_interchange_sample_set_t), intent(in) :: reference, candidate
        real(dp), intent(in) :: coordinate_tolerance
        real(dp), intent(out) :: absolute_error, relative_error, maximum_error
        integer, intent(out) :: status

        complex(dp), allocatable :: difference(:, :)
        real(dp) :: reference_norm, squared_error
        integer :: validation_status

        absolute_error = huge(1.0_dp)
        relative_error = huge(1.0_dp)
        maximum_error = huge(1.0_dp)
        status = 1
        if (.not. validate_complex_interchange_samples(reference, validation_status) .or. &
            .not. validate_complex_interchange_samples(candidate, validation_status)) return
        if (.not. ieee_is_finite(coordinate_tolerance) .or. coordinate_tolerance < 0.0_dp .or. &
            reference%spatial_dimension /= candidate%spatial_dimension .or. &
            reference%component_count /= candidate%component_count .or. &
            reference%sample_count /= candidate%sample_count) return
        if (maxval(abs(reference%coordinates - candidate%coordinates)) > &
            coordinate_tolerance .or. maxval(abs(reference%weights - candidate%weights)) > &
            coordinate_tolerance) return

        allocate(difference(reference%component_count, reference%sample_count))
        difference = candidate%values - reference%values
        squared_error = sum(spread(reference%weights, 1, reference%component_count)* &
            abs(difference)**2)
        reference_norm = sqrt(sum(spread(reference%weights, 1, reference%component_count)* &
            abs(reference%values)**2))
        absolute_error = sqrt(squared_error)
        relative_error = absolute_error/max(reference_norm, tiny(1.0_dp))
        maximum_error = maxval(abs(difference))
        status = 0
    end subroutine compare_complex_interchange_samples

    subroutine assign_complex_interchange_samples(lhs, rhs)
        class(complex_interchange_sample_set_t), intent(out) :: lhs
        type(complex_interchange_sample_set_t), intent(in) :: rhs

        call clear_complex_interchange_samples(lhs)
        lhs%schema_version = rhs%schema_version
        lhs%producer = rhs%producer
        lhs%provenance = rhs%provenance
        lhs%spatial_dimension = rhs%spatial_dimension
        lhs%component_count = rhs%component_count
        lhs%sample_count = rhs%sample_count
        if (allocated(rhs%coordinates)) allocate(lhs%coordinates, source=rhs%coordinates)
        if (allocated(rhs%values)) allocate(lhs%values, source=rhs%values)
        if (allocated(rhs%weights)) allocate(lhs%weights, source=rhs%weights)
    end subroutine assign_complex_interchange_samples

    subroutine clear_complex_interchange_samples(samples)
        type(complex_interchange_sample_set_t), intent(inout) :: samples

        if (allocated(samples%coordinates)) deallocate(samples%coordinates)
        if (allocated(samples%values)) deallocate(samples%values)
        if (allocated(samples%weights)) deallocate(samples%weights)
        samples%schema_version = "fortfem-complex-samples-1"
        samples%producer = ""
        samples%provenance = ""
        samples%spatial_dimension = 0
        samples%component_count = 0
        samples%sample_count = 0
    end subroutine clear_complex_interchange_samples

    logical function finite_complex(value) result(valid)
        complex(dp), intent(in) :: value(:, :)

        valid = all(ieee_is_finite(real(value, dp))) .and. &
            all(ieee_is_finite(aimag(value)))
    end function finite_complex

end module fortfem_complex_interchange_samples
