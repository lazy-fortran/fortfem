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
    public :: compare_complex_interchange_samples_jvp
    public :: compare_complex_interchange_samples_vjp
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

    subroutine compare_complex_interchange_samples_jvp( &
            reference, candidate, coordinate_tolerance, reference_values_dot, &
            candidate_values_dot, weights_dot, absolute_error_dot, &
            relative_error_dot, status)
        !! Differentiate complex weighted L2 and relative errors.
        !!
        !! The derivative is real-valued and uses
        !! `Re(conjg(value)*value_dot)`.  Coordinates are fixed and reference
        !! and candidate share one validated weight vector.  The maximum-norm
        !! diagnostic is intentionally excluded because its active component
        !! is nonsmooth at ties.
        type(complex_interchange_sample_set_t), intent(in) :: reference, candidate
        real(dp), intent(in) :: coordinate_tolerance
        complex(dp), intent(in) :: reference_values_dot(:, :)
        complex(dp), intent(in) :: candidate_values_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: absolute_error_dot, relative_error_dot
        integer, intent(out) :: status

        real(dp) :: absolute_error, relative_error, maximum_error
        real(dp) :: reference_norm, squared_error_dot, reference_norm_dot
        complex(dp), allocatable :: difference(:, :), difference_dot(:, :)

        absolute_error_dot = 0.0_dp
        relative_error_dot = 0.0_dp
        status = 1
        call compare_complex_interchange_samples( &
            reference, candidate, coordinate_tolerance, absolute_error, &
            relative_error, maximum_error, status)
        if (status /= 0) return
        status = 1
        if (any(shape(reference_values_dot) /= shape(reference%values)) .or. &
            any(shape(candidate_values_dot) /= shape(candidate%values)) .or. &
            size(weights_dot) /= reference%sample_count) return
        if (.not. finite_complex(reference_values_dot) .or. &
            .not. finite_complex(candidate_values_dot) .or. &
            .not. all(ieee_is_finite(weights_dot))) return
        if (absolute_error <= tiny(1.0_dp)) return
        allocate(difference(reference%component_count, reference%sample_count), &
            difference_dot(reference%component_count, reference%sample_count))
        difference = candidate%values - reference%values
        difference_dot = candidate_values_dot - reference_values_dot
        squared_error_dot = sum( &
            spread(weights_dot, 1, reference%component_count)*abs(difference)**2 + &
            2.0_dp*spread(reference%weights, 1, reference%component_count)* &
            real(conjg(difference)*difference_dot, dp))
        absolute_error_dot = squared_error_dot/(2.0_dp*absolute_error)

        reference_norm = sqrt(sum( &
            spread(reference%weights, 1, reference%component_count)* &
            abs(reference%values)**2))
        if (reference_norm <= tiny(1.0_dp)) return
        reference_norm_dot = sum( &
            spread(weights_dot, 1, reference%component_count)* &
            abs(reference%values)**2 + 2.0_dp* &
            spread(reference%weights, 1, reference%component_count)* &
            real(conjg(reference%values)*reference_values_dot, dp))/ &
            (2.0_dp*reference_norm)
        relative_error_dot = absolute_error_dot/reference_norm - &
            absolute_error*reference_norm_dot/reference_norm**2
        status = 0
    end subroutine compare_complex_interchange_samples_jvp

    subroutine compare_complex_interchange_samples_vjp( &
            reference, candidate, coordinate_tolerance, absolute_error_bar, &
            relative_error_bar, reference_values_bar, candidate_values_bar, &
            weights_bar, status)
        !! Apply the real-part complex VJP of weighted error metrics.
        type(complex_interchange_sample_set_t), intent(in) :: reference, candidate
        real(dp), intent(in) :: coordinate_tolerance
        real(dp), intent(in) :: absolute_error_bar, relative_error_bar
        complex(dp), intent(out) :: reference_values_bar(:, :)
        complex(dp), intent(out) :: candidate_values_bar(:, :)
        real(dp), intent(out) :: weights_bar(:)
        integer, intent(out) :: status

        real(dp) :: absolute_error, relative_error, maximum_error
        real(dp) :: reference_norm, absolute_partial, reference_partial
        complex(dp), allocatable :: difference(:, :)
        real(dp), allocatable :: weights_matrix(:, :)

        reference_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        candidate_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        weights_bar = 0.0_dp
        status = 1
        call compare_complex_interchange_samples( &
            reference, candidate, coordinate_tolerance, absolute_error, &
            relative_error, maximum_error, status)
        if (status /= 0) return
        status = 1
        if (any(shape(reference_values_bar) /= shape(reference%values)) .or. &
            any(shape(candidate_values_bar) /= shape(candidate%values)) .or. &
            size(weights_bar) /= reference%sample_count) return
        if (.not. ieee_is_finite(absolute_error_bar) .or. &
            .not. ieee_is_finite(relative_error_bar)) return
        if (absolute_error <= tiny(1.0_dp)) return
        allocate(difference(reference%component_count, reference%sample_count), &
            weights_matrix(reference%component_count, reference%sample_count))
        difference = candidate%values - reference%values
        weights_matrix = spread(reference%weights, 1, reference%component_count)
        reference_norm = sqrt(sum(weights_matrix*abs(reference%values)**2))
        if (reference_norm <= tiny(1.0_dp)) return
        absolute_partial = absolute_error_bar + relative_error_bar/reference_norm
        reference_partial = -relative_error_bar*absolute_error/reference_norm**2
        candidate_values_bar = absolute_partial*weights_matrix*difference/absolute_error
        reference_values_bar = -candidate_values_bar + reference_partial* &
            weights_matrix*reference%values/reference_norm
        weights_bar = absolute_partial*0.5_dp*sum(abs(difference)**2, dim=1)/ &
            absolute_error
        weights_bar = weights_bar + reference_partial*0.5_dp* &
            sum(abs(reference%values)**2, dim=1)/reference_norm
        status = 0
    end subroutine compare_complex_interchange_samples_vjp

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
