module fortfem_equilibrium_sample_adapter
    !! Select already-sampled equilibrium coefficients for code comparisons.
    !!
    !! This is deliberately not an interpolator or an equilibrium reader.  An
    !! external adapter owns resampling and coordinate conventions; this
    !! routine only projects the physical samples stored in the neutral
    !! equilibrium interchange record onto the common weighted sample-set
    !! schema.
    use fortfem_equilibrium_interchange, only: &
        equilibrium_interchange_t, validate_equilibrium_interchange
    use fortfem_interchange_samples, only: &
        initialize_interchange_samples, interchange_sample_set_t
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: build_equilibrium_interchange_sample_set

contains

    subroutine build_equilibrium_interchange_sample_set( &
            interchange, component_indices, weights, producer, provenance, &
            samples, status)
        !! Project selected physical coefficient components onto a common grid.
        !!
        !! `component_indices` uses the component ordering in
        !! `interchange%coefficient_values`, not field-local component names.
        !! The caller supplies positive quadrature or comparison weights and
        !! explicit provenance for the resulting sample set.
        type(equilibrium_interchange_t), intent(in) :: interchange
        integer, intent(in) :: component_indices(:)
        real(dp), intent(in) :: weights(:)
        character(*), intent(in) :: producer, provenance
        type(interchange_sample_set_t), intent(inout) :: samples
        integer, intent(out) :: status

        integer :: validation_status, component
        real(dp), allocatable :: selected_values(:, :)

        status = 1
        if (.not. validate_equilibrium_interchange( &
            interchange, validation_status)) return
        if (size(component_indices) < 1 .or. &
            size(weights) /= interchange%sample_count) return
        if (any(component_indices < 1) .or. &
            any(component_indices > interchange%coefficient_component_count)) return
        do component = 1, size(component_indices) - 1
            if (any(component_indices(component + 1:) == &
                component_indices(component))) return
        end do

        allocate(selected_values(size(component_indices), interchange%sample_count))
        do component = 1, size(component_indices)
            selected_values(component, :) = &
                interchange%coefficient_values(component_indices(component), :)
        end do
        call initialize_interchange_samples( &
            samples, interchange%physical_coordinates, selected_values, weights, &
            producer, provenance, status)
    end subroutine build_equilibrium_interchange_sample_set

end module fortfem_equilibrium_sample_adapter
