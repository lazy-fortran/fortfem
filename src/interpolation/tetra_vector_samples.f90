module fortfem_tetra_vector_samples
    !! Topology-aware samples of a physical vector field on one tetrahedron.
    use fortfem_kinds, only: dp
    implicit none
    private

    type :: tetra_vector_samples_t
        real(dp), allocatable :: edge(:, :, :)
        real(dp), allocatable :: face(:, :, :)
        real(dp), allocatable :: volume(:, :)
    end type tetra_vector_samples_t

    type :: tetra_vector_sample_gradients_t
        !! First index is field component; second is physical coordinate.
        real(dp), allocatable :: edge(:, :, :, :)
        real(dp), allocatable :: face(:, :, :, :)
        real(dp), allocatable :: volume(:, :, :)
    end type tetra_vector_sample_gradients_t

    interface assignment(=)
        module procedure assign_tetra_vector_samples
        module procedure assign_tetra_vector_sample_gradients
    end interface

    public :: assignment(=)
    public :: tetra_vector_sample_gradients_t
    public :: tetra_vector_samples_t
    public :: zero_tetra_vector_samples_like

contains

    subroutine zero_tetra_vector_samples_like(template, samples)
        type(tetra_vector_samples_t), intent(in) :: template
        type(tetra_vector_samples_t), intent(out) :: samples

        if (allocated(template%edge)) then
            allocate(samples%edge, mold=template%edge)
            samples%edge = 0.0_dp
        end if
        if (allocated(template%face)) then
            allocate(samples%face, mold=template%face)
            samples%face = 0.0_dp
        end if
        if (allocated(template%volume)) then
            allocate(samples%volume, mold=template%volume)
            samples%volume = 0.0_dp
        end if
    end subroutine zero_tetra_vector_samples_like

    subroutine assign_tetra_vector_samples(left, right)
        type(tetra_vector_samples_t), intent(out) :: left
        type(tetra_vector_samples_t), intent(in) :: right

        if (allocated(right%edge)) allocate(left%edge, source=right%edge)
        if (allocated(right%face)) allocate(left%face, source=right%face)
        if (allocated(right%volume)) allocate(left%volume, source=right%volume)
    end subroutine assign_tetra_vector_samples

    subroutine assign_tetra_vector_sample_gradients(left, right)
        type(tetra_vector_sample_gradients_t), intent(out) :: left
        type(tetra_vector_sample_gradients_t), intent(in) :: right

        if (allocated(right%edge)) allocate(left%edge, source=right%edge)
        if (allocated(right%face)) allocate(left%face, source=right%face)
        if (allocated(right%volume)) then
            allocate(left%volume, source=right%volume)
        end if
    end subroutine assign_tetra_vector_sample_gradients

end module fortfem_tetra_vector_samples
