module fortfem_fci_support_geometry
    !! Geometry-side factors for the FCI support-operator volumes.
    !!
    !! Field-line tracing supplies forward and backward toroidal flux
    !! expansion integrals.  This module combines them with a plane-cell area
    !! and the local toroidal magnetic field, leaving tracing and mesh lookup
    !! in their own services.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: compute_fci_staggered_flux_box_volumes

contains

    subroutine compute_fci_staggered_flux_box_volumes( &
            forward_flux_expansion, backward_flux_expansion, base_cell_area, &
            toroidal_field, staggered_volumes, status)
        !! Compute positive staggered flux-box volumes.
        !!
        !! The algebraic contract is
        !!
        !!   omega = (v_forward + v_backward) A_plane B_toroidal.
        !!
        !! The factors are normally obtained by integrating a field-line
        !! tracer on either side of a staggered plane.  All arrays are
        !! pointwise and must have the same positive, finite length.
        real(dp), intent(in) :: forward_flux_expansion(:)
        real(dp), intent(in) :: backward_flux_expansion(:)
        real(dp), intent(in) :: base_cell_area(:)
        real(dp), intent(in) :: toroidal_field(:)
        real(dp), intent(out) :: staggered_volumes(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: point_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI support-volume factors have incompatible arrays")
        staggered_volumes = 0.0_dp
        point_count = size(forward_flux_expansion)
        if (point_count < 1) return
        if (size(backward_flux_expansion) /= point_count .or. &
            size(base_cell_area) /= point_count .or. &
            size(toroidal_field) /= point_count .or. &
            size(staggered_volumes) /= point_count) return
        if (any(.not. ieee_is_finite(forward_flux_expansion)) .or. &
            any(.not. ieee_is_finite(backward_flux_expansion)) .or. &
            any(.not. ieee_is_finite(base_cell_area)) .or. &
            any(.not. ieee_is_finite(toroidal_field))) return
        if (any(forward_flux_expansion <= 0.0_dp) .or. &
            any(backward_flux_expansion <= 0.0_dp) .or. &
            any(base_cell_area <= 0.0_dp) .or. &
            any(toroidal_field <= 0.0_dp)) return
        staggered_volumes = (forward_flux_expansion + &
            backward_flux_expansion)*base_cell_area*toroidal_field
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_fci_staggered_flux_box_volumes

end module fortfem_fci_support_geometry
