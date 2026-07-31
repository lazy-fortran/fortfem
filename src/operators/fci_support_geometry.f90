module fortfem_fci_support_geometry
    !! Geometry-side factors for the FCI support-operator volumes.
    !!
    !! Field-line tracing supplies forward and backward toroidal flux
    !! expansion integrals.  This module combines them with a plane-cell area
    !! and the local toroidal magnetic field, leaving tracing and mesh lookup
    !! in their own services.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_generated_fci_support_volume, only: &
        generated_fci_staggered_flux_box_volume
    use fortfem_generated_fci_support_volume_jvp, only: &
        generated_fci_staggered_flux_box_volume_jvp
    use fortfem_generated_fci_support_volume_vjp, only: &
        generated_fci_staggered_flux_box_volume_vjp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: compute_fci_staggered_flux_box_volumes
    public :: compute_fci_staggered_flux_box_volumes_jvp
    public :: compute_fci_staggered_flux_box_volumes_vjp

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

        integer :: point_count, point

        call validate_support_volume_factors( &
            forward_flux_expansion, backward_flux_expansion, base_cell_area, &
            toroidal_field, point_count, status)
        if (status%code /= FORTSPARSE_OK) then
            staggered_volumes = 0.0_dp
            return
        end if
        if (size(staggered_volumes) /= point_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI support volumes have an incompatible output")
            staggered_volumes = 0.0_dp
            return
        end if
        do point = 1, point_count
            call generated_fci_staggered_flux_box_volume( &
                forward_flux_expansion(point), backward_flux_expansion(point), &
                base_cell_area(point), toroidal_field(point), &
                staggered_volumes(point))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_fci_staggered_flux_box_volumes

    subroutine compute_fci_staggered_flux_box_volumes_jvp( &
            forward_flux_expansion, backward_flux_expansion, base_cell_area, &
            toroidal_field, forward_flux_expansion_dot, &
            backward_flux_expansion_dot, base_cell_area_dot, toroidal_field_dot, &
            staggered_volumes_dot, status)
        !! Apply the FortSym-generated JVP of the support-volume product.
        real(dp), intent(in) :: forward_flux_expansion(:)
        real(dp), intent(in) :: backward_flux_expansion(:)
        real(dp), intent(in) :: base_cell_area(:)
        real(dp), intent(in) :: toroidal_field(:)
        real(dp), intent(in) :: forward_flux_expansion_dot(:)
        real(dp), intent(in) :: backward_flux_expansion_dot(:)
        real(dp), intent(in) :: base_cell_area_dot(:)
        real(dp), intent(in) :: toroidal_field_dot(:)
        real(dp), intent(out) :: staggered_volumes_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: point_count, point

        call validate_support_volume_factors( &
            forward_flux_expansion, backward_flux_expansion, base_cell_area, &
            toroidal_field, point_count, status)
        staggered_volumes_dot = 0.0_dp
        if (status%code /= FORTSPARSE_OK) return
        if (size(staggered_volumes_dot) /= point_count .or. &
            size(forward_flux_expansion_dot) /= point_count .or. &
            size(backward_flux_expansion_dot) /= point_count .or. &
            size(base_cell_area_dot) /= point_count .or. &
            size(toroidal_field_dot) /= point_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI support-volume JVP has incompatible tangents")
            return
        end if
        if (any(.not. ieee_is_finite(forward_flux_expansion_dot)) .or. &
            any(.not. ieee_is_finite(backward_flux_expansion_dot)) .or. &
            any(.not. ieee_is_finite(base_cell_area_dot)) .or. &
            any(.not. ieee_is_finite(toroidal_field_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI support-volume JVP has non-finite tangents")
            return
        end if
        do point = 1, point_count
            call generated_fci_staggered_flux_box_volume_jvp( &
                forward_flux_expansion(point), backward_flux_expansion(point), &
                base_cell_area(point), toroidal_field(point), &
                forward_flux_expansion_dot(point), &
                backward_flux_expansion_dot(point), base_cell_area_dot(point), &
                toroidal_field_dot(point), staggered_volumes_dot(point))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_fci_staggered_flux_box_volumes_jvp

    subroutine compute_fci_staggered_flux_box_volumes_vjp( &
            forward_flux_expansion, backward_flux_expansion, base_cell_area, &
            toroidal_field, staggered_volume_bar, forward_flux_expansion_bar, &
            backward_flux_expansion_bar, base_cell_area_bar, toroidal_field_bar, &
            status)
        !! Apply the FortSym-generated VJP of the support-volume product.
        real(dp), intent(in) :: forward_flux_expansion(:)
        real(dp), intent(in) :: backward_flux_expansion(:)
        real(dp), intent(in) :: base_cell_area(:)
        real(dp), intent(in) :: toroidal_field(:)
        real(dp), intent(in) :: staggered_volume_bar(:)
        real(dp), intent(out) :: forward_flux_expansion_bar(:)
        real(dp), intent(out) :: backward_flux_expansion_bar(:)
        real(dp), intent(out) :: base_cell_area_bar(:)
        real(dp), intent(out) :: toroidal_field_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: point_count, point

        forward_flux_expansion_bar = 0.0_dp
        backward_flux_expansion_bar = 0.0_dp
        base_cell_area_bar = 0.0_dp
        toroidal_field_bar = 0.0_dp
        call validate_support_volume_factors( &
            forward_flux_expansion, backward_flux_expansion, base_cell_area, &
            toroidal_field, point_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(staggered_volume_bar) /= point_count .or. &
            size(forward_flux_expansion_bar) /= point_count .or. &
            size(backward_flux_expansion_bar) /= point_count .or. &
            size(base_cell_area_bar) /= point_count .or. &
            size(toroidal_field_bar) /= point_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI support-volume VJP has incompatible cotangents")
            return
        end if
        if (any(.not. ieee_is_finite(staggered_volume_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI support-volume VJP has non-finite cotangents")
            return
        end if
        do point = 1, point_count
            call generated_fci_staggered_flux_box_volume_vjp( &
                forward_flux_expansion(point), backward_flux_expansion(point), &
                base_cell_area(point), toroidal_field(point), &
                staggered_volume_bar(point), forward_flux_expansion_bar(point), &
                backward_flux_expansion_bar(point), base_cell_area_bar(point), &
                toroidal_field_bar(point))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_fci_staggered_flux_box_volumes_vjp

    subroutine validate_support_volume_factors( &
            forward_flux_expansion, backward_flux_expansion, base_cell_area, &
            toroidal_field, point_count, status)
        real(dp), intent(in) :: forward_flux_expansion(:)
        real(dp), intent(in) :: backward_flux_expansion(:)
        real(dp), intent(in) :: base_cell_area(:)
        real(dp), intent(in) :: toroidal_field(:)
        integer, intent(out) :: point_count
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI support-volume factors have incompatible arrays")
        point_count = size(forward_flux_expansion)
        if (point_count < 1) return
        if (size(backward_flux_expansion) /= point_count .or. &
            size(base_cell_area) /= point_count .or. &
            size(toroidal_field) /= point_count) return
        if (any(.not. ieee_is_finite(forward_flux_expansion)) .or. &
            any(.not. ieee_is_finite(backward_flux_expansion)) .or. &
            any(.not. ieee_is_finite(base_cell_area)) .or. &
            any(.not. ieee_is_finite(toroidal_field))) return
        if (any(forward_flux_expansion <= 0.0_dp) .or. &
            any(backward_flux_expansion <= 0.0_dp) .or. &
            any(base_cell_area <= 0.0_dp) .or. &
            any(toroidal_field <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_support_volume_factors

end module fortfem_fci_support_geometry
