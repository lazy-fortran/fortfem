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
    use fortfem_generated_fci_quadrilateral_area, only: &
        generated_fci_quadrilateral_cell_area
    use fortfem_generated_fci_quadrilateral_area_jvp, only: &
        generated_fci_quadrilateral_cell_area_jvp
    use fortfem_generated_fci_quadrilateral_area_vjp, only: &
        generated_fci_quadrilateral_cell_area_vjp
    use fortfem_generated_fci_curved_quadrilateral_area, only: &
        generated_fci_curved_quadrilateral_cell_area
    use fortfem_generated_fci_curved_quadrilateral_area_jvp, only: &
        generated_fci_curved_quadrilateral_cell_area_jvp
    use fortfem_generated_fci_curved_quadrilateral_area_vjp, only: &
        generated_fci_curved_quadrilateral_cell_area_vjp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: compute_fci_staggered_flux_box_volumes
    public :: compute_fci_staggered_flux_box_volumes_jvp
    public :: compute_fci_staggered_flux_box_volumes_vjp
    public :: compute_fci_quadrilateral_cell_areas_2d
    public :: compute_fci_quadrilateral_cell_areas_2d_jvp
    public :: compute_fci_quadrilateral_cell_areas_2d_vjp
    public :: compute_fci_curved_quadrilateral_cell_areas_2d
    public :: compute_fci_curved_quadrilateral_cell_areas_2d_jvp
    public :: compute_fci_curved_quadrilateral_cell_areas_2d_vjp

contains

    subroutine compute_fci_quadrilateral_cell_areas_2d( &
            cell_vertices, areas, status)
        !! Compute signed-positive areas for counter-clockwise quadrilaterals.
        !!
        !! `cell_vertices(:, :, cell)` stores the four vertices in boundary
        !! order.  Clockwise, self-intersecting, and degenerate cells are
        !! rejected so the returned measure is positive and differentiable on
        !! a fixed orientation/topology.
        real(dp), intent(in) :: cell_vertices(:, :, :)
        real(dp), intent(out) :: areas(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: cell_count, cell

        areas = 0.0_dp
        call validate_quadrilateral_vertices(cell_vertices, cell_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(areas) /= cell_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadrilateral areas have an incompatible output")
            return
        end if
        do cell = 1, cell_count
            call generated_fci_quadrilateral_cell_area( &
                cell_vertices(1, 1, cell), cell_vertices(2, 1, cell), &
                cell_vertices(1, 2, cell), cell_vertices(2, 2, cell), &
                cell_vertices(1, 3, cell), cell_vertices(2, 3, cell), &
                cell_vertices(1, 4, cell), cell_vertices(2, 4, cell), &
                areas(cell))
            if (.not. ieee_is_finite(areas(cell)) .or. areas(cell) <= 0.0_dp) then
                areas = 0.0_dp
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI quadrilateral cells must be counter-clockwise and "// &
                    "nondegenerate")
                return
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_fci_quadrilateral_cell_areas_2d

    subroutine compute_fci_quadrilateral_cell_areas_2d_jvp( &
            cell_vertices, cell_vertices_dot, areas_dot, status)
        !! Apply the fixed-topology JVP of quadrilateral cell areas.
        real(dp), intent(in) :: cell_vertices(:, :, :)
        real(dp), intent(in) :: cell_vertices_dot(:, :, :)
        real(dp), intent(out) :: areas_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: cell_count, cell
        real(dp) :: area

        areas_dot = 0.0_dp
        call validate_quadrilateral_vertices(cell_vertices, cell_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(cell_vertices_dot, 1) /= 2 .or. &
            size(cell_vertices_dot, 2) /= 4 .or. &
            size(cell_vertices_dot, 3) /= cell_count .or. &
            size(areas_dot) /= cell_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadrilateral area JVP has incompatible arrays")
            return
        end if
        if (any(.not. ieee_is_finite(cell_vertices_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadrilateral area JVP has non-finite tangents")
            return
        end if
        do cell = 1, cell_count
            call generated_fci_quadrilateral_cell_area( &
                cell_vertices(1, 1, cell), cell_vertices(2, 1, cell), &
                cell_vertices(1, 2, cell), cell_vertices(2, 2, cell), &
                cell_vertices(1, 3, cell), cell_vertices(2, 3, cell), &
                cell_vertices(1, 4, cell), cell_vertices(2, 4, cell), area)
            if (.not. ieee_is_finite(area) .or. area <= 0.0_dp) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI quadrilateral area JVP crosses an orientation event")
                return
            end if
            call generated_fci_quadrilateral_cell_area_jvp( &
                cell_vertices(1, 1, cell), cell_vertices(2, 1, cell), &
                cell_vertices(1, 2, cell), cell_vertices(2, 2, cell), &
                cell_vertices(1, 3, cell), cell_vertices(2, 3, cell), &
                cell_vertices(1, 4, cell), cell_vertices(2, 4, cell), &
                cell_vertices_dot(1, 1, cell), cell_vertices_dot(2, 1, cell), &
                cell_vertices_dot(1, 2, cell), cell_vertices_dot(2, 2, cell), &
                cell_vertices_dot(1, 3, cell), cell_vertices_dot(2, 3, cell), &
                cell_vertices_dot(1, 4, cell), cell_vertices_dot(2, 4, cell), &
                areas_dot(cell))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_fci_quadrilateral_cell_areas_2d_jvp

    subroutine compute_fci_quadrilateral_cell_areas_2d_vjp( &
            cell_vertices, areas_bar, cell_vertices_bar, status)
        !! Apply the real VJP of fixed-topology quadrilateral areas.
        real(dp), intent(in) :: cell_vertices(:, :, :)
        real(dp), intent(in) :: areas_bar(:)
        real(dp), intent(out) :: cell_vertices_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: cell_count, cell
        real(dp) :: area

        cell_vertices_bar = 0.0_dp
        call validate_quadrilateral_vertices(cell_vertices, cell_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(areas_bar) /= cell_count .or. &
            size(cell_vertices_bar, 1) /= 2 .or. &
            size(cell_vertices_bar, 2) /= 4 .or. &
            size(cell_vertices_bar, 3) /= cell_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadrilateral area VJP has incompatible arrays")
            return
        end if
        if (any(.not. ieee_is_finite(areas_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadrilateral area VJP has non-finite cotangents")
            return
        end if
        do cell = 1, cell_count
            call generated_fci_quadrilateral_cell_area( &
                cell_vertices(1, 1, cell), cell_vertices(2, 1, cell), &
                cell_vertices(1, 2, cell), cell_vertices(2, 2, cell), &
                cell_vertices(1, 3, cell), cell_vertices(2, 3, cell), &
                cell_vertices(1, 4, cell), cell_vertices(2, 4, cell), area)
            if (.not. ieee_is_finite(area) .or. area <= 0.0_dp) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI quadrilateral area VJP crosses an orientation event")
                return
            end if
            call generated_fci_quadrilateral_cell_area_vjp( &
                cell_vertices(1, 1, cell), cell_vertices(2, 1, cell), &
                cell_vertices(1, 2, cell), cell_vertices(2, 2, cell), &
                cell_vertices(1, 3, cell), cell_vertices(2, 3, cell), &
                cell_vertices(1, 4, cell), cell_vertices(2, 4, cell), &
                areas_bar(cell), cell_vertices_bar(1, 1, cell), &
                cell_vertices_bar(2, 1, cell), cell_vertices_bar(1, 2, cell), &
                cell_vertices_bar(2, 2, cell), cell_vertices_bar(1, 3, cell), &
                cell_vertices_bar(2, 3, cell), cell_vertices_bar(1, 4, cell), &
                cell_vertices_bar(2, 4, cell))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_fci_quadrilateral_cell_areas_2d_vjp

    subroutine compute_fci_curved_quadrilateral_cell_areas_2d( &
            cell_vertices, edge_controls, areas, status)
        !! Compute positive Green areas for quadratic Bezier-edge cells.
        !!
        !! `edge_controls(:, edge, cell)` is the control point for the
        !! boundary edge from vertex `edge` to the next cyclic vertex.  The
        !! area is differentiable for fixed orientation and connectivity.
        real(dp), intent(in) :: cell_vertices(:, :, :)
        real(dp), intent(in) :: edge_controls(:, :, :)
        real(dp), intent(out) :: areas(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: cell_count, cell

        areas = 0.0_dp
        call validate_curved_quadrilateral_vertices( &
            cell_vertices, edge_controls, cell_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(areas) /= cell_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI curved quadrilateral areas have incompatible output")
            return
        end if
        do cell = 1, cell_count
            call generated_fci_curved_quadrilateral_cell_area( &
                cell_vertices(1, 1, cell), cell_vertices(2, 1, cell), &
                cell_vertices(1, 2, cell), cell_vertices(2, 2, cell), &
                cell_vertices(1, 3, cell), cell_vertices(2, 3, cell), &
                cell_vertices(1, 4, cell), cell_vertices(2, 4, cell), &
                edge_controls(1, 1, cell), edge_controls(2, 1, cell), &
                edge_controls(1, 2, cell), edge_controls(2, 2, cell), &
                edge_controls(1, 3, cell), edge_controls(2, 3, cell), &
                edge_controls(1, 4, cell), edge_controls(2, 4, cell), &
                areas(cell))
            if (.not. ieee_is_finite(areas(cell)) .or. areas(cell) <= 0.0_dp) then
                areas = 0.0_dp
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI curved quadrilateral cells must have positive area")
                return
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_fci_curved_quadrilateral_cell_areas_2d

    subroutine compute_fci_curved_quadrilateral_cell_areas_2d_jvp( &
            cell_vertices, edge_controls, cell_vertices_dot, &
            edge_controls_dot, areas_dot, status)
        !! Apply the fixed-topology JVP of quadratic Bezier-edge cell areas.
        real(dp), intent(in) :: cell_vertices(:, :, :)
        real(dp), intent(in) :: edge_controls(:, :, :)
        real(dp), intent(in) :: cell_vertices_dot(:, :, :)
        real(dp), intent(in) :: edge_controls_dot(:, :, :)
        real(dp), intent(out) :: areas_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: cell_count, cell
        real(dp) :: area

        areas_dot = 0.0_dp
        call validate_curved_quadrilateral_vertices( &
            cell_vertices, edge_controls, cell_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(cell_vertices_dot, 1) /= 2 .or. &
            size(cell_vertices_dot, 2) /= 4 .or. &
            size(cell_vertices_dot, 3) /= cell_count .or. &
            size(edge_controls_dot, 1) /= 2 .or. &
            size(edge_controls_dot, 2) /= 4 .or. &
            size(edge_controls_dot, 3) /= cell_count .or. &
            size(areas_dot) /= cell_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI curved quadrilateral area JVP has incompatible arrays")
            return
        end if
        if (any(.not. ieee_is_finite(cell_vertices_dot)) .or. &
            any(.not. ieee_is_finite(edge_controls_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI curved quadrilateral area JVP has non-finite tangents")
            return
        end if
        do cell = 1, cell_count
            call generated_fci_curved_quadrilateral_cell_area( &
                cell_vertices(1, 1, cell), cell_vertices(2, 1, cell), &
                cell_vertices(1, 2, cell), cell_vertices(2, 2, cell), &
                cell_vertices(1, 3, cell), cell_vertices(2, 3, cell), &
                cell_vertices(1, 4, cell), cell_vertices(2, 4, cell), &
                edge_controls(1, 1, cell), edge_controls(2, 1, cell), &
                edge_controls(1, 2, cell), edge_controls(2, 2, cell), &
                edge_controls(1, 3, cell), edge_controls(2, 3, cell), &
                edge_controls(1, 4, cell), edge_controls(2, 4, cell), area)
            if (.not. ieee_is_finite(area) .or. area <= 0.0_dp) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI curved quadrilateral JVP crosses an area event")
                return
            end if
            call generated_fci_curved_quadrilateral_cell_area_jvp( &
                cell_vertices(1, 1, cell), cell_vertices(2, 1, cell), &
                cell_vertices(1, 2, cell), cell_vertices(2, 2, cell), &
                cell_vertices(1, 3, cell), cell_vertices(2, 3, cell), &
                cell_vertices(1, 4, cell), cell_vertices(2, 4, cell), &
                edge_controls(1, 1, cell), edge_controls(2, 1, cell), &
                edge_controls(1, 2, cell), edge_controls(2, 2, cell), &
                edge_controls(1, 3, cell), edge_controls(2, 3, cell), &
                edge_controls(1, 4, cell), edge_controls(2, 4, cell), &
                cell_vertices_dot(1, 1, cell), cell_vertices_dot(2, 1, cell), &
                cell_vertices_dot(1, 2, cell), cell_vertices_dot(2, 2, cell), &
                cell_vertices_dot(1, 3, cell), cell_vertices_dot(2, 3, cell), &
                cell_vertices_dot(1, 4, cell), cell_vertices_dot(2, 4, cell), &
                edge_controls_dot(1, 1, cell), edge_controls_dot(2, 1, cell), &
                edge_controls_dot(1, 2, cell), edge_controls_dot(2, 2, cell), &
                edge_controls_dot(1, 3, cell), edge_controls_dot(2, 3, cell), &
                edge_controls_dot(1, 4, cell), edge_controls_dot(2, 4, cell), &
                areas_dot(cell))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_fci_curved_quadrilateral_cell_areas_2d_jvp

    subroutine compute_fci_curved_quadrilateral_cell_areas_2d_vjp( &
            cell_vertices, edge_controls, areas_bar, cell_vertices_bar, &
            edge_controls_bar, status)
        !! Apply the real VJP of fixed-topology Bezier-edge cell areas.
        real(dp), intent(in) :: cell_vertices(:, :, :)
        real(dp), intent(in) :: edge_controls(:, :, :)
        real(dp), intent(in) :: areas_bar(:)
        real(dp), intent(out) :: cell_vertices_bar(:, :, :)
        real(dp), intent(out) :: edge_controls_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: cell_count, cell
        real(dp) :: area

        cell_vertices_bar = 0.0_dp
        edge_controls_bar = 0.0_dp
        call validate_curved_quadrilateral_vertices( &
            cell_vertices, edge_controls, cell_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(areas_bar) /= cell_count .or. &
            size(cell_vertices_bar, 1) /= 2 .or. &
            size(cell_vertices_bar, 2) /= 4 .or. &
            size(cell_vertices_bar, 3) /= cell_count .or. &
            size(edge_controls_bar, 1) /= 2 .or. &
            size(edge_controls_bar, 2) /= 4 .or. &
            size(edge_controls_bar, 3) /= cell_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI curved quadrilateral area VJP has incompatible arrays")
            return
        end if
        if (any(.not. ieee_is_finite(areas_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI curved quadrilateral area VJP has non-finite cotangents")
            return
        end if
        do cell = 1, cell_count
            call generated_fci_curved_quadrilateral_cell_area( &
                cell_vertices(1, 1, cell), cell_vertices(2, 1, cell), &
                cell_vertices(1, 2, cell), cell_vertices(2, 2, cell), &
                cell_vertices(1, 3, cell), cell_vertices(2, 3, cell), &
                cell_vertices(1, 4, cell), cell_vertices(2, 4, cell), &
                edge_controls(1, 1, cell), edge_controls(2, 1, cell), &
                edge_controls(1, 2, cell), edge_controls(2, 2, cell), &
                edge_controls(1, 3, cell), edge_controls(2, 3, cell), &
                edge_controls(1, 4, cell), edge_controls(2, 4, cell), area)
            if (.not. ieee_is_finite(area) .or. area <= 0.0_dp) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI curved quadrilateral VJP crosses an area event")
                return
            end if
            call generated_fci_curved_quadrilateral_cell_area_vjp( &
                cell_vertices(1, 1, cell), cell_vertices(2, 1, cell), &
                cell_vertices(1, 2, cell), cell_vertices(2, 2, cell), &
                cell_vertices(1, 3, cell), cell_vertices(2, 3, cell), &
                cell_vertices(1, 4, cell), cell_vertices(2, 4, cell), &
                edge_controls(1, 1, cell), edge_controls(2, 1, cell), &
                edge_controls(1, 2, cell), edge_controls(2, 2, cell), &
                edge_controls(1, 3, cell), edge_controls(2, 3, cell), &
                edge_controls(1, 4, cell), edge_controls(2, 4, cell), &
                areas_bar(cell), cell_vertices_bar(1, 1, cell), &
                cell_vertices_bar(2, 1, cell), cell_vertices_bar(1, 2, cell), &
                cell_vertices_bar(2, 2, cell), cell_vertices_bar(1, 3, cell), &
                cell_vertices_bar(2, 3, cell), cell_vertices_bar(1, 4, cell), &
                cell_vertices_bar(2, 4, cell), edge_controls_bar(1, 1, cell), &
                edge_controls_bar(2, 1, cell), edge_controls_bar(1, 2, cell), &
                edge_controls_bar(2, 2, cell), edge_controls_bar(1, 3, cell), &
                edge_controls_bar(2, 3, cell), edge_controls_bar(1, 4, cell), &
                edge_controls_bar(2, 4, cell))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compute_fci_curved_quadrilateral_cell_areas_2d_vjp

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

    subroutine validate_quadrilateral_vertices( &
            cell_vertices, cell_count, status)
        real(dp), intent(in) :: cell_vertices(:, :, :)
        integer, intent(out) :: cell_count
        type(fortsparse_status_t), intent(out) :: status
        integer :: cell, vertex, other

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI quadrilateral vertices have incompatible shape")
        cell_count = size(cell_vertices, 3)
        if (size(cell_vertices, 1) /= 2 .or. &
            size(cell_vertices, 2) /= 4 .or. cell_count < 1) return
        if (any(.not. ieee_is_finite(cell_vertices))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadrilateral vertices contain non-finite values")
            return
        end if
        do cell = 1, cell_count
            do vertex = 1, 4
                do other = vertex + 1, 4
                    if (all(cell_vertices(:, vertex, cell) == &
                        cell_vertices(:, other, cell))) then
                        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            "FCI quadrilateral cells contain repeated vertices")
                        return
                    end if
                end do
            end do
            if (segments_intersect(cell_vertices(1, 1, cell), &
                cell_vertices(2, 1, cell), cell_vertices(1, 2, cell), &
                cell_vertices(2, 2, cell), cell_vertices(1, 3, cell), &
                cell_vertices(2, 3, cell), cell_vertices(1, 4, cell), &
                cell_vertices(2, 4, cell)) .or. &
                segments_intersect(cell_vertices(1, 2, cell), &
                cell_vertices(2, 2, cell), cell_vertices(1, 3, cell), &
                cell_vertices(2, 3, cell), cell_vertices(1, 4, cell), &
                cell_vertices(2, 4, cell), cell_vertices(1, 1, cell), &
                cell_vertices(2, 1, cell))) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI quadrilateral cells are self-intersecting")
                return
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_quadrilateral_vertices

    subroutine validate_curved_quadrilateral_vertices( &
            cell_vertices, edge_controls, cell_count, status)
        real(dp), intent(in) :: cell_vertices(:, :, :)
        real(dp), intent(in) :: edge_controls(:, :, :)
        integer, intent(out) :: cell_count
        type(fortsparse_status_t), intent(out) :: status

        call validate_quadrilateral_vertices(cell_vertices, cell_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(edge_controls, 1) /= 2 .or. &
            size(edge_controls, 2) /= 4 .or. &
            size(edge_controls, 3) /= cell_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI curved quadrilateral controls have incompatible shape")
            return
        end if
        if (any(.not. ieee_is_finite(edge_controls))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI curved quadrilateral controls contain non-finite values")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_curved_quadrilateral_vertices

    pure logical function segments_intersect( &
            first_x, first_y, second_x, second_y, third_x, third_y, &
            fourth_x, fourth_y)
        real(dp), intent(in) :: first_x, first_y, second_x, second_y
        real(dp), intent(in) :: third_x, third_y, fourth_x, fourth_y
        real(dp) :: first_cross, second_cross, third_cross, fourth_cross

        first_cross = orientation(first_x, first_y, second_x, second_y, &
            third_x, third_y)
        second_cross = orientation(first_x, first_y, second_x, second_y, &
            fourth_x, fourth_y)
        third_cross = orientation(third_x, third_y, fourth_x, fourth_y, &
            first_x, first_y)
        fourth_cross = orientation(third_x, third_y, fourth_x, fourth_y, &
            second_x, second_y)
        segments_intersect = .false.
        if (((first_cross > 0.0_dp .and. second_cross < 0.0_dp) .or. &
            (first_cross < 0.0_dp .and. second_cross > 0.0_dp)) .and. &
            ((third_cross > 0.0_dp .and. fourth_cross < 0.0_dp) .or. &
            (third_cross < 0.0_dp .and. fourth_cross > 0.0_dp))) then
            segments_intersect = .true.
            return
        end if
        if (first_cross == 0.0_dp .and. point_on_segment( &
            third_x, third_y, first_x, first_y, second_x, second_y)) &
            segments_intersect = .true.
        if (second_cross == 0.0_dp .and. point_on_segment( &
            fourth_x, fourth_y, first_x, first_y, second_x, second_y)) &
            segments_intersect = .true.
        if (third_cross == 0.0_dp .and. point_on_segment( &
            first_x, first_y, third_x, third_y, fourth_x, fourth_y)) &
            segments_intersect = .true.
        if (fourth_cross == 0.0_dp .and. point_on_segment( &
            second_x, second_y, third_x, third_y, fourth_x, fourth_y)) &
            segments_intersect = .true.
    end function segments_intersect

    pure real(dp) function orientation(first_x, first_y, second_x, second_y, &
            third_x, third_y) result(cross)
        real(dp), intent(in) :: first_x, first_y, second_x, second_y
        real(dp), intent(in) :: third_x, third_y

        cross = (second_x - first_x)*(third_y - first_y) - &
            (second_y - first_y)*(third_x - first_x)
    end function orientation

    pure logical function point_on_segment(point_x, point_y, first_x, first_y, &
            second_x, second_y)
        real(dp), intent(in) :: point_x, point_y, first_x, first_y
        real(dp), intent(in) :: second_x, second_y

        point_on_segment = point_x >= min(first_x, second_x) .and. &
            point_x <= max(first_x, second_x) .and. &
            point_y >= min(first_y, second_y) .and. &
            point_y <= max(first_y, second_y)
    end function point_on_segment

end module fortfem_fci_support_geometry
