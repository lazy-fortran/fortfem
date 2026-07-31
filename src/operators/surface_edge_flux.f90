module fortfem_surface_edge_flux
    !! Geometry-to-edge contraction for tangential surface flux.
    !!
    !! For an edge quadrature row q and oriented surface conormal c, the
    !! integrated edge flux is
    !!
    !!   f_e = sum_q w_qe dot(c_qe, K_qe).
    !!
    !! The conormal includes the chosen edge orientation; its magnitude may
    !! carry a metric factor.  This module does not impose a surface law or
    !! require the conormal to be unit length.  The resulting edge flux can
    !! be passed to the topology-only surface-edge balance primitive.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_surface_edge_flux
    public :: assemble_surface_edge_flux_jvp
    public :: assemble_surface_edge_flux_vjp

contains

    subroutine assemble_surface_edge_flux( &
            edge_conormal, edge_weights, surface_current, edge_flux, status)
        !! Contract oriented conormal quadrature against a surface current.
        real(dp), intent(in) :: edge_conormal(:, :, :)
        real(dp), intent(in) :: edge_weights(:, :)
        real(dp), intent(in) :: surface_current(:, :, :)
        real(dp), intent(out) :: edge_flux(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, edge_count, quadrature, edge

        edge_flux = 0.0_dp
        call validate_surface_edge_flux_inputs( &
            edge_conormal, edge_weights, surface_current, &
            quadrature_count, edge_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(edge_flux) /= edge_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface edge flux output has incompatible size")
            return
        end if

        do quadrature = 1, quadrature_count
            do edge = 1, edge_count
                edge_flux(edge) = edge_flux(edge) + &
                    edge_weights(quadrature, edge)*dot_product( &
                    edge_conormal(quadrature, edge, :), &
                    surface_current(quadrature, edge, :))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_edge_flux

    subroutine assemble_surface_edge_flux_jvp( &
            edge_conormal, edge_weights, surface_current, &
            edge_conormal_dot, edge_weights_dot, surface_current_dot, &
            edge_flux_dot, status)
        !! Apply the product-rule directional derivative of edge fluxes.
        real(dp), intent(in) :: edge_conormal(:, :, :)
        real(dp), intent(in) :: edge_weights(:, :)
        real(dp), intent(in) :: surface_current(:, :, :)
        real(dp), intent(in) :: edge_conormal_dot(:, :, :)
        real(dp), intent(in) :: edge_weights_dot(:, :)
        real(dp), intent(in) :: surface_current_dot(:, :, :)
        real(dp), intent(out) :: edge_flux_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, edge_count, quadrature, edge

        edge_flux_dot = 0.0_dp
        call validate_surface_edge_flux_inputs( &
            edge_conormal, edge_weights, surface_current, &
            quadrature_count, edge_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(edge_flux_dot) /= edge_count .or. &
            size(edge_conormal_dot, 1) /= quadrature_count .or. &
            size(edge_conormal_dot, 2) /= edge_count .or. &
            size(edge_conormal_dot, 3) /= 3 .or. &
            size(edge_weights_dot, 1) /= quadrature_count .or. &
            size(edge_weights_dot, 2) /= edge_count .or. &
            size(surface_current_dot, 1) /= quadrature_count .or. &
            size(surface_current_dot, 2) /= edge_count .or. &
            size(surface_current_dot, 3) /= 3 .or. &
            any(.not. ieee_is_finite(edge_conormal_dot)) .or. &
            any(.not. ieee_is_finite(edge_weights_dot)) .or. &
            any(.not. ieee_is_finite(surface_current_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface edge flux JVP received incompatible increments")
            return
        end if

        do quadrature = 1, quadrature_count
            do edge = 1, edge_count
                edge_flux_dot(edge) = edge_flux_dot(edge) + &
                    edge_weights_dot(quadrature, edge)*dot_product( &
                    edge_conormal(quadrature, edge, :), &
                    surface_current(quadrature, edge, :)) + &
                    edge_weights(quadrature, edge)*dot_product( &
                    edge_conormal_dot(quadrature, edge, :), &
                    surface_current(quadrature, edge, :)) + &
                    edge_weights(quadrature, edge)*dot_product( &
                    edge_conormal(quadrature, edge, :), &
                    surface_current_dot(quadrature, edge, :))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_edge_flux_jvp

    subroutine assemble_surface_edge_flux_vjp( &
            edge_conormal, edge_weights, surface_current, edge_flux_bar, &
            edge_conormal_bar, edge_weights_bar, surface_current_bar, status)
        !! Apply the real transpose of the edge-flux contraction.
        real(dp), intent(in) :: edge_conormal(:, :, :)
        real(dp), intent(in) :: edge_weights(:, :)
        real(dp), intent(in) :: surface_current(:, :, :)
        real(dp), intent(in) :: edge_flux_bar(:)
        real(dp), intent(out) :: edge_conormal_bar(:, :, :)
        real(dp), intent(out) :: edge_weights_bar(:, :)
        real(dp), intent(out) :: surface_current_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, edge_count, quadrature, edge

        edge_conormal_bar = 0.0_dp
        edge_weights_bar = 0.0_dp
        surface_current_bar = 0.0_dp
        call validate_surface_edge_flux_inputs( &
            edge_conormal, edge_weights, surface_current, &
            quadrature_count, edge_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(edge_flux_bar) /= edge_count .or. &
            size(edge_conormal_bar, 1) /= quadrature_count .or. &
            size(edge_conormal_bar, 2) /= edge_count .or. &
            size(edge_conormal_bar, 3) /= 3 .or. &
            size(edge_weights_bar, 1) /= quadrature_count .or. &
            size(edge_weights_bar, 2) /= edge_count .or. &
            size(surface_current_bar, 1) /= quadrature_count .or. &
            size(surface_current_bar, 2) /= edge_count .or. &
            size(surface_current_bar, 3) /= 3 .or. &
            any(.not. ieee_is_finite(edge_flux_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface edge flux VJP received incompatible cotangents")
            return
        end if

        do quadrature = 1, quadrature_count
            do edge = 1, edge_count
                edge_conormal_bar(quadrature, edge, :) = &
                    edge_conormal_bar(quadrature, edge, :) + &
                    edge_weights(quadrature, edge)* &
                    surface_current(quadrature, edge, :)*edge_flux_bar(edge)
                surface_current_bar(quadrature, edge, :) = &
                    surface_current_bar(quadrature, edge, :) + &
                    edge_weights(quadrature, edge)* &
                    edge_conormal(quadrature, edge, :)*edge_flux_bar(edge)
                edge_weights_bar(quadrature, edge) = &
                    edge_weights_bar(quadrature, edge) + edge_flux_bar(edge)* &
                    dot_product(edge_conormal(quadrature, edge, :), &
                    surface_current(quadrature, edge, :))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_edge_flux_vjp

    subroutine validate_surface_edge_flux_inputs( &
            edge_conormal, edge_weights, surface_current, &
            quadrature_count, edge_count, status)
        real(dp), intent(in) :: edge_conormal(:, :, :)
        real(dp), intent(in) :: edge_weights(:, :)
        real(dp), intent(in) :: surface_current(:, :, :)
        integer, intent(out) :: quadrature_count, edge_count
        type(fortsparse_status_t), intent(out) :: status

        quadrature_count = size(edge_conormal, 1)
        edge_count = size(edge_conormal, 2)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface edge flux received incompatible arrays")
        if (quadrature_count < 1 .or. edge_count < 1) return
        if (size(edge_conormal, 3) /= 3 .or. &
            size(edge_weights, 1) /= quadrature_count .or. &
            size(edge_weights, 2) /= edge_count .or. &
            size(surface_current, 1) /= quadrature_count .or. &
            size(surface_current, 2) /= edge_count .or. &
            size(surface_current, 3) /= 3) return
        if (any(.not. ieee_is_finite(edge_conormal)) .or. &
            any(.not. ieee_is_finite(edge_weights)) .or. &
            any(.not. ieee_is_finite(surface_current)) .or. &
            any(edge_weights <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_surface_edge_flux_inputs

end module fortfem_surface_edge_flux
