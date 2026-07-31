module fortfem_surface_edge_balance
    !! Topology-only discrete divergence ledger for integrated surface flux.
    !!
    !! A caller supplies an oriented vertex-edge incidence matrix and the
    !! already integrated scalar flux through each oriented surface edge. The
    !! resulting incidence transpose is deliberately geometry-agnostic: edge
    !! conormals, quadrature, and the tangential surface-current law remain in
    !! the caller's trace space.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_surface_edge_flux_balance
    public :: assemble_surface_edge_flux_balance_jvp
    public :: assemble_surface_edge_flux_balance_vjp

contains

    subroutine assemble_surface_edge_flux_balance( &
            edge_boundary, edge_flux, vertex_balance, global_balance, status)
        integer, intent(in) :: edge_boundary(:, :)
        real(dp), intent(in) :: edge_flux(:)
        real(dp), intent(out) :: vertex_balance(:)
        real(dp), intent(out) :: global_balance
        type(fortsparse_status_t), intent(out) :: status
        integer :: vertex, edge

        vertex_balance = 0.0_dp
        global_balance = 0.0_dp
        call validate_edge_balance_inputs( &
            edge_boundary, edge_flux, vertex_balance, status)
        if (status%code /= FORTSPARSE_OK) return

        do edge = 1, size(edge_boundary, 2)
            do vertex = 1, size(edge_boundary, 1)
                vertex_balance(vertex) = vertex_balance(vertex) + &
                    real(edge_boundary(vertex, edge), dp)*edge_flux(edge)
            end do
        end do
        global_balance = sum(vertex_balance)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_edge_flux_balance

    subroutine assemble_surface_edge_flux_balance_jvp( &
            edge_boundary, edge_flux, edge_flux_dot, vertex_balance_dot, &
            global_balance_dot, status)
        integer, intent(in) :: edge_boundary(:, :)
        real(dp), intent(in) :: edge_flux(:), edge_flux_dot(:)
        real(dp), intent(out) :: vertex_balance_dot(:)
        real(dp), intent(out) :: global_balance_dot
        type(fortsparse_status_t), intent(out) :: status
        integer :: vertex, edge

        vertex_balance_dot = 0.0_dp
        global_balance_dot = 0.0_dp
        call validate_edge_balance_inputs( &
            edge_boundary, edge_flux, vertex_balance_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(edge_flux_dot) /= size(edge_boundary, 2) .or. &
            any(.not. ieee_is_finite(edge_flux_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface edge-balance JVP has incompatible increments")
            return
        end if

        do edge = 1, size(edge_boundary, 2)
            do vertex = 1, size(edge_boundary, 1)
                vertex_balance_dot(vertex) = vertex_balance_dot(vertex) + &
                    real(edge_boundary(vertex, edge), dp)*edge_flux_dot(edge)
            end do
        end do
        global_balance_dot = sum(vertex_balance_dot)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_edge_flux_balance_jvp

    subroutine assemble_surface_edge_flux_balance_vjp( &
            edge_boundary, edge_flux, vertex_balance_bar, global_balance_bar, &
            edge_flux_bar, status)
        integer, intent(in) :: edge_boundary(:, :)
        real(dp), intent(in) :: edge_flux(:), vertex_balance_bar(:)
        real(dp), intent(in) :: global_balance_bar
        real(dp), intent(out) :: edge_flux_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: vertex, edge

        edge_flux_bar = 0.0_dp
        call validate_edge_balance_inputs( &
            edge_boundary, edge_flux, vertex_balance_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(edge_flux_bar) /= size(edge_boundary, 2) .or. &
            any(.not. ieee_is_finite(vertex_balance_bar)) .or. &
            .not. ieee_is_finite(global_balance_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface edge-balance VJP has incompatible cotangents")
            return
        end if

        do edge = 1, size(edge_boundary, 2)
            do vertex = 1, size(edge_boundary, 1)
                edge_flux_bar(edge) = edge_flux_bar(edge) + &
                    real(edge_boundary(vertex, edge), dp) * &
                    (vertex_balance_bar(vertex) + global_balance_bar)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_surface_edge_flux_balance_vjp

    subroutine validate_edge_balance_inputs( &
            edge_boundary, edge_flux, vertex_balance, status)
        integer, intent(in) :: edge_boundary(:, :)
        real(dp), intent(in) :: edge_flux(:), vertex_balance(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: edge

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface edge balance has incompatible arrays")
        if (size(edge_flux) /= size(edge_boundary, 2) .or. &
            size(vertex_balance) /= size(edge_boundary, 1)) return
        if (any(.not. ieee_is_finite(edge_flux)) .or. &
            any(.not. ieee_is_finite(vertex_balance))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface edge balance received non-finite data")
            return
        end if
        do edge = 1, size(edge_boundary, 2)
            if (sum(edge_boundary(:, edge)) /= 0) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_edge_balance_inputs

end module fortfem_surface_edge_balance
