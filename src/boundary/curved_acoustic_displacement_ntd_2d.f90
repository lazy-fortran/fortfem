module fortfem_curved_acoustic_displacement_ntd_2d
    !! Acoustic displacement-to-pressure map on a closed polygonal interface.
    !!
    !! For the time convention exp(-i omega t), the fluid momentum equation
    !! gives q = d_n p = rho omega^2 u_n.  The outgoing exterior Calderon
    !! equation
    !!
    !!     (1/2 M - K) p + V q = 0
    !!
    !! is discretized with continuous P1 pressure and panelwise P0 normal
    !! displacement.  Panel orientation must be counter-clockwise, so the
    !! right normal points from the bounded FEM domain into the exterior.
    use fortfem_helmholtz_boundary_operators_2d, only: &
        assemble_helmholtz_double_layer_mixed_linear, &
        assemble_helmholtz_hypersingular_linear, &
        assemble_helmholtz_single_layer_constant
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    implicit none

    private

    public :: apply_curved_acoustic_displacement_ntd_2d

contains

    subroutine apply_curved_acoustic_displacement_ntd_2d( &
            panel_start, panel_end, panel_nodes, wavenumber, density, &
            angular_frequency, quadrature_order, normal_displacement, &
            pressure, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(in) :: wavenumber, density, angular_frequency
        integer, intent(in) :: quadrature_order
        complex(dp), intent(in) :: normal_displacement(:)
        complex(dp), intent(out) :: pressure(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: double_layer(:, :), hypersingular(:, :)
        complex(dp), allocatable :: matrix(:, :)
        complex(dp), allocatable :: neumann_trace(:), right_hand_side(:)
        complex(dp), allocatable :: single_layer(:, :)
        real(dp), allocatable :: mass(:, :)
        integer :: node_count, operator_status, panel_count

        pressure = (0.0_dp, 0.0_dp)
        status = 1
        panel_count = size(panel_start, 2)
        node_count = size(pressure)
        if (panel_count < 1 .or. node_count /= panel_count) return
        if (size(panel_start, 1) /= 2 .or. size(panel_end, 1) /= 2) return
        if (size(panel_end, 2) /= panel_count) return
        if (size(panel_nodes, 1) /= 2 .or. &
            size(panel_nodes, 2) /= panel_count) return
        if (any(panel_nodes < 1) .or. any(panel_nodes > node_count)) return
        if (size(normal_displacement) /= panel_count) return
        if (wavenumber <= 0.0_dp .or. density <= 0.0_dp .or. &
            angular_frequency <= 0.0_dp .or. quadrature_order < 1) return

        allocate(mass(panel_count, node_count))
        call assemble_boundary_mass( &
            panel_start, panel_end, panel_nodes, mass, operator_status)
        if (operator_status /= 0) return
        allocate(double_layer(panel_count, node_count))
        call assemble_helmholtz_double_layer_mixed_linear( &
            panel_start, panel_end, panel_nodes, node_count, wavenumber, &
            quadrature_order, double_layer, operator_status)
        if (operator_status /= 0) return
        allocate(single_layer(panel_count, panel_count))
        call assemble_helmholtz_single_layer_constant( &
            panel_start, panel_end, wavenumber, quadrature_order, &
            single_layer, operator_status)
        if (operator_status /= 0) return
        allocate(hypersingular(node_count, node_count))
        call assemble_helmholtz_hypersingular_linear( &
            panel_start, panel_end, panel_nodes, node_count, wavenumber, &
            quadrature_order, hypersingular, operator_status)
        if (operator_status /= 0) return

        allocate(matrix(panel_count, node_count))
        allocate(neumann_trace(panel_count), right_hand_side(panel_count))
        neumann_trace = density*angular_frequency**2*normal_displacement
        ! Burton--Miller combination of both exterior Calderon equations.
        ! Multiplication of the first equation by i*k gives it the same
        ! physical scaling as the hypersingular second equation and removes
        ! the fictitious interior resonances of either equation alone.
        matrix = hypersingular + cmplx(0.0_dp, wavenumber, dp)*( &
            0.5_dp*mass - double_layer)
        right_hand_side = -matmul( &
            0.5_dp*transpose(mass) + transpose(double_layer), &
            neumann_trace) - cmplx(0.0_dp, wavenumber, dp)* &
            matmul(single_layer, neumann_trace)
        call dense_solve(matrix, right_hand_side, pressure, operator_status)
        if (operator_status /= 0) then
            status = 2
            pressure = (0.0_dp, 0.0_dp)
            return
        end if
        status = 0
    end subroutine apply_curved_acoustic_displacement_ntd_2d

    pure subroutine assemble_boundary_mass( &
            panel_start, panel_end, panel_nodes, mass, status)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        integer, intent(in) :: panel_nodes(:, :)
        real(dp), intent(out) :: mass(:, :)
        integer, intent(out) :: status

        real(dp) :: length
        integer :: endpoint, panel

        mass = 0.0_dp
        status = 1
        do panel = 1, size(panel_start, 2)
            length = norm2(panel_end(:, panel) - panel_start(:, panel))
            if (length <= 0.0_dp) return
            do endpoint = 1, 2
                mass(panel, panel_nodes(endpoint, panel)) = &
                    mass(panel, panel_nodes(endpoint, panel)) + 0.5_dp*length
            end do
        end do
        status = 0
    end subroutine assemble_boundary_mass

end module fortfem_curved_acoustic_displacement_ntd_2d
