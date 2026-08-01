module fortfem_period_constraints
    !! Fixed-topology complex cycle-period constraints.
    !!
    !! A real cycle basis contracts a complex H(curl) edge field into its
    !! periods.  The block is intentionally neutral: cycle orientation,
    !! physical units, and gauge selection are caller-owned.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: assemble_period_constraints
    public :: assemble_period_constraints_jvp
    public :: assemble_period_constraints_vjp

contains

    subroutine assemble_period_constraints( &
            cycles, edge_values, target_periods, residual, status)
        real(dp), intent(in) :: cycles(:, :)
        complex(dp), intent(in) :: edge_values(:), target_periods(:)
        complex(dp), intent(out) :: residual(:)
        integer, intent(out) :: status

        integer :: edge, period

        residual = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_shapes(cycles, edge_values, target_periods, residual) .or. &
            .not. finite_data(cycles, edge_values, target_periods)) return
        do period = 1, size(cycles, 2)
            residual(period) = -target_periods(period)
            do edge = 1, size(cycles, 1)
                residual(period) = residual(period) + &
                    cycles(edge, period)*edge_values(edge)
            end do
        end do
        status = 0
    end subroutine assemble_period_constraints

    subroutine assemble_period_constraints_jvp( &
            cycles, edge_values, target_periods, cycles_dot, edge_values_dot, &
            target_periods_dot, residual_dot, status)
        real(dp), intent(in) :: cycles(:, :), cycles_dot(:, :)
        complex(dp), intent(in) :: edge_values(:), target_periods(:)
        complex(dp), intent(in) :: edge_values_dot(:), target_periods_dot(:)
        complex(dp), intent(out) :: residual_dot(:)
        integer, intent(out) :: status

        integer :: edge, period

        residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_shapes(cycles, edge_values, target_periods, residual_dot) .or. &
            any(shape(cycles_dot) /= shape(cycles)) .or. &
            size(edge_values_dot) /= size(edge_values) .or. &
            size(target_periods_dot) /= size(target_periods) .or. &
            .not. finite_data(cycles, edge_values, target_periods) .or. &
            .not. all(ieee_is_finite(cycles_dot)) .or. &
            .not. finite_complex(edge_values_dot) .or. &
            .not. finite_complex(target_periods_dot)) return
        do period = 1, size(cycles, 2)
            residual_dot(period) = -target_periods_dot(period)
            do edge = 1, size(cycles, 1)
                residual_dot(period) = residual_dot(period) + &
                    cycles_dot(edge, period)*edge_values(edge) + &
                    cycles(edge, period)*edge_values_dot(edge)
            end do
        end do
        status = 0
    end subroutine assemble_period_constraints_jvp

    subroutine assemble_period_constraints_vjp( &
            cycles, edge_values, target_periods, residual_bar, cycles_bar, &
            edge_values_bar, target_periods_bar, status)
        real(dp), intent(in) :: cycles(:, :)
        complex(dp), intent(in) :: edge_values(:), target_periods(:)
        complex(dp), intent(in) :: residual_bar(:)
        real(dp), intent(out) :: cycles_bar(:, :)
        complex(dp), intent(out) :: edge_values_bar(:), target_periods_bar(:)
        integer, intent(out) :: status

        integer :: edge, period

        cycles_bar = 0.0_dp
        edge_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        target_periods_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. valid_shapes(cycles, edge_values, target_periods, residual_bar) .or. &
            any(shape(cycles_bar) /= shape(cycles)) .or. &
            size(edge_values_bar) /= size(edge_values) .or. &
            size(target_periods_bar) /= size(target_periods) .or. &
            .not. finite_data(cycles, edge_values, target_periods) .or. &
            .not. finite_complex(residual_bar)) return
        do period = 1, size(cycles, 2)
            target_periods_bar(period) = -residual_bar(period)
            do edge = 1, size(cycles, 1)
                cycles_bar(edge, period) = real(conjg(residual_bar(period))* &
                    edge_values(edge), dp)
                edge_values_bar(edge) = edge_values_bar(edge) + &
                    cycles(edge, period)*residual_bar(period)
            end do
        end do
        status = 0
    end subroutine assemble_period_constraints_vjp

    logical function valid_shapes(cycles, edge_values, target_periods, residual) &
            result(valid)
        real(dp), intent(in) :: cycles(:, :)
        complex(dp), intent(in) :: edge_values(:), target_periods(:), residual(:)

        valid = size(cycles, 1) > 0 .and. size(cycles, 2) > 0 .and. &
            size(edge_values) == size(cycles, 1) .and. &
            size(target_periods) == size(cycles, 2) .and. &
            size(residual) == size(cycles, 2)
    end function valid_shapes

    logical function finite_data(cycles, edge_values, target_periods) result(valid)
        real(dp), intent(in) :: cycles(:, :)
        complex(dp), intent(in) :: edge_values(:), target_periods(:)

        valid = all(ieee_is_finite(cycles)) .and. finite_complex(edge_values) .and. &
            finite_complex(target_periods)
    end function finite_data

    logical function finite_complex(values) result(valid)
        complex(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex

end module fortfem_period_constraints
