module fortfem_fourier_lines
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_line_tables, only: tensor_line_tables_t, LINE_FOURIER
    use fortfem_fourier_mode_registry, only: fourier_mode_registry_t, &
        initialize_fourier_mode_registry, evaluate_fourier_mode
    use fortsparse, only: fortsparse_status_t
    implicit none
    private
    public :: build_fourier_line, build_fourier_line_pair
contains
    subroutine build_fourier_line(m, line, status, nq, nfp)
        integer, intent(in) :: m
        type(tensor_line_tables_t), intent(out) :: line
        integer, intent(out) :: status
        integer, optional, intent(in) :: nq, nfp
        type(fourier_mode_registry_t) :: registry
        type(fortsparse_status_t) :: result
        complex(dp) :: value, dr, dt, dx
        integer :: order, periods, i, q, mode(1)
        integer, parameter :: zero(1) = [0]

        status = 1
        if (m < 0) return
        periods = 1
        if (present(nfp)) periods = nfp
        order = 2*m + 1
        if (present(nq)) order = nq
        if (periods < 1 .or. order < 2*m + 1) return
        line%family = LINE_FOURIER
        line%degree = m
        line%nelem = 1
        line%nloc = 2*m + 1
        line%nq = order
        line%ndof = line%nloc
        line%periodic = .true.
        line%x0 = 0.0_dp
        line%x1 = 2*acos(-1.0_dp)/real(periods, dp)
        allocate(line%xq(order, 1), line%wq(order, 1), &
            line%b(order, line%nloc, 1), line%d(order, line%nloc, 1), &
            line%dof(line%nloc, 1))
        do i = 1, line%nloc
            line%dof(i, 1) = i
        end do
        line%wq = line%x1/real(order, dp)
        line%b(:, 1, 1) = 1.0_dp
        line%d(:, 1, 1) = 0.0_dp
        do q = 1, order
            line%xq(q, 1) = real(q - 1, dp)*line%wq(q, 1)
        end do
        ! A single-mode registry keeps its validation independent of M.
        do i = 1, m
            mode(1) = i
            call initialize_fourier_mode_registry(registry, zero, mode, &
                periods, 0.0_dp, 0.0_dp, .false., radial_powers=zero, status=result)
            if (result%code /= 0) return
            do q = 1, order
                call evaluate_fourier_mode(registry, 1, 1.0_dp, 0.0_dp, &
                    line%xq(q, 1), value, dr, dt, dx, result)
                if (result%code /= 0) return
                line%b(q, 2*i, 1) = real(value, dp)
                line%b(q, 2*i + 1, 1) = aimag(value)
                line%d(q, 2*i, 1) = real(dx, dp)
                line%d(q, 2*i + 1, 1) = aimag(dx)
            end do
        end do
        status = 0
    end subroutine build_fourier_line

    subroutine build_fourier_line_pair(m, v0, v1, incidence, status, nq, nfp)
        integer, intent(in) :: m
        type(tensor_line_tables_t), intent(out) :: v0, v1
        real(dp), allocatable, intent(out) :: incidence(:, :)
        integer, intent(out) :: status
        integer, optional, intent(in) :: nq, nfp
        integer :: i, periods
        real(dp) :: frequency

        call build_fourier_line(m, v0, status, nq, nfp)
        if (status /= 0) return
        ! The derivative image is the zero-mean span, with no constant mode.
        v1%family = LINE_FOURIER
        v1%degree = m
        v1%nelem = 1
        v1%nloc = 2*m
        v1%nq = v0%nq
        v1%ndof = 2*m
        v1%periodic = .true.
        v1%x0 = v0%x0
        v1%x1 = v0%x1
        v1%xq = v0%xq
        v1%wq = v0%wq
        v1%b = v0%b(:, 2:, :)
        v1%d = v0%d(:, 2:, :)
        allocate(v1%dof(2*m, 1), incidence(2*m, 2*m + 1))
        do i = 1, 2*m
            v1%dof(i, 1) = i
        end do
        incidence = 0.0_dp
        periods = 1
        if (present(nfp)) periods = nfp
        do i = 1, m
            frequency = real(i, dp)*real(periods, dp)
            incidence(2*i - 1, 2*i + 1) = frequency
            incidence(2*i, 2*i) = -frequency
        end do
    end subroutine build_fourier_line_pair
end module fortfem_fourier_lines
