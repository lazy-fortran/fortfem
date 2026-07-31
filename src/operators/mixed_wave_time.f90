module fortfem_mixed_wave_time
    !! Structure-preserving time stepping for a mixed first-order wave state.
    !!
    !! The semidiscrete equations are
    !!
    !!   M_q q_dot + C^T v = 0,
    !!   M_v v_dot - C q = 0.
    !!
    !! Implicit midpoint is the Cayley map of this port-Hamiltonian block. For
    !! symmetric positive mass blocks it preserves the quadratic energy and is
    !! exactly reversible under a signed time-step reversal.
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, FORTSPARSE_SINGULAR
    implicit none
    private

    public :: advance_mixed_wave_midpoint
    public :: advance_mixed_wave_symplectic_euler

contains

    subroutine advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, time_step, q, v, status)
        !! Advance one implicit-midpoint step of the mixed wave state.
        real(dp), intent(in) :: mass_q(:, :)
        real(dp), intent(in) :: mass_v(:, :)
        real(dp), intent(in) :: coupling(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(inout) :: q(:)
        real(dp), intent(inout) :: v(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: nq, nv, total_size, info
        real(dp), allocatable :: midpoint_matrix(:, :), right_hand_side(:)
        real(dp), allocatable :: state_next(:)

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed wave midpoint received incompatible blocks")
        nq = size(mass_q, 1)
        nv = size(mass_v, 1)
        if (nq < 1 .or. nv < 1) return
        if (size(mass_q, 2) /= nq) return
        if (size(mass_v, 2) /= nv) return
        if (size(coupling, 1) /= nv) return
        if (size(coupling, 2) /= nq) return
        if (size(q) /= nq) return
        if (size(v) /= nv) return

        total_size = nq + nv
        allocate( &
            midpoint_matrix(total_size, total_size), &
            right_hand_side(total_size), state_next(total_size))
        midpoint_matrix = 0.0_dp
        midpoint_matrix(:nq, :nq) = mass_q
        midpoint_matrix(:nq, nq + 1:total_size) = &
            0.5_dp*time_step*transpose(coupling)
        midpoint_matrix(nq + 1:total_size, :nq) = &
            -0.5_dp*time_step*coupling
        midpoint_matrix(nq + 1:total_size, nq + 1:total_size) = mass_v

        right_hand_side(:nq) = matmul(mass_q, q) - &
            0.5_dp*time_step*matmul(transpose(coupling), v)
        right_hand_side(nq + 1:total_size) = matmul(mass_v, v) + &
            0.5_dp*time_step*matmul(coupling, q)

        call dense_solve(midpoint_matrix, right_hand_side, state_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave midpoint block is singular")
            return
        end if
        q = state_next(:nq)
        v = state_next(nq + 1:total_size)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_mixed_wave_midpoint

    subroutine advance_mixed_wave_symplectic_euler( &
            mass_q, mass_v, coupling, time_step, q, v, status)
        !! Advance one partitioned symplectic-Euler wave step.
        !!
        !! The velocity/flux block is updated first,
        !!
        !!   M_v v_{n+1} = M_v v_n + h C q_n,
        !!
        !! followed by the coordinate/pressure block,
        !!
        !!   M_q q_{n+1} = M_q q_n - h C^T v_{n+1}.
        !!
        !! For compatible canonical blocks this is a symplectic, generally
        !! energy-non-preserving first-order method.  The state is only
        !! committed after both mass solves succeed.
        real(dp), intent(in) :: mass_q(:, :)
        real(dp), intent(in) :: mass_v(:, :)
        real(dp), intent(in) :: coupling(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(inout) :: q(:)
        real(dp), intent(inout) :: v(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: nq, nv, info
        real(dp), allocatable :: q_rhs(:), v_rhs(:), q_next(:), v_next(:)

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed wave symplectic Euler received incompatible blocks")
        nq = size(mass_q, 1)
        nv = size(mass_v, 1)
        if (nq < 1 .or. nv < 1) return
        if (size(mass_q, 2) /= nq) return
        if (size(mass_v, 2) /= nv) return
        if (size(coupling, 1) /= nv) return
        if (size(coupling, 2) /= nq) return
        if (size(q) /= nq) return
        if (size(v) /= nv) return

        allocate(q_rhs(nq), v_rhs(nv), q_next(nq), v_next(nv))
        v_rhs = matmul(mass_v, v) + time_step*matmul(coupling, q)
        call dense_solve(mass_v, v_rhs, v_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave symplectic Euler velocity mass is singular")
            return
        end if
        q_rhs = matmul(mass_q, q) - &
            time_step*matmul(transpose(coupling), v_next)
        call dense_solve(mass_q, q_rhs, q_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave symplectic Euler coordinate mass is singular")
            return
        end if
        q = q_next
        v = v_next
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_mixed_wave_symplectic_euler

end module fortfem_mixed_wave_time
