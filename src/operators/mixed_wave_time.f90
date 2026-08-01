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
    public :: advance_mixed_wave_midpoint_jvp
    public :: advance_mixed_wave_midpoint_vjp
    public :: advance_mixed_wave_symplectic_euler
    public :: advance_mixed_wave_symplectic_euler_jvp
    public :: advance_mixed_wave_symplectic_euler_vjp
    public :: advance_mixed_wave_strang
    public :: advance_mixed_wave_strang_jvp
    public :: advance_mixed_wave_strang_vjp

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

    subroutine advance_mixed_wave_midpoint_jvp( &
            mass_q, mass_v, coupling, time_step, q, v, time_step_dot, q_dot, &
            v_dot, q_next_dot, v_next_dot, status)
        !! Apply the tangent of one implicit-midpoint wave step.
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), coupling(:, :)
        real(dp), intent(in) :: time_step, q(:), v(:), time_step_dot
        real(dp), intent(in) :: q_dot(:), v_dot(:)
        real(dp), intent(out) :: q_next_dot(:), v_next_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: nq, nv, total_size, info
        real(dp), allocatable :: matrix_a(:, :), rhs(:), state_next(:)
        real(dp), allocatable :: rhs_dot(:), state_next_dot(:)

        q_next_dot = 0.0_dp
        v_next_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed wave midpoint JVP received incompatible blocks")
        nq = size(mass_q, 1)
        nv = size(mass_v, 1)
        if (nq < 1 .or. nv < 1) return
        if (size(mass_q, 2) /= nq .or. size(mass_v, 2) /= nv) return
        if (size(coupling, 1) /= nv .or. size(coupling, 2) /= nq) return
        if (size(q) /= nq .or. size(v) /= nv) return
        if (size(q_dot) /= nq .or. size(v_dot) /= nv) return
        if (size(q_next_dot) /= nq .or. size(v_next_dot) /= nv) return

        total_size = nq + nv
        allocate(matrix_a(total_size, total_size), rhs(total_size), &
            state_next(total_size), rhs_dot(total_size), &
            state_next_dot(total_size))
        call assemble_midpoint_system( &
            mass_q, mass_v, coupling, time_step, q, v, matrix_a, rhs)
        call dense_solve(matrix_a, rhs, state_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave midpoint JVP block is singular")
            return
        end if
        rhs_dot(:nq) = matmul(mass_q, q_dot) - &
            0.5_dp*time_step*matmul(transpose(coupling), v_dot) - &
            0.5_dp*time_step_dot*matmul(transpose(coupling), v) - &
            0.5_dp*time_step_dot*matmul(transpose(coupling), &
            state_next(nq + 1:total_size))
        rhs_dot(nq + 1:total_size) = matmul(mass_v, v_dot) + &
            0.5_dp*time_step*matmul(coupling, q_dot) + &
            0.5_dp*time_step_dot*matmul(coupling, q) + &
            0.5_dp*time_step_dot*matmul(coupling, state_next(:nq))
        call dense_solve(matrix_a, rhs_dot, state_next_dot, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave midpoint JVP block is singular")
            return
        end if
        q_next_dot = state_next_dot(:nq)
        v_next_dot = state_next_dot(nq + 1:total_size)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_mixed_wave_midpoint_jvp

    subroutine advance_mixed_wave_midpoint_vjp( &
            mass_q, mass_v, coupling, time_step, q, v, q_next_bar, v_next_bar, &
            q_bar, v_bar, time_step_bar, status)
        !! Apply the real adjoint of one implicit-midpoint wave step.
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), coupling(:, :)
        real(dp), intent(in) :: time_step, q(:), v(:)
        real(dp), intent(in) :: q_next_bar(:), v_next_bar(:)
        real(dp), intent(out) :: q_bar(:), v_bar(:), time_step_bar
        type(fortsparse_status_t), intent(out) :: status

        integer :: nq, nv, total_size, info
        real(dp), allocatable :: matrix_a(:, :), rhs(:), state_next(:)
        real(dp), allocatable :: state_next_bar(:), state_bar(:)

        q_bar = 0.0_dp
        v_bar = 0.0_dp
        time_step_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed wave midpoint VJP received incompatible blocks")
        nq = size(mass_q, 1)
        nv = size(mass_v, 1)
        if (nq < 1 .or. nv < 1) return
        if (size(mass_q, 2) /= nq .or. size(mass_v, 2) /= nv) return
        if (size(coupling, 1) /= nv .or. size(coupling, 2) /= nq) return
        if (size(q) /= nq .or. size(v) /= nv) return
        if (size(q_next_bar) /= nq .or. size(v_next_bar) /= nv) return
        if (size(q_bar) /= nq .or. size(v_bar) /= nv) return

        total_size = nq + nv
        allocate(matrix_a(total_size, total_size), rhs(total_size), &
            state_next(total_size), state_next_bar(total_size), &
            state_bar(total_size))
        call assemble_midpoint_system( &
            mass_q, mass_v, coupling, time_step, q, v, matrix_a, rhs)
        call dense_solve(matrix_a, rhs, state_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave midpoint VJP block is singular")
            return
        end if
        state_next_bar(:nq) = q_next_bar
        state_next_bar(nq + 1:total_size) = v_next_bar
        call dense_solve(transpose(matrix_a), state_next_bar, state_bar, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave midpoint VJP transpose block is singular")
            return
        end if
        q_bar = matmul(transpose(mass_q), state_bar(:nq)) + &
            0.5_dp*time_step*matmul(transpose(coupling), &
            state_bar(nq + 1:total_size))
        v_bar = matmul(transpose(mass_v), state_bar(nq + 1:total_size)) - &
            0.5_dp*time_step*matmul(coupling, state_bar(:nq))
        time_step_bar = -0.5_dp*dot_product(state_bar(:nq), &
            matmul(transpose(coupling), v + state_next(nq + 1:total_size))) + &
            0.5_dp*dot_product(state_bar(nq + 1:total_size), &
            matmul(coupling, q + state_next(:nq)))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_mixed_wave_midpoint_vjp

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

    subroutine advance_mixed_wave_symplectic_euler_jvp( &
            mass_q, mass_v, coupling, time_step, q, v, time_step_dot, q_dot, &
            v_dot, q_next_dot, v_next_dot, status)
        !! Apply the tangent of one partitioned symplectic-Euler step.
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), coupling(:, :)
        real(dp), intent(in) :: time_step, q(:), v(:), time_step_dot
        real(dp), intent(in) :: q_dot(:), v_dot(:)
        real(dp), intent(out) :: q_next_dot(:), v_next_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: nq, nv, info
        real(dp), allocatable :: v_rhs(:), v_next(:), v_rhs_dot(:)
        real(dp), allocatable :: q_rhs_dot(:)

        q_next_dot = 0.0_dp
        v_next_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed wave symplectic Euler JVP received incompatible blocks")
        nq = size(mass_q, 1)
        nv = size(mass_v, 1)
        if (nq < 1 .or. nv < 1) return
        if (size(mass_q, 2) /= nq .or. size(mass_v, 2) /= nv) return
        if (size(coupling, 1) /= nv .or. size(coupling, 2) /= nq) return
        if (size(q) /= nq .or. size(v) /= nv) return
        if (size(q_dot) /= nq .or. size(v_dot) /= nv) return
        if (size(q_next_dot) /= nq .or. size(v_next_dot) /= nv) return

        allocate(v_rhs(nv), v_next(nv), v_rhs_dot(nv), q_rhs_dot(nq))
        v_rhs = matmul(mass_v, v) + time_step*matmul(coupling, q)
        call dense_solve(mass_v, v_rhs, v_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave symplectic Euler JVP velocity mass is singular")
            return
        end if
        v_rhs_dot = matmul(mass_v, v_dot) + &
            time_step*matmul(coupling, q_dot) + &
            time_step_dot*matmul(coupling, q)
        call dense_solve(mass_v, v_rhs_dot, v_next_dot, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave symplectic Euler JVP velocity mass is singular")
            return
        end if
        q_rhs_dot = matmul(mass_q, q_dot) - &
            time_step*matmul(transpose(coupling), v_next_dot) - &
            time_step_dot*matmul(transpose(coupling), v_next)
        call dense_solve(mass_q, q_rhs_dot, q_next_dot, info)
        if (info /= 0) then
            q_next_dot = 0.0_dp
            v_next_dot = 0.0_dp
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave symplectic Euler JVP coordinate mass is singular")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_mixed_wave_symplectic_euler_jvp

    subroutine advance_mixed_wave_symplectic_euler_vjp( &
            mass_q, mass_v, coupling, time_step, q, v, q_next_bar, v_next_bar, &
            q_bar, v_bar, time_step_bar, status)
        !! Apply the real adjoint of one partitioned symplectic-Euler step.
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), coupling(:, :)
        real(dp), intent(in) :: time_step, q(:), v(:)
        real(dp), intent(in) :: q_next_bar(:), v_next_bar(:)
        real(dp), intent(out) :: q_bar(:), v_bar(:), time_step_bar
        type(fortsparse_status_t), intent(out) :: status

        integer :: nq, nv, info
        real(dp), allocatable :: v_rhs(:), v_next(:), q_rhs_bar(:)
        real(dp), allocatable :: v_next_bar_work(:), v_rhs_bar(:)

        q_bar = 0.0_dp
        v_bar = 0.0_dp
        time_step_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed wave symplectic Euler VJP received incompatible blocks")
        nq = size(mass_q, 1)
        nv = size(mass_v, 1)
        if (nq < 1 .or. nv < 1) return
        if (size(mass_q, 2) /= nq .or. size(mass_v, 2) /= nv) return
        if (size(coupling, 1) /= nv .or. size(coupling, 2) /= nq) return
        if (size(q) /= nq .or. size(v) /= nv) return
        if (size(q_next_bar) /= nq .or. size(v_next_bar) /= nv) return
        if (size(q_bar) /= nq .or. size(v_bar) /= nv) return

        allocate(v_rhs(nv), v_next(nv), q_rhs_bar(nq), &
            v_next_bar_work(nv), v_rhs_bar(nv))
        v_rhs = matmul(mass_v, v) + time_step*matmul(coupling, q)
        call dense_solve(mass_v, v_rhs, v_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave symplectic Euler VJP velocity mass is singular")
            return
        end if
        call dense_solve(transpose(mass_q), q_next_bar, q_rhs_bar, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave symplectic Euler VJP coordinate mass is singular")
            return
        end if
        q_bar = matmul(transpose(mass_q), q_rhs_bar)
        v_next_bar_work = v_next_bar - &
            time_step*matmul(coupling, q_rhs_bar)
        time_step_bar = -dot_product(q_rhs_bar, matmul(transpose(coupling), v_next))

        call dense_solve(transpose(mass_v), v_next_bar_work, v_rhs_bar, info)
        if (info /= 0) then
            q_bar = 0.0_dp
            time_step_bar = 0.0_dp
            call status_set(status, FORTSPARSE_SINGULAR, &
                "mixed wave symplectic Euler VJP velocity mass is singular")
            return
        end if
        q_bar = q_bar + time_step*matmul(transpose(coupling), v_rhs_bar)
        v_bar = matmul(transpose(mass_v), v_rhs_bar)
        time_step_bar = time_step_bar + &
            dot_product(v_rhs_bar, matmul(coupling, q))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_mixed_wave_symplectic_euler_vjp

    subroutine advance_mixed_wave_strang( &
            mass_q, mass_v, coupling_a, coupling_b, time_step, q, v, status)
        !! Advance a symmetric A(h/2)-B(h)-A(h/2) Cayley split.
        !!
        !! Each substep is an implicit midpoint map for a separately
        !! caller-owned mixed Hamiltonian block.  The composition is symmetric
        !! and reversible for signed time steps; dissipative operators must be
        !! composed outside this routine.
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :)
        real(dp), intent(in) :: coupling_a(:, :), coupling_b(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(inout) :: q(:), v(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: q_one(:), v_one(:), q_two(:), v_two(:)
        integer :: nq, nv

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed wave Strang split received incompatible blocks")
        nq = size(mass_q, 1)
        nv = size(mass_v, 1)
        if (nq < 1 .or. nv < 1) return
        if (size(mass_q, 2) /= nq .or. size(mass_v, 2) /= nv) return
        if (size(coupling_a, 1) /= nv .or. size(coupling_a, 2) /= nq) return
        if (size(coupling_b, 1) /= nv .or. size(coupling_b, 2) /= nq) return
        if (size(q) /= nq .or. size(v) /= nv) return

        allocate(q_one(nq), v_one(nv), q_two(nq), v_two(nv))
        q_one = q
        v_one = v
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling_a, 0.5_dp*time_step, q_one, v_one, status)
        if (status%code /= FORTSPARSE_OK) return
        q_two = q_one
        v_two = v_one
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling_b, time_step, q_two, v_two, status)
        if (status%code /= FORTSPARSE_OK) return
        q_one = q_two
        v_one = v_two
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling_a, 0.5_dp*time_step, q_one, v_one, status)
        if (status%code /= FORTSPARSE_OK) return
        q = q_one
        v = v_one
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_mixed_wave_strang

    subroutine advance_mixed_wave_strang_jvp( &
            mass_q, mass_v, coupling_a, coupling_b, time_step, q, v, &
            time_step_dot, q_dot, v_dot, q_next_dot, v_next_dot, status)
        !! Apply the tangent of the symmetric mixed-wave Strang split.
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :)
        real(dp), intent(in) :: coupling_a(:, :), coupling_b(:, :)
        real(dp), intent(in) :: time_step, q(:), v(:), time_step_dot
        real(dp), intent(in) :: q_dot(:), v_dot(:)
        real(dp), intent(out) :: q_next_dot(:), v_next_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: q_one(:), v_one(:), q_two(:), v_two(:)
        real(dp), allocatable :: q_one_dot(:), v_one_dot(:), q_two_dot(:), v_two_dot(:)
        integer :: nq, nv

        q_next_dot = 0.0_dp
        v_next_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed wave Strang JVP received incompatible blocks")
        nq = size(mass_q, 1)
        nv = size(mass_v, 1)
        if (nq < 1 .or. nv < 1) return
        if (size(mass_q, 2) /= nq .or. size(mass_v, 2) /= nv) return
        if (size(coupling_a, 1) /= nv .or. size(coupling_a, 2) /= nq) return
        if (size(coupling_b, 1) /= nv .or. size(coupling_b, 2) /= nq) return
        if (size(q) /= nq .or. size(v) /= nv) return
        if (size(q_dot) /= nq .or. size(v_dot) /= nv) return
        if (size(q_next_dot) /= nq .or. size(v_next_dot) /= nv) return

        allocate(q_one(nq), v_one(nv), q_two(nq), v_two(nv), &
            q_one_dot(nq), v_one_dot(nv), q_two_dot(nq), v_two_dot(nv))
        q_one = q
        v_one = v
        q_one_dot = q_dot
        v_one_dot = v_dot
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling_a, 0.5_dp*time_step, q_one, v_one, status)
        if (status%code /= FORTSPARSE_OK) return
        call advance_mixed_wave_midpoint_jvp( &
            mass_q, mass_v, coupling_a, 0.5_dp*time_step, q, v, &
            0.5_dp*time_step_dot, q_dot, v_dot, q_one_dot, v_one_dot, status)
        if (status%code /= FORTSPARSE_OK) return

        q_two = q_one
        v_two = v_one
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling_b, time_step, q_two, v_two, status)
        if (status%code /= FORTSPARSE_OK) return
        call advance_mixed_wave_midpoint_jvp( &
            mass_q, mass_v, coupling_b, time_step, q_one, v_one, time_step_dot, &
            q_one_dot, v_one_dot, q_two_dot, v_two_dot, status)
        if (status%code /= FORTSPARSE_OK) return

        q_one = q_two
        v_one = v_two
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling_a, 0.5_dp*time_step, q_one, v_one, status)
        if (status%code /= FORTSPARSE_OK) return
        call advance_mixed_wave_midpoint_jvp( &
            mass_q, mass_v, coupling_a, 0.5_dp*time_step, q_two, v_two, &
            0.5_dp*time_step_dot, q_two_dot, v_two_dot, q_next_dot, &
            v_next_dot, status)
    end subroutine advance_mixed_wave_strang_jvp

    subroutine advance_mixed_wave_strang_vjp( &
            mass_q, mass_v, coupling_a, coupling_b, time_step, q, v, &
            q_next_bar, v_next_bar, q_bar, v_bar, time_step_bar, status)
        !! Apply the real adjoint of the symmetric mixed-wave Strang split.
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :)
        real(dp), intent(in) :: coupling_a(:, :), coupling_b(:, :)
        real(dp), intent(in) :: time_step, q(:), v(:)
        real(dp), intent(in) :: q_next_bar(:), v_next_bar(:)
        real(dp), intent(out) :: q_bar(:), v_bar(:), time_step_bar
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: q_one(:), v_one(:), q_two(:), v_two(:)
        real(dp), allocatable :: q_bar_stage_three(:), v_bar_stage_three(:)
        real(dp), allocatable :: q_bar_stage_two(:), v_bar_stage_two(:)
        real(dp) :: time_step_bar_one, time_step_bar_two, time_step_bar_three
        integer :: nq, nv

        q_bar = 0.0_dp
        v_bar = 0.0_dp
        time_step_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "mixed wave Strang VJP received incompatible blocks")
        nq = size(mass_q, 1)
        nv = size(mass_v, 1)
        if (nq < 1 .or. nv < 1) return
        if (size(mass_q, 2) /= nq .or. size(mass_v, 2) /= nv) return
        if (size(coupling_a, 1) /= nv .or. size(coupling_a, 2) /= nq) return
        if (size(coupling_b, 1) /= nv .or. size(coupling_b, 2) /= nq) return
        if (size(q) /= nq .or. size(v) /= nv) return
        if (size(q_next_bar) /= nq .or. size(v_next_bar) /= nv) return
        if (size(q_bar) /= nq .or. size(v_bar) /= nv) return

        allocate(q_one(nq), v_one(nv), q_two(nq), v_two(nv), &
            q_bar_stage_three(nq), v_bar_stage_three(nv), &
            q_bar_stage_two(nq), v_bar_stage_two(nv))
        q_one = q
        v_one = v
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling_a, 0.5_dp*time_step, q_one, v_one, status)
        if (status%code /= FORTSPARSE_OK) return
        q_two = q_one
        v_two = v_one
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling_b, time_step, q_two, v_two, status)
        if (status%code /= FORTSPARSE_OK) return

        call advance_mixed_wave_midpoint_vjp( &
            mass_q, mass_v, coupling_a, 0.5_dp*time_step, q_two, v_two, &
            q_next_bar, v_next_bar, q_bar_stage_three, v_bar_stage_three, &
            time_step_bar_three, status)
        if (status%code /= FORTSPARSE_OK) return
        call advance_mixed_wave_midpoint_vjp( &
            mass_q, mass_v, coupling_b, time_step, q_one, v_one, &
            q_bar_stage_three, v_bar_stage_three, q_bar_stage_two, &
            v_bar_stage_two, &
            time_step_bar_two, status)
        if (status%code /= FORTSPARSE_OK) return
        call advance_mixed_wave_midpoint_vjp( &
            mass_q, mass_v, coupling_a, 0.5_dp*time_step, q, v, &
            q_bar_stage_two, v_bar_stage_two, q_bar, v_bar, &
            time_step_bar_one, status)
        if (status%code /= FORTSPARSE_OK) return
        time_step_bar = 0.5_dp*time_step_bar_one + time_step_bar_two + &
            0.5_dp*time_step_bar_three
    end subroutine advance_mixed_wave_strang_vjp

    subroutine assemble_midpoint_system( &
            mass_q, mass_v, coupling, time_step, q, v, matrix_a, rhs)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), coupling(:, :)
        real(dp), intent(in) :: time_step, q(:), v(:)
        real(dp), intent(out) :: matrix_a(:, :), rhs(:)

        integer :: nq, nv, total_size

        nq = size(mass_q, 1)
        nv = size(mass_v, 1)
        total_size = nq + nv
        matrix_a = 0.0_dp
        matrix_a(:nq, :nq) = mass_q
        matrix_a(:nq, nq + 1:total_size) = &
            0.5_dp*time_step*transpose(coupling)
        matrix_a(nq + 1:total_size, :nq) = &
            -0.5_dp*time_step*coupling
        matrix_a(nq + 1:total_size, nq + 1:total_size) = mass_v
        rhs(:nq) = matmul(mass_q, q) - &
            0.5_dp*time_step*matmul(transpose(coupling), v)
        rhs(nq + 1:total_size) = matmul(mass_v, v) + &
            0.5_dp*time_step*matmul(coupling, q)
    end subroutine assemble_midpoint_system

end module fortfem_mixed_wave_time
