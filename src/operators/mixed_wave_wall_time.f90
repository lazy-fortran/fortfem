module fortfem_mixed_wave_wall_time
    !! Structure-preserving midpoint coupling of a mixed wave port to a wall RL block.
    !!
    !! The neutral semidiscrete equations are
    !!
    !!   Mq qdot + C^T v = 0,
    !!   Mv vdot - C q - P^T i = 0,
    !!   L  idot + R i + P v = 0.
    !!
    !! Implicit midpoint preserves the skew wave/port power exchange and
    !! dissipates only i^T R i.  Geometry, port orientation, and constitutive
    !! normalization remain caller-owned.
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, FORTSPARSE_SINGULAR
    implicit none
    private

    public :: advance_mixed_wave_wall_midpoint
    public :: advance_mixed_wave_wall_midpoint_jvp
    public :: advance_mixed_wave_wall_midpoint_vjp
    public :: evaluate_mixed_wave_wall_energy_balance

contains

    subroutine advance_mixed_wave_wall_midpoint( &
            mass_q, mass_v, inductance, resistance, coupling, port, time_step, &
            q, v, current, status)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), inductance(:, :)
        real(dp), intent(in) :: resistance(:, :), coupling(:, :), port(:, :)
        real(dp), intent(in) :: time_step
        real(dp), intent(inout) :: q(:), v(:), current(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: matrix(:, :), rhs(:), state_next(:)
        integer :: nq, nv, ni, total_size, info

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "coupled wave-wall midpoint received incompatible blocks")
        if (.not. valid_inputs(mass_q, mass_v, inductance, resistance, coupling, &
            port, time_step, q, v, current)) return
        nq = size(q)
        nv = size(v)
        ni = size(current)
        total_size = nq + nv + ni
        allocate(matrix(total_size, total_size), rhs(total_size), state_next(total_size))
        call assemble_system( &
            mass_q, mass_v, inductance, resistance, coupling, port, time_step, &
            q, v, current, matrix, rhs)
        call dense_solve(matrix, rhs, state_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "coupled wave-wall midpoint block is singular")
            return
        end if
        q = state_next(:nq)
        v = state_next(nq + 1:nq + nv)
        current = state_next(nq + nv + 1:total_size)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_mixed_wave_wall_midpoint

    subroutine advance_mixed_wave_wall_midpoint_jvp( &
            mass_q, mass_v, inductance, resistance, coupling, port, time_step, &
            q, v, current, mass_q_dot, mass_v_dot, inductance_dot, resistance_dot, &
            coupling_dot, port_dot, time_step_dot, q_dot, v_dot, current_dot, &
            q_next_dot, v_next_dot, current_next_dot, status)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), inductance(:, :)
        real(dp), intent(in) :: resistance(:, :), coupling(:, :), port(:, :)
        real(dp), intent(in) :: time_step, q(:), v(:), current(:)
        real(dp), intent(in) :: mass_q_dot(:, :), mass_v_dot(:, :), inductance_dot(:, :)
        real(dp), intent(in) :: resistance_dot(:, :), coupling_dot(:, :), port_dot(:, :)
        real(dp), intent(in) :: time_step_dot, q_dot(:), v_dot(:), current_dot(:)
        real(dp), intent(out) :: q_next_dot(:), v_next_dot(:), current_next_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: matrix(:, :), rhs(:), state_next(:)
        real(dp), allocatable :: matrix_dot(:, :), rhs_dot(:), state_next_dot(:)
        integer :: nq, nv, ni, total_size, info

        q_next_dot = 0.0_dp
        v_next_dot = 0.0_dp
        current_next_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "coupled wave-wall midpoint JVP received incompatible blocks")
        if (.not. valid_inputs(mass_q, mass_v, inductance, resistance, coupling, &
            port, time_step, q, v, current)) return
        if (.not. same_shapes(mass_q_dot, mass_q, mass_v_dot, mass_v, &
            inductance_dot, inductance, resistance_dot, resistance, coupling_dot, coupling, &
            port_dot, port) .or. size(q_dot) /= size(q) .or. size(v_dot) /= size(v) .or. &
            size(current_dot) /= size(current) .or. size(q_next_dot) /= size(q) .or. &
            size(v_next_dot) /= size(v) .or. size(current_next_dot) /= size(current)) return
        nq = size(q)
        nv = size(v)
        ni = size(current)
        total_size = nq + nv + ni
        allocate(matrix(total_size, total_size), rhs(total_size), state_next(total_size), &
            matrix_dot(total_size, total_size), rhs_dot(total_size), &
            state_next_dot(total_size))
        call assemble_system( &
            mass_q, mass_v, inductance, resistance, coupling, port, time_step, &
            q, v, current, matrix, rhs)
        call assemble_system_dot( &
            mass_q, mass_v, inductance, resistance, coupling, port, time_step, q, v, current, &
            mass_q_dot, mass_v_dot, inductance_dot, resistance_dot, coupling_dot, port_dot, &
            time_step_dot, q_dot, v_dot, current_dot, matrix_dot, rhs_dot)
        call dense_solve(matrix, rhs, state_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "coupled wave-wall midpoint JVP block is singular")
            return
        end if
        call dense_solve(matrix, rhs_dot - matmul(matrix_dot, state_next), state_next_dot, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "coupled wave-wall midpoint JVP tangent block is singular")
            return
        end if
        q_next_dot = state_next_dot(:nq)
        v_next_dot = state_next_dot(nq + 1:nq + nv)
        current_next_dot = state_next_dot(nq + nv + 1:total_size)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_mixed_wave_wall_midpoint_jvp

    subroutine advance_mixed_wave_wall_midpoint_vjp( &
            mass_q, mass_v, inductance, resistance, coupling, port, time_step, &
            q, v, current, q_next_bar, v_next_bar, current_next_bar, mass_q_bar, &
            mass_v_bar, inductance_bar, resistance_bar, coupling_bar, port_bar, &
            q_bar, v_bar, current_bar, time_step_bar, status)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), inductance(:, :)
        real(dp), intent(in) :: resistance(:, :), coupling(:, :), port(:, :)
        real(dp), intent(in) :: time_step, q(:), v(:), current(:)
        real(dp), intent(in) :: q_next_bar(:), v_next_bar(:), current_next_bar(:)
        real(dp), intent(out) :: mass_q_bar(:, :), mass_v_bar(:, :), inductance_bar(:, :)
        real(dp), intent(out) :: resistance_bar(:, :), coupling_bar(:, :), port_bar(:, :)
        real(dp), intent(out) :: q_bar(:), v_bar(:), current_bar(:), time_step_bar
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: matrix(:, :), rhs(:), state_next(:), lambda(:), matrix_bar(:, :)
        integer :: nq, nv, ni, total_size, info
        real(dp) :: bq_bar(size(q)), bv_bar(size(v)), bi_bar(size(current))

        mass_q_bar = 0.0_dp
        mass_v_bar = 0.0_dp
        inductance_bar = 0.0_dp
        resistance_bar = 0.0_dp
        coupling_bar = 0.0_dp
        port_bar = 0.0_dp
        q_bar = 0.0_dp
        v_bar = 0.0_dp
        current_bar = 0.0_dp
        time_step_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "coupled wave-wall midpoint VJP received incompatible blocks")
        if (.not. valid_inputs(mass_q, mass_v, inductance, resistance, coupling, &
            port, time_step, q, v, current)) return
        if (size(q_next_bar) /= size(q) .or. size(v_next_bar) /= size(v) .or. &
            size(current_next_bar) /= size(current) .or. any(shape(mass_q_bar) /= shape(mass_q)) .or. &
            any(shape(mass_v_bar) /= shape(mass_v)) .or. any(shape(inductance_bar) /= shape(inductance)) .or. &
            any(shape(resistance_bar) /= shape(resistance)) .or. any(shape(coupling_bar) /= shape(coupling)) .or. &
            any(shape(port_bar) /= shape(port)) .or. size(q_bar) /= size(q) .or. &
            size(v_bar) /= size(v) .or. size(current_bar) /= size(current)) return
        nq = size(q)
        nv = size(v)
        ni = size(current)
        total_size = nq + nv + ni
        allocate(matrix(total_size, total_size), rhs(total_size), state_next(total_size), &
            lambda(total_size), matrix_bar(total_size, total_size))
        call assemble_system( &
            mass_q, mass_v, inductance, resistance, coupling, port, time_step, &
            q, v, current, matrix, rhs)
        call dense_solve(matrix, rhs, state_next, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "coupled wave-wall midpoint VJP block is singular")
            return
        end if
        rhs(:nq) = q_next_bar
        rhs(nq + 1:nq + nv) = v_next_bar
        rhs(nq + nv + 1:total_size) = current_next_bar
        call dense_solve(transpose(matrix), rhs, lambda, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "coupled wave-wall midpoint VJP transpose block is singular")
            return
        end if
        matrix_bar = -outer_product(lambda, state_next)
        bq_bar = lambda(:nq)
        bv_bar = lambda(nq + 1:nq + nv)
        bi_bar = lambda(nq + nv + 1:total_size)
        mass_q_bar = matrix_bar(:nq, :nq) + outer_product(bq_bar, q)
        mass_v_bar = matrix_bar(nq + 1:nq + nv, nq + 1:nq + nv) + outer_product(bv_bar, v)
        inductance_bar = matrix_bar(nq + nv + 1:total_size, nq + nv + 1:total_size) + &
            outer_product(bi_bar, current)
        resistance_bar = 0.5_dp*time_step*matrix_bar(nq + nv + 1:total_size, nq + nv + 1:total_size) - &
            0.5_dp*time_step*outer_product(bi_bar, current)
        coupling_bar = 0.5_dp*time_step*transpose(matrix_bar(:nq, nq + 1:nq + nv)) - &
            0.5_dp*time_step*matrix_bar(nq + 1:nq + nv, :nq) - &
            0.5_dp*time_step*outer_product(v, bq_bar) + &
            0.5_dp*time_step*outer_product(bv_bar, q)
        port_bar = -0.5_dp*time_step*transpose(matrix_bar(nq + 1:nq + nv, nq + nv + 1:total_size)) + &
            0.5_dp*time_step*matrix_bar(nq + nv + 1:total_size, nq + 1:nq + nv) + &
            0.5_dp*time_step*outer_product(current, bv_bar) - &
            0.5_dp*time_step*outer_product(bi_bar, v)
        q_bar = matmul(transpose(mass_q), bq_bar) + &
            0.5_dp*time_step*matmul(transpose(coupling), bv_bar)
        v_bar = -0.5_dp*time_step*matmul(coupling, bq_bar) + &
            matmul(transpose(mass_v), bv_bar) - 0.5_dp*time_step*matmul(transpose(port), bi_bar)
        current_bar = matmul(transpose(inductance), bi_bar) - &
            0.5_dp*time_step*matmul(transpose(resistance), bi_bar) + &
            0.5_dp*time_step*matmul(port, bv_bar)
        time_step_bar = 0.5_dp*sum(matrix_bar(:nq, nq + 1:nq + nv)*transpose(coupling)) - &
            0.5_dp*sum(matrix_bar(nq + 1:nq + nv, :nq)*coupling) - &
            0.5_dp*sum(matrix_bar(nq + 1:nq + nv, nq + nv + 1:total_size)*transpose(port)) + &
            0.5_dp*sum(matrix_bar(nq + nv + 1:total_size, nq + 1:nq + nv)*port) + &
            0.5_dp*sum(matrix_bar(nq + nv + 1:total_size, nq + nv + 1:total_size)*resistance) - &
            0.5_dp*dot_product(bq_bar, matmul(transpose(coupling), v)) + &
            0.5_dp*dot_product(bv_bar, matmul(coupling, q) + matmul(transpose(port), current)) - &
            0.5_dp*dot_product(bi_bar, matmul(resistance, current) + matmul(port, v))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine advance_mixed_wave_wall_midpoint_vjp

    subroutine evaluate_mixed_wave_wall_energy_balance( &
            mass_q, mass_v, inductance, resistance, coupling, port, time_step, &
            q, v, current, q_next, v_next, current_next, energy_n, energy_next, &
            dissipation, balance, status)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), inductance(:, :)
        real(dp), intent(in) :: resistance(:, :), coupling(:, :), port(:, :)
        real(dp), intent(in) :: time_step, q(:), v(:), current(:), q_next(:), v_next(:), current_next(:)
        real(dp), intent(out) :: energy_n, energy_next, dissipation, balance
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: current_mid(size(current))

        energy_n = 0.0_dp
        energy_next = 0.0_dp
        dissipation = 0.0_dp
        balance = huge(1.0_dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "coupled wave-wall energy ledger received incompatible blocks")
        if (.not. valid_inputs(mass_q, mass_v, inductance, resistance, coupling, &
            port, time_step, q, v, current) .or. size(q_next) /= size(q) .or. &
            size(v_next) /= size(v) .or. size(current_next) /= size(current)) return
        current_mid = 0.5_dp*(current + current_next)
        energy_n = 0.5_dp*(dot_product(q, matmul(mass_q, q)) + &
            dot_product(v, matmul(mass_v, v)) + dot_product(current, matmul(inductance, current)))
        energy_next = 0.5_dp*(dot_product(q_next, matmul(mass_q, q_next)) + &
            dot_product(v_next, matmul(mass_v, v_next)) + &
            dot_product(current_next, matmul(inductance, current_next)))
        dissipation = time_step*dot_product(current_mid, matmul(resistance, current_mid))
        balance = energy_next - energy_n + dissipation
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_mixed_wave_wall_energy_balance

    subroutine assemble_system( &
            mass_q, mass_v, inductance, resistance, coupling, port, time_step, &
            q, v, current, matrix, rhs)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), inductance(:, :)
        real(dp), intent(in) :: resistance(:, :), coupling(:, :), port(:, :), time_step
        real(dp), intent(in) :: q(:), v(:), current(:)
        real(dp), intent(out) :: matrix(:, :), rhs(:)
        integer :: nq, nv, ni, total_size

        nq = size(q)
        nv = size(v)
        ni = size(current)
        total_size = nq + nv + ni
        matrix = 0.0_dp
        matrix(:nq, :nq) = mass_q
        matrix(:nq, nq + 1:nq + nv) = 0.5_dp*time_step*transpose(coupling)
        matrix(nq + 1:nq + nv, :nq) = -0.5_dp*time_step*coupling
        matrix(nq + 1:nq + nv, nq + 1:nq + nv) = mass_v
        matrix(nq + 1:nq + nv, nq + nv + 1:total_size) = -0.5_dp*time_step*transpose(port)
        matrix(nq + nv + 1:total_size, nq + 1:nq + nv) = 0.5_dp*time_step*port
        matrix(nq + nv + 1:total_size, nq + nv + 1:total_size) = inductance + &
            0.5_dp*time_step*resistance
        rhs(:nq) = matmul(mass_q, q) - 0.5_dp*time_step*matmul(transpose(coupling), v)
        rhs(nq + 1:nq + nv) = matmul(mass_v, v) + 0.5_dp*time_step*matmul(coupling, q) + &
            0.5_dp*time_step*matmul(transpose(port), current)
        rhs(nq + nv + 1:total_size) = matmul(inductance - 0.5_dp*time_step*resistance, current) - &
            0.5_dp*time_step*matmul(port, v)
    end subroutine assemble_system

    subroutine assemble_system_dot( &
            mass_q, mass_v, inductance, resistance, coupling, port, time_step, q, v, current, &
            mass_q_dot, mass_v_dot, inductance_dot, resistance_dot, coupling_dot, port_dot, &
            time_step_dot, q_dot, v_dot, current_dot, matrix_dot, rhs_dot)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), inductance(:, :), resistance(:, :)
        real(dp), intent(in) :: coupling(:, :), port(:, :), time_step, q(:), v(:), current(:)
        real(dp), intent(in) :: mass_q_dot(:, :), mass_v_dot(:, :), inductance_dot(:, :)
        real(dp), intent(in) :: resistance_dot(:, :), coupling_dot(:, :), port_dot(:, :)
        real(dp), intent(in) :: time_step_dot, q_dot(:), v_dot(:), current_dot(:)
        real(dp), intent(out) :: matrix_dot(:, :), rhs_dot(:)
        integer :: nq, nv, ni, total_size

        nq = size(q)
        nv = size(v)
        ni = size(current)
        total_size = nq + nv + ni
        matrix_dot = 0.0_dp
        matrix_dot(:nq, :nq) = mass_q_dot
        matrix_dot(:nq, nq + 1:nq + nv) = 0.5_dp*(time_step_dot*transpose(coupling) + &
            time_step*transpose(coupling_dot))
        matrix_dot(nq + 1:nq + nv, :nq) = -0.5_dp*(time_step_dot*coupling + time_step*coupling_dot)
        matrix_dot(nq + 1:nq + nv, nq + 1:nq + nv) = mass_v_dot
        matrix_dot(nq + 1:nq + nv, nq + nv + 1:total_size) = -0.5_dp*(time_step_dot*transpose(port) + &
            time_step*transpose(port_dot))
        matrix_dot(nq + nv + 1:total_size, nq + 1:nq + nv) = 0.5_dp*(time_step_dot*port + time_step*port_dot)
        matrix_dot(nq + nv + 1:total_size, nq + nv + 1:total_size) = inductance_dot + &
            0.5_dp*(time_step_dot*resistance + time_step*resistance_dot)
        rhs_dot(:nq) = matmul(mass_q_dot, q) + matmul(mass_q, q_dot) - &
            0.5_dp*(time_step_dot*matmul(transpose(coupling), v) + &
            time_step*matmul(transpose(coupling_dot), v) + time_step*matmul(transpose(coupling), v_dot))
        rhs_dot(nq + 1:nq + nv) = matmul(mass_v_dot, v) + matmul(mass_v, v_dot) + &
            0.5_dp*(time_step_dot*(matmul(coupling, q) + matmul(transpose(port), current)) + &
            time_step*(matmul(coupling_dot, q) + matmul(coupling, q_dot) + &
            matmul(transpose(port_dot), current) + matmul(transpose(port), current_dot)))
        rhs_dot(nq + nv + 1:total_size) = matmul(inductance_dot - 0.5_dp*(time_step_dot*resistance + &
            time_step*resistance_dot), current) + matmul(inductance - 0.5_dp*time_step*resistance, current_dot) - &
            0.5_dp*(time_step_dot*matmul(port, v) + time_step*matmul(port_dot, v) + &
            time_step*matmul(port, v_dot))
    end subroutine assemble_system_dot

    logical function valid_inputs( &
            mass_q, mass_v, inductance, resistance, coupling, port, time_step, q, v, current) result(valid)
        real(dp), intent(in) :: mass_q(:, :), mass_v(:, :), inductance(:, :), resistance(:, :)
        real(dp), intent(in) :: coupling(:, :), port(:, :), time_step, q(:), v(:), current(:)
        integer :: nq, nv, ni

        nq = size(q)
        nv = size(v)
        ni = size(current)
        valid = time_step > 0.0_dp .and. nq > 0 .and. nv > 0 .and. ni > 0 .and. &
            all(shape(mass_q) == [nq, nq]) .and. all(shape(mass_v) == [nv, nv]) .and. &
            all(shape(inductance) == [ni, ni]) .and. all(shape(resistance) == [ni, ni]) .and. &
            all(shape(coupling) == [nv, nq]) .and. all(shape(port) == [ni, nv])
    end function valid_inputs

    logical function same_shapes(a, a_ref, b, b_ref, c, c_ref, d, d_ref, e, e_ref, f, f_ref) result(valid)
        real(dp), intent(in) :: a(:, :), a_ref(:, :), b(:, :), b_ref(:, :), c(:, :), c_ref(:, :)
        real(dp), intent(in) :: d(:, :), d_ref(:, :), e(:, :), e_ref(:, :), f(:, :), f_ref(:, :)

        valid = all(shape(a) == shape(a_ref)) .and. all(shape(b) == shape(b_ref)) .and. &
            all(shape(c) == shape(c_ref)) .and. all(shape(d) == shape(d_ref)) .and. &
            all(shape(e) == shape(e_ref)) .and. all(shape(f) == shape(f_ref))
    end function same_shapes

    pure function outer_product(first, second) result(product)
        real(dp), intent(in) :: first(:), second(:)
        real(dp) :: product(size(first), size(second))

        product = spread(first, 2, size(second))*spread(second, 1, size(first))
    end function outer_product

end module fortfem_mixed_wave_wall_time
