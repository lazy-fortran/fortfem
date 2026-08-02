program test_mixed_wave_pressure_displacement_gallery
    !! Independent oracle for the mixed pressure/displacement gallery.
    use check, only: check_condition, check_summary
    use fortfem_time, only: advance_mixed_wave_midpoint
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: n = 3, steps = 12
    real(dp), parameter :: h = 0.031_dp
    real(dp), parameter :: mass_q(n, n) = reshape([ &
        1.40_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.90_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.70_dp], [n, n])
    real(dp), parameter :: mass_v(n, n) = reshape([ &
        0.80_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.30_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.10_dp], [n, n])
    real(dp), parameter :: coupling(n, n) = reshape([ &
        0.95_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.35_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.75_dp], [n, n])
    real(dp), parameter :: q_initial(n) = [0.65_dp, -0.35_dp, 0.22_dp]
    real(dp), parameter :: v_initial(n) = [-0.18_dp, 0.42_dp, 0.31_dp]
    real(dp) :: q(n), v(n), q_exact(n), v_exact(n)
    real(dp) :: q_reverse(n), v_reverse(n)
    real(dp) :: map(2*n, 2*n), oracle(2*n, 2*n)
    real(dp) :: energy_metric(2*n, 2*n), symplectic_form(2*n, 2*n)
    real(dp) :: energy_initial, energy, maximum_error, maximum_drift
    real(dp) :: reversal_error, map_error, symplectic_error
    type(fortsparse_status_t) :: status
    integer :: step

    q = q_initial
    v = v_initial
    energy_initial = wave_energy(q, v)
    maximum_error = 0.0_dp
    maximum_drift = 0.0_dp
    do step = 1, steps
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, h, q, v, status)
        call check_condition(status%code == 0, &
            "mixed pressure/displacement midpoint accepts the state blocks")
        call independent_exact_state(real(step, dp)*h, q_exact, v_exact)
        maximum_error = max(maximum_error, maxval(abs(q - q_exact)))
        energy = wave_energy(q, v)
        maximum_drift = max(maximum_drift, &
            abs(energy/energy_initial - 1.0_dp))
    end do
    call check_condition(maximum_error < 4.0e-4_dp, &
        "midpoint state agrees with the independent modal solution")
    call check_condition(maximum_drift < 3.0e-13_dp, &
        "midpoint preserves the mixed pressure/displacement energy")

    q_reverse = q
    v_reverse = v
    do step = 1, steps
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, -h, q_reverse, v_reverse, status)
    end do
    reversal_error = max(maxval(abs(q_reverse - q_initial)), &
        maxval(abs(v_reverse - v_initial)))
    call check_condition(reversal_error < 4.0e-13_dp, &
        "mixed midpoint is reversible for pressure and displacement states")

    call build_api_map(map, status)
    call build_independent_cayley_map(oracle)
    map_error = maxval(abs(map - oracle))
    call check_condition(map_error < 4.0e-13_dp, &
        "gallery midpoint map matches the independent Cayley oracle")
    call build_energy_metric(energy_metric)
    call build_symplectic_form(symplectic_form)
    symplectic_error = maxval(abs(matmul(transpose(map), &
        matmul(symplectic_form, map)) - symplectic_form))
    call check_condition(symplectic_error < 5.0e-12_dp, &
        "mixed pressure/displacement map preserves its symplectic form")

    call check_summary("mixed pressure/displacement gallery oracle")
    if (maximum_error >= 4.0e-4_dp .or. maximum_drift >= 3.0e-13_dp .or. &
            reversal_error >= 4.0e-13_dp .or. map_error >= 4.0e-13_dp .or. &
            symplectic_error >= 5.0e-12_dp) error stop 1

contains

    subroutine independent_exact_state(current_time, q_value, v_value)
        real(dp), intent(in) :: current_time
        real(dp), intent(out) :: q_value(:), v_value(:)
        real(dp) :: frequency, scale_q, scale_v
        integer :: mode

        do mode = 1, n
            frequency = coupling(mode, mode)/sqrt( &
                mass_q(mode, mode)*mass_v(mode, mode))
            scale_q = sqrt(mass_v(mode, mode)/mass_q(mode, mode))
            scale_v = sqrt(mass_q(mode, mode)/mass_v(mode, mode))
            q_value(mode) = q_initial(mode)*cos(frequency*current_time) - &
                scale_q*v_initial(mode)*sin(frequency*current_time)
            v_value(mode) = v_initial(mode)*cos(frequency*current_time) + &
                scale_v*q_initial(mode)*sin(frequency*current_time)
        end do
    end subroutine independent_exact_state

    subroutine build_api_map(map_value, map_status)
        real(dp), intent(out) :: map_value(:, :)
        type(fortsparse_status_t), intent(out) :: map_status
        real(dp) :: q_work(n), v_work(n)
        integer :: column

        map_value = 0.0_dp
        do column = 1, 2*n
            q_work = 0.0_dp
            v_work = 0.0_dp
            if (column <= n) then
                q_work(column) = 1.0_dp
            else
                v_work(column - n) = 1.0_dp
            end if
            call advance_mixed_wave_midpoint( &
                mass_q, mass_v, coupling, h, q_work, v_work, map_status)
            map_value(:, column) = [q_work, v_work]
            if (map_status%code /= 0) return
        end do
    end subroutine build_api_map

    subroutine build_independent_cayley_map(map_value)
        real(dp), intent(out) :: map_value(:, :)
        real(dp) :: left(2*n, 2*n), right(2*n, 2*n)
        real(dp) :: identity(2*n, 2*n)

        left = 0.0_dp
        right = 0.0_dp
        left(:n, :n) = mass_q
        left(:n, n + 1:2*n) = 0.5_dp*h*transpose(coupling)
        left(n + 1:2*n, :n) = -0.5_dp*h*coupling
        left(n + 1:2*n, n + 1:2*n) = mass_v
        right(:n, :n) = mass_q
        right(:n, n + 1:2*n) = -0.5_dp*h*transpose(coupling)
        right(n + 1:2*n, :n) = 0.5_dp*h*coupling
        right(n + 1:2*n, n + 1:2*n) = mass_v
        identity = eye(2*n)
        call solve_matrix(left, right, map_value)
    end subroutine build_independent_cayley_map

    subroutine build_energy_metric(metric)
        real(dp), intent(out) :: metric(:, :)

        metric = 0.0_dp
        metric(:n, :n) = mass_q
        metric(n + 1:2*n, n + 1:2*n) = mass_v
    end subroutine build_energy_metric

    subroutine build_symplectic_form(form)
        real(dp), intent(out) :: form(:, :)
        integer :: mode

        form = 0.0_dp
        do mode = 1, n
            form(mode, n + mode) = mass_q(mode, mode)*mass_v(mode, mode)/ &
                coupling(mode, mode)
            form(n + mode, mode) = -form(mode, n + mode)
        end do
    end subroutine build_symplectic_form

    subroutine solve_matrix(matrix, rhs, solution)
        real(dp), intent(in) :: matrix(:, :), rhs(:, :)
        real(dp), intent(out) :: solution(:, :)
        real(dp) :: augmented(size(matrix, 1), size(matrix, 1) + size(rhs, 2))
        real(dp) :: row_work(size(augmented, 2)), pivot, factor
        integer :: dimension, row, column, pivot_row

        dimension = size(matrix, 1)
        augmented(:, :dimension) = matrix
        augmented(:, dimension + 1:) = rhs
        do column = 1, dimension
            pivot_row = column - 1 + maxloc( &
                abs(augmented(column:dimension, column)), dim=1)
            pivot = augmented(pivot_row, column)
            if (abs(pivot) < 1.0e-13_dp) error stop "oracle matrix singular"
            if (pivot_row /= column) then
                row_work = augmented(column, :)
                augmented(column, :) = augmented(pivot_row, :)
                augmented(pivot_row, :) = row_work
            end if
            augmented(column, :) = augmented(column, :)/pivot
            do row = 1, dimension
                if (row == column) cycle
                factor = augmented(row, column)
                augmented(row, :) = augmented(row, :) - &
                    factor*augmented(column, :)
            end do
        end do
        solution = augmented(:, dimension + 1:)
    end subroutine solve_matrix

    pure function eye(dimension) result(identity)
        integer, intent(in) :: dimension
        real(dp) :: identity(dimension, dimension)
        integer :: diagonal

        identity = 0.0_dp
        do diagonal = 1, dimension
            identity(diagonal, diagonal) = 1.0_dp
        end do
    end function eye

    pure function wave_energy(q_value, v_value) result(value)
        real(dp), intent(in) :: q_value(:), v_value(:)
        real(dp) :: value

        value = 0.5_dp*(dot_product(q_value, matmul(mass_q, q_value)) + &
            dot_product(v_value, matmul(mass_v, v_value)))
    end function wave_energy

end program test_mixed_wave_pressure_displacement_gallery
