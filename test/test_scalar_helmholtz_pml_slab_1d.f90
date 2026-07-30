program test_scalar_helmholtz_pml_slab_1d
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_scalar_helmholtz_pml_slab_1d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wave_number = 3.0_dp
    real(dp), parameter :: physical_end = 1.0_dp
    real(dp), parameter :: outer_end = 2.0_dp
    real(dp), parameter :: sigma_max = 4.0_dp
    integer, parameter :: counts(2) = [41, 81]
    complex(dp), allocatable :: exact(:), solution(:)
    complex(dp) :: outgoing, incoming, amplitude_in, amplitude_out
    real(dp), allocatable :: nodes(:)
    real(dp) :: errors(2), reflection, expected_reflection
    integer :: i, level, status
    logical :: all_passed

    all_passed = .true.
    do level = 1, 2
        allocate(nodes(counts(level)))
        do i = 1, size(nodes)
            nodes(i) = outer_end*real(i - 1, dp)/real(size(nodes) - 1, dp)
        end do
        call solve_scalar_helmholtz_pml_slab_1d( &
            nodes, physical_end, wave_number, sigma_max, 2, &
            cmplx(1.0_dp, 0.0_dp, dp), solution, status)
        if (status /= 0) error stop "scalar Helmholtz slab PML solve failed"
        call exact_stretched_solution(nodes, exact)
        errors(level) = sqrt(sum(abs(solution - exact)**2))/ &
            sqrt(sum(abs(exact)**2))
        if (level == 2) then
            outgoing = exp(cmplx(0.0_dp, wave_number*nodes(2), dp))
            incoming = exp(cmplx(0.0_dp, -wave_number*nodes(2), dp))
            amplitude_in = (solution(2) - outgoing)/(incoming - outgoing)
            amplitude_out = cmplx(1.0_dp, 0.0_dp, dp) - amplitude_in
            reflection = abs(amplitude_in/amplitude_out)
            expected_reflection = exp( &
                -2.0_dp*sigma_max*(outer_end - physical_end)/3.0_dp)
            call record_condition( &
                abs(log(reflection/expected_reflection)) < 8.0e-2_dp, &
                "PML reflection matches the polynomial-stretch prediction")
        end if
        deallocate(nodes, solution, exact)
    end do

    call record_condition(errors(2) < 0.3_dp*errors(1), &
        "slab PML solution converges under uniform refinement")
    call record_condition(errors(2) < 2.0e-3_dp, &
        "slab PML matches the exact stretched-coordinate solution")
    call check_summary("Scalar Helmholtz slab PML")
    if (.not. all_passed) error stop 1

contains

    subroutine exact_stretched_solution(x, values)
        real(dp), intent(in) :: x(:)
        complex(dp), allocatable, intent(out) :: values(:)

        complex(dp) :: stretched_outer
        complex(dp) :: stretched_x
        integer :: node

        allocate(values(size(x)))
        stretched_outer = cmplx(outer_end, &
            sigma_max*(outer_end - physical_end)/ &
            (3.0_dp*wave_number), dp)
        do node = 1, size(x)
            if (x(node) <= physical_end) then
                stretched_x = cmplx(x(node), 0.0_dp, dp)
            else
                stretched_x = cmplx(x(node), &
                    sigma_max*(x(node) - physical_end)**3/ &
                    (3.0_dp*wave_number*(outer_end - physical_end)**2), dp)
            end if
            values(node) = sin(wave_number*(stretched_outer - stretched_x))/ &
                sin(wave_number*stretched_outer)
        end do
    end subroutine exact_stretched_solution

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_scalar_helmholtz_pml_slab_1d
