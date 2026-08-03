program test_circular_dtn
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        apply_circular_helmholtz_dtn, &
        circular_helmholtz_dtn_eigenvalue
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), parameter :: reference(4) = [ &
        (-0.22566209326700434_dp, 0.53649229362815990_dp), &
        (-0.33345841256565933_dp, 0.39593835603309263_dp), &
        (-1.9255971453720786_dp, 1.3723160047139202_dp), &
        (-4.6087886367050120_dp, 0.06400083907630437_dp)]
    integer, parameter :: modes(4) = [0, 1, 3, 8]
    real(dp), parameter :: radii(4) = [2.0_dp, 2.0_dp, 0.75_dp, 1.25_dp]
    real(dp), parameter :: arguments(4) = [1.0_dp, 1.0_dp, 2.5_dp, 5.0_dp]

    complex(dp) :: eigenvalue, negative_eigenvalue
    real(dp) :: relative_error
    integer :: case_id, status
    logical :: all_passed

    all_passed = .true.
    do case_id = 1, size(modes)
        call circular_helmholtz_dtn_eigenvalue( &
            modes(case_id), arguments(case_id) / radii(case_id), &
            radii(case_id), eigenvalue, status)
        relative_error = abs(eigenvalue - reference(case_id)) / &
            abs(reference(case_id))
        call record_condition(status == 0 .and. relative_error < 2.0e-13_dp, &
            "Circular DtN eigenvalue matches SciPy 1.18")
        call record_condition(aimag(eigenvalue) > 0.0_dp, &
            "Outgoing circular DtN mode has positive radiated flux")

        call circular_helmholtz_dtn_eigenvalue( &
            -modes(case_id), arguments(case_id) / radii(case_id), &
            radii(case_id), negative_eigenvalue, status)
        call record_condition(abs(negative_eigenvalue - eigenvalue) < &
            2.0e-13_dp * abs(eigenvalue), &
            "Positive and negative Fourier modes share one eigenvalue")
    end do

    call circular_helmholtz_dtn_eigenvalue( &
        0, 0.0_dp, 1.0_dp, eigenvalue, status)
    call record_condition(status /= 0, &
        "Circular DtN rejects zero wavenumber")
    call circular_helmholtz_dtn_eigenvalue( &
        0, 1.0_dp, 0.0_dp, eigenvalue, status)
    call record_condition(status /= 0, &
        "Circular DtN rejects zero radius")
    call test_modal_application()

    call check_summary("Circular Helmholtz DtN")
    if (.not. all_passed) error stop 1

contains

    subroutine test_modal_application()
        integer, parameter :: point_count = 13
        complex(dp), parameter :: lambda2 = &
            (-0.9338381358297664_dp, 1.2589265099477043_dp)
        complex(dp), parameter :: lambda5 = &
            (-4.4432429649354090_dp, 0.00644848493650646_dp)
        complex(dp) :: derivative(point_count), expected(point_count)
        complex(dp) :: mode2, mode5, trace(point_count)
        real(dp) :: angle, discarded_relative_norm
        integer :: point

        do point = 1, point_count
            angle = 2.0_dp * acos(-1.0_dp) * &
                real(point - 1, dp) / real(point_count, dp)
            mode2 = exp(cmplx(0.0_dp, 2.0_dp * angle, dp))
            mode5 = exp(cmplx(0.0_dp, 5.0_dp * angle, dp))
            trace(point) = 3.0_dp * mode2 + 4.0_dp * mode5
            expected(point) = 3.0_dp * lambda2 * mode2
        end do

        call apply_circular_helmholtz_dtn( &
            trace, 2.0_dp, 1.0_dp, 2, derivative, &
            discarded_relative_norm, status)
        call record_condition(status == 0 .and. &
            maxval(abs(derivative - expected)) < 2.0e-12_dp, &
            "Truncated circular DtN retains the exact low mode")
        call record_condition( &
            abs(discarded_relative_norm - 0.8_dp) < 2.0e-14_dp, &
            "Circular DtN reports the exact discarded trace norm")

        do point = 1, point_count
            angle = 2.0_dp * acos(-1.0_dp) * &
                real(point - 1, dp) / real(point_count, dp)
            mode2 = exp(cmplx(0.0_dp, 2.0_dp * angle, dp))
            mode5 = exp(cmplx(0.0_dp, 5.0_dp * angle, dp))
            expected(point) = &
                3.0_dp * lambda2 * mode2 + 4.0_dp * lambda5 * mode5
        end do
        call apply_circular_helmholtz_dtn( &
            trace, 2.0_dp, 1.0_dp, 6, derivative, &
            discarded_relative_norm, status)
        call record_condition(status == 0 .and. &
            maxval(abs(derivative - expected)) < 5.0e-12_dp, &
            "Untruncated circular DtN applies all resolvable modes")
        call record_condition(discarded_relative_norm < 2.0e-15_dp, &
            "Untruncated circular DtN reports zero discarded norm")
    end subroutine test_modal_application

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_circular_dtn
