program test_direct_fourier_transform_ad
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: &
        direct_fourier_forward_jvp, direct_fourier_forward_vjp, &
        direct_fourier_plan_t, initialize_direct_fourier_plan
    use fortfem_kinds, only: dp, pi
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: n_sample = 5, n_mode = 4
    real(dp), parameter :: angles(n_sample) = &
        [0.13_dp, 0.91_dp, 2.02_dp, 3.77_dp, 5.61_dp]
    real(dp), parameter :: angle_tangent(n_sample) = &
        [0.21_dp, -0.17_dp, 0.08_dp, 0.31_dp, -0.12_dp]
    integer, parameter :: modes(n_mode) = [-3, -1, 0, 2]
    complex(dp), parameter :: samples(n_sample) = [ &
        cmplx(0.9_dp, -0.4_dp, dp), cmplx(-1.2_dp, 0.7_dp, dp), &
        cmplx(0.3_dp, 1.1_dp, dp), cmplx(1.7_dp, -0.2_dp, dp), &
        cmplx(-0.6_dp, -0.9_dp, dp)]
    complex(dp), parameter :: sample_tangent(n_sample) = [ &
        cmplx(-0.2_dp, 0.5_dp, dp), cmplx(0.4_dp, -0.3_dp, dp), &
        cmplx(0.8_dp, 0.1_dp, dp), cmplx(-0.7_dp, 0.6_dp, dp), &
        cmplx(0.2_dp, -0.4_dp, dp)]
    complex(dp), parameter :: coefficient_cotangent(n_mode) = [ &
        cmplx(0.6_dp, -0.8_dp, dp), cmplx(-0.3_dp, 0.2_dp, dp), &
        cmplx(1.1_dp, 0.4_dp, dp), cmplx(-0.5_dp, -0.7_dp, dp)]
    real(dp), parameter :: step = 2.0e-7_dp
    logical :: all_passed
    real(dp) :: lhs, rhs
    real(dp) :: angle_cotangent(n_sample)
    complex(dp) :: tangent(n_mode), finite_difference(n_mode)
    complex(dp) :: sample_cotangent(n_sample)
    type(direct_fourier_plan_t) :: plan, plus_plan, minus_plan, bad_plan
    type(fortsparse_status_t) :: status

    all_passed = .true.
    call initialize_direct_fourier_plan(plan, angles, modes, 3, status=status)
    call record(status%code == 0, &
        "nonuniform signed-mode transform plan initializes")

    call direct_fourier_forward_jvp(plan, samples, sample_tangent, &
        angle_tangent, tangent, status)
    call record(status%code == 0, &
        "Fourier JVP accepts finite sample and angle tangents")
    call independent_central_difference(finite_difference)
    call record(maxval(abs(tangent - finite_difference)) < 2.0e-8_dp, &
        "Fourier JVP agrees with an independent central-difference oracle")

    call direct_fourier_forward_vjp(plan, samples, coefficient_cotangent, &
        sample_cotangent, angle_cotangent, status)
    call record(status%code == 0, &
        "Fourier VJP accepts a finite coefficient cotangent")
    lhs = real(sum(conjg(coefficient_cotangent)*tangent), dp)
    rhs = real(sum(conjg(sample_cotangent)*sample_tangent), dp) + &
        sum(angle_cotangent*angle_tangent)
    call record(abs(lhs - rhs) < 2.0e-13_dp, &
        "Fourier JVP and VJP satisfy the real-part complex adjoint identity")

    call direct_fourier_forward_jvp(plan, samples, sample_tangent(1:4), &
        angle_tangent, tangent, status)
    call record(status%code /= 0 .and. all(tangent == cmplx(0.0_dp, 0.0_dp, dp)), &
        "Fourier JVP rejects an invalid sample-tangent shape and clears output")
    call direct_fourier_forward_vjp(plan, samples, coefficient_cotangent(1:3), &
        sample_cotangent, angle_cotangent, status)
    call record(status%code /= 0 .and. &
        all(sample_cotangent == cmplx(0.0_dp, 0.0_dp, dp)) .and. &
        all(angle_cotangent == 0.0_dp), &
        "Fourier VJP rejects an invalid cotangent shape and clears outputs")

    call initialize_direct_fourier_plan(bad_plan, &
        [angles(1), angles(2), angles(1) + 2.0_dp*pi, angles(4), angles(5)], &
        modes, 3, status=status)
    call record(status%code /= 0, &
        "Fourier plan rejects duplicate sample locations modulo one period")

    call check_summary("direct Fourier transform derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine independent_central_difference(result)
        complex(dp), intent(out) :: result(n_mode)

        integer :: mode, sample
        complex(dp) :: plus_value, minus_value
        real(dp) :: plus_phase, minus_phase, normalization

        call initialize_direct_fourier_plan(plus_plan, angles + step*angle_tangent, &
            modes, 3, status=status)
        call initialize_direct_fourier_plan(minus_plan, angles - step*angle_tangent, &
            modes, 3, status=status)
        call record(status%code == 0, &
            "central-difference perturbations preserve transform metadata")
        normalization = 1.0_dp/sqrt(real(n_sample, dp))
        do mode = 1, n_mode
            plus_value = cmplx(0.0_dp, 0.0_dp, dp)
            minus_value = cmplx(0.0_dp, 0.0_dp, dp)
            do sample = 1, n_sample
                plus_phase = real(modes(mode), dp)*( &
                    angles(sample) + step*angle_tangent(sample))
                minus_phase = real(modes(mode), dp)*( &
                    angles(sample) - step*angle_tangent(sample))
                plus_value = plus_value + normalization*cmplx( &
                    cos(plus_phase), -sin(plus_phase), dp)*( &
                    samples(sample) + step*sample_tangent(sample))
                minus_value = minus_value + normalization*cmplx( &
                    cos(minus_phase), -sin(minus_phase), dp)*( &
                    samples(sample) - step*sample_tangent(sample))
            end do
            result(mode) = (plus_value - minus_value)/(2.0_dp*step)
        end do
    end subroutine independent_central_difference

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_direct_fourier_transform_ad
