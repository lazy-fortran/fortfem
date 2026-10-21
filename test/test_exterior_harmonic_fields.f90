program test_exterior_harmonic_fields
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    use fortfem_exterior_harmonic_fields, only: evaluate_current_loops, &
                                                prescribe_exterior_periods
    use fortfem_generated_current_segment, only: current_segment
    use fortfem_generated_current_segment_jvp, only: current_segment_jvp
    use fortfem_generated_current_segment_vjp, only: current_segment_vjp
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: nv = 512, nc = 512
    real(dp) :: wire(3, nv), x(3, nc), a(3, nc, 1), b(3, nc, 1)
    real(dp) :: periods(2, 1), target(2), coeff(1), t, theta, dl(3), circ
    real(dp) :: near_field, far_field, expected
    integer :: i, status
    logical :: passed
    passed = .true.
    do i = 1, nv
        t = 2*pi*real(i - 1, dp)/nv
        wire(:, i) = [2*cos(t), 2*sin(t), 0.0_dp]
    end do
    x = 0
    x(3, :) = 1
    call evaluate_current_loops(wire, [1, nv + 1], x, a, b, status)
    expected = 4.0_dp/(2*5.0_dp**1.5_dp)
    call record(status == 0, 'closed-loop exterior field construction')
    call record(abs(b(3, 1, 1)/expected - 1) < 2.e-5_dp, &
                'independent circular-loop axial field oracle')
    call record(maxval(abs(b(1:2, 1, 1))) < 1.e-14_dp, 'axial transverse symmetry')
    ! Meridian oriented by increasing theta has linking number -1 for +phi.
    do i = 1, nc
        theta = 2*pi*(real(i, dp) - 0.5_dp)/nc
        x(:, i) = [2 + 0.6_dp*cos(theta), 0.0_dp, 0.6_dp*sin(theta)]
    end do
    call evaluate_current_loops(wire, [1, nv + 1], x, a, b, status)
    circ = 0
    do i = 1, nc
        theta = 2*pi*(real(i, dp) - 0.5_dp)/nc
        dl = [-0.6_dp*sin(theta), 0.0_dp, 0.6_dp*cos(theta)]*(2*pi/nc)
        circ = circ + dot_product(b(:, i, 1), dl)
    end do
    periods(1, 1) = circ
    call record(abs(circ + 1) < 2.e-10_dp, 'independent poloidal linking period')
    ! Unlinked longitude has zero period. Both fields and cycles are external.
    do i = 1, nc
        theta = 2*pi*(real(i, dp) - 0.5_dp)/nc
        x(:, i) = [2.6_dp*cos(theta), 2.6_dp*sin(theta), 0.0_dp]
    end do
    call evaluate_current_loops(wire, [1, nv + 1], x, a, b, status)
    circ = 0
    do i = 1, nc
        theta = 2*pi*(real(i, dp) - 0.5_dp)/nc
        dl = [-2.6_dp*sin(theta), 2.6_dp*cos(theta), 0.0_dp]*(2*pi/nc)
        circ = circ + dot_product(b(:, i, 1), dl)
    end do
    periods(2, 1) = circ
    call record(abs(circ) < 1.e-13_dp, 'unlinked toroidal period is zero')
    target = [3.0_dp, 0.0_dp]
    call prescribe_exterior_periods(periods, target, [1], 1.e-8_dp, coeff, status)
    call record(status == 0 .and. abs(coeff(1) + 3) < 1.e-8_dp, &
                'prescribed compatible poloidal and toroidal data')
    target(2) = 1
    call prescribe_exterior_periods(periods, target, [1], 1.e-8_dp, coeff, status)
    call record(status == 3, 'incompatible exterior periods are rejected')
    x = 0; x(3, :) = 100
    call evaluate_current_loops(wire, [1, nv + 1], x, a, b, status)
    near_field = b(3, 1, 1)
    x(3, :) = 200
    call evaluate_current_loops(wire, [1, nv + 1], x, a, b, status)
    far_field = b(3, 1, 1)
    call record(abs(near_field/far_field - 8) < 0.004_dp, 'dipole decay at infinity')
    x(:, 1) = wire(:, 1)
    call evaluate_current_loops(wire, [1, nv + 1], x, a, b, status)
    call record(status /= 0, 'wire singularity returns status before evaluation')
    call derivatives()
    call multiple_handles()
    call check_summary('Exterior loop harmonic fields')
    if (.not. passed) error stop 1
contains
    subroutine derivatives()
        real(dp) :: v(6), d(6), jd(6), bars(6), z(6), yp(6), ym(6), vp(6), vm(6)
        real(dp), parameter :: step = 1.e-5_dp
        v = [1.2_dp, -0.3_dp, 0.6_dp, 0.8_dp, -0.5_dp, 0.4_dp]
        d = [0.1_dp, -0.2_dp, 0.3_dp, 0.2_dp, 0.1_dp, -0.1_dp]
        z = [0.2_dp, 0.3_dp, -0.1_dp, 0.4_dp, -0.2_dp, 0.5_dp]
        call current_segment_jvp(v(1), v(2), v(3), v(4), v(5), v(6), &
           d(1), d(2), d(3), d(4), d(5), d(6), jd(1), jd(2), jd(3), jd(4), jd(5), jd(6))
        call current_segment_vjp(v(1), v(2), v(3), v(4), v(5), v(6), &
                                 z(1), z(2), z(3), z(4), z(5), z(6), &
                                 bars(1), bars(2), bars(3), bars(4), bars(5), bars(6))
        call record(abs(dot_product(z, jd) - dot_product(bars, d)) < 1.e-12_dp, &
                    'segment JVP/VJP duality')
        vp = v + step*d; vm = v - step*d
        call current_segment(vp(1), vp(2), vp(3), vp(4), vp(5), vp(6), &
                             yp(1), yp(2), yp(3), yp(4), yp(5), yp(6))
        call current_segment(vm(1), vm(2), vm(3), vm(4), vm(5), vm(6), &
                             ym(1), ym(2), ym(3), ym(4), ym(5), ym(6))
        call record(maxval(abs((yp - ym)/(2*step) - jd)) < 1.e-8_dp, &
                    'segment derivatives match independent differences')
    end subroutine
    subroutine multiple_handles()
        real(dp) :: wires(3, 2*nv), points(3, nc), ap(3, nc, 2), bp(3, nc, 2)
        real(dp) :: matrix(2, 2), c(2), angle, delta(3)
        integer :: handle, j, h
        wires(:, :nv) = wire; wires(:, nv + 1:) = wire
        wires(1, nv + 1:) = wires(1, nv + 1:) + 8
        do handle = 1, 2
            do j = 1, nc
                angle = 2*pi*(real(j, dp) - 0.5_dp)/nc
                points(:, j) = [2 + 8*real(handle - 1, dp) + 0.6_dp*cos(angle), &
                                0.0_dp, 0.6_dp*sin(angle)]
            end do
       call evaluate_current_loops(wires, [1, nv + 1, 2*nv + 1], points, ap, bp, status)
            if (status /= 0) error stop 'two handle evaluation'
            matrix(handle, :) = 0
            do j = 1, nc
                angle = 2*pi*(real(j, dp) - 0.5_dp)/nc
                delta = [-0.6_dp*sin(angle), 0.0_dp, 0.6_dp*cos(angle)]*(2*pi/nc)
                do h = 1, 2
                 matrix(handle, h) = matrix(handle, h) + dot_product(bp(:, j, h), delta)
                end do
            end do
        end do
        call prescribe_exterior_periods(matrix, [2.0_dp, -3.0_dp], [1, 2], &
                                        1.e-8_dp, c, status)
        call record(status == 0 .and. maxval(abs(c - [-2.0_dp, 3.0_dp])) < 1.e-8_dp, &
                    'independent linking cycles prescribe two exterior handle periods')
    end subroutine
    subroutine record(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok, label)
    end subroutine
end program
