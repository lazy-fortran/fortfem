program test_constrained_newton
    !! fortfem_constrained_newton: KKT Newton by LDL^T with inertia.
    !!
    !! 1. An energy with a continuous symmetry (points in the plane, invariant
    !!    under a common rotation): without a constraint the Hessian at the
    !!    minimum has a zero eigenvalue; with the moment constraint orthogonal
    !!    to the rotation at the start the KKT matrix is regular, the
    !!    multiplier vanishes, the angle of the start is kept and the inertia
    !!    on the constraint has neither negative nor zero eigenvalues. The same
    !!    with the constraint following the state.
    !! 2. A quadratic saddle: the inertia on the tangent space of a
    !!    constraint is the KKT count minus the number of constraints, for a
    !!    constraint that removes the negative direction and for one that
    !!    does not; an inactive (zero) row does not count.
    !! 3. Fixed slots stay; a start off the constraint is moved onto it.
    !! 4. A curved stiff valley (a penalised circle with a weak pull along
    !!    it): steps of the energy test alone against the energy test with
    !!    relaxed steps, for a growing stiffness.
    !! 5. Invalid input is refused.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_constrained_newton, only: newton_constraints_t, newton_options_t, &
        newton_report_t, newton_constraints_init, newton_constraints_add, &
        constrained_newton, NEWTON_ACCEPT_ENERGY, NEWTON_ACCEPT_GRADIENT, &
        NEWTON_ACCEPT_ENERGY_RELAXED, NEWTON_CONVERGED, NEWTON_BAD_INPUT
    use fortsparse_csr, only: csr_t
    implicit none
    integer, parameter :: np = 4, n = 2*np
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: stiffness = 1.0_dp, saddle(4) = [1.0_dp, -2.0_dp, 3.0_dp, 0.5_dp]
    integer :: which = 1

    call symmetry()
    call inertia_on_constraints()
    call fixed_and_offset()
    call stiff_valley()
    call invalid()
    call check_summary("Constrained Newton")

contains

    subroutine dense_upper(h, hessian)
        !! Upper triangle of the dense symmetric h by rows, all entries.
        real(dp), intent(in) :: h(:, :)
        type(csr_t), intent(out) :: hessian
        integer :: i, j, k, m

        m = size(h, 1)
        hessian%nrow = m
        hessian%ncol = m
        allocate (hessian%row_ptr(m + 1), hessian%col(m*(m + 1)/2), &
                  hessian%val(m*(m + 1)/2))
        k = 0
        do i = 1, m
            hessian%row_ptr(i) = k + 1
            do j = i, m
                k = k + 1
                hessian%col(k) = j
                hessian%val(k) = h(i, j)
            end do
        end do
        hessian%row_ptr(m + 1) = k + 1
    end subroutine dense_upper

    subroutine energy_system(z, energy, gradient, hessian)
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: energy
        real(dp), intent(out), optional :: gradient(:)
        type(csr_t), intent(out), optional :: hessian
        real(dp) :: g(size(z)), h(size(z), size(z)), r2, d(2)
        integer :: i, j, a, b

        g = 0
        h = 0
        energy = 0
        select case (which)
        case (1)
            ! points u_i = z(2i-1:2i): sum (|u_i|^2 - 1)^2/4 + 0.15 sum_(i<j)
            ! i |u_i - u_j|^2, invariant under a common rotation of all points
            do i = 1, np
                a = 2*i - 1
                r2 = z(a)**2 + z(a + 1)**2
                energy = energy + 0.25_dp*(r2 - 1)**2
                g(a:a + 1) = g(a:a + 1) + (r2 - 1)*z(a:a + 1)
                h(a, a) = h(a, a) + r2 - 1 + 2*z(a)**2
                h(a + 1, a + 1) = h(a + 1, a + 1) + r2 - 1 + 2*z(a + 1)**2
                h(a, a + 1) = h(a, a + 1) + 2*z(a)*z(a + 1)
                h(a + 1, a) = h(a + 1, a) + 2*z(a)*z(a + 1)
                do j = i + 1, np
                    b = 2*j - 1
                    d = z(a:a + 1) - z(b:b + 1)
                    energy = energy + 0.15_dp*i*sum(d**2)
                    g(a:a + 1) = g(a:a + 1) + 0.3_dp*i*d
                    g(b:b + 1) = g(b:b + 1) - 0.3_dp*i*d
                    h(a, a) = h(a, a) + 0.3_dp*i
                    h(a + 1, a + 1) = h(a + 1, a + 1) + 0.3_dp*i
                    h(b, b) = h(b, b) + 0.3_dp*i
                    h(b + 1, b + 1) = h(b + 1, b + 1) + 0.3_dp*i
                    h(a, b) = h(a, b) - 0.3_dp*i
                    h(b, a) = h(b, a) - 0.3_dp*i
                    h(a + 1, b + 1) = h(a + 1, b + 1) - 0.3_dp*i
                    h(b + 1, a + 1) = h(b + 1, a + 1) - 0.3_dp*i
                end do
            end do
        case (2)
            ! quadratic saddle
            energy = 0.5_dp*sum(saddle*z**2)
            g = saddle*z
            do i = 1, size(z)
                h(i, i) = saddle(i)
            end do
        case (3)
            ! stiffness (|u|^2 - 1)^2/2 - x in the plane: a valley along the
            ! unit circle with its minimum on the positive x axis
            r2 = z(1)**2 + z(2)**2
            energy = 0.5_dp*stiffness*(r2 - 1)**2 - z(1)
            g = 2*stiffness*(r2 - 1)*z
            g(1) = g(1) - 1
            do i = 1, 2
                do j = 1, 2
                    h(i, j) = 4*stiffness*z(i)*z(j)
                end do
                h(i, i) = h(i, i) + 2*stiffness*(r2 - 1)
            end do
        end select
        if (present(gradient)) gradient = g
        if (present(hessian)) call dense_upper(h, hessian)
    end subroutine energy_system

    subroutine rotation_rows(z, c)
        !! The moment of the rotation generator (-y_i, x_i) at z, as one
        !! constraint through z.
        real(dp), intent(in) :: z(:)
        type(newton_constraints_t), intent(inout) :: c
        real(dp) :: row(n)
        integer :: i

        do i = 1, np
            row(2*i - 1) = -z(2*i)
            row(2*i) = z(2*i - 1)
        end do
        call newton_constraints_init(c)
        call newton_constraints_add(c, [(i, i=1, n)], row, dot_product(row, z))
    end subroutine rotation_rows

    subroutine symmetry()
        type(newton_constraints_t) :: c, none
        type(newton_options_t) :: options
        type(newton_report_t) :: report
        real(dp) :: z0(n), z(n), angle0, angle
        integer :: status, i, fixed(0), pass

        which = 1
        do i = 1, np
            z0(2*i - 1) = (0.7_dp + 0.1_dp*i)*cos(0.4_dp + 0.2_dp*i)
            z0(2*i) = (0.7_dp + 0.1_dp*i)*sin(0.4_dp + 0.2_dp*i)
        end do
        call newton_constraints_init(none)
        z = z0
        call constrained_newton(energy_system, z, fixed, none, options, report, status)
        call check_condition(status == NEWTON_CONVERGED, &
            "symmetric energy: converges without a constraint (shift)")
        call check_condition(report%zero_eigenvalues == 1 .or. &
            .not. report%inertia_reliable, &
            "symmetric energy: the free Hessian is singular along the orbit")
        do pass = 1, 2
            z = z0
            call rotation_rows(z0, c)
            if (pass == 1) then
                call constrained_newton(energy_system, z, fixed, c, options, report, status)
            else
                call constrained_newton(energy_system, z, fixed, c, options, report, status, &
                                        rotation_rows)
            end if
            call check_condition(status == NEWTON_CONVERGED .and. report%converged, &
                "symmetric energy: converges with the moment constraint")
            call check_condition(report%negative_eigenvalues == 0 .and. &
                report%zero_eigenvalues == 0 .and. report%inertia_reliable, &
                "symmetric energy: regular and positive on the constraint")
            call check_condition(report%multiplier_norm < 1.0e-9_dp, &
                "symmetric energy: the multiplier of a symmetry vanishes")
            ! all points coincide on the unit circle at the minimum
            call check_condition(abs(norm2(z(1:2)) - 1) < 1.0e-8_dp .and. &
                maxval(abs(z(3:n:2) - z(1))) < 1.0e-8_dp, &
                "symmetric energy: the minimum is reached")
            if (pass == 1) then
                ! the moment of the start fixes the angle
                angle0 = 0
                do i = 1, np
                    angle0 = angle0 + (-z0(2*i)*z0(2*i - 1) + z0(2*i - 1)*z0(2*i))
                end do
                angle = 0
                do i = 1, np
                    angle = angle + (-z0(2*i)*z(2*i - 1) + z0(2*i - 1)*z(2*i))
                end do
                call check_condition(abs(angle - angle0) < 1.0e-10_dp, &
                    "symmetric energy: the constraint holds at the solution")
            end if
        end do
    end subroutine symmetry

    subroutine inertia_on_constraints()
        type(newton_constraints_t) :: c
        type(newton_options_t) :: options
        type(newton_report_t) :: report
        real(dp) :: z(4)
        integer :: status, fixed(0)

        which = 2
        options%acceptance = NEWTON_ACCEPT_GRADIENT
        call newton_constraints_init(c)
        z = [0.3_dp, -0.2_dp, 0.5_dp, 0.1_dp]
        call constrained_newton(energy_system, z, fixed, c, options, report, status)
        call check_condition(status == NEWTON_CONVERGED .and. maxval(abs(z)) < 1.0e-12_dp &
            .and. report%negative_eigenvalues == 1 .and. report%steps == 1, &
            "saddle: Newton reaches it in one step and counts one negative eigenvalue")
        ! z_2 = 0 removes the negative direction
        call newton_constraints_add(c, [2], [1.0_dp], 0.0_dp)
        z = [0.3_dp, 0.0_dp, 0.5_dp, 0.1_dp]
        call constrained_newton(energy_system, z, fixed, c, options, report, status)
        call check_condition(status == NEWTON_CONVERGED .and. &
            report%negative_eigenvalues == 0 .and. report%zero_eigenvalues == 0, &
            "saddle: positive on a constraint across the negative direction")
        ! an inactive row and a constraint that keeps the negative direction
        call newton_constraints_init(c)
        call newton_constraints_add(c, [1, 3], [0.0_dp, 0.0_dp], 0.0_dp)
        call newton_constraints_add(c, [1, 3], [1.0_dp, -1.0_dp], 0.0_dp)
        z = [0.3_dp, -0.2_dp, 0.3_dp, 0.1_dp]
        call constrained_newton(energy_system, z, fixed, c, options, report, status)
        call check_condition(status == NEWTON_CONVERGED .and. &
            report%negative_eigenvalues == 1 .and. report%zero_eigenvalues == 0, &
            "saddle: one negative eigenvalue on the other constraint, zero row inactive")
    end subroutine inertia_on_constraints

    subroutine fixed_and_offset()
        type(newton_constraints_t) :: c
        type(newton_options_t) :: options
        type(newton_report_t) :: report
        real(dp) :: z(4)
        integer :: status

        which = 2
        saddle = [1.0_dp, 2.0_dp, 3.0_dp, 0.5_dp]
        call newton_constraints_init(c)
        ! z_1 + z_3 = 1 from a start off it; z_4 fixed
        call newton_constraints_add(c, [1, 3], [1.0_dp, 1.0_dp], 1.0_dp)
        z = [0.3_dp, -0.2_dp, 0.5_dp, 0.7_dp]
        call constrained_newton(energy_system, z, [4], c, options, report, status)
        ! minimum of z_1^2/2 + 3 z_3^2/2 on z_1 + z_3 = 1: z_1 = 3/4, z_3 = 1/4
        call check_condition(status == NEWTON_CONVERGED .and. &
            abs(z(1) - 0.75_dp) < 1.0e-8_dp .and. abs(z(3) - 0.25_dp) < 1.0e-8_dp .and. &
            abs(z(2)) < 1.0e-8_dp .and. z(4) == 0.7_dp, &
            "quadratic: the constrained minimum, the fixed slot kept")
        call check_condition(abs(report%multiplier_norm - 0.75_dp) < 1.0e-8_dp, &
            "quadratic: the multiplier is the constraint force")
        saddle = [1.0_dp, -2.0_dp, 3.0_dp, 0.5_dp]
    end subroutine fixed_and_offset

    subroutine stiff_valley()
        type(newton_constraints_t) :: c
        type(newton_options_t) :: options
        type(newton_report_t) :: report
        real(dp) :: z(2)
        integer :: status, fixed(0), k, steps(3, 2), fallback(3)

        which = 3
        call newton_constraints_init(c)
        options%max_steps = 3000
        do k = 1, 3
            stiffness = 10.0_dp**(2*k)
            options%acceptance = NEWTON_ACCEPT_ENERGY
            z = [cos(0.6_dp*pi), sin(0.6_dp*pi)]
            call constrained_newton(energy_system, z, fixed, c, options, report, status)
            call check_condition(status == NEWTON_CONVERGED .and. at_minimum(z), &
                "stiff valley: the energy test reaches the minimum")
            steps(k, 1) = report%steps
            options%acceptance = NEWTON_ACCEPT_ENERGY_RELAXED
            z = [cos(0.6_dp*pi), sin(0.6_dp*pi)]
            call constrained_newton(energy_system, z, fixed, c, options, report, status)
            call check_condition(status == NEWTON_CONVERGED .and. at_minimum(z), &
                "stiff valley: the relaxed test reaches the minimum")
            call check_condition(report%relaxed_steps > 0, &
                "stiff valley: steps are taken by the relaxed test")
            steps(k, 2) = report%steps
            fallback(k) = report%relaxed_steps
        end do
        print '(a,3i5,a,3i5,a,3i5)', 'stiff valley, stiffness 1e2, 1e4, 1e6: steps energy', &
            steps(:, 1), '  relaxed', steps(:, 2), '  of them relaxed pairs', fallback
        call check_condition(all(2*steps(2:, 2) < steps(2:, 1)), &
            "stiff valley: the relaxed test needs less than half the steps")

    end subroutine stiff_valley

    logical function at_minimum(x)
        !! y = 0 and 2 stiffness (x^2 - 1) x = 1, x > 0
        real(dp), intent(in) :: x(2)

        at_minimum = abs(x(2)) < 1.0e-6_dp .and. x(1) > 0 .and. &
            abs(2*stiffness*(x(1)**2 - 1)*x(1) - 1) < 1.0e-6_dp
    end function at_minimum

    subroutine invalid()
        type(newton_constraints_t) :: c
        type(newton_options_t) :: options
        type(newton_report_t) :: report
        real(dp) :: z(4)
        integer :: status, fixed(0)

        which = 2
        z = 0.1_dp
        call newton_constraints_init(c)
        call newton_constraints_add(c, [1, 1], [1.0_dp, 1.0_dp], 0.0_dp)
        call constrained_newton(energy_system, z, fixed, c, options, report, status)
        call check_condition(status == NEWTON_BAD_INPUT, "a slot twice in a row is refused")
        call newton_constraints_init(c)
        call constrained_newton(energy_system, z, [5], c, options, report, status)
        call check_condition(status == NEWTON_BAD_INPUT, "a fixed slot out of range is refused")
        options%acceptance = 7
        call constrained_newton(energy_system, z, fixed, c, options, report, status)
        call check_condition(status == NEWTON_BAD_INPUT, "an unknown acceptance is refused")
    end subroutine invalid

end program test_constrained_newton
