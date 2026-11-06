program test_meridian_bspline_n
    !! Per-mode B-spline curl complex with FEM-BEM exterior closure.
    !!  1. curl_n grad_n = 0 and div_n curl_n = 0 as matrix identities, and
    !!     the V0 mass integrates R exactly.
    !!  2. The load of a curl-free field, volume plus surface term, vanishes.
    !!  3. A compactly supported divergence-free B (n = 2; n = 0 with net
    !!     toroidal flux) is recovered, zero outside its support, with no
    !!     exterior field and a vanishing gauge multiplier.
    !!  4. The n = 0 field of a circular filament matches its elliptic
    !!     integral potential.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_meridian_bspline_n, only: meridian_bspline_t, &
        meridian_bspline_build, meridian_bspline_load_curl, &
        meridian_bspline_load_value
    use fortfem_meridian_feec_n, only: meridian_op_complex, &
        meridian_zcsc_matmul
    use fortsparse, only: csc_z_t
    use meridian_bspline_fields, only: centre, gauge_gradient, &
        disc_quadrature, circle_quadrature, cross, manufactured_solve, &
        loop_solve
    implicit none

    logical :: all_passed

    all_passed = .true.
    call test_sequence()
    call test_curl_free_load()
    call test_manufactured()
    call test_loop()
    call check_summary("Meridian per-n B-spline complex")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    subroutine test_sequence()
        type(meridian_bspline_t) :: sp
        type(csc_z_t) :: g, c, d, cg, dc
        integer :: status, n, st

        call meridian_bspline_build(1.5_dp, 2.5_dp, -0.4_dp, 0.6_dp, 5, 4, 3, &
            sp, status)
        call record(status == 0, "build")
        do n = 0, 3, 3
            call meridian_op_complex(sp%grad, n, g, st)
            call meridian_op_complex(sp%curl, n, c, st)
            call meridian_op_complex(sp%div, n, d, st)
            call meridian_zcsc_matmul(c, g, cg, st)
            call meridian_zcsc_matmul(d, c, dc, st)
            call record(maxval(abs(cg%val)) < 1.0e-12_dp*maxval(abs(g%val)) &
                *maxval(abs(c%val)), "curl_n grad_n = 0")
            call record(maxval(abs(dc%val)) < 1.0e-12_dp*maxval(abs(d%val)) &
                *maxval(abs(c%val)), "div_n curl_n = 0")
        end do
        call record(abs(sum(sp%mass(0)%val) - 0.5_dp*(2.5_dp**2 - 1.5_dp**2)) &
            < 1.0e-13_dp, "V0 mass integrates R")
    end subroutine test_sequence

    subroutine test_curl_free_load()
        !! (grad_n chi, curl_n A')_V + <n x grad_n chi, A'>_dV = 0.
        type(meridian_bspline_t) :: sp
        real(dp), allocatable :: x(:, :), w(:), nrm(:, :)
        complex(dp), allocatable :: v(:, :), f(:), fb(:)
        complex(dp) :: g(3), nv(3)
        integer :: status, k

        call meridian_bspline_build(2.0_dp, 3.0_dp, -0.5_dp, 0.5_dp, 8, 8, 3, &
            sp, status)
        allocate (f(sp%n1), fb(sp%n1))
        f = (0.0_dp, 0.0_dp)
        fb = f
        call disc_quadrature(centre, 0.4_dp, 64, 128, x, w)
        allocate (v(3, size(w)))
        do k = 1, size(w)
            call gauge_gradient(x(:, k), 2, g)
            v(:, k) = g
        end do
        call meridian_bspline_load_curl(sp, 2, x, w, v, f, status)
        call meridian_bspline_load_curl(sp, 2, x, w, v, fb, status)
        call circle_quadrature(centre, 0.4_dp, 1024, x, w, nrm)
        deallocate (v)
        allocate (v(3, size(w)))
        do k = 1, size(w)
            call gauge_gradient(x(:, k), 2, g)
            nv = [cmplx(nrm(1, k), 0.0_dp, dp), (0.0_dp, 0.0_dp), &
                cmplx(nrm(2, k), 0.0_dp, dp)]
            v(:, k) = cross(nv, g)
        end do
        call meridian_bspline_load_value(sp, 2, x, w, v, f, status)
        call record(status == 0, "loads assembled")
        print '(a,es10.2)', "  curl-free load / volume term:", &
            maxval(abs(f))/maxval(abs(fb))
        ! limited by the polar rule across the knot lines
        call record(maxval(abs(f)) < 1.0e-4_dp*maxval(abs(fb)), &
            "curl-free source has no load")
    end subroutine test_curl_free_load

    subroutine test_manufactured()
        real(dp) :: err_b, err_out, ext, pmax
        integer :: status, n

        do n = 0, 2, 2
            call manufactured_solve(16, 3, n, .true., err_b, err_out, ext, &
                pmax, status)
            call record(status == 0, "manufactured solve")
            print '(a,i0,a,4es10.2)', "  n = ", n, &
                ": |B_h - B|, |B_h| outside, W_ext/W_int, |p|:", err_b, &
                err_out, ext, pmax
            call record(err_b < 2.5e-2_dp, "B recovered")
            call record(err_out < 3.0e-3_dp, "zero outside the support")
            call record(ext < 1.0e-9_dp, "no exterior field")
            ! the caller quadrature of the curl-free part limits |p|
            call record(pmax < 1.0e-4_dp, "gauge multiplier vanishes")
        end do
    end subroutine test_manufactured

    subroutine test_loop()
        real(dp) :: err
        integer :: status

        call loop_solve(16, 2, 2.5_dp, err, status)
        call record(status == 0, "filament solve")
        print '(a,es10.2)', "  filament: max |A_phi,h - A_phi| / max |A_phi|", &
            err
        call record(err < 3.0e-3_dp, "filament potential")
    end subroutine test_loop

end program test_meridian_bspline_n
