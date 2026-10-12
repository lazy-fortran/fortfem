program test_meridian_feec_exact_sequence
    !! Per-n twisted FEEC complex on curved meridional meshes:
    !!  - curl_n grad_n = 0 and div_n curl_n = 0 as matrix identities for
    !!    n = 0, 1, 3, 6, orders 1-3 and geometry orders 1-3,
    !!  - the operators represent the physical grad, curl and div with
    !!    their i n and 1/R terms (L2-projected analytic fields: the
    !!    discrete operator applied to the projection converges to the
    !!    analytic operator at the approximation order).
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_op_complex, meridian_zcsc_matmul, meridian_error_norm, &
        meridian_stiffness, meridian_load, meridian_real_to_complex
    use fortsparse, only: csc_conjugate_transpose, csc_matvec, csc_z_t, &
        fortsparse_status_t
    use meridian_test_fields, only: mode_n, scalar_f, grad_f, vec_a, curl_a, &
        vec_b, div_b, solve_z, zadd
    use meridian_test_meshes, only: make_annulus_mesh, &
        make_warped_rectangle_mesh
    implicit none

    integer, parameter :: ns(4) = [0, 1, 3, 6]
    logical :: all_passed

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    call identities()
    call physical_operators()
    call check_summary("Meridian per-n FEEC complex: exact sequence")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    function max_abs(z) result(m)
        type(csc_z_t), intent(in) :: z
        real(dp) :: m
        m = 0.0_dp
        if (z%nnz > 0) m = maxval(abs(z%val))
    end function max_abs

    subroutine identities()
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(csc_z_t) :: g, c, d, cg, dc
        character(len=80) :: tag
        integer :: geo, k, i, status, shape_id
        real(dp) :: e1, e2, s1, s2

        do shape_id = 1, 2
            do geo = 1, 3
                if (shape_id == 1) then
                    call make_annulus_mesh(2, 8, 0.4_dp, 1.0_dp, 3.0_dp, geo, mm)
                else
                    call make_warped_rectangle_mesh(3, geo, mm)
                end if
                do k = 1, 3
                    call meridian_feec_build(mm, k, fs, status)
                    write (tag, '(a,i0,a,i0,a,i0)') "shape ", shape_id, &
                        " geo ", geo, " order ", k
                    call record(status == 0, "build "//trim(tag))
                    if (status /= 0) cycle
                    e1 = 0.0_dp
                    e2 = 0.0_dp
                    s1 = 0.0_dp
                    s2 = 0.0_dp
                    do i = 1, size(ns)
                        call meridian_op_complex(fs%grad, ns(i), g, status)
                        call meridian_op_complex(fs%curl, ns(i), c, status)
                        call meridian_op_complex(fs%div, ns(i), d, status)
                        call meridian_zcsc_matmul(c, g, cg, status)
                        call meridian_zcsc_matmul(d, c, dc, status)
                        e1 = max(e1, max_abs(cg))
                        e2 = max(e2, max_abs(dc))
                        s1 = max(s1, max_abs(g))
                        s2 = max(s2, max_abs(c))
                    end do
                    call record(e1 < 1.0e-12_dp*max(1.0_dp, s1*s2), &
                        "curl_n grad_n = 0 for n = 0,1,3,6, "//trim(tag))
                    call record(e2 < 1.0e-12_dp*max(1.0_dp, s2), &
                        "div_n curl_n = 0 for n = 0,1,3,6, "//trim(tag))
                    call record(s1 > 0.5_dp .and. s2 > 0.5_dp, &
                        "operators are nontrivial, "//trim(tag))
                    print '(a,a,2es10.2)', "  residuals ", trim(tag), e1, e2
                end do
            end do
        end do
    end subroutine identities

    subroutine op_error(fs, space, op_sel, fun_in, fun_out, n, err)
        !! || op_n P_space(fun_in) - fun_out ||_{L2, space + 1}.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: space, op_sel, n
        interface
            subroutine fun_in(x, v)
                import :: dp
                real(dp), intent(in) :: x(2)
                complex(dp), intent(out) :: v(3)
            end subroutine fun_in
            subroutine fun_out(x, v)
                import :: dp
                real(dp), intent(in) :: x(2)
                complex(dp), intent(out) :: v(3)
            end subroutine fun_out
        end interface
        real(dp), intent(out) :: err
        type(csc_z_t) :: op, oph, kz, mz, lhs
        complex(dp), allocatable :: u(:), w(:), b1(:), b2(:), rhs(:)
        type(fortsparse_status_t) :: fst
        logical :: ok
        integer :: st

        mode_n = n
        select case (op_sel)
        case (0)
            call meridian_op_complex(fs%grad, n, op, st)
        case (1)
            call meridian_op_complex(fs%curl, n, op, st)
        case default
            call meridian_op_complex(fs%div, n, op, st)
        end select
        ! energy projection: (op^H M op + M) u = b + op^H b_out
        call meridian_stiffness(fs, op_sel, n, kz, st)
        call meridian_real_to_complex(fs%mass(space), mz)
        call zadd(kz, mz, lhs)
        call meridian_load(fs, space, fun_in, b1, st)
        call meridian_load(fs, space + 1, fun_out, b2, st)
        call csc_conjugate_transpose(op, oph, fst)
        rhs = b1 + csc_matvec(oph, b2)
        allocate (u(size(rhs)))
        call solve_z(lhs, rhs, u, ok)
        w = csc_matvec(op, u)
        call meridian_error_norm(fs, space + 1, w, fun_out, err, st)
    end subroutine op_error

    subroutine physical_operators()
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        real(dp) :: eg(2), ec(2), ed(2), rg, rc, rd
        integer :: k, i, lvl, status, nx(2)
        character(len=60) :: tag

        nx = [3, 6]
        do k = 1, 2
            do i = 1, size(ns)
                do lvl = 1, 2
                    call make_warped_rectangle_mesh(nx(lvl), 3, mm)
                    call meridian_feec_build(mm, k, fs, status)
                    if (status /= 0) error stop "build failed"
                    call op_error(fs, 0, 0, scalar_f, grad_f, ns(i), eg(lvl))
                    call op_error(fs, 1, 1, vec_a, curl_a, ns(i), ec(lvl))
                    call op_error(fs, 2, 2, vec_b, div_b, ns(i), ed(lvl))
                end do
                rg = log(eg(1)/eg(2))/log(2.0_dp)
                rc = log(ec(1)/ec(2))/log(2.0_dp)
                rd = log(ed(1)/ed(2))/log(2.0_dp)
                write (tag, '(a,i0,a,i0)') "order ", k, " n ", ns(i)
                print '(a,a,3f7.2,a,3es10.2)', "  rates grad/curl/div ", &
                    trim(tag), rg, rc, rd, "  fine errors", eg(2), ec(2), ed(2)
                call record(rg > real(k, dp) - 0.4_dp, &
                    "grad_n represents (grad f, i n f/R) at order, "//trim(tag))
                call record(rc > real(k, dp) - 0.4_dp, &
                    "curl_n represents the physical curl, "//trim(tag))
                call record(rd > real(k, dp) - 0.4_dp, &
                    "div_n represents (1/R) d_R (R B_R) + i n B_phi/R + d_Z B_Z, "// &
                    trim(tag))
            end do
        end do
    end subroutine physical_operators

end program test_meridian_feec_exact_sequence
