program test_meridian_feec_manufactured
    !! Manufactured solutions of the per-n Dirichlet problems on a curved
    !! annulus (geometry order = element order):
    !!   V0  -div_n(R grad_n f)/R + f = F        (resistive diffusion)
    !!   V1  curl_n curl_n A + A = F             (ideal/resistive operator)
    !!   V2  -grad_n div_n B + B = F             (grad-div on H(div_n))
    !! with exact solutions vanishing on the boundary.  The right-hand sides
    !! are formed by finite differences of the analytic fields, independent
    !! of the discrete operators.  Expected: L2 order k + 1 for V0 and k for
    !! V1, V2; energy (operator) order k.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_op_complex, meridian_error_norm, meridian_stiffness, &
        meridian_load, meridian_real_to_complex
    use fortsparse, only: csc_z_t, csc_matvec
    use meridian_test_fields, only: mode_n, solve_z, zadd, zero_rows, &
        fun_astar, fun_bstar, fun_f1, fun_bs, fun_dbs, fun_f2, fun_fs, &
        fun_gfs, fun_f0, ann_c, ann_r1, ann_r2
    use meridian_test_meshes, only: make_annulus_mesh
    implicit none

    logical :: all_passed

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    ann_c = 3.0_dp
    ann_r1 = 0.4_dp
    ann_r2 = 1.0_dp
    call run_space(0, 1, [0, 1, 3, 6], 3)
    call run_space(0, 2, [0, 3, 6], 3)
    call run_space(0, 3, [3], 2)
    call run_space(1, 1, [0, 1, 3, 6], 3)
    call run_space(1, 2, [0, 3, 6], 3)
    call run_space(1, 3, [3], 2)
    call run_space(2, 1, [0, 1, 3, 6], 3)
    call run_space(2, 2, [0, 3, 6], 3)
    call run_space(2, 3, [3], 2)
    call check_summary("Meridian per-n FEEC: manufactured solutions")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    subroutine solve_problem(fs, space, n, e_l2, e_op)
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: space, n
        real(dp), intent(out) :: e_l2, e_op
        type(csc_z_t) :: kz, mz, lhs, op
        complex(dp), allocatable :: rhs(:), u(:), w(:)
        logical :: ok
        integer :: st

        mode_n = n
        call meridian_stiffness(fs, space, n, kz, st)
        call meridian_real_to_complex(fs%mass(space), mz)
        call zadd(kz, mz, lhs)
        select case (space)
        case (0)
            call meridian_load(fs, 0, fun_f0, rhs, st)
        case (1)
            call meridian_load(fs, 1, fun_f1, rhs, st)
        case default
            call meridian_load(fs, 2, fun_f2, rhs, st)
        end select
        where (fs%sp(space)%bnd) rhs = (0.0_dp, 0.0_dp)
        call zero_rows(lhs, fs%sp(space)%bnd)
        allocate (u(size(rhs)))
        call solve_z(lhs, rhs, u, ok)
        if (.not. ok) error stop "solve failed"
        select case (space)
        case (0)
            call meridian_error_norm(fs, 0, u, fun_fs, e_l2, st)
            call meridian_op_complex(fs%grad, n, op, st)
            w = csc_matvec(op, u)
            call meridian_error_norm(fs, 1, w, fun_gfs, e_op, st)
        case (1)
            call meridian_error_norm(fs, 1, u, fun_astar, e_l2, st)
            call meridian_op_complex(fs%curl, n, op, st)
            w = csc_matvec(op, u)
            call meridian_error_norm(fs, 2, w, fun_bstar, e_op, st)
        case default
            call meridian_error_norm(fs, 2, u, fun_bs, e_l2, st)
            call meridian_op_complex(fs%div, n, op, st)
            w = csc_matvec(op, u)
            call meridian_error_norm(fs, 3, w, fun_dbs, e_op, st)
        end select
    end subroutine solve_problem

    subroutine run_space(space, k, ns, nlev)
        integer, intent(in) :: space, k, ns(:), nlev
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        real(dp) :: el(3), eo(3), r_l2, r_op, want_l2
        integer :: i, lev, status, nr, nth
        character(len=70) :: tag

        do i = 1, size(ns)
            do lev = 1, nlev
                nr = 2**(lev - 1)*2
                nth = 2**(lev - 1)*8
                call make_annulus_mesh(nr, nth, ann_r1, ann_r2, ann_c, k, mm)
                call meridian_feec_build(mm, k, fs, status)
                if (status /= 0) error stop "build failed"
                call solve_problem(fs, space, ns(i), el(lev), eo(lev))
            end do
            r_l2 = log(el(nlev - 1)/el(nlev))/log(2.0_dp)
            r_op = log(eo(nlev - 1)/eo(nlev))/log(2.0_dp)
            want_l2 = real(k, dp)
            if (space == 0) want_l2 = real(k + 1, dp)
            write (tag, '(a,i0,a,i0,a,i0)') "space V", space, " order ", k, &
                " n ", ns(i)
            print '(a,a,a,2f7.2,a,2es10.2)', "  ", trim(tag), " rates L2/op", &
                r_l2, r_op, "  fine errors", el(nlev), eo(nlev)
            call record(r_l2 > want_l2 - 0.35_dp, &
                "L2 convergence order, "//trim(tag))
            call record(r_op > real(k, dp) - 0.35_dp, &
                "operator (energy) convergence order, "//trim(tag))
        end do
    end subroutine run_space

end program test_meridian_feec_manufactured
