program test_meridian_feec_cohomology
    !! Dimensions of the cohomology of the per-n complex, from the ranks of
    !! the assembled operators.  The cohomology of the axisymmetric domain
    !! times a circle is carried by n = 0; every n /= 0 complex is acyclic
    !! (twisted by a nontrivial character).
    !!   annulus (hole): natural b = (1, 2, 1, 0), relative (0, 1, 2, 1)
    !!   topological disk: natural (1, 1, 0, 0), relative (0, 0, 1, 1)
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_op_complex
    use fortsparse, only: csc_z_t
    use meridian_test_meshes, only: make_annulus_mesh, &
        make_warped_rectangle_mesh
    implicit none

    interface
        subroutine zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
                work, lwork, rwork, info)
            import :: dp
            character(len=1), intent(in) :: jobu, jobvt
            integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
            complex(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: s(*)
            complex(dp), intent(out) :: u(ldu, *), vt(ldvt, *), work(*)
            real(dp), intent(out) :: rwork(*)
            integer, intent(out) :: info
        end subroutine zgesvd
    end interface

    logical :: all_passed
    integer :: k, geo

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    do k = 1, 3
        do geo = 1, 3, 2
            call run_domain(.true., k, geo)
            call run_domain(.false., k, geo)
        end do
    end do
    call check_summary("Meridian per-n FEEC: harmonic field dimensions")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    function dense_restricted(z, keep_row, keep_col) result(a)
        type(csc_z_t), intent(in) :: z
        logical, intent(in) :: keep_row(:), keep_col(:)
        complex(dp), allocatable :: a(:, :)
        integer, allocatable :: ri(:), ci(:)
        integer :: i, j, p, nr, nc

        allocate (ri(z%nrow), ci(z%ncol))
        ri = 0
        ci = 0
        nr = 0
        nc = 0
        do i = 1, z%nrow
            if (keep_row(i)) then
                nr = nr + 1
                ri(i) = nr
            end if
        end do
        do j = 1, z%ncol
            if (keep_col(j)) then
                nc = nc + 1
                ci(j) = nc
            end if
        end do
        allocate (a(nr, nc))
        a = (0.0_dp, 0.0_dp)
        do j = 1, z%ncol
            do p = z%col_ptr(j), z%col_ptr(j + 1) - 1
                if (ci(j) > 0 .and. ri(z%row_idx(p)) > 0) &
                    a(ri(z%row_idx(p)), ci(j)) = z%val(p)
            end do
        end do
    end function dense_restricted

    function leaked(z, keep_row, keep_col) result(m)
        !! Largest entry mapping retained columns into dropped rows.
        type(csc_z_t), intent(in) :: z
        logical, intent(in) :: keep_row(:), keep_col(:)
        real(dp) :: m
        integer :: j, p

        m = 0.0_dp
        do j = 1, z%ncol
            if (.not. keep_col(j)) cycle
            do p = z%col_ptr(j), z%col_ptr(j + 1) - 1
                if (.not. keep_row(z%row_idx(p))) m = max(m, abs(z%val(p)))
            end do
        end do
    end function leaked

    function matrix_rank(a) result(r)
        complex(dp), intent(inout) :: a(:, :)
        integer :: r
        real(dp), allocatable :: s(:), rwork(:)
        complex(dp), allocatable :: work(:), u(:, :), vt(:, :)
        complex(dp) :: wq(1)
        integer :: m, n, info

        m = size(a, 1)
        n = size(a, 2)
        r = 0
        if (m == 0 .or. n == 0) return
        allocate (s(min(m, n)), rwork(5*min(m, n)), u(1, 1), vt(1, 1))
        call zgesvd('N', 'N', m, n, a, m, s, u, 1, vt, 1, wq, -1, rwork, info)
        allocate (work(max(1, int(real(wq(1))))))
        call zgesvd('N', 'N', m, n, a, m, s, u, 1, vt, 1, work, size(work), &
            rwork, info)
        r = count(s > 1.0e-8_dp*s(1))
    end function matrix_rank

    subroutine run_domain(annulus, k, geo)
        logical, intent(in) :: annulus
        integer, intent(in) :: k, geo
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(csc_z_t) :: d0, d1, d2
        complex(dp), allocatable :: a0(:, :), a1(:, :), a2(:, :)
        logical, allocatable :: all0(:), all1(:), all2(:), all3(:)
        logical, allocatable :: in0(:), in1(:), in2(:), in3(:)
        integer, parameter :: ns(4) = [0, 1, 3, 6]
        integer :: n(4), expect(4), i, rel, status, b(0:3), r(0:2), idx
        character(len=90) :: tag
        real(dp) :: leak

        if (annulus) then
            call make_annulus_mesh(2, 8, 0.4_dp, 1.0_dp, 3.0_dp, geo, mm)
        else
            call make_warped_rectangle_mesh(3, geo, mm)
        end if
        call meridian_feec_build(mm, k, fs, status)
        if (status /= 0) error stop "build failed"
        n = [fs%sp(0)%ndof, fs%sp(1)%ndof, fs%sp(2)%ndof, fs%sp(3)%ndof]
        allocate (all0(n(1)), all1(n(2)), all2(n(3)), all3(n(4)))
        all0 = .true.
        all1 = .true.
        all2 = .true.
        all3 = .true.
        in0 = .not. fs%sp(0)%bnd
        in1 = .not. fs%sp(1)%bnd
        in2 = .not. fs%sp(2)%bnd
        in3 = .not. fs%sp(3)%bnd
        do rel = 0, 1
            do i = 1, size(ns)
                call meridian_op_complex(fs%grad, ns(i), d0, status)
                call meridian_op_complex(fs%curl, ns(i), d1, status)
                call meridian_op_complex(fs%div, ns(i), d2, status)
                if (rel == 0) then
                    a0 = dense_restricted(d0, all1, all0)
                    a1 = dense_restricted(d1, all2, all1)
                    a2 = dense_restricted(d2, all3, all2)
                    leak = 0.0_dp
                else
                    a0 = dense_restricted(d0, in1, in0)
                    a1 = dense_restricted(d1, in2, in1)
                    a2 = dense_restricted(d2, in3, in2)
                    leak = max(leaked(d0, in1, in0), leaked(d1, in2, in1), &
                        leaked(d2, in3, in2))
                end if
                r(0) = matrix_rank(a0)
                r(1) = matrix_rank(a1)
                r(2) = matrix_rank(a2)
                idx = 0
                if (rel == 0) then
                    b(0) = n(1) - r(0)
                    b(1) = n(2) - r(1) - r(0)
                    b(2) = n(3) - r(2) - r(1)
                    b(3) = n(4) - r(2)
                else
                    b(0) = count(in0) - r(0)
                    b(1) = count(in1) - r(1) - r(0)
                    b(2) = count(in2) - r(2) - r(1)
                    b(3) = count(in3) - r(2)
                end if
                expect = 0
                if (ns(i) == 0) then
                    if (annulus .and. rel == 0) expect = [1, 2, 1, 0]
                    if (annulus .and. rel == 1) expect = [0, 1, 2, 1]
                    if (.not. annulus .and. rel == 0) expect = [1, 1, 0, 0]
                    if (.not. annulus .and. rel == 1) expect = [0, 0, 1, 1]
                end if
                write (tag, '(a,l1,a,i0,a,i0,a,i0,a,i0)') "annulus ", annulus, &
                    " order ", k, " geo ", geo, " n ", ns(i), " relative ", rel
                call record(all(b == expect), "Betti numbers "// &
                    trim(tag))
                call record(leak < 1.0e-12_dp, "boundary conditions preserved "// &
                    trim(tag))
            end do
        end do
    end subroutine run_domain

end program test_meridian_feec_cohomology
