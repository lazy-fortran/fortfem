program test_spline_space_c1
    !! C1 (Argyris AS-G1) space on the polygonal 5-block O-grid: the
    !! dimension formula of Kapl, Sangalli, Takacs, full rank of the
    !! extraction, and continuity of values and gradients across every block
    !! side and at the valence-3 vertices for random coefficients.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, SPLINE_C1
    use fortfem_spline_eval, only: spline_evaluate_at
    implicit none

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

    type(block_layout_t) :: layout
    type(spline_space_t) :: space
    real(dp), allocatable :: coef(:, :)
    integer :: status, degree, cells, expected

    call layout_ogrid([0.0_dp, 0.0_dp], sqrt(2.0_dp), 0.5_dp, layout, status, &
        polygon=.true.)
    do degree = 3, 5
        do cells = 4, 6, 2
            call c1_space(degree, cells, space, status)
            expected = 5*(2*(cells - 1) + degree + 1 - 4)**2 + &
                12*(2*((cells - 1) + degree + 1) - 11) + 6*8
            call check_condition(status == 0 .and. space%dof_count == expected, &
                "C1 dimension equals the Kapl-Sangalli-Takacs count")
            if (status /= 0) print *, "status", status
            if (space%dof_count /= expected) print *, degree, cells, &
                space%dof_count, expected
        end do
    end do

    call c1_space(3, 4, space, status)
    call check_condition(status == 0 .and. extraction_rank_ok(space), &
        "C1 extraction has full row rank")
    allocate(coef(1, space%dof_count))
    call random_number(coef)
    call check_condition(max_jump(space, layout, coef) < 1.0e-11_dp, &
        "random C1 fields have continuous values and gradients across sides")
    call check_condition(vertex_gradient_spread(space, layout, coef) < 1.0e-11_dp, &
        "gradients agree at the valence-3 vertices from all blocks")
    call c0_space(3, 4, space, status)
    deallocate(coef)
    allocate(coef(1, space%dof_count))
    call random_number(coef)
    call check_condition(status == 0 .and. max_jump(space, layout, coef) > 1.0e-2_dp, &
        "the same check detects the gradient jumps of a C0 space")
    call check_summary("spline_space_c1")

contains

    subroutine c1_space(degree, cells, space, status)
        integer, intent(in) :: degree, cells
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status

        type(spline_space_options_t) :: options

        options%degree = degree
        options%regularity = SPLINE_C1
        allocate(options%knots(2, 5))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
    end subroutine c1_space

    subroutine c0_space(degree, cells, space, status)
        integer, intent(in) :: degree, cells
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status

        type(spline_space_options_t) :: options

        options%degree = degree
        allocate(options%knots(2, 5))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
    end subroutine c0_space

    logical function extraction_rank_ok(space) result(ok)
        type(spline_space_t), intent(in) :: space

        real(dp), allocatable :: e(:, :), gram(:, :), w(:), work(:)
        integer :: j, q, n, info

        n = space%dof_count
        allocate(e(n, space%broken_count))
        e = 0.0_dp
        do j = 1, space%broken_count
            do q = space%broken_ptr(j), space%broken_ptr(j + 1) - 1
                e(space%broken_dof(q), j) = e(space%broken_dof(q), j) + &
                    space%broken_weight(q)
            end do
        end do
        gram = matmul(e, transpose(e))
        allocate(w(n), work(10*n))
        call dsyev('N', 'U', n, gram, n, w, work, 10*n, info)
        ok = info == 0 .and. w(1) > 1.0e-10_dp*w(n)
    end function extraction_rank_ok

    real(dp) function max_jump(space, layout, coef) result(jump)
        type(spline_space_t), intent(in) :: space
        type(block_layout_t), intent(in) :: layout
        real(dp), intent(in) :: coef(:, :)

        real(dp) :: mine(3), theirs(3), t
        integer :: k, s, k2, s2, q, status

        jump = 0.0_dp
        do k = 1, layout%block_count
            do s = 1, 4
                k2 = layout%neighbour(s, k)
                if (k2 == 0) cycle
                s2 = layout%neighbour_side(s, k)
                do q = 0, 8
                    t = real(q, dp)/8.0_dp
                    call spline_evaluate_at(space, coef, k, side_xi(s, t), 1, &
                        mine, status)
                    if (layout%reversed(s, k)) then
                        call spline_evaluate_at(space, coef, k2, &
                            side_xi(s2, 1.0_dp - t), 1, theirs, status)
                    else
                        call spline_evaluate_at(space, coef, k2, side_xi(s2, t), 1, &
                            theirs, status)
                    end if
                    jump = max(jump, maxval(abs(mine - theirs))/ &
                        max(1.0_dp, maxval(abs(mine))))
                end do
            end do
        end do
    end function max_jump

    real(dp) function vertex_gradient_spread(space, layout, coef) result(spread)
        type(spline_space_t), intent(in) :: space
        type(block_layout_t), intent(in) :: layout
        real(dp), intent(in) :: coef(:, :)

        real(dp) :: jet(3), first(3)
        real(dp), parameter :: corner_xi(2, 4) = reshape([0.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [2, 4])
        integer :: v, k, c, status
        logical :: seen

        spread = 0.0_dp
        do v = 1, 4
            seen = .false.
            do k = 1, layout%block_count
                do c = 1, 4
                    if (layout%corners(c, k) /= v) cycle
                    call spline_evaluate_at(space, coef, k, corner_xi(:, c), 1, &
                        jet, status)
                    if (.not. seen) then
                        first = jet
                        seen = .true.
                    end if
                    spread = max(spread, maxval(abs(jet - first)))
                end do
            end do
        end do
    end function vertex_gradient_spread

    pure function side_xi(side, t) result(xi)
        integer, intent(in) :: side
        real(dp), intent(in) :: t
        real(dp) :: xi(2)

        select case (side)
        case (1)
            xi = [t, 0.0_dp]
        case (2)
            xi = [1.0_dp, t]
        case (3)
            xi = [t, 1.0_dp]
        case default
            xi = [0.0_dp, t]
        end select
    end function side_xi

end program test_spline_space_c1
