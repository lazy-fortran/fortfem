module spline_singular_kernel
    !! W = |grad u|^2/2 + u^2/2 for the pointwise-path checks.
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    implicit none
    private

    public :: plain_kernel_t

    type, extends(pointwise_energy_t) :: plain_kernel_t
    contains
        procedure :: evaluate
    end type plain_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(plain_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status
        integer :: q

        associate (unused => cell, u2 => weight, u3 => first_point, u4 => self)
        end associate
        do q = 1, npoint
            w(q) = 0.5_dp*sum(jet(1:3, q)**2)
            if (need >= 1) gradient(:, q) = jet(1:3, q)
            if (need >= 2) hessian(:, q) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
        end do
        status = 0
    end subroutine evaluate

end module spline_singular_kernel

program test_spline_singular
    !! Polar (Duffy) quadrature at a point singularity: on the cells at X the
    !! integral of 1/|s - X| becomes smooth, so the total error of the
    !! square around X drops by more than 20 (it stays O(h), set by the
    !! nearly singular neighbour cells). With 8 x 8 points per triangle the
    !! pointwise cells reproduce the tensor energy of a polynomial field,
    !! and Hessians and evaluation agree.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_quadrature_point_singular, spline_evaluate, spline_evaluate_at
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr, &
        spline_energy_hessian_apply
    use spline_singular_kernel, only: plain_kernel_t
    implicit none

    real(dp), parameter :: x0(2) = [0.3_dp, -0.2_dp]
    type(block_layout_t) :: layout
    type(spline_space_t) :: space
    type(spline_quadrature_t) :: quad, adapted
    real(dp) :: exact, a, standard(3), singular(3)
    integer :: status, level

    call layout_ogrid(x0, 1.0_dp, 0.5_dp, layout, status, polygon=.true.)
    a = 1.0_dp/sqrt(2.0_dp)
    exact = 8.0_dp*a*log(1.0_dp + sqrt(2.0_dp))
    do level = 1, 3
        call build(2**(level + 1))
        call spline_quadrature_build(space, 4, quad, status)
        call spline_quadrature_point_singular(space, quad, x0, adapted, status)
        if (status /= 0) error stop "singular"
        standard(level) = abs(inverse_distance(quad) - exact)
        singular(level) = abs(inverse_distance(adapted) - exact)
    end do
    print '(a,3es10.2,a,3es10.2)', "1/r errors, Gauss", standard, "  polar", singular
    call check_condition(all(singular < standard/20.0_dp) .and. &
        singular(3) < singular(1)/3.5_dp, &
        "polar rule integrates 1/|s - X| 20 times better, both errors O(h)")
    call build(4)
    call spline_quadrature_build(space, 16, quad, status)
    call spline_quadrature_point_singular(space, quad, x0, adapted, status)
    call check_condition(status == 0 .and. count(adapted%pointwise) == 4 .and. &
        pointwise_ok(), "pointwise cells reproduce energy, Hessians and evaluation")
    call spline_quadrature_build(space, 3, quad, status)
    call spline_quadrature_point_singular(space, quad, x0, adapted, status)
    call check_condition(status == 1, "odd point counts are rejected")
    call check_summary("spline_singular")

contains

    subroutine build(cells)
        integer, intent(in) :: cells
        type(spline_space_options_t) :: options
        integer :: status

        options%degree = 3
        allocate(options%knots(2, 5))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
        if (status /= 0) error stop "space"
    end subroutine build

    real(dp) function inverse_distance(q) result(total)
        type(spline_quadrature_t), intent(in) :: q
        integer :: i

        total = 0.0_dp
        do i = 1, q%point_count
            total = total + q%weight(i)/norm2(q%s(:, i) - x0)
        end do
    end function inverse_distance

    logical function pointwise_ok() result(ok)
        type(plain_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), v(:, :), hv(:, :), g1(:, :), g2(:, :)
        real(dp), allocatable :: jet(:, :), cv(:)
        real(dp) :: e1, e2, one(1)
        integer :: status, i, k, q

        allocate(c(1, space%dof_count), v(1, space%dof_count), hv(1, space%dof_count), &
            g1(1, space%dof_count), g2(1, space%dof_count), cv(space%dof_count))
        call random_number(c)
        call random_number(v)
        call spline_energy(space, quad, kernel, c, e1, status, g1)
        call spline_energy(space, adapted, kernel, c, e2, status, g2)
        ok = abs(e1 - e2) < 1.0e-12_dp*abs(e1) .and. &
            maxval(abs(g1 - g2)) < 1.0e-12_dp*maxval(abs(g1))
        call spline_energy_hessian_csr(space, adapted, kernel, c, matrix, status)
        call spline_energy_hessian_apply(space, adapted, kernel, c, v, hv, status)
        do i = 1, matrix%n
            cv(i) = 0.0_dp
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                cv(i) = cv(i) + matrix%val(k)*v(1, matrix%colidx(k))
            end do
        end do
        ok = ok .and. maxval(abs(cv - hv(1, :))) < 1.0e-11_dp*maxval(abs(cv))
        allocate(jet(1, adapted%point_count))
        call spline_evaluate(space, adapted, c, 0, jet, status)
        do q = 1, adapted%point_count, 7
            call spline_evaluate_at(space, c, adapted%block(q), adapted%xi(:, q), 0, &
                one, status)
            ok = ok .and. abs(one(1) - jet(1, q)) < 1.0e-12_dp
        end do
    end function pointwise_ok

end program test_spline_singular
