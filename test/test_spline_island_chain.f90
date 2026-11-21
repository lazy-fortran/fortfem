module spline_island_chain_problem
    !! Island chain of psi = (r - r0)^2/2 - eps cos(m theta) in polar
    !! coordinates (r, theta) about the axis (r_axis, 0): O-points at
    !! theta = 2 pi i/m, X-points at theta = (2i - 1) pi/m, separatrix
    !! r = r0 -+ 2 sqrt(eps) cos(m (theta - theta_O)/2) (inner, outer
    !! branch), elliptic wall about the axis.
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    use fortfem_spline_separatrix, only: island_chain_curves_t
    implicit none
    private

    public :: projection_kernel_t, chain, smooth

    real(dp), parameter, public :: r_axis = 1.5_dp, r0 = 0.5_dp
    real(dp), parameter, public :: wall_a = 0.9_dp, wall_b = 1.0_dp

    type, extends(pointwise_energy_t) :: projection_kernel_t
        !! w = u^2/2 - f u: the Hessian is the mass matrix.
        real(dp), allocatable :: f(:)
    contains
        procedure :: evaluate
    end type projection_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(projection_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status

        real(dp) :: f
        integer :: q

        associate (unused => cell, unused_weight => weight)
        end associate
        do q = 1, npoint
            f = self%f(first_point + q - 1)
            w(q) = 0.5_dp*jet(1, q)**2 - f*jet(1, q)
            if (need >= 1) gradient(:, q) = [jet(1, q) - f, 0.0_dp, 0.0_dp]
            if (need >= 2) hessian(:, q) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        end do
        status = 0
    end subroutine evaluate

    pure real(dp) function smooth(s) result(u)
        real(dp), intent(in) :: s(2)

        u = cos(2.0_dp*s(1))*sin(1.7_dp*s(2)) + s(1)*s(2)
    end function smooth

    function chain(m, half_width) result(curves)
        !! half_width = 2 sqrt(eps), the island half width at the O-point.
        integer, intent(in) :: m
        real(dp), intent(in) :: half_width
        type(island_chain_curves_t) :: curves

        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: theta, theta_o, r
        integer :: i, j, k
        integer, parameter :: n = 201

        curves%axis = [r_axis, 0.0_dp]
        allocate(curves%xpoint(2, m), curves%opoint(2, m), curves%inner(m), &
            curves%outer(m))
        do i = 1, m
            curves%xpoint(:, i) = point(r0, (2*i - 1)*pi/m)
            curves%opoint(:, i) = point(r0, 2*i*pi/m)
            theta_o = 2*i*pi/m
            allocate(curves%inner(i)%points(2, n), curves%outer(i)%points(2, n))
            do j = 1, n
                theta = theta_o + pi/m*(2.0_dp*(j - 1)/(n - 1) - 1.0_dp)
                do k = 1, 2
                    r = r0 + (2*k - 3)*half_width*cos(0.5_dp*m*(theta - theta_o))
                    if (k == 1) then
                        curves%inner(i)%points(:, j) = point(r, theta)
                    else
                        curves%outer(i)%points(:, j) = point(r, theta)
                    end if
                end do
            end do
        end do
        allocate(curves%wall(2, 400))
        do j = 1, 400
            theta = 2.0_dp*pi*(j - 1)/400
            curves%wall(:, j) = [r_axis + wall_a*cos(theta), wall_b*sin(theta)]
        end do

    contains

        pure function point(radius, angle) result(s)
            real(dp), intent(in) :: radius, angle
            real(dp) :: s(2)

            s = [r_axis + radius*cos(angle), radius*sin(angle)]
        end function point

    end function chain

end module spline_island_chain_problem

program test_spline_island_chain
    !! Island chain layouts (doc/design/multipatch_spline_spaces.md,
    !! separatrix layouts) for m = 2, 3, 5: block counts, X-points of
    !! valence 8, the axis as a vertex of valence m (m >= 3), a regular core
    !! quad (m = 2), order p + 1 for a smooth L2 projection, a jump across
    !! island 1 resolved only in the cut space, conforming graded knots and
    !! a moved chain on the same space.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, LAYOUT_SEPARATRIX
    use fortfem_spline_separatrix, only: layout_island_chain, island_chain_curves_t
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_space_update_geometry, spline_grade_knots
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr
    use spline_island_chain_problem, only: projection_kernel_t, chain, smooth, &
        r_axis, wall_a, wall_b
    implicit none

    interface
        subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine dposv
    end interface

    real(dp), parameter :: pi = acos(-1.0_dp), width = 0.12_dp
    integer, parameter :: ms(3) = [2, 3, 5]
    integer :: im

    do im = 1, 3
        call check_chain(ms(im))
    end do
    call check_summary("spline_island_chain")

contains

    subroutine check_chain(m)
        integer, intent(in) :: m

        type(block_layout_t) :: layout, moved
        type(spline_space_t) :: space
        real(dp) :: smooth_error(2), cut_error(2), whole_error(2), area, rate
        real(dp) :: cut_rate, error, angle(4)
        character(len=8) :: label
        integer :: status, level, ncore, dofs, axis_vertex
        logical :: valence_ok

        write (label, '(a,i0)') "m = ", m
        call layout_island_chain(chain(m, width), 0.5_dp, layout, status)
        if (status /= 0) then
            call check_condition(.false., trim(label)//": layout")
            return
        end if
        valence_ok = all(xpoint_valence(layout, chain(m, width)) == 8)
        ncore = 2*m + merge(1, m, m == 2)
        call check_condition(layout%block_count == merge(19, 10*m, m == 2) .and. &
            valence_ok .and. count(layout%side_tag == LAYOUT_SEPARATRIX) == 8*m, &
            trim(label)//": block count, X-points of valence 8, 4m branch pieces")
        if (m == 2) then
            angle = corner_angles(layout, ncore)
            call check_condition(minval(angle) > 60.0_dp .and. &
                maxval(angle) < 120.0_dp, &
                trim(label)//": the core quad is regular")
        else
            axis_vertex = minloc(norm2(layout%vertices - spread([r_axis, 0.0_dp], &
                2, layout%vertex_count), dim=1), 1)
            call check_condition(norm2(layout%vertices(:, axis_vertex) - &
                [r_axis, 0.0_dp]) < 1.0e-14_dp .and. &
                layout%valence(axis_vertex) == m, &
                trim(label)//": the axis is a core vertex of valence m")
        end if

        do level = 1, 2
            call build(layout, 4*level, .false., space, status)
            if (status /= 0) error stop "space"
            call project(space, ncore, .false., smooth_error(level), area)
            call project(space, ncore, .true., whole_error(level), area)
            call build(layout, 4*level, .true., space, status)
            if (status /= 0) error stop "cut space"
            call project(space, ncore, .true., cut_error(level), area)
        end do
        rate = log(smooth_error(1)/smooth_error(2))/log(2.0_dp)
        cut_rate = log(cut_error(1)/cut_error(2))/log(2.0_dp)
        print '(a,2es11.3,a,f6.2,a,2es11.3,a,f6.2,a,2es11.3,a,es10.2)', &
            trim(label)//" smooth", smooth_error, " rate", rate, " jump cut", &
            cut_error, " rate", cut_rate, " jump C0", whole_error, " area error", &
            abs(area - pi*wall_a*wall_b)
        call check_condition(rate > 3.5_dp .and. &
            abs(area - pi*wall_a*wall_b) < 1.0e-4_dp, &
            trim(label)//": smooth projection at order p + 1")
        call check_condition(cut_rate > 3.5_dp .and. &
            whole_error(2) > 1.0e2_dp*cut_error(2), &
            trim(label)//": jump across island 1 at order p + 1 in the cut space only")

        call build(layout, 8, .true., space, status, grade=.true.)
        call check_condition(status == 0, &
            trim(label)//": knots graded towards the separatrix conform")
        call build(layout, 8, .true., space, status)
        dofs = space%dof_count
        call layout_island_chain(chain(m, 0.9_dp*width), 0.5_dp, moved, status)
        call spline_space_update_geometry(space, moved, status)
        call project(space, ncore, .false., error, area)
        call check_condition(status == 0 .and. space%dof_count == dofs .and. &
            error < 2.0_dp*smooth_error(2), &
            trim(label)//": a moved chain keeps the space and its accuracy")
    end subroutine check_chain

    function xpoint_valence(layout, curves) result(valence)
        type(block_layout_t), intent(in) :: layout
        type(island_chain_curves_t), intent(in) :: curves
        integer :: valence(size(curves%xpoint, 2))

        integer :: i, v

        do i = 1, size(curves%xpoint, 2)
            v = minloc(norm2(layout%vertices - spread(curves%xpoint(:, i), 2, &
                layout%vertex_count), dim=1), 1)
            valence(i) = layout%valence(v)
            if (norm2(layout%vertices(:, v) - curves%xpoint(:, i)) > 1.0e-14_dp) &
                valence(i) = 0
        end do
    end function xpoint_valence

    function corner_angles(layout, block) result(angle)
        !! Interior angles in degrees of a block with straight sides.
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: block
        real(dp) :: angle(4)

        real(dp) :: a(2), b(2), c(2)
        integer :: j

        do j = 1, 4
            a = layout%vertices(:, layout%corners(modulo(j - 2, 4) + 1, block))
            b = layout%vertices(:, layout%corners(j, block))
            c = layout%vertices(:, layout%corners(modulo(j, 4) + 1, block))
            angle(j) = acos(dot_product(a - b, c - b)/(norm2(a - b)*norm2(c - b)))* &
                180.0_dp/pi
        end do
    end function corner_angles

    subroutine build(layout, cells, cut, space, status, grade)
        type(block_layout_t), intent(in) :: layout
        integer, intent(in) :: cells
        logical, intent(in) :: cut
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status
        logical, intent(in), optional :: grade

        type(spline_space_options_t) :: options

        options%degree = 3
        allocate(options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = cells
        if (cut) options%cut_tags = [LAYOUT_SEPARATRIX]
        if (present(grade)) then
            call spline_grade_knots(layout, LAYOUT_SEPARATRIX, 1.0_dp, &
                options%knots, status, power=4.0_dp)
            if (status /= 0) return
        end if
        call spline_space_build(layout, options, space, status)
    end subroutine build

    subroutine project(space, ncore, jump, error, area)
        !! L2 projection of smooth(s), plus 1 in island 1 when jump.
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: ncore
        logical, intent(in) :: jump
        real(dp), intent(out) :: error, area

        type(spline_quadrature_t) :: quad
        type(projection_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), g(:, :), jet(:, :), a(:, :), rhs(:, :)
        real(dp) :: energy
        integer :: status, q, row, k

        call spline_quadrature_build(space, 5, quad, status)
        if (status /= 0) error stop "quadrature"
        allocate(kernel%f(quad%point_count))
        do q = 1, quad%point_count
            kernel%f(q) = smooth(quad%s(:, q))
            if (jump .and. quad%block(q) > ncore .and. quad%block(q) <= ncore + 5) &
                kernel%f(q) = kernel%f(q) + 1.0_dp
        end do
        allocate(c(1, space%dof_count), g(1, space%dof_count))
        c = 0.0_dp
        call spline_energy(space, quad, kernel, c, energy, status, g)
        call spline_energy_hessian_csr(space, quad, kernel, c, matrix, status)
        allocate(a(matrix%n, matrix%n), rhs(matrix%n, 1))
        a = 0.0_dp
        do row = 1, matrix%n
            do k = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                a(row, matrix%colidx(k)) = matrix%val(k)
            end do
        end do
        rhs(:, 1) = -g(1, :)
        call dposv("L", matrix%n, 1, a, matrix%n, rhs, matrix%n, status)
        if (status /= 0) error stop "dposv"
        c(1, :) = rhs(:, 1)
        allocate(jet(1, quad%point_count))
        call spline_evaluate(space, quad, c, 0, jet, status)
        error = sqrt(sum(quad%weight*(jet(1, :) - kernel%f)**2))
        area = sum(quad%weight)
    end subroutine project

end program test_spline_island_chain
