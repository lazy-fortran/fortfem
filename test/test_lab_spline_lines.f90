program test_lab_spline_lines
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use check, only: check_condition, check_summary
    use fortfem_tensor_line_tables, only: tensor_line_tables_t
    use fortfem_bspline_lines, only: build_bspline_line, build_bspline_line_pair
    use fortfem_fourier_lines, only: build_fourier_line, build_fourier_line_pair
    use fortfem_bspline_polar, only: evaluate_periodic_bspline_basis
    use fortfem_assembly_bspline_2d, only: assemble_bspline_h1_operator_csc
    use fortsparse, only: csc_t, fortsparse_status_t, csc_matvec
    implicit none
    type(tensor_line_tables_t) :: v0, v1
    real(dp), allocatable :: incidence(:, :), mass(:, :), stiff(:, :)
    real(dp), allocatable :: values(:), deriv(:), actual(:), expected(:)
    real(dp), parameter :: breaks(*) = [2.0_dp, 2.1_dp, 2.7_dp, 3.0_dp]
    real(dp), parameter :: open_knots(*) = &
        [2.0_dp, 2.0_dp, 2.0_dp, 2.1_dp, 2.7_dp, 3.0_dp, 3.0_dp, 3.0_dp]
    real(dp) :: pi, diagonal, frequency
    integer :: p, k, periodic, status, e, q, i, r, mode

    do periodic = 0, 1
        do p = 1, 5
            do k = 0, p - 1
                call build_bspline_line_pair(p, v0, v1, incidence, status, &
                    nelem=7, continuity=k, periodic=periodic == 1, nq=p + 2)
                call check_condition(status == 0, 'Spline pair builds')
                if (status /= 0) error stop 'Cannot test invalid pair'
                call check_condition(maxval(abs(sum(v0%b, dim=2) - 1)) &
                    < 3e-15_dp, 'Partition of unity at every point')
                call check_condition(maxval(abs(sum(v0%d, dim=2))) &
                    < 2e-13_dp, 'Derivative annihilates constants')
                call commuting(v0, v1, incidence)
                call check_condition(v1%ndof == v0%ndof - 1 + periodic, &
                    'Open/periodic exact sequence dimensions')
                do r = 0, 2*p + 1
                    call check_condition(abs(sum(v0%wq*v0%xq**r) &
                        - 1.0_dp/real(r + 1, dp)) < 2e-15_dp, &
                        'Gauss polynomial moments')
                end do
            end do
        end do
    end do
    call build_bspline_line_pair(2, v0, v1, incidence, status, &
        knots=open_knots, nq=4)
    call check_condition(status == 0, 'Nonuniform open knots')
    call commuting(v0, v1, incidence)
    call check_condition(abs(sum(v0%wq*v0%xq**5) &
        - (3.0_dp**6 - 2.0_dp**6)/6) < 1e-13_dp, &
        'Nonunit interval polynomial integral')
    call build_bspline_line_pair(2, v0, v1, incidence, status, &
        knots=[0.0_dp, 0.0_dp, 0.0_dp, 0.2_dp, 0.2_dp, &
        0.8_dp, 0.8_dp, 1.0_dp, 1.0_dp, 1.0_dp], continuity=0)
    call check_condition(status == 0, 'Nonuniform open C0 pair')
    call commuting(v0, v1, incidence)
    call build_bspline_line_pair(3, v0, v1, incidence, status, &
        knots=breaks, continuity=1, periodic=.true., nq=5)
    call check_condition(status == 0, 'Nonuniform periodic breaks')
    call check_condition(maxval(abs(sum(v0%b, dim=2) - 1)) < 3e-15_dp, &
        'Nonuniform periodic partition of unity')
    call commuting(v0, v1, incidence)
    call check_condition(any(v0%dof(:, v0%nelem) < v0%dof(1, v0%nelem)), &
        'Periodic dof wrap occurs')
    call build_bspline_line(3, v0, status, nelem=7, periodic=.true., nq=4)
    allocate(actual(7))
    do e = 1, v0%nelem
        do q = 1, v0%nq
            call evaluate_periodic_bspline_basis(7, 3, v0%xq(q, e), &
                values, deriv, status)
            actual = 0
            do i = 1, v0%nloc
                ! The extended-knot convention shifts the cardinal index by p.
                r = modulo(v0%dof(i, e) - 1 - 3, 7) + 1
                actual(r) = actual(r) + v0%b(q, i, e)
            end do
            call check_condition(maxval(abs(actual - values)) < 3e-15_dp, &
                'Periodic values match independent cardinal evaluator')
        end do
    end do
    ! A single cell wraps several local functions to the same dof.
    call build_bspline_line_pair(4, v0, v1, incidence, status, &
        nelem=1, periodic=.true.)
    call check_condition(status == 0, 'Small periodic space permits repeated dofs')
    call commuting(v0, v1, incidence)
    call build_bspline_line(0, v0, status, nelem=3)
    call check_condition(status == 0, 'Piecewise constants build')
    call build_bspline_line(2, v0, status, nelem=0)
    call check_condition(status /= 0, 'Reject empty mesh')
    call build_bspline_line(2, v0, status, nelem=2, continuity=2)
    call check_condition(status /= 0, 'Reject invalid continuity')
    call build_bspline_line(2, v0, status, knots=[0.0_dp, 1.0_dp])
    call check_condition(status /= 0, 'Reject short open knot vector')
    call build_bspline_line(2, v0, status, nelem=2, knots=open_knots)
    call check_condition(status /= 0, 'Reject ambiguous mesh input')
    call build_bspline_line_pair(1, v0, v1, incidence, status, &
        knots=[0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp, 1.0_dp])
    call check_condition(status /= 0, 'Reject excess endpoint multiplicity')
    pi = acos(-1.0_dp)
    do p = 0, 5
        call build_fourier_line_pair(p, v0, v1, incidence, status, &
            nq=2*p + 1, nfp=3)
        call check_condition(status == 0, 'Fourier pair at minimum quadrature')
        call commuting(v0, v1, incidence)
        call matrices(v0, mass, stiff)
        allocate(expected(v0%ndof))
        do i = 1, v0%ndof
            expected = 0
            diagonal = pi/3
            if (i == 1) diagonal = 2*pi/3
            expected(i) = diagonal
            call check_condition(maxval(abs(mass(:, i) - expected)) &
                < 1e-13_dp, 'Fourier mass diagonal and trig products exact')
            mode = i/2
            frequency = real(3*mode, dp)
            expected(i) = diagonal*frequency**2
            call check_condition(maxval(abs(stiff(:, i) - expected)) &
                < 2e-12_dp, 'Fourier stiffness has (k*nfp)^2 eigenvalues')
        end do
        deallocate(expected)
    end do
    call build_fourier_line(3, v0, status, nq=6)
    call check_condition(status /= 0, 'Reject Fourier aliasing quadrature')
    call build_fourier_line(1, v0, status, nfp=0)
    call check_condition(status /= 0, 'Reject invalid field periods')
    call assembly_oracles()
    call check_summary('lab spline and Fourier lines')
contains
    subroutine assembly_oracles()
        type(tensor_line_tables_t) :: x, y
        type(csc_t) :: assembled
        type(fortsparse_status_t) :: result
        real(dp), allocatable :: mx(:, :), sx(:, :), my(:, :), sy(:, :)
        real(dp), allocatable :: cp(:, :, :), weights(:, :), unit(:), column(:)
        real(dp), allocatable :: marginal(:, :)
        real(dp), parameter :: ty(*) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
        real(dp) :: gx, gy, target
        integer :: i, j, ix, iy, jx, jy, nx, ny, mode, status

        call build_bspline_line(2, x, status, knots=open_knots, nq=4)
        call check_condition(status == 0, 'Assembly oracle x line builds')
        call build_bspline_line(1, y, status, knots=ty, nq=4)
        call check_condition(status == 0, 'Assembly oracle y line builds')
        call matrices(x, mx, sx)
        call matrices(y, my, sy)
        nx = x%ndof
        ny = y%ndof
        allocate(cp(2, nx, ny), weights(nx, ny), unit(nx*ny), marginal(nx, nx))
        weights = 1
        do iy = 1, ny
            gy = ty(iy + 1)
            do ix = 1, nx
                gx = sum(open_knots(ix + 1:ix + 2))/2
                cp(1, ix, iy) = gx
                cp(2, ix, iy) = gy
            end do
        end do
        do mode = 0, 1
            call assemble_bspline_h1_operator_csc(open_knots, ty, 2, 1, &
                cp, weights, 4, assembled, result, &
                mass_coefficient=real(1 - mode, dp), &
                stiffness_coefficient=real(mode, dp))
            call check_condition(result%code == 0, 'Existing assembly succeeds')
            if (result%code /= 0) error stop 'Invalid assembly oracle'
            marginal = 0
            do jy = 1, ny
                do jx = 1, nx
                    j = jx + nx*(jy - 1)
                    unit = 0
                    unit(j) = 1
                    column = csc_matvec(assembled, unit)
                    do iy = 1, ny
                        do ix = 1, nx
                            i = ix + nx*(iy - 1)
                            if (mode == 0) then
                                target = mx(ix, jx)*my(iy, jy)
                            else
                                target = sx(ix, jx)*my(iy, jy) &
                                    + mx(ix, jx)*sy(iy, jy)
                            end if
                            call check_condition(abs(column(i) - target) &
                                < 1e-13_dp, '2D tensor equals B-spline assembly')
                            marginal(ix, jx) = marginal(ix, jx) + column(i)
                        end do
                    end do
                end do
            end do
            if (mode == 0) then
                call check_condition(maxval(abs(marginal - mx)) < 1e-13_dp, &
                    '1D mass equals constant-transverse assembly restriction')
            else
                call check_condition(maxval(abs(marginal - sx)) < 1e-13_dp, &
                    '1D stiffness equals constant-transverse restriction')
            end if
        end do
    end subroutine assembly_oracles

    subroutine commuting(a, b, incidence)
        type(tensor_line_tables_t), intent(in) :: a, b
        real(dp), intent(in) :: incidence(:, :)
        real(dp) :: db(a%ndof), vb(b%ndof), reconstructed(a%ndof)
        integer :: e, q, j
        do e = 1, a%nelem
            do q = 1, a%nq
                db = 0
                vb = 0
                do j = 1, a%nloc
                    db(a%dof(j, e)) = db(a%dof(j, e)) + a%d(q, j, e)
                end do
                do j = 1, b%nloc
                    vb(b%dof(j, e)) = vb(b%dof(j, e)) + b%b(q, j, e)
                end do
                do j = 1, a%ndof
                    reconstructed(j) = dot_product(vb, incidence(:, j))
                end do
                call check_condition(maxval(abs(db - reconstructed)) &
                    < 3e-13_dp, 'Every basis derivative commutes with incidence')
            end do
        end do
    end subroutine commuting

    subroutine matrices(a, m, s)
        type(tensor_line_tables_t), intent(in) :: a
        real(dp), allocatable, intent(out) :: m(:, :), s(:, :)
        integer :: e, i, j, ii, jj
        allocate(m(a%ndof, a%ndof), s(a%ndof, a%ndof))
        m = 0
        s = 0
        do e = 1, a%nelem
            do j = 1, a%nloc
                jj = a%dof(j, e)
                do i = 1, a%nloc
                    ii = a%dof(i, e)
                    m(ii, jj) = m(ii, jj) &
                        + sum(a%wq(:, e)*a%b(:, i, e)*a%b(:, j, e))
                    s(ii, jj) = s(ii, jj) &
                        + sum(a%wq(:, e)*a%d(:, i, e)*a%d(:, j, e))
                end do
            end do
        end do
    end subroutine matrices
end program test_lab_spline_lines
