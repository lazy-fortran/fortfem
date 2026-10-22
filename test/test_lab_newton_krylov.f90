module test_lab_newton_krylov_problem
    ! Semilinear MMS on the unit-disk O-grid through the axis:
    ! -lap u + u^3 = f, u = cos(a r^2), a = pi/2, u = 0 at r = 1.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_newton_krylov, only: nonlinear_system_t
    use fortfem_tensor_operator, only: tensor_operator_t, apply_tensor_operator
    use fortfem_tensor_solver, only: tensor_multigrid_t, tensor_vcycle
    use fortfem_tensor_pointwise, only: tensor_quadrature_values, &
                                        tensor_weighted_mass, tensor_weighted_load
    implicit none
    private
    public :: cubic_system_t, exact, identity_coefficient
    real(dp), parameter :: a = acos(-1.0_dp)/2
    type, extends(nonlinear_system_t) :: cubic_system_t
        type(tensor_operator_t) :: op
        type(tensor_multigrid_t) :: mg
        real(dp), allocatable :: uq(:, :), source(:, :), dsource(:, :)
    contains
        procedure :: residual, linearise, jacobian, precondition
    end type
contains
    function exact(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = cos(a*sum(x(1:2)**2))
    end function

    subroutine identity_coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(k, 1)
            k(i, i) = 1
        end do
        c = 0
        if (size(x) < 0) c = 1
    end subroutine

    subroutine pointwise(this, u)
        class(cubic_system_t), intent(inout) :: this
        real(dp), intent(in) :: u(:)
        real(dp) :: r2, f
        integer :: e, q
        if (.not. allocated(this%uq)) then
            allocate (this%uq(product(this%op%no), this%op%mesh%nelem))
            allocate (this%source, this%dsource, mold=this%uq)
        end if
        call tensor_quadrature_values(this%op, u, this%uq)
        do e = 1, this%op%mesh%nelem
            do q = 1, size(this%uq, 1)
                r2 = sum(this%op%points(1:2, q, e)**2)
                f = 4*a*sin(a*r2) + 4*a**2*r2*cos(a*r2) + cos(a*r2)**3
                this%source(q, e) = f - this%uq(q, e)**3
                this%dsource(q, e) = 3*this%uq(q, e)**2
            end do
        end do
    end subroutine

    subroutine residual(this, u, r, status)
        class(cubic_system_t), intent(inout) :: this
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: r(:)
        integer, intent(out) :: status
        real(dp), allocatable :: load(:)
        allocate (load(size(u)))
        call pointwise(this, u)
        call tensor_weighted_load(this%op, this%source, load)
        call apply_tensor_operator(this%op, u, r)
        r = r - load
        status = 0
    end subroutine

    subroutine linearise(this, u, status)
        class(cubic_system_t), intent(inout) :: this
        real(dp), intent(in) :: u(:)
        integer, intent(out) :: status
        call pointwise(this, u)
        status = 0
    end subroutine

    subroutine jacobian(this, v, w, status)
        class(cubic_system_t), intent(inout) :: this
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: w(:)
        integer, intent(out) :: status
        real(dp), allocatable :: mv(:)
        allocate (mv(size(v)))
        call apply_tensor_operator(this%op, v, w)
        call tensor_weighted_mass(this%op, this%dsource, v, mv)
        w = w + mv
        status = 0
    end subroutine

    subroutine precondition(this, v, w, status)
        class(cubic_system_t), intent(inout) :: this
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: w(:)
        integer, intent(out) :: status
        call tensor_vcycle(this%mg, 1, v, w, status)
    end subroutine
end module

program test_lab_newton_krylov
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use test_lab_newton_krylov_problem, only: cubic_system_t, exact, &
                                              identity_coefficient
    use fortfem_newton_krylov, only: newton_options_t, newton_info_t, &
                                     newton_krylov_solve
    use fortfem_tensor_mesh, only: tensor_mesh_t, ogrid_tensor_mesh
    use fortfem_tensor_operator, only: build_tensor_operator, tensor_l2_error
    use fortfem_tensor_solver, only: build_tensor_multigrid
    use fortfem_tensor_points, only: tensor_locator_t, build_tensor_locator, &
                                     locate_tensor_point, evaluate_tensor_point, &
                                     tensor_field_extremum
    implicit none
    integer :: failures, p, k
    real(dp) :: errors(3)
    failures = 0
    call check_points()
    do k = 1, 3
        p = 2*k
        call solve_case(p, errors(k))
    end do
    ! Spectral p-convergence on the curved through-axis O-grid.
    call expect(errors(2) < 0.05_dp*errors(1), 'p=4 error well below p=2')
    call expect(errors(3) < 0.05_dp*errors(2), 'p=6 error well below p=4')
    call expect(errors(3) < 3.e-7_dp, 'p=6 relative L2 error < 3e-7')
    if (failures > 0) error stop 'test_lab_newton_krylov failed'
    print '(a)', 'test_lab_newton_krylov passed'
contains
    subroutine expect(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', label
        end if
    end subroutine

    subroutine solve_case(p, relative)
        integer, intent(in) :: p
        real(dp), intent(out) :: relative
        type(tensor_mesh_t) :: mesh
        type(cubic_system_t) :: system
        type(newton_options_t) :: options
        type(newton_info_t) :: info
        real(dp), allocatable :: u(:), v(:), r0(:), r1(:), jv(:), zero(:)
        real(dp) :: h, fd_error
        integer :: status, i
        call ogrid_tensor_mesh(2, 2, p, mesh, status)
        call expect(status == 0, 'O-grid mesh')
        call build_tensor_operator(mesh, identity_coefficient, system%op, status)
        call expect(status == 0, 'operator')
        call build_tensor_multigrid(mesh, identity_coefficient, system%mg, status)
        call expect(status == 0, 'p-multigrid')
        allocate (u(mesh%ndof), v(mesh%ndof), r0(mesh%ndof), r1(mesh%ndof), &
                  jv(mesh%ndof), zero(mesh%ndof))
        ! Jacobian action against a central difference of the residual.
        do i = 1, mesh%ndof
            u(i) = 0.3_dp*sin(real(i, dp))
            v(i) = cos(real(3*i, dp))
        end do
        where (mesh%boundary)
            u = 0
            v = 0
        end where
        call system%linearise(u, status)
        call system%jacobian(v, jv, status)
        h = 1.e-6_dp
        call system%residual(u + h*v, r1, status)
        call system%residual(u - h*v, r0, status)
        fd_error = norm2((r1 - r0)/(2*h) - jv)/norm2(jv)
        call expect(fd_error < 1.e-7_dp, 'Jacobian matches finite difference')
        u = 0
        options%rtol = 1.e-11_dp
        call newton_krylov_solve(system, u, options, info)
        call expect(info%status == 0, 'Newton-Krylov converged')
        call expect(info%newton <= 8, 'at most eight Newton steps')
        call expect(maxval(info%krylov_history(1:info%newton)) <= 25, &
                    'p-MG preconditioned Krylov iterations bounded')
        zero = 0
        relative = tensor_l2_error(system%op, u, exact)/ &
                   tensor_l2_error(system%op, zero, exact)
        print '(a,i0,a,i0,a,i0,a,i0,a,es10.3,a,es10.3)', 'p=', p, ' dofs=', &
            mesh%ndof, ' newton=', info%newton, ' krylov=', info%krylov, &
            ' fd=', fd_error, ' rel_l2=', relative
    end subroutine

    subroutine check_points()
        ! A quadratic field is exact on the affine core block for p >= 2.
        type(tensor_mesh_t) :: mesh
        type(tensor_locator_t) :: locator
        real(dp), allocatable :: u(:)
        real(dp) :: x(2), xi(2), value, gradient(2), extremum(2), peak
        integer :: status, i, cell
        call ogrid_tensor_mesh(2, 2, 3, mesh, status)
        call build_tensor_locator(mesh, locator, status)
        call expect(status == 0, 'locator')
        allocate (u(mesh%ndof))
        do i = 1, mesh%ndof
            u(i) = quadratic(mesh%coordinates(:, i))
        end do
        cell = 0
        x = [0.07_dp, -0.11_dp]
        call locate_tensor_point(mesh, locator, x, cell, xi, status)
        call expect(status == 0, 'locate core point')
        call evaluate_tensor_point(mesh, locator, u, cell, xi, value, gradient, &
                                   status)
        call expect(abs(value - quadratic(x)) < 1.e-13_dp, 'point value')
        call expect(norm2(gradient - [-2*(x(1) - 0.1_dp), &
                                      -4*(x(2) + 0.05_dp)]) < 1.e-12_dp, &
                    'point gradient')
        x = [0.5_dp, 0.6_dp]
        cell = 0
        call locate_tensor_point(mesh, locator, x, cell, xi, status)
        call expect(status == 0, 'locate curved outer-block point')
        cell = 0
        call tensor_field_extremum(mesh, locator, u, cell, xi, extremum, peak, &
                                   status, maximum=.true.)
        call expect(status == 0, 'extremum search')
        call expect(norm2(extremum - [0.1_dp, -0.05_dp]) < 1.e-8_dp, &
                    'extremum position')
        call expect(abs(peak - 0.25_dp) < 1.e-13_dp, 'extremum value')
    end subroutine

    pure function quadratic(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = 0.25_dp - (x(1) - 0.1_dp)**2 - 2*(x(2) + 0.05_dp)**2
    end function
end program
