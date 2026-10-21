module sipg_test_problems
    !! Diffusion coefficient A = 1 + amplitude x y and manufactured data for the
    !! SIPG equilibration test (homogeneous Dirichlet).
    use fortfem_dg_curved_2d, only: dg_coefficients_t, dg_data_t
    use fortfem_kinds, only: dp
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    type, extends(dg_coefficients_t), public :: diffusion_t
        real(dp) :: amplitude = 0.0_dp   !< A = 1 + amplitude x y
    contains
        procedure :: evaluate => diffusion_evaluate
    end type diffusion_t

    type, extends(dg_data_t), public :: problem_t
        logical :: polynomial = .false.
    contains
        procedure :: source => problem_source
        procedure :: boundary => problem_boundary
    end type problem_t

contains

    subroutine diffusion_evaluate(self, t, x, kk, b, c, status)
        class(diffusion_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        complex(dp), intent(out) :: kk(:, :, :)
        real(dp), intent(out) :: b(:, :)
        complex(dp), intent(out) :: c(:)
        integer, intent(out) :: status
        integer :: q

        associate (unused => t)
        end associate
        do q = 1, size(x, 2)
            kk(:, :, q) = (0.0_dp, 0.0_dp)
            kk(1, 1, q) = cmplx(1.0_dp + self%amplitude*x(1, q)*x(2, q), 0.0_dp, dp)
            kk(2, 2, q) = kk(1, 1, q)
        end do
        b = 0.0_dp
        c = (0.0_dp, 0.0_dp)
        status = 0
    end subroutine diffusion_evaluate

    subroutine problem_source(self, t, x, f, status)
        class(problem_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        complex(dp), intent(out) :: f(:)
        integer, intent(out) :: status
        real(dp) :: u, aa, value
        integer :: q

        associate (unused => t)
        end associate
        do q = 1, size(x, 2)
            if (self%polynomial) then
                value = 2.0_dp*(x(1, q)*(1.0_dp - x(1, q)) &
                    + x(2, q)*(1.0_dp - x(2, q)))
            else
                u = sin(pi*x(1, q))*sin(pi*x(2, q))
                aa = 1.0_dp + x(1, q)*x(2, q)
                ! -div(A grad u) = -(A_x u_x + A_y u_y) + 2 pi^2 A u.
                value = -pi*(x(2, q)*cos(pi*x(1, q))*sin(pi*x(2, q)) &
                    + x(1, q)*sin(pi*x(1, q))*cos(pi*x(2, q))) &
                    + 2.0_dp*pi*pi*aa*u
            end if
            f(q) = cmplx(value, 0.0_dp, dp)
        end do
        status = 0
    end subroutine problem_source

    subroutine problem_boundary(self, t, x, f, status)
        class(problem_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        complex(dp), intent(out) :: f(:)
        integer, intent(out) :: status

        associate (unused => t, unused2 => self, unused3 => x)
        end associate
        f = (0.0_dp, 0.0_dp)
        status = 0
    end subroutine problem_boundary

end module sipg_test_problems

program test_equilibrated_sipg_2d
    !! Equilibrated-residual bound of the SIPG solver (dg_curved_assemble) on
    !! affine triangle meshes of the unit square, homogeneous Dirichlet data.
    !! Oracles: (1) k = 4, A = 1, u = x(1-x)y(1-y) is in the broken space, so
    !! the SIPG solution is exact and every term of the bound vanishes (this
    !! exercises the patch equilibration, lifting and Oswald average);
    !! (2) u = sin(pi x) sin(pi y), A = 1 + x y, k = 1..3 on three meshes: the
    !! bound is >= the exact broken energy error ||A^(1/2) grad_h (u - u_h)||
    !! (quadrature of degree 2k + 10), the DG-norm bound is >= the DG-norm
    !! error, the efficiency index stays below a measured cap and the bound
    !! converges at the energy rate k.
    use check, only: check_condition, check_summary
    use fortfem_core, only: mesh_t, rectangle_mesh
    use fortfem_dg_curved_2d, only: DG_FORM_CONSERVATIVE, &
        dg_curved_assemble, dg_curved_rhs, dg_curved_space_init, &
        dg_curved_space_t
    use fortfem_equilibrated_nonconforming_2d, only: &
        nonconforming_estimate_t, NONCONFORMING_BAD_INPUT, NONCONFORMING_OK, &
        sipg_energy_bound
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_straight, &
        meridian_mesh_t
    use sipg_test_problems, only: diffusion_t, problem_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        evaluate_triangle_lagrange_basis, initialize_triangle_lagrange_basis, &
        triangle_lagrange_basis_t
    use fortfem_sparse_direct, only: sparse_direct_solve_once
    use fortsparse, only: csc_z_t, fortsparse_status_t, &
        status_ok
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: friedrichs = 1.0_dp/(pi*sqrt(2.0_dp))
    !> Measured indices 1.19-1.42 (k = 1..3); caps leave headroom.
    real(dp), parameter :: efficiency_cap(3) = [1.8_dp, 1.8_dp, 1.8_dp]

    type(nonconforming_estimate_t) :: est
    real(dp) :: energy_error, dg_error, bound(3), index, rate
    integer :: k, level, vertex_count

    call run_case(5, 4, .true., est, energy_error, dg_error)
    call check_condition(est%bound < 1.0e-9_dp, &
        "exact polynomial solution: SIPG bound vanishes")
    call check_condition(est%nonconformity < 1.0e-9_dp .and. &
        est%jump_seminorm < 1.0e-9_dp, "exact solution is conforming")
    print '(a, 4es10.2)', "polynomial case bound, error, nu, J ", est%bound, &
        energy_error, est%nonconformity, est%jump_seminorm

    do k = 1, 3
        do level = 1, 3
            vertex_count = 1 + 2**(level + 1)
            call run_case(vertex_count, k, .false., est, energy_error, dg_error)
            bound(level) = est%bound
            index = est%bound/energy_error
            call check_condition(est%bound >= energy_error, &
                "SIPG bound is above the broken energy error")
            call check_condition(est%potential_bound >= dg_error, &
                "SIPG DG-norm bound is above the DG-norm error")
            call check_condition(index <= efficiency_cap(k), &
                "SIPG efficiency index below the cap")
            ! The DG lifting makes every interior patch problem compatible; the
            ! remainder is the solver's quadrature of f (without the lifting the
            ! defect is 1e-2 to 1e-6 relative to the bound).
            call check_condition(est%divergence_defect < 1.0e-5_dp*est%bound, &
                "lifted patch loads leave no divergence defect")
            print '(a, i2, a, i3, 7(a, es10.2))', "k", k, " n", vertex_count - 1, &
                " err", energy_error, " bound", est%bound, " index", index, &
                " rho", est%flux_defect, " nu", est%nonconformity, " J", &
                est%jump_seminorm, " div", est%divergence_defect
        end do
        rate = log(bound(2)/bound(3))/log(2.0_dp)
        call check_condition(rate > real(k, dp) - 0.3_dp, &
            "SIPG bound converges at the energy rate k")
        print '(a, i2, a, f6.2)', "k", k, " bound rate", rate
    end do
    call check_rejections()
    call check_summary("equilibrated_sipg_2d")

contains

    subroutine check_rejections()
        !! Curved cells and the inconsistent k = 0 SIPG space are rejected.
        type(mesh_t) :: mesh
        type(meridian_mesh_t) :: mm
        type(diffusion_t) :: coef
        type(problem_t) :: data
        type(dg_curved_space_t) :: space
        type(nonconforming_estimate_t) :: estimate
        complex(dp), allocatable :: u(:)
        integer :: st, geo, degree

        mesh = rectangle_mesh(3, 3, [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
        do geo = 1, 2
            degree = 2 - mod(geo, 2)*2
            call meridian_mesh_straight(mesh%data, geo, mm, st)
            call dg_curved_space_init(mm, degree, space, st)
            allocate (u(space%ndof))
            u = (0.0_dp, 0.0_dp)
            call sipg_energy_bound(space, coef, data, u, 2*max(degree, 1) + 8, &
                [(1.0_dp, st = 1, space%nt)], friedrichs, estimate, st)
            call check_condition(st == NONCONFORMING_BAD_INPUT, &
                "curved geometry or k = 0 is rejected")
            deallocate (u)
        end do
    end subroutine check_rejections

    subroutine run_case(vertex_count, degree, polynomial, estimate, energy_err, &
            dg_err)
        integer, intent(in) :: vertex_count, degree
        logical, intent(in) :: polynomial
        type(nonconforming_estimate_t), intent(out) :: estimate
        real(dp), intent(out) :: energy_err, dg_err
        type(mesh_t) :: mesh
        type(diffusion_t) :: coef
        type(problem_t) :: data
        type(dg_curved_space_t) :: space
        type(csc_z_t) :: a
        type(fortsparse_status_t) :: fst
        complex(dp), allocatable :: rhs(:), u(:)
        real(dp), allocatable :: lower(:)
        integer :: st
        type(meridian_mesh_t) :: mm

        mesh = rectangle_mesh(vertex_count, vertex_count, &
            [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
        call meridian_mesh_straight(mesh%data, 1, mm, st)
        if (st /= 0) error stop "mesh"
        call dg_curved_space_init(mm, degree, space, st)
        if (st /= 0) error stop "space"
        coef%amplitude = merge(0.0_dp, 1.0_dp, polynomial)
        data%polynomial = polynomial
        allocate (rhs(space%ndof), u(space%ndof), lower(space%nt))
        call dg_curved_assemble(space, coef, DG_FORM_CONSERVATIVE, a, st)
        if (st /= 0) error stop "assemble"
        call dg_curved_rhs(space, coef, data, DG_FORM_CONSERVATIVE, rhs, st)
        if (st /= 0) error stop "rhs"
        call sparse_direct_solve_once(a, rhs, u, fst)
        if (.not. status_ok(fst)) error stop "solve"
        lower = 1.0_dp
        call sipg_energy_bound(space, coef, data, u, 2*degree + 8, lower, &
            friedrichs, estimate, st)
        call check_condition(st == NONCONFORMING_OK, "SIPG estimate succeeds")
        if (st /= NONCONFORMING_OK) then
            print '(a, i0)', "SIPG estimate status ", st
            error stop "SIPG estimate failed"
        end if
        energy_err = broken_energy_error(space, u, polynomial, coef%amplitude)
        dg_err = sqrt(energy_err**2 + estimate%jump_seminorm**2)
    end subroutine run_case

    function broken_energy_error(space, u, polynomial, amplitude) result(err)
        !! ||A^(1/2) grad_h (u - u_h)||, exact solution, rule of degree 2k + 10.
        type(dg_curved_space_t), intent(in) :: space
        complex(dp), intent(in) :: u(:)
        logical, intent(in) :: polynomial
        real(dp), intent(in) :: amplitude
        real(dp) :: err
        type(triangle_lagrange_basis_t) :: basis
        real(dp), allocatable :: xi(:), eta(:), w(:), val(:), grad(:, :)
        real(dp) :: v(2, 3), jac(2, 2), det, jit(2, 2), x(2), gu(2), gh(2), aa
        integer :: nl, t, q, st

        nl = space%nl
        call initialize_triangle_lagrange_basis(space%k, basis, st)
        call triangle_duffy_quadrature(2*space%k + 10, xi, eta, w, st)
        allocate (val(nl), grad(2, nl))
        err = 0.0_dp
        do t = 1, space%nt
            v = space%mm%topology%vertices(:, space%mm%topology%triangles(1:3, t))
            jac(:, 1) = v(:, 2) - v(:, 1)
            jac(:, 2) = v(:, 3) - v(:, 1)
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            jit(1, :) = [jac(2, 2), -jac(2, 1)]/det
            jit(2, :) = [-jac(1, 2), jac(1, 1)]/det
            do q = 1, size(w)
                x = v(:, 1) + xi(q)*jac(:, 1) + eta(q)*jac(:, 2)
                call evaluate_triangle_lagrange_basis(basis, xi(q), eta(q), val, &
                    grad, st)
                gh = matmul(jit, matmul(grad, real(u((t - 1)*nl + 1:t*nl), dp)))
                gu = exact_gradient(x, polynomial)
                aa = 1.0_dp + amplitude*x(1)*x(2)
                err = err + w(q)*det*aa*sum((gu - gh)**2)
            end do
        end do
        err = sqrt(err)
    end function broken_energy_error

    function exact_gradient(x, polynomial) result(g)
        real(dp), intent(in) :: x(2)
        logical, intent(in) :: polynomial
        real(dp) :: g(2)

        if (polynomial) then
            g = [(1.0_dp - 2.0_dp*x(1))*x(2)*(1.0_dp - x(2)), &
                x(1)*(1.0_dp - x(1))*(1.0_dp - 2.0_dp*x(2))]
        else
            g = pi*[cos(pi*x(1))*sin(pi*x(2)), sin(pi*x(1))*cos(pi*x(2))]
        end if
    end function exact_gradient

end program test_equilibrated_sipg_2d
