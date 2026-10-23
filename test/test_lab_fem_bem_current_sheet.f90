program test_lab_fem_bem_current_sheet
    !! Vacuum-field free-boundary model: a prescribed divergence-free surface
    !! current K = grad(Phi) x n on the torus boundary Gamma (mu0 = 1). Its
    !! field is B = -grad(psi) off Gamma with psi the double-layer potential
    !! of Phi, i.e. the transmission problem u - w = -Phi, d_n u = d_n w.
    !! The FEM/BEM field (interior FEM gradient on the axis core cells,
    !! exterior Green representation) is compared with an independent
    !! Biot-Savart surface integral of K on the same discrete surface.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t, ogrid_tensor_mesh, build_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator
    use fortfem_tensor_solver, only: tensor_solve_info_t
    use fortfem_tensor_fem_bem, only: tensor_fem_bem_t, build_tensor_fem_bem, &
                                      solve_tensor_fem_bem, tensor_fem_bem_rhs, &
                                      tensor_fem_bem_exterior_flux, &
                                      free_tensor_fem_bem, &
                                      tensor_h1_seminorm_error
    use fortfem_parametric_bem_tensor, only: bem_tensor_exterior_field
    use fem_bem_oracle, only: unit_coefficient
    implicit none
    integer, parameter :: nout = 10
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: outside(3, nout), e(2, 3)
    type(tensor_fem_bem_t), allocatable, target :: sheet
    integer :: failures = 0, its(3), values(3), k
    character(32) :: arg
    outside = reshape([0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.3_dp, &
                       0.0_dp, -1.0_dp, -0.4_dp, 5.6_dp, 0.0_dp, 0.0_dp, &
                       0.0_dp, -5.6_dp, 0.5_dp, -4.2_dp, 3.9_dp, -0.6_dp, &
                       3.0_dp, 0.0_dp, 2.3_dp, 0.0_dp, -3.0_dp, -2.3_dp, &
                       -2.2_dp, -2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], [3, nout])
    print '(a)', 'FB1_SHEET n,nphi,p,dofs,trace,its,ext_rel_max,core_rel_l2'
    if (command_argument_count() >= 3) then
        do k = 1, 3
            call get_command_argument(k, arg)
            read (arg, *) values(k)
        end do
        call case(values(1), values(2), values(3), e(:, 1), its(1))
        stop
    end if
    call case(1, 6, 2, e(:, 1), its(1))
    call case(1, 6, 4, e(:, 2), its(2))
    call case(2, 12, 2, e(:, 3), its(3))
    call check(all(e(:, 2) < e(:, 1)), 'p-refinement reduces exterior and core B error')
    call check(all(e(:, 3) < e(:, 1)), 'h-refinement reduces exterior and core B error')
    call check(maxval(e(:, 2)) < 5.e-3_dp, 'p=4 B errors below 5e-3')
    if (failures /= 0) error stop 'test_lab_fem_bem_current_sheet failed'
    print '(a)', 'test_lab_fem_bem_current_sheet passed'
contains
    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (condition) then
            print '(a,a)', '[PASS] ', label
        else
            print '(a,a)', '[FAIL] ', label
            failures = failures + 1
        end if
    end subroutine

    function phi(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = x(3) + x(1)*x(2)*x(3)/9
    end function

    function minus_phi(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = -phi(x)
    end function

    pure function grad_phi(x) result(g)
        real(dp), intent(in) :: x(3)
        real(dp) :: g(3)
        g = [x(2)*x(3)/9, x(1)*x(3)/9, 1 + x(1)*x(2)/9]
    end function

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function

    subroutine biot_savart(fb, x, nq, b)
        !! B(x) = 1/(4 pi) int K(y) x (x - y)/|x - y|^3 dA_y, tensor Gauss.
        type(tensor_fem_bem_t), intent(in) :: fb
        real(dp), intent(in) :: x(3)
        integer, intent(in) :: nq
        real(dp), intent(out) :: b(3)
        real(dp) :: t(nq), w(nq), y(3), a1(3), a2(3), n(3), d(3)
        integer :: panel, i, j, status
        call gauss01(nq, t, w)
        b = 0
        do panel = 1, fb%trace%npanel
            do j = 1, nq
                do i = 1, nq
                    call fb%trace%surface%evaluate(panel, [t(i), t(j)], y, a1, a2, &
                                                   status)
                    n = cross(a1, a2)
                    d = x - y
                    ! |a1 x a2| dA cancels the unit normal normalisation.
                    b = b + w(i)*w(j)*cross(cross(grad_phi(y), n), d)/norm2(d)**3
                end do
            end do
        end do
        b = b/(4*pi)
    end subroutine

    subroutine gauss01(n, t, w)
        integer, intent(in) :: n
        real(dp), intent(out) :: t(n), w(n)
        real(dp) :: x, p0, p1, p2, dp1
        integer :: i, k, it
        do i = 1, n
            x = cos(pi*(i - 0.25_dp)/(n + 0.5_dp))
            do it = 1, 100
                p0 = 1
                p1 = x
                do k = 2, n
                    p2 = ((2*k - 1)*x*p1 - (k - 1)*p0)/k
                    p0 = p1
                    p1 = p2
                end do
                dp1 = n*(x*p1 - p0)/(x*x - 1)
                x = x - p1/dp1
                if (abs(p1/dp1) < 1.e-15_dp) exit
            end do
            t(i) = (1 - x)/2
            w(i) = 1/((1 - x*x)*dp1*dp1)
        end do
    end subroutine

    subroutine minus_b(x, g)
        !! grad u = -B on the interior core, from the Biot-Savart reference.
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: g(:)
        real(dp) :: b(3)
        call biot_savart(sheet, x(1:3), 48, b)
        g(1:3) = -b
    end subroutine

    subroutine case(n, nphi, p, errors, iterations)
        !! errors(1): max over exterior points |B_h - B| / max |B|;
        !! errors(2): relative L2 error of B on the axis core cells.
        integer, intent(in) :: n, nphi, p
        real(dp), intent(out) :: errors(2)
        integer, intent(out) :: iterations
        type(tensor_mesh_t) :: mesh, core
        type(tensor_operator_t) :: core_op
        type(tensor_solve_info_t) :: info
        real(dp), allocatable :: load(:), rhs(:), x(:), trace(:), flux(:)
        real(dp), allocatable :: potential(:), gradient(:, :), zero(:)
        real(dp) :: reference(3, nout), check_ref(3), scale, bnorm
        integer :: status, i, iz, cells(nphi*n*n)
        allocate (sheet)
        call ogrid_tensor_mesh(3, n, p, mesh, status, nphi=nphi)
        if (status /= 0) error stop 'mesh'
        mesh%boundary = .false.
        call build_tensor_fem_bem(mesh, unit_coefficient, sheet, status)
        if (status /= 0) error stop 'build'
        allocate (load(mesh%ndof), x(mesh%ndof), zero(mesh%ndof))
        load = 0
        x = 0
        zero = 0
        call tensor_fem_bem_rhs(sheet, load, rhs, status, dirichlet=minus_phi)
        if (status /= 0) error stop 'rhs'
        call solve_tensor_fem_bem(sheet, rhs, 1.e-11_dp, x, info)
        iterations = info%iterations
        ! Exterior trace w = u + Phi and its flux, then -grad of the
        ! representation at points well away from Gamma.
        allocate (trace(sheet%trace%ntrace), potential(nout), gradient(3, nout))
        do i = 1, sheet%trace%ntrace
            trace(i) = x(sheet%trace%trace_dof(i)) + &
                       phi(mesh%coordinates(:, sheet%trace%trace_dof(i)))
        end do
        call tensor_fem_bem_exterior_flux(sheet, trace, flux, status)
        if (status /= 0) error stop 'flux'
        call bem_tensor_exterior_field(sheet%trace%surface, sheet%space, trace, &
                                       flux, outside, potential, gradient, status, 40)
        if (status /= 0) error stop 'field'
        do i = 1, nout
            call biot_savart(sheet, outside(:, i), 48, reference(:, i))
        end do
        call biot_savart(sheet, outside(:, 2), 64, check_ref)
        if (norm2(check_ref - reference(:, 2)) > 1.e-9_dp*norm2(check_ref)) &
            error stop 'Biot-Savart reference not converged'
        scale = maxval(norm2(reference, dim=1))
        errors(1) = maxval(norm2(-gradient - reference, dim=1))/scale
        ! Axis core: the n*n block-0 cells of every toroidal layer.
        do iz = 0, nphi - 1
            do i = 1, n*n
                cells(iz*n*n + i) = iz*5*n*n + i
            end do
        end do
        call build_tensor_mesh(3, p, mesh%geometry_degree, mesh%vertices(:, cells), &
                               mesh%geometry(:, :, cells), core, status, &
                               gather=mesh%dof(:, cells))
        if (status /= 0) error stop 'core mesh'
        core%boundary = .false.
        call build_tensor_operator(core, unit_coefficient, core_op, status, nq=p + 2)
        if (status /= 0) error stop 'core operator'
        bnorm = tensor_h1_seminorm_error(core_op, zero, minus_b)
        errors(2) = tensor_h1_seminorm_error(core_op, x, minus_b)/bnorm
        print '(a,6(i0,a),2(es10.3,a))', 'FB1_SHEET ', n, ',', nphi, ',', p, ',', &
            mesh%ndof, ',', sheet%trace%ntrace, ',', iterations, ',', errors(1), &
            ',', errors(2), ''
        call free_tensor_fem_bem(sheet)
        deallocate (sheet)
    end subroutine
end program
