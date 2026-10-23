module fortfem_tensor_fem_bem
    !! Matrix-free tensor-product FEM interior coupled to the exterior Laplace
    !! problem through the symmetric BEM Dirichlet-to-Neumann Schur complement.
    !!
    !! Interior: -div(K grad u) + c u = f in the meshed region, boundary Gamma
    !! with outward normal n. Exterior: harmonic w, decaying at infinity.
    !! Transmission: u - w = g_D and K grad u . n - grad w . n = g_N on Gamma.
    !! With the exterior energy map S w = -d_n w (fortfem_parametric_laplace_bem_3d)
    !! the coupled weak form is
    !!     a(u, v) + <S R u, R v> = (f, v) + <g_N, R v> + <S g_D, R v>,
    !! where R restricts volume GLL face nodes to the Q_p BEM trace. The
    !! operator A + R^T S R is SPD because S is positive on constants, so the
    !! coupled solve is PCG. The preconditioner is the symmetric block
    !! factorisation with the Dirichlet p-multigrid V-cycle on the interior
    !! block and schur_scale*S (Cholesky) for the boundary Schur complement
    !! S + S_int: both are H^(1/2) maps, so S_int is replaced by a multiple of
    !! S. The default 3 minimised PCG iterations on the O-grid torus scans.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_get_wtime
    use fortfem_tensor_mesh, only: tensor_mesh_t
    use fortfem_tensor_operator, only: tensor_operator_t, apply_tensor_operator, &
                                       diffusion_coefficient, scalar_function, metric_at
    use fortfem_tensor_contractions, only: tensor_contract
    use fortfem_tensor_solver, only: tensor_multigrid_t, build_tensor_multigrid, &
                                     tensor_vcycle, free_tensor_multigrid, &
                                     tensor_solve_info_t
    use fortfem_gll_line, only: gll_rule
    use fortfem_tensor_boundary_trace, only: tensor_boundary_trace_t, &
                                             build_tensor_boundary_trace
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t, bem_exterior_dtn_t, &
                                                 build_exterior_dtn
    use fortfem_parametric_bem_tensor, only: bem_tensor_options_t, bem_tensor_space_t, &
                                             build_bem_tensor_space, &
                                             assemble_parametric_bem_tensor, &
                                             bem_tensor_boundary_load, boundary_function
    implicit none
    private
    public :: tensor_fem_bem_options_t, tensor_fem_bem_t, build_tensor_fem_bem
    public :: apply_tensor_fem_bem, solve_tensor_fem_bem, tensor_fem_bem_rhs
    public :: tensor_fem_bem_exterior_flux, free_tensor_fem_bem
    public :: tensor_h1_seminorm_error, gradient_function

    type :: tensor_fem_bem_options_t
        !! flux_degree < 0 selects max(trace degree - 1, 0).
        integer :: flux_degree = -1
        type(bem_tensor_options_t) :: quadrature
        real(dp) :: schur_scale = 3.0_dp
    end type

    type :: tensor_fem_bem_t
        type(tensor_boundary_trace_t) :: trace
        type(bem_tensor_space_t) :: space
        type(bem_operators_t) :: bem
        type(bem_exterior_dtn_t) :: dtn
        real(dp), allocatable :: steklov(:, :), schur(:, :)
        logical, allocatable :: on_trace(:)
        type(tensor_multigrid_t) :: mg
        !! Wall times in seconds: BEM assembly, DtN build, p-MG setup, Cholesky.
        real(dp) :: time_bem = 0, time_dtn = 0, time_mg = 0, time_schur = 0
    end type

    abstract interface
        subroutine gradient_function(x, g)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: g(:)
        end subroutine
    end interface

    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine
        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine
    end interface
contains
    subroutine build_tensor_fem_bem(mesh, coefficient, fb, status, options)
        !! Status 1 rejects the mesh, 2 the boundary topology or geometry,
        !! 3 the BEM/DtN build, 4 the multigrid, 5 a non-SPD Schur block.
        type(tensor_mesh_t), intent(in) :: mesh
        procedure(diffusion_coefficient) :: coefficient
        type(tensor_fem_bem_t), intent(out) :: fb
        integer, intent(out) :: status
        type(tensor_fem_bem_options_t), intent(in), optional :: options
        type(tensor_fem_bem_options_t) :: opt
        type(tensor_mesh_t) :: interior
        real(dp), allocatable :: x(:), w(:)
        real(dp) :: start
        integer :: flux_degree, n, info
        if (present(options)) opt = options
        call build_tensor_boundary_trace(mesh, fb%trace, status)
        if (status /= 0) return
        call gll_rule(mesh%degree, x, w, status)
        if (status /= 0) return
        flux_degree = opt%flux_degree
        if (flux_degree < 0) flux_degree = max(mesh%degree - 1, 0)
        call build_bem_tensor_space((x + 1)/2, flux_degree, fb%trace%panel_dof, &
                                    fb%space, status)
        if (status /= 0) then
            status = 1
            return
        end if
        start = omp_get_wtime()
        call assemble_parametric_bem_tensor(fb%trace%surface, fb%space, &
                                            opt%quadrature, fb%bem, status)
        fb%time_bem = omp_get_wtime() - start
        if (status /= 0) then
            status = 3
            return
        end if
        fb%bem%genus = fb%trace%genus
        start = omp_get_wtime()
        call build_exterior_dtn(fb%bem, fb%dtn, status)
        if (status == 0) call fb%dtn%matrix(fb%steklov, status)
        fb%time_dtn = omp_get_wtime() - start
        if (status /= 0) then
            status = 3
            return
        end if
        ! Interior block: the same operator with Dirichlet rows on the trace.
        start = omp_get_wtime()
        interior = mesh
        interior%boundary = .false.
        interior%boundary(fb%trace%trace_dof) = .true.
        fb%on_trace = interior%boundary
        call build_tensor_multigrid(interior, coefficient, fb%mg, status)
        fb%time_mg = omp_get_wtime() - start
        if (status /= 0) then
            status = 4
            return
        end if
        start = omp_get_wtime()
        n = fb%trace%ntrace
        fb%schur = opt%schur_scale*fb%steklov
        call dpotrf('L', n, fb%schur, n, info)
        fb%time_schur = omp_get_wtime() - start
        status = merge(0, 5, info == 0)
    end subroutine

    subroutine symmetric_matvec(a, x, y)
        real(dp), intent(in) :: a(:, :), x(:)
        real(dp), intent(out) :: y(:)
        integer :: i
        !$omp parallel do default(shared) schedule(static)
        do i = 1, size(x)
            y(i) = dot_product(a(:, i), x)
        end do
        !$omp end parallel do
    end subroutine

    subroutine apply_tensor_fem_bem(fb, x, y)
        !! y = (A + R^T S R) x.
        type(tensor_fem_bem_t), intent(in) :: fb
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: s(fb%trace%ntrace)
        call apply_tensor_operator(fb%mg%level(1)%op, x, y, eliminate=.false.)
        call symmetric_matvec(fb%steklov, x(fb%trace%trace_dof), s)
        y(fb%trace%trace_dof) = y(fb%trace%trace_dof) + s
    end subroutine

    subroutine precondition(fb, r, z, status)
        type(tensor_fem_bem_t), intent(inout) :: fb
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        real(dp), allocatable :: ri(:), y(:), t(:), zb(:, :)
        integer :: n
        n = fb%trace%ntrace
        allocate (ri(size(r)), y(size(r)), t(size(r)), zb(n, 1))
        ri = r
        where (fb%on_trace) ri = 0
        call tensor_vcycle(fb%mg, 1, ri, y, status)
        if (status /= 0) return
        call apply_tensor_operator(fb%mg%level(1)%op, y, t, eliminate=.false.)
        zb(:, 1) = r(fb%trace%trace_dof) - t(fb%trace%trace_dof)
        call dpotrs('L', n, 1, fb%schur, n, zb, n, status)
        if (status /= 0) return
        t = 0
        t(fb%trace%trace_dof) = zb(:, 1)
        call apply_tensor_operator(fb%mg%level(1)%op, t, ri, eliminate=.false.)
        where (fb%on_trace) ri = 0
        call tensor_vcycle(fb%mg, 1, ri, z, status)
        if (status /= 0) return
        z = y - z
        z(fb%trace%trace_dof) = zb(:, 1)
    end subroutine

    subroutine solve_tensor_fem_bem(fb, rhs, tol, x, info, maxiter)
        !! PCG on the coupled SPD system, starting from the supplied x.
        type(tensor_fem_bem_t), intent(inout) :: fb
        real(dp), intent(in) :: rhs(:), tol
        real(dp), intent(inout) :: x(:)
        type(tensor_solve_info_t), intent(out) :: info
        integer, intent(in), optional :: maxiter
        real(dp), allocatable :: r(:), z(:), p(:), ap(:)
        real(dp) :: rho, previous, alpha, denom, initial
        integer :: limit, it, status, n
        info%status = 1
        n = size(x)
        if (.not. allocated(fb%on_trace)) return
        if (size(rhs) /= n .or. size(fb%on_trace) /= n .or. tol <= 0) return
        limit = 300
        if (present(maxiter)) limit = maxiter
        allocate (r(n), z(n), p(n), ap(n))
        call apply_tensor_fem_bem(fb, x, ap)
        r = rhs - ap
        initial = max(norm2(rhs), tiny(1.0_dp))
        info%residual = norm2(r)/initial
        if (info%residual <= tol) then
            info%status = 0
            return
        end if
        previous = 0
        do it = 1, limit
            call precondition(fb, r, z, status)
            if (status /= 0) exit
            rho = dot_product(r, z)
            if (rho <= 0) exit
            if (it == 1) then
                p = z
            else
                p = z + (rho/previous)*p
            end if
            call apply_tensor_fem_bem(fb, p, ap)
            denom = dot_product(p, ap)
            if (denom <= 0) exit
            alpha = rho/denom
            x = x + alpha*p
            r = r - alpha*ap
            info%iterations = it
            info%residual = norm2(r)/initial
            if (info%residual <= tol) then
                call apply_tensor_fem_bem(fb, x, ap)
                r = rhs - ap
                info%residual = norm2(r)/initial
                if (info%residual <= tol) then
                    info%status = 0
                    exit
                end if
            end if
            previous = rho
        end do
    end subroutine

    subroutine tensor_fem_bem_rhs(fb, volume_load, rhs, status, neumann, dirichlet)
        !! rhs = volume_load + R^T <g_N, phi> + R^T S I(g_D). volume_load is a
        !! Neumann load (no boundary rows removed), e.g. tensor_load on an
        !! operator whose mesh has no Dirichlet nodes. g_D is interpolated at
        !! the trace nodes.
        type(tensor_fem_bem_t), intent(in) :: fb
        real(dp), intent(in) :: volume_load(:)
        real(dp), allocatable, intent(out) :: rhs(:)
        integer, intent(out) :: status
        procedure(boundary_function), optional :: neumann
        procedure(scalar_function), optional :: dirichlet
        real(dp), allocatable :: load(:), gd(:), s(:)
        integer :: i
        status = 1
        if (size(volume_load) /= size(fb%on_trace)) return
        rhs = volume_load
        if (present(neumann)) then
            call bem_tensor_boundary_load(fb%trace%surface, fb%space, neumann, &
                                          load, status)
            if (status /= 0) return
            rhs(fb%trace%trace_dof) = rhs(fb%trace%trace_dof) + load
        end if
        if (present(dirichlet)) then
            allocate (gd(fb%trace%ntrace), s(fb%trace%ntrace))
            do i = 1, fb%trace%ntrace
                gd(i) = dirichlet(fb%mg%level(1)%op%mesh%coordinates(:, &
                                  fb%trace%trace_dof(i)))
            end do
            call symmetric_matvec(fb%steklov, gd, s)
            rhs(fb%trace%trace_dof) = rhs(fb%trace%trace_dof) + s
        end if
        status = 0
    end subroutine

    subroutine tensor_fem_bem_exterior_flux(fb, exterior_trace, flux, status)
        !! Physical d_n w of the decaying exterior field with trace w, in the
        !! discontinuous Q_pf flux space of fb%space.
        type(tensor_fem_bem_t), intent(in) :: fb
        real(dp), intent(in) :: exterior_trace(:)
        real(dp), allocatable, intent(out) :: flux(:)
        integer, intent(out) :: status
        real(dp), allocatable :: load(:)
        allocate (flux(fb%space%nflux), load(fb%trace%ntrace))
        call fb%dtn%apply(exterior_trace, flux, load, status)
    end subroutine

    subroutine free_tensor_fem_bem(fb)
        type(tensor_fem_bem_t), intent(inout) :: fb
        call free_tensor_multigrid(fb%mg)
    end subroutine

    function tensor_h1_seminorm_error(op, x, exact) result(error)
        !! |u_h - u|_{H^1} with the operator's quadrature. The reference
        !! gradient error e = grad_ref u_h - J^T grad u is measured in the
        !! generated identity-coefficient metric w det J^{-1} J^{-T}.
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: x(:)
        procedure(gradient_function) :: exact
        real(dp) :: error
        real(dp), allocatable :: local(:), grad(:, :), t1(:), t2(:)
        real(dp) :: g(3, 3), m, det, identity(3, 3), e(3), gx(3)
        integer :: el, q, dir, ws, dim, mq, c(3)
        dim = op%mesh%dim
        mq = product(op%no)
        ws = maxval([op%ni, op%no])**3
        allocate (local(product(op%ni)), grad(mq, 3), t1(ws), t2(ws))
        identity = 0
        identity(1, 1) = 1; identity(2, 2) = 1; identity(3, 3) = 1
        error = 0
        grad = 0
        do el = 1, op%mesh%nelem
            local = x(op%mesh%dof(:, el))
            c = op%cell(:, el)
            do dir = 1, dim
                select case (dir)
                case (1)
                    call tensor_contract(op%ni, op%no, op%line(1)%d(:, :, c(1)), &
                                         op%line(2)%b(:, :, c(2)), &
                                         op%line(3)%b(:, :, c(3)), local, &
                                         grad(:, dir), t1, t2)
                case (2)
                    call tensor_contract(op%ni, op%no, op%line(1)%b(:, :, c(1)), &
                                         op%line(2)%d(:, :, c(2)), &
                                         op%line(3)%b(:, :, c(3)), local, &
                                         grad(:, dir), t1, t2)
                case default
                    call tensor_contract(op%ni, op%no, op%line(1)%b(:, :, c(1)), &
                                         op%line(2)%b(:, :, c(2)), &
                                         op%line(3)%d(:, :, c(3)), local, &
                                         grad(:, dir), t1, t2)
                end select
            end do
            do q = 1, mq
                gx = 0
                call exact(op%points(:, q, el), gx(1:dim))
                e = grad(q, :) - matmul(transpose(op%jac(:, :, q, el)), gx)
                call metric_at(op%jac(:, :, q, el), identity, op%weight(q, el), 0.0_dp, &
                               g, m, det)
                error = error + dot_product(e(1:dim), matmul(g(1:dim, 1:dim), e(1:dim)))
            end do
        end do
        error = sqrt(max(error, 0.0_dp))
    end function
end module
