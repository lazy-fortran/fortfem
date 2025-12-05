module fortfem_krylov_solvers
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: gmres_impl, bicgstab_impl

contains

    subroutine gmres_impl(A, b, x, tol, max_iter, restart, tol_type, verbosity, &
        converged, iterations, restarts_count, final_resid)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: tol
        integer, intent(in) :: max_iter, restart, verbosity
        character(len=*), intent(in) :: tol_type
        logical, intent(out) :: converged
        integer, intent(out) :: iterations, restarts_count
        real(dp), intent(out) :: final_resid

        real(dp), allocatable :: V(:,:), H(:,:), r(:), w(:), y(:), c(:), s(:), g(:)
        real(dp) :: beta, residual_norm, initial_norm, tolerance
        integer :: n, m, m_used, total_iter

        n = size(x)
        m = min(restart, n)

        allocate(V(n, m+1), H(m+1, m), r(n), w(n), y(m), c(m), s(m), g(m+1))

        r = b - matmul(A, x)
        initial_norm = sqrt(dot_product(r, r))
        residual_norm = initial_norm

        if (trim(tol_type) == "absolute") then
            tolerance = tol
        else
            tolerance = min(tol, tol * initial_norm)
        end if

        total_iter = 0
        restarts_count = 0
        converged = .false.

        do while (total_iter < max_iter .and. .not. converged)
            beta = sqrt(dot_product(r, r))
            residual_norm = beta

            if (beta <= tolerance) then
                converged = .true.
                exit
            end if

            V(:, 1) = r / beta
            g = 0.0_dp
            g(1) = beta
            H = 0.0_dp

            call gmres_arnoldi_cycle(A, V, H, c, s, g, m, total_iter, max_iter, &
                tolerance, verbosity, converged, residual_norm, m_used)

            call gmres_solve_update(H, V, g, y, x, m_used)

            if (.not. converged) then
                r = b - matmul(A, x)
                residual_norm = sqrt(dot_product(r, r))
                restarts_count = restarts_count + 1
            end if
        end do

        iterations = total_iter
        final_resid = residual_norm

        deallocate(V, H, r, w, y, c, s, g)
    end subroutine gmres_impl

    subroutine gmres_arnoldi_cycle(A, V, H, c, s, g, m, total_iter, max_iter, &
        tolerance, verbosity, converged, residual_norm, m_used)
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(inout) :: V(:,:), H(:,:), c(:), s(:), g(:)
        integer, intent(in) :: m, max_iter, verbosity
        integer, intent(inout) :: total_iter
        real(dp), intent(in) :: tolerance
        logical, intent(out) :: converged
        real(dp), intent(inout) :: residual_norm
        integer, intent(out) :: m_used

        real(dp), allocatable :: w(:)
        integer :: n, j, i

        n = size(V, 1)
        allocate(w(n))
        converged = .false.
        m_used = 0

        do j = 1, m
            total_iter = total_iter + 1
            if (total_iter > max_iter) exit

            w = matmul(A, V(:, j))

            do i = 1, j
                H(i, j) = dot_product(w, V(:, i))
                w = w - H(i, j) * V(:, i)
            end do

            H(j+1, j) = sqrt(dot_product(w, w))
            m_used = j

            if (H(j+1, j) < 1.0e-14_dp) exit

            V(:, j+1) = w / H(j+1, j)

            call apply_givens_rotations(H, c, s, g, j, residual_norm)

            if (verbosity > 0) then
                write(*,'(A,I4,A,E12.5)') "GMRES iter ", total_iter, &
                    " residual: ", residual_norm
            end if

            if (residual_norm <= tolerance) then
                converged = .true.
                exit
            end if
        end do

        deallocate(w)
    end subroutine gmres_arnoldi_cycle

    subroutine apply_givens_rotations(H, c, s, g, j, residual_norm)
        real(dp), intent(inout) :: H(:,:), c(:), s(:), g(:)
        integer, intent(in) :: j
        real(dp), intent(out) :: residual_norm

        real(dp) :: temp
        integer :: k

        do k = 1, j-1
            temp = c(k) * H(k, j) + s(k) * H(k+1, j)
            H(k+1, j) = -s(k) * H(k, j) + c(k) * H(k+1, j)
            H(k, j) = temp
        end do

        if (abs(H(j+1, j)) < 1.0e-14_dp) then
            c(j) = 1.0_dp
            s(j) = 0.0_dp
        else
            temp = sqrt(H(j, j)**2 + H(j+1, j)**2)
            c(j) = H(j, j) / temp
            s(j) = H(j+1, j) / temp
        end if

        H(j, j) = c(j) * H(j, j) + s(j) * H(j+1, j)
        H(j+1, j) = 0.0_dp

        g(j+1) = -s(j) * g(j)
        g(j) = c(j) * g(j)

        residual_norm = abs(g(j+1))
    end subroutine apply_givens_rotations

    subroutine gmres_solve_update(H, V, g, y, x, m)
        real(dp), intent(in) :: H(:,:), V(:,:), g(:)
        real(dp), intent(inout) :: y(:), x(:)
        integer, intent(in) :: m

        integer :: i, j

        do i = m, 1, -1
            y(i) = g(i)
            do j = i+1, m
                y(i) = y(i) - H(i, j) * y(j)
            end do
            y(i) = y(i) / H(i, i)
        end do

        do j = 1, m
            x = x + y(j) * V(:, j)
        end do
    end subroutine gmres_solve_update

    subroutine bicgstab_impl(A, b, x, precond_diag, precond_L, precond_U, &
        use_precond, tol, max_iter, tol_type, verbosity, &
        converged, iterations, final_resid)
        real(dp), intent(in) :: A(:,:), b(:)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in), optional :: precond_diag(:), precond_L(:,:), precond_U(:,:)
        logical, intent(in) :: use_precond
        real(dp), intent(in) :: tol
        integer, intent(in) :: max_iter, verbosity
        character(len=*), intent(in) :: tol_type
        logical, intent(out) :: converged
        integer, intent(out) :: iterations
        real(dp), intent(out) :: final_resid

        real(dp), allocatable :: r(:), r0(:), p(:), v(:), s(:), t(:), z(:), y(:)
        real(dp) :: rho, rho_old, alpha, omega, beta
        real(dp) :: residual_norm, initial_norm, tolerance
        integer :: n, iter

        n = size(x)
        allocate(r(n), r0(n), p(n), v(n), s(n), t(n), z(n), y(n))

        r = b - matmul(A, x)
        r0 = r
        initial_norm = sqrt(dot_product(r, r))

        if (trim(tol_type) == "absolute") then
            tolerance = tol
        else
            tolerance = tol * initial_norm
        end if

        p = r
        rho = 1.0_dp
        alpha = 1.0_dp
        omega = 1.0_dp

        converged = .false.

        do iter = 1, max_iter
            rho_old = rho
            rho = dot_product(r0, r)

            if (abs(rho) < 1.0e-14_dp) exit

            beta = (rho / rho_old) * (alpha / omega)
            p = r + beta * (p - omega * v)

            call bicgstab_matvec_step(A, p, v, z, use_precond, precond_diag, &
                precond_L, precond_U)

            alpha = rho / dot_product(r0, v)
            s = r - alpha * v

            residual_norm = sqrt(dot_product(s, s))
            if (residual_norm <= tolerance) then
                call bicgstab_update_early(x, z, p, alpha, use_precond)
                converged = .true.
                exit
            end if

            call bicgstab_matvec_step(A, s, t, y, use_precond, precond_diag, &
                precond_L, precond_U)

            omega = dot_product(t, s) / dot_product(t, t)

            call bicgstab_update_full(x, z, y, p, s, alpha, omega, use_precond)

            r = s - omega * t
            residual_norm = sqrt(dot_product(r, r))

            if (verbosity > 0) then
                write(*,'(A,I4,A,E12.5)') "BiCGSTAB iter ", iter, &
                    " residual: ", residual_norm
            end if

            if (residual_norm <= tolerance) then
                converged = .true.
                exit
            end if

            if (abs(omega) < 1.0e-14_dp) exit
        end do

        iterations = iter
        final_resid = residual_norm

        deallocate(r, r0, p, v, s, t, z, y)
    end subroutine bicgstab_impl

    subroutine bicgstab_matvec_step(A, input, output, precond_out, use_precond, &
        precond_diag, precond_L, precond_U)
        real(dp), intent(in) :: A(:,:), input(:)
        real(dp), intent(out) :: output(:), precond_out(:)
        logical, intent(in) :: use_precond
        real(dp), intent(in), optional :: precond_diag(:), precond_L(:,:), precond_U(:,:)

        if (use_precond .and. present(precond_L) .and. present(precond_U)) then
            call apply_ilu(precond_L, precond_U, input, precond_out)
            output = matmul(A, precond_out)
        else
            output = matmul(A, input)
            precond_out = input
        end if
    end subroutine bicgstab_matvec_step

    subroutine bicgstab_update_early(x, z, p, alpha, use_precond)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: z(:), p(:), alpha
        logical, intent(in) :: use_precond

        if (use_precond) then
            x = x + alpha * z
        else
            x = x + alpha * p
        end if
    end subroutine bicgstab_update_early

    subroutine bicgstab_update_full(x, z, y, p, s, alpha, omega, use_precond)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: z(:), y(:), p(:), s(:), alpha, omega
        logical, intent(in) :: use_precond

        if (use_precond) then
            x = x + alpha * z + omega * y
        else
            x = x + alpha * p + omega * s
        end if
    end subroutine bicgstab_update_full

    subroutine apply_ilu(L, U, r, z)
        real(dp), intent(in) :: L(:,:), U(:,:), r(:)
        real(dp), intent(out) :: z(:)

        real(dp), allocatable :: y(:)
        integer :: n, i, j

        n = size(r)
        allocate(y(n))

        do i = 1, n
            y(i) = r(i)
            do j = 1, i-1
                y(i) = y(i) - L(i, j) * y(j)
            end do
        end do

        do i = n, 1, -1
            z(i) = y(i)
            do j = i+1, n
                z(i) = z(i) - U(i, j) * z(j)
            end do
            z(i) = z(i) / U(i, i)
        end do

        deallocate(y)
    end subroutine apply_ilu

end module fortfem_krylov_solvers
