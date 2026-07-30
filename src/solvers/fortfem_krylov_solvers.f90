module fortfem_krylov_solvers
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: bicgstab_impl

contains

    subroutine bicgstab_impl(A, b, x, precond_diag, precond_L, precond_U, &
            use_precond, tol, max_iter, tol_type, verbosity, &
            converged, iterations, final_resid)
        real(dp), intent(in) :: A(:, :), b(:)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in), optional :: precond_diag(:), precond_L(:, :), &
            precond_U(:, :)
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
        allocate (r(n), r0(n), p(n), v(n), s(n), t(n), z(n), y(n))

        r = b - matmul(A, x)
        r0 = r
        initial_norm = sqrt(dot_product(r, r))

        if (trim(tol_type) == "absolute") then
            tolerance = tol
        else
            tolerance = tol*initial_norm
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

            beta = (rho/rho_old)*(alpha/omega)
            p = r + beta*(p - omega*v)

            call bicgstab_matvec_step(A, p, v, z, use_precond, precond_diag, &
                precond_L, precond_U)

            alpha = rho/dot_product(r0, v)
            s = r - alpha*v

            residual_norm = sqrt(dot_product(s, s))
            if (residual_norm <= tolerance) then
                call bicgstab_update_early(x, z, p, alpha, use_precond)
                converged = .true.
                exit
            end if

            call bicgstab_matvec_step(A, s, t, y, use_precond, precond_diag, &
                precond_L, precond_U)

            omega = dot_product(t, s)/dot_product(t, t)

            call bicgstab_update_full(x, z, y, p, s, alpha, omega, use_precond)

            r = s - omega*t
            residual_norm = sqrt(dot_product(r, r))

            if (verbosity > 0) then
                write (*, '(A,I4,A,E12.5)') "BiCGSTAB iter ", iter, &
                    " residual: ", residual_norm
            end if

            if (residual_norm <= tolerance) then
                converged = .true.
                exit
            end if

            if (abs(omega) < 1.0e-14_dp) exit
        end do

        if (.not. converged) then
            if (residual_norm <= tolerance) then
                converged = .true.
            end if
        end if

        iterations = iter
        final_resid = residual_norm

        deallocate (r, r0, p, v, s, t, z, y)
    end subroutine bicgstab_impl

    subroutine bicgstab_matvec_step(A, input, output, precond_out, use_precond, &
            precond_diag, precond_L, precond_U)
        real(dp), intent(in) :: A(:, :), input(:)
        real(dp), intent(out) :: output(:), precond_out(:)
        logical, intent(in) :: use_precond
        real(dp), intent(in), optional :: precond_diag(:), precond_L(:, :), &
            precond_U(:, :)

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
            x = x + alpha*z
        else
            x = x + alpha*p
        end if
    end subroutine bicgstab_update_early

    subroutine bicgstab_update_full(x, z, y, p, s, alpha, omega, use_precond)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: z(:), y(:), p(:), s(:), alpha, omega
        logical, intent(in) :: use_precond

        if (use_precond) then
            x = x + alpha*z + omega*y
        else
            x = x + alpha*p + omega*s
        end if
    end subroutine bicgstab_update_full

    subroutine apply_ilu(L, U, r, z)
        real(dp), intent(in) :: L(:, :), U(:, :), r(:)
        real(dp), intent(out) :: z(:)

        real(dp), allocatable :: y(:)
        integer :: n, i, j

        n = size(r)
        allocate (y(n))

        do i = 1, n
            y(i) = r(i)
            do j = 1, i - 1
                y(i) = y(i) - L(i, j)*y(j)
            end do
        end do

        do i = n, 1, -1
            z(i) = y(i)
            do j = i + 1, n
                z(i) = z(i) - U(i, j)*z(j)
            end do
            z(i) = z(i)/U(i, i)
        end do

        deallocate (y)
    end subroutine apply_ilu

end module fortfem_krylov_solvers
