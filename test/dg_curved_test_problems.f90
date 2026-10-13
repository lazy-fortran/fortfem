module dg_curved_test_problems
    !! Manufactured complex solution and coefficients for the curved DG tests:
    !! u = (1 + i/2) sin(1.3 x + 0.4) cos(1.1 y) + i x^2 y / 4 with
    !! b = (1 + 0.2 x, 0.5 - 0.1 y) (div b = 0.1), a constant complex tensor
    !! K and a constant complex reaction c.
    use fortfem_kinds, only: dp
    use fortfem_dg_curved_2d, only: DG_FORM_CONSERVATIVE, dg_coefficients_t, &
        dg_curved_qpoint, dg_curved_space_t, dg_data_t
    implicit none
    private

    complex(dp), parameter :: ci = (0.0_dp, 1.0_dp)

    type, extends(dg_coefficients_t), public :: mms_coef_t
        complex(dp) :: kk(2, 2) = (0.0_dp, 0.0_dp)
        complex(dp) :: c = (0.0_dp, 0.0_dp)
        real(dp) :: bscale = 1.0_dp
    contains
        procedure :: evaluate => mms_coefficients
    end type mms_coef_t

    type, extends(dg_data_t), public :: mms_data_t
        type(mms_coef_t) :: coef
        integer :: form = DG_FORM_CONSERVATIVE
        logical :: unit = .false.   !! u = 1 boundary data, f = 0
    contains
        procedure :: source => mms_source
        procedure :: boundary => mms_boundary
    end type mms_data_t

    public :: exact_u, anisotropic_tensor, l2_error

contains

    pure subroutine velocity(x, b)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: b(2)

        b = [1.0_dp + 0.2_dp*x(1), 0.5_dp - 0.1_dp*x(2)]
    end subroutine velocity

    pure subroutine exact_u(x, u, g, h)
        !! Value, gradient and Hessian of the manufactured solution.
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: u, g(2), h(2, 2)

        complex(dp), parameter :: a = (1.0_dp, 0.5_dp)
        real(dp) :: s, c, sy, cy

        s = sin(1.3_dp*x(1) + 0.4_dp)
        c = cos(1.3_dp*x(1) + 0.4_dp)
        sy = sin(1.1_dp*x(2))
        cy = cos(1.1_dp*x(2))
        u = a*s*cy + ci*0.25_dp*x(1)**2*x(2)
        g(1) = a*1.3_dp*c*cy + ci*0.5_dp*x(1)*x(2)
        g(2) = -a*1.1_dp*s*sy + ci*0.25_dp*x(1)**2
        h(1, 1) = -a*1.69_dp*s*cy + ci*0.5_dp*x(2)
        h(2, 2) = -a*1.21_dp*s*cy
        h(1, 2) = -a*1.43_dp*c*sy + ci*0.5_dp*x(1)
        h(2, 1) = h(1, 2)
    end subroutine exact_u

    pure function anisotropic_tensor(scale, theta, eps) result(kk)
        !! scale * R(theta) diag(1, eps) R(theta)^T.
        complex(dp), intent(in) :: scale
        real(dp), intent(in) :: theta, eps
        complex(dp) :: kk(2, 2)

        real(dp) :: r(2, 2), d(2, 2)

        r = reshape([cos(theta), sin(theta), -sin(theta), cos(theta)], [2, 2])
        d = reshape([1.0_dp, 0.0_dp, 0.0_dp, eps], [2, 2])
        kk = scale*matmul(r, matmul(d, transpose(r)))
    end function anisotropic_tensor

    subroutine mms_coefficients(self, t, x, kk, b, c, status)
        class(mms_coef_t), intent(in) :: self
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
            kk(:, :, q) = self%kk
            call velocity(x(:, q), b(:, q))
            b(:, q) = self%bscale*b(:, q)
            c(q) = self%c
        end do
        status = 0
    end subroutine mms_coefficients

    subroutine mms_source(self, t, x, f, status)
        class(mms_data_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        complex(dp), intent(out) :: f(:)
        integer, intent(out) :: status

        complex(dp) :: u, g(2), h(2, 2)
        real(dp) :: b(2)
        integer :: q

        associate (unused => t)
        end associate
        status = 0
        if (self%unit) then
            f = (0.0_dp, 0.0_dp)
            return
        end if
        do q = 1, size(x, 2)
            call exact_u(x(:, q), u, g, h)
            call velocity(x(:, q), b)
            b = self%coef%bscale*b
            f(q) = b(1)*g(1) + b(2)*g(2) + self%coef%c*u - sum(self%coef%kk*h)
            if (self%form == DG_FORM_CONSERVATIVE) then
                f(q) = f(q) + 0.1_dp*self%coef%bscale*u
            end if
        end do
    end subroutine mms_source

    subroutine mms_boundary(self, t, x, f, status)
        class(mms_data_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        complex(dp), intent(out) :: f(:)
        integer, intent(out) :: status

        complex(dp) :: g(2), h(2, 2)
        integer :: q

        associate (unused => t)
        end associate
        status = 0
        if (self%unit) then
            f = (1.0_dp, 0.0_dp)
            return
        end if
        do q = 1, size(x, 2)
            call exact_u(x(:, q), f(q), g, h)
        end do
    end subroutine mms_boundary

    function l2_error(space, u) result(err)
        !! || u_h - u ||_L2 over the curved mesh.
        type(dg_curved_space_t), intent(in) :: space
        complex(dp), intent(in) :: u(:)
        real(dp) :: err

        real(dp) :: x(2), w, val(space%nl), grad(2, space%nl)
        complex(dp) :: ue, g(2), h(2, 2), uh
        integer :: t, q, st, base

        err = 0.0_dp
        do t = 1, space%nt
            base = (t - 1)*space%nl
            do q = 1, space%nq
                call dg_curved_qpoint(space, t, q, x, w, val, grad, st)
                if (st /= 0) error stop "l2_error: bad geometry"
                call exact_u(x, ue, g, h)
                uh = sum(val*u(base + 1:base + space%nl))
                err = err + w*abs(uh - ue)**2
            end do
        end do
        err = sqrt(err)
    end function l2_error

end module dg_curved_test_problems
