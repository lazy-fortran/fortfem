module fem_bem_oracle
    !! Exact transmission data on the through-axis O-grid torus: an interior
    !! field from point charges outside the solid torus (in the hole, above it
    !! and outside it) plus a non-harmonic quadratic, and a decaying exterior
    !! field from point charges inside the solid torus. Both fields are exact
    !! on any isoparametric approximation of the torus that keeps the charges
    !! on the stated side, so the discrete problem has no geometry error.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: interior_u, interior_gradient, interior_forcing, exterior_u
    public :: exterior_gradient, jump_dirichlet, jump_neumann, exterior_normal
    public :: unit_coefficient, zero_function
    real(dp), parameter :: outer(3, 3) = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 2.5_dp, 5.5_dp, 0.0_dp, 0.5_dp], [3, 3])
    real(dp), parameter :: outer_q(3) = [1.0_dp, -0.7_dp, 2.0_dp]
    real(dp), parameter :: inner(3, 3) = reshape([3.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, -3.0_dp, 0.2_dp, -3.0_dp, 0.0_dp, -0.1_dp], [3, 3])
    real(dp), parameter :: inner_q(3) = [1.0_dp, 0.5_dp, -0.8_dp]
contains
    subroutine unit_coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        k = 0
        k(1, 1) = 1; k(2, 2) = 1
        if (size(x) == 3) k(3, 3) = 1
        c = 0
    end subroutine

    function zero_function(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = 0*x(1)
    end function

    function charges(x, y, q) result(value)
        real(dp), intent(in) :: x(:), y(:, :), q(:)
        real(dp) :: value
        integer :: k
        value = 0
        do k = 1, size(q)
            value = value + q(k)/norm2(x(1:3) - y(:, k))
        end do
    end function

    subroutine charge_gradient(x, y, q, g)
        real(dp), intent(in) :: x(:), y(:, :), q(:)
        real(dp), intent(out) :: g(3)
        integer :: k
        g = 0
        do k = 1, size(q)
            g = g - q(k)*(x(1:3) - y(:, k))/norm2(x(1:3) - y(:, k))**3
        end do
    end subroutine

    function interior_u(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = charges(x, outer, outer_q) + 0.1_dp*x(3)**2 + &
                0.05_dp*(x(1)**2 - x(2)**2)
    end function

    subroutine interior_gradient(x, g)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: g(:)
        real(dp) :: c(3)
        call charge_gradient(x, outer, outer_q, c)
        g(1:3) = c + [0.1_dp*x(1), -0.1_dp*x(2), 0.2_dp*x(3)]
    end subroutine

    function interior_forcing(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = -0.2_dp + 0*x(1)
    end function

    function exterior_u(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = charges(x, inner, inner_q)
    end function

    subroutine exterior_gradient(x, g)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: g(:)
        call charge_gradient(x, inner, inner_q, g(1:3))
    end subroutine

    function jump_dirichlet(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = interior_u(x) - exterior_u(x)
    end function

    function jump_neumann(x, normal) result(value)
        real(dp), intent(in) :: x(3), normal(3)
        real(dp) :: value, gi(3), ge(3)
        call interior_gradient(x, gi)
        call exterior_gradient(x, ge)
        value = dot_product(gi - ge, normal)
    end function

    function exterior_normal(x, normal) result(value)
        real(dp), intent(in) :: x(3), normal(3)
        real(dp) :: value, ge(3)
        call exterior_gradient(x, ge)
        value = dot_product(ge, normal)
    end function
end module
