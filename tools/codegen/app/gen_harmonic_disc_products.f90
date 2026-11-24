program gen_harmonic_disc_products
    !! Local kernels of fortfem_harmonic_disc (doc/HARMONIC_DISC_MAP.md).
    !!
    !! Complex arithmetic is written in real and imaginary parts.
    !!   cauchy_terms: q = w/(y - z), p = w/(y - z)^2, the summands of the
    !!     barycentric Cauchy formula for v and v'.
    !!   nystrom_offdiag: w/(y_j - y_k)/(2 pi i), the double-layer entry.
    !!   nystrom_self: y''/(2 y') h/(2 pi i), its limit on the diagonal.
    !!   interior_limit_self: d mu h/(2 pi i), the diagonal term of the
    !!     interior limit of the Cauchy integral of the density.
    !!   boundary_weight: |y'| (1 + (scale kappa)^2)^(1/3), the density of
    !!     the boundary correspondence in the curve parameter, kappa =
    !!     Im(conj(y') y'')/|y'|^3 the curvature.
    !!   spline_jacobian_determinant: det of the 2x2 geometry Jacobian.
    !! The certificate of a spline space needs, per cell and direction, the
    !! Bernstein coefficients of a polynomial of degree m - 1, m = 2 p, from
    !! its values at the m Chebyshev nodes x_i = (1 - cos(pi (i - 1/2)/m))/2.
    !! The inverse of B_ij = binom(m-1, j-1) x_i^(j-1) (1 - x_i)^(m-j) is
    !! computed in exact rational arithmetic and rounded once. The nodes
    !! come from a 100-digit cosine series. Above degree five the inverse
    !! uses the rounded binary64 nodes to keep the rational solve in range.
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/), operator(**), pi_expr
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: generated_path, fortsym_revision
    use fortfem_codegen_rational, only: q_t, q, q_text, q_zero, q_literal, &
        q_real, q_binomial, operator(+), operator(-), operator(*), &
        operator(/)
    implicit none

    integer, parameter :: max_degree = 7
    integer :: unit

    open(newunit=unit, file=generated_path("fortfem_harmonic_disc_products.f90"), &
        status="replace", action="write")
    call emit_products(unit)
    call emit_bernstein(unit)
    close(unit)
    write(*, "(a)") "Harmonic disc kernels written"

contains

    subroutine emit_products(unit)
        integer, intent(in) :: unit

        type(expr_t) :: yr, yi, zr, zi, wr, wi, rr, ri, d, qr, qi, pr, pim
        type(expr_t) :: yjr, yji, ykr, yki, dyr, dyi, ddyr, ddyi, h, dmu
        type(expr_t) :: ar, ai, br, bi, two_pi, scale, speed2, cross

        yr = sym("yr"); yi = sym("yi"); zr = sym("zr"); zi = sym("zi")
        wr = sym("wr"); wi = sym("wi")
        two_pi = num(2)*pi_expr()
        ! q = w/(y - z) and p = q/(y - z)
        rr = yr - zr
        ri = yi - zi
        d = rr*rr + ri*ri
        qr = (wr*rr + wi*ri)/d
        qi = (wi*rr - wr*ri)/d
        pr = (qr*rr + qi*ri)/d
        pim = (qi*rr - qr*ri)/d
        call emit("harmonic_disc_cauchy_terms", &
            [str("yr"), str("yi"), str("zr"), str("zi"), str("wr"), str("wi")], &
            [qr, qi, pr, pim], [str("qr"), str("qi"), str("pr"), str("pim")])
        ! w/(y_j - y_k)/(2 pi i) = -i q/(2 pi)
        yjr = sym("yjr"); yji = sym("yji"); ykr = sym("ykr"); yki = sym("yki")
        rr = yjr - ykr
        ri = yji - yki
        d = rr*rr + ri*ri
        qr = (wr*rr + wi*ri)/d
        qi = (wi*rr - wr*ri)/d
        call emit("harmonic_disc_nystrom_offdiag", &
            [str("wr"), str("wi"), str("yjr"), str("yji"), str("ykr"), str("yki")], &
            [qi/two_pi, -qr/two_pi], [str("cr"), str("ci")])
        ! ddy/(2 dy) h/(2 pi i)
        ddyr = sym("ddyr"); ddyi = sym("ddyi"); dyr = sym("dyr"); dyi = sym("dyi")
        h = sym("h")
        d = num(2)*(dyr*dyr + dyi*dyi)
        ar = (ddyr*dyr + ddyi*dyi)/d*h
        ai = (ddyi*dyr - ddyr*dyi)/d*h
        call emit("harmonic_disc_nystrom_self", &
            [str("ddyr"), str("ddyi"), str("dyr"), str("dyi"), str("h")], &
            [ai/two_pi, -ar/two_pi], [str("cr"), str("ci")])
        ! d mu h/(2 pi i)
        dmu = sym("dmu")
        br = num(0)*dmu
        bi = -dmu*h/two_pi
        call emit("harmonic_disc_interior_limit_self", [str("dmu"), str("h")], &
            [br, bi], [str("er"), str("ei")])
        ! |y'| (1 + (scale kappa)^2)^(1/3), kappa^2 = cross^2/|y'|^6
        scale = sym("scale")
        speed2 = dyr*dyr + dyi*dyi
        cross = dyr*ddyi - dyi*ddyr
        call emit("harmonic_disc_boundary_weight", &
            [str("dyr"), str("dyi"), str("ddyr"), str("ddyi"), str("scale")], &
            [speed2**(num(1)/num(2))*(num(1) + scale*scale*cross*cross/ &
            (speed2*speed2*speed2))**(num(1)/num(3))], [str("weight")])
        ! jet(:, 2) = (dx/du, dy/du), jet(:, 3) = (dx/dv, dy/dv)
        call emit("spline_jacobian_determinant", &
            [str("dxdu"), str("dydu"), str("dxdv"), str("dydv")], &
            [sym("dxdu")*sym("dydv") - sym("dxdv")*sym("dydu")], [str("det")])
    end subroutine emit_products

    subroutine emit(name, args, roots, outputs)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)
        type(expr_t), intent(in) :: roots(:)

        type(kernel_spec_t) :: spec
        character(:), allocatable :: code

        spec%name = str(name)
        spec%module_name = str("fortfem_generated_"//name)
        spec%generator = str("gen_harmonic_disc_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%args = args
        spec%outputs = outputs
        code = chars(emit_kernel(roots, spec))
        write(unit, "(a)") code(:len(code) - 1)
    end subroutine emit

    subroutine emit_bernstein(unit)
        integer, intent(in) :: unit

        type(q_t), allocatable :: nodes(:), inverse(:, :)
        integer :: p, m

        write(unit, "(a)") "! Generated by fortsym. Do not edit."
        write(unit, "(a)") "! Generator: gen_harmonic_disc_products"
        write(unit, "(a)") "! Generator revision: "//fortsym_revision()
        write(unit, "(a)") "! Regenerate with: cd tools/codegen && ./generate.sh"
        write(unit, "(a)") "! Every value is an exact rational rounded once to binary64."
        write(unit, "(a)") ""
        write(unit, "(a)") "module fortfem_generated_spline_bernstein_interpolation"
        write(unit, "(a)") "    use, intrinsic :: iso_fortran_env, only: dp => real64"
        write(unit, "(a)") "    implicit none"
        write(unit, "(a)") "    private"
        write(unit, "(a)") "    public :: spline_bernstein_interpolation"
        write(unit, "(a)") "    public :: spline_bernstein_max_degree"
        write(unit, "(a,i0)") "    integer, parameter :: spline_bernstein_max_degree = ", &
            max_degree
        write(unit, "(a)") "contains"
        write(unit, "(a)") ""
        write(unit, "(a)") "    subroutine spline_bernstein_interpolation(degree, nodes, inverse, status)"
        write(unit, "(a)") "        !! Chebyshev nodes on [0, 1] and the inverse of the Bernstein matrix"
        write(unit, "(a)") "        !! B(i, j) = binom(m-1, j-1) x_i^(j-1) (1-x_i)^(m-j), m = 2*degree:"
        write(unit, "(a)") "        !! Bernstein coefficients = inverse * values at the nodes."
        write(unit, "(a)") "        !! status 1: degree outside 1..spline_bernstein_max_degree."
        write(unit, "(a)") "        integer, intent(in) :: degree"
        write(unit, "(a)") "        real(dp), allocatable, intent(out) :: nodes(:), inverse(:, :)"
        write(unit, "(a)") "        integer, intent(out) :: status"
        write(unit, "(a)") ""
        write(unit, "(a)") "        status = 1"
        write(unit, "(a)") "        select case (degree)"
        do p = 1, max_degree
            m = 2*p
            call bernstein_inverse(m, nodes, inverse)
            write(unit, "(a,i0,a)") "        case (", p, ")"
            write(unit, "(a,i0,a)") "            allocate(nodes(", m, "))"
            write(unit, "(a,i0,a,i0,a)") "            allocate(inverse(", m, ", ", m, "))"
            call array_assignment("nodes", nodes)
            call array_assignment("inverse", reshape(inverse, [m*m]), m)
        end do
        write(unit, "(a)") "        case default"
        write(unit, "(a)") "            return"
        write(unit, "(a)") "        end select"
        write(unit, "(a)") "        status = 0"
        write(unit, "(a)") "    end subroutine spline_bernstein_interpolation"
        write(unit, "(a)") ""
        write(unit, "(a)") "end module fortfem_generated_spline_bernstein_interpolation"
    end subroutine emit_bernstein

    subroutine array_assignment(lhs, values, square)
        !! Column-major array constructor wrapped at 88 columns.
        character(*), intent(in) :: lhs
        type(q_t), intent(in) :: values(:)
        integer, intent(in), optional :: square

        character(:), allocatable :: current, item
        integer :: i

        if (present(square)) then
            write(unit, "(a)") "            "//lhs//" = reshape([ &"
        else
            write(unit, "(a)") "            "//lhs//" = [ &"
        end if
        current = "                "
        do i = 1, size(values)
            item = q_literal(values(i))
            if (i < size(values)) then
                item = item//","
            else if (present(square)) then
                item = item//"], [ &"
            else
                item = item//"]"
            end if
            if (len(current) + len(item) + 3 > 88) then
                write(unit, "(a)") current//" &"
                current = "                "
            end if
            if (len(current) > 16) current = current//" "
            current = current//item
        end do
        write(unit, "(a)") current
        if (present(square)) then
            write(unit, "(a,i0,a,i0,a)") "                ", square, ", ", square, "])"
        end if
    end subroutine array_assignment

    subroutine bernstein_inverse(m, nodes, inverse)
        integer, intent(in) :: m
        type(q_t), allocatable, intent(out) :: nodes(:), inverse(:, :)

        type(q_t), allocatable :: work(:, :)
        type(q_t) :: pi_q, theta, factor, pivot_value, swap
        integer :: i, j, k, r, pivot

        pi_q = q_text("31415926535897932384626433832795028841971693993751058209749445923078"// &
            "164062862089986280348253421170679/1"//repeat("0", 100))
        if (q_real(pi_q) /= acos(-1.0d0)) error stop "pi constant"
        allocate(nodes(m), inverse(m, m), work(m, 2*m))
        do i = 1, m
            theta = pi_q*q(2*i - 1)/q(2*m)
            nodes(i) = (q(1) - cos_series(theta))/q(2)
            if (m > 10) nodes(i) = binary64_node(q_real(nodes(i)))
        end do
        do i = 1, m
            do j = 1, m
                work(i, j) = q_binomial(m - 1, j - 1)*power(nodes(i), j - 1)* &
                    power(q(1) - nodes(i), m - j)
                work(i, m + j) = q(0)
            end do
            work(i, m + i) = q(1)
        end do
        do k = 1, m
            pivot = 0
            do r = k, m
                if (.not. q_zero(work(r, k))) then
                    pivot = r
                    exit
                end if
            end do
            if (pivot == 0) error stop "singular Bernstein matrix"
            do j = 1, 2*m
                swap = work(k, j)
                work(k, j) = work(pivot, j)
                work(pivot, j) = swap
            end do
            pivot_value = work(k, k)
            do j = 1, 2*m
                work(k, j) = work(k, j)/pivot_value
            end do
            do r = 1, m
                if (r == k) cycle
                factor = work(r, k)
                if (q_zero(factor)) cycle
                do j = 1, 2*m
                    work(r, j) = work(r, j) - factor*work(k, j)
                end do
            end do
        end do
        inverse = work(:, m + 1:2*m)
    end subroutine bernstein_inverse

    function binary64_node(value) result(node)
        real(real64), intent(in) :: value
        type(q_t) :: node
        character(32) :: significand
        integer :: k, bits

        bits = digits(value)
        write(significand, '(i0)') int(scale(fraction(value), bits), int64)
        node = q_text(trim(significand))
        do k = 1, bits - exponent(value)
            node = node/q(2)
        end do
    end function binary64_node

    function cos_series(x) result(c)
        type(q_t), intent(in) :: x
        type(q_t) :: c

        type(q_t) :: term
        integer :: k

        c = q(1)
        term = q(1)
        do k = 1, 60
            term = term*x*x/q((2*k - 1)*(2*k))
            if (mod(k, 2) == 1) then
                c = c - term
            else
                c = c + term
            end if
        end do
    end function cos_series

    function power(x, n) result(r)
        type(q_t), intent(in) :: x
        integer, intent(in) :: n
        type(q_t) :: r

        integer :: i

        r = q(1)
        do i = 1, n
            r = r*x
        end do
    end function power

end program gen_harmonic_disc_products
