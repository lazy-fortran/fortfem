module fortfem_codegen_rational
    !! Exact rational scalars for reference-table generators. Every value is a
    !! canonical FLINT fmpq string owned by FortSym; there is no floating
    !! arithmetic. A value is rounded to double only once, by MPFR, when it is
    !! emitted as a literal.
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use fortsym_exact, only: exact_add, exact_div, exact_mul, exact_normalize, &
        exact_sub, exact_to_real
    use fortsym_string, only: chars, str, str_t
    implicit none

    private

    public :: q_t
    public :: q, q_frac, q_text, q_zero, q_literal, q_real
    public :: q_factorial, q_binomial
    public :: operator(+), operator(-), operator(*), operator(/)
    public :: operator(==)

    type :: q_t
        type(str_t) :: value
    end type q_t

    interface operator(+)
        module procedure add_qq
    end interface
    interface operator(-)
        module procedure sub_qq, neg_q
    end interface
    interface operator(*)
        module procedure mul_qq, mul_iq
    end interface
    interface operator(/)
        module procedure div_qq
    end interface
    interface operator(==)
        module procedure eq_qq
    end interface

contains

    function q(value) result(r)
        integer, intent(in) :: value
        type(q_t) :: r

        character(32) :: buffer

        write(buffer, "(i0)") value
        r = q_text(trim(buffer))
    end function q

    function q_frac(numerator, denominator) result(r)
        integer, intent(in) :: numerator, denominator
        type(q_t) :: r

        r = q(numerator) / q(denominator)
    end function q_frac

    function q_text(text) result(r)
        character(*), intent(in) :: text
        type(q_t) :: r

        logical :: ok

        r%value = exact_normalize(text, ok)
        if (.not. ok) error stop "invalid exact rational text"
    end function q_text

    pure logical function q_zero(a)
        type(q_t), intent(in) :: a

        q_zero = chars(a%value) == "0"
    end function q_zero

    function q_real(a) result(value)
        !! Correctly rounded (round-to-nearest) binary64 value of a.
        type(q_t), intent(in) :: a
        real(real64) :: value

        logical :: ok

        value = exact_to_real(chars(a%value), ok)
        if (.not. ok) error stop "exact rational is outside binary64"
    end function q_real

    function q_literal(a) result(text)
        !! Shortest decimal that reads back as the correctly rounded binary64
        !! value, so the compiler reproduces it without a second rounding.
        type(q_t), intent(in) :: a
        character(:), allocatable :: text

        character(48) :: buffer, form
        character(:), allocatable :: mantissa
        real(real64) :: value, parsed
        integer :: digits, position, exponent

        if (q_zero(a)) then
            text = "0.0_dp"
            return
        end if
        value = q_real(a)
        do digits = 1, 17
            write(form, "(a,i0,a,i0,a)") "(es", digits + 10, ".", &
                max(digits - 1, 1), "e3)"
            write(buffer, form) value
            read(buffer, *) parsed
            if (parsed == value) exit
        end do
        if (parsed /= value) error stop "binary64 literal does not round-trip"
        buffer = adjustl(buffer)
        position = index(buffer, "E")
        mantissa = buffer(:position - 1)
        read(buffer(position + 1:), *) exponent
        do while (len(mantissa) > 3)
            if (mantissa(len(mantissa):len(mantissa)) /= "0") exit
            if (mantissa(len(mantissa) - 1:len(mantissa) - 1) == ".") exit
            mantissa = mantissa(:len(mantissa) - 1)
        end do
        if (exponent == 0) then
            text = mantissa//"_dp"
        else
            write(buffer, "(i0)") exponent
            text = mantissa//"e"//trim(buffer)//"_dp"
        end if
    end function q_literal

    function q_factorial(n) result(r)
        integer, intent(in) :: n
        type(q_t) :: r

        integer, parameter :: cached = 64
        type(q_t), save :: table(0:cached)
        logical, save :: ready = .false.
        integer :: i

        if (.not. ready) then
            table(0) = q(1)
            do i = 1, cached
                table(i) = table(i - 1) * q(i)
            end do
            ready = .true.
        end if
        if (n < 0 .or. n > cached) error stop "factorial argument out of range"
        r = table(n)
    end function q_factorial

    function q_binomial(n, k) result(r)
        integer, intent(in) :: n, k
        type(q_t) :: r

        if (k < 0 .or. k > n) then
            r = q(0)
            return
        end if
        r = q_factorial(n) / (q_factorial(k) * q_factorial(n - k))
    end function q_binomial

    function add_qq(a, b) result(r)
        type(q_t), intent(in) :: a, b
        type(q_t) :: r

        logical :: ok

        r%value = exact_add(chars(a%value), chars(b%value), ok)
        if (.not. ok) error stop "exact rational addition failed"
    end function add_qq

    function sub_qq(a, b) result(r)
        type(q_t), intent(in) :: a, b
        type(q_t) :: r

        logical :: ok

        r%value = exact_sub(chars(a%value), chars(b%value), ok)
        if (.not. ok) error stop "exact rational subtraction failed"
    end function sub_qq

    function neg_q(a) result(r)
        type(q_t), intent(in) :: a
        type(q_t) :: r

        r = q(0) - a
    end function neg_q

    function mul_qq(a, b) result(r)
        type(q_t), intent(in) :: a, b
        type(q_t) :: r

        logical :: ok

        r%value = exact_mul(chars(a%value), chars(b%value), ok)
        if (.not. ok) error stop "exact rational multiplication failed"
    end function mul_qq

    function mul_iq(i, b) result(r)
        integer, intent(in) :: i
        type(q_t), intent(in) :: b
        type(q_t) :: r

        r = q(i) * b
    end function mul_iq

    function div_qq(a, b) result(r)
        type(q_t), intent(in) :: a, b
        type(q_t) :: r

        logical :: ok

        if (q_zero(b)) error stop "exact rational division by zero"
        r%value = exact_div(chars(a%value), chars(b%value), ok)
        if (.not. ok) error stop "exact rational division failed"
    end function div_qq

    pure logical function eq_qq(a, b)
        type(q_t), intent(in) :: a, b

        eq_qq = chars(a%value) == chars(b%value)
    end function eq_qq

end module fortfem_codegen_rational
