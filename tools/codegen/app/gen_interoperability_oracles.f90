program gen_interoperability_oracles
    use fortsym_arena, only: arena_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(+), operator(-), &
        operator(**), pi_expr, sin, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: ampere_curl, ampere_field(2), pi, poisson, roots(7)
    type(expr_t) :: x, y
    type(engine_result_t) :: result
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code
    integer :: ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    x = sym(arena, "x")
    y = sym(arena, "y")
    pi = pi_expr(arena)

    poisson = sin(pi*x)*sin(pi*y)
    ampere_field = [sin(pi*y), sin(pi*x)]
    ampere_curl = diff(ampere_field(2), x) - diff(ampere_field(1), y)
    roots(1) = poisson
    roots(2) = -(diff(diff(poisson, x), x) + diff(diff(poisson, y), y))
    roots(3:4) = ampere_field
    roots(5) = ampere_curl
    roots(6) = diff(ampere_curl, y) + ampere_field(1)
    roots(7) = -diff(ampere_curl, x) + ampere_field(2)
    do root = 1, size(roots)
        result = engine%simplify(roots(root))
        if (result%ok) roots(root) = result%value
    end do

    spec%name = str("generated_interoperability_oracles")
    spec%module_name = str("fortfem_generated_interoperability_oracles")
    spec%mode = KERNEL_SUBROUTINE
    spec%generator = str("gen_interoperability_oracles")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    spec%pure_procedure = .true.
    allocate(spec%args(2), spec%outputs(1), spec%output_shapes(1))
    allocate(spec%output_references(7))
    spec%args = [str("x"), str("y")]
    spec%outputs = [str("values")]
    spec%output_shapes = [str("(7)")]
    do root = 1, size(roots)
        spec%output_references(root) = str("values("//integer_text(root)//")")
    end do

    code = chars(emit_kernel(roots, spec))
    open (newunit=unit, &
        file=generated_path("fortfem_interoperability_oracles.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write interoperability oracle kernel"
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)

contains

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text

        character(32) :: buffer

        write (buffer, "(i0)") value
        text = trim(buffer)
    end function integer_text

end program gen_interoperability_oracles
