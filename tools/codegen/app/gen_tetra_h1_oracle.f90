program gen_tetra_h1_oracle
    use fortsym_arena, only: arena_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(+), operator(-), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: bubble, roots(5), x, y, z
    type(engine_result_t) :: result
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code
    integer :: ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    x = sym(arena, "x")
    y = sym(arena, "y")
    z = sym(arena, "z")

    bubble = x*y*z*(1 - x - y - z)
    roots(1) = bubble
    roots(2) = diff(bubble, x)
    roots(3) = diff(bubble, y)
    roots(4) = diff(bubble, z)
    roots(5) = -( &
        diff(diff(bubble, x), x) + diff(diff(bubble, y), y) + &
        diff(diff(bubble, z), z))
    do root = 1, size(roots)
        result = engine%simplify(roots(root))
        if (result%ok) roots(root) = result%value
    end do

    spec%name = str("generated_tetra_h1_oracle")
    spec%module_name = str("fortfem_generated_tetra_h1_oracle")
    spec%mode = KERNEL_SUBROUTINE
    spec%generator = str("gen_tetra_h1_oracle")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    spec%pure_procedure = .true.
    allocate(spec%args(3), spec%outputs(1), spec%output_shapes(1))
    allocate(spec%output_references(5))
    spec%args = [str("x"), str("y"), str("z")]
    spec%outputs = [str("values")]
    spec%output_shapes = [str("(5)")]
    do root = 1, size(roots)
        spec%output_references(root) = str("values("//integer_text(root)//")")
    end do

    code = chars(emit_kernel(roots, spec))
    open (newunit=unit, &
        file=generated_path("fortfem_tetra_h1_oracle.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write tetrahedral H1 oracle kernel"
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

end program gen_tetra_h1_oracle
