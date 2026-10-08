program gen_tetra_whitney_products
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
    type(engine_result_t) :: result
    type(expr_t) :: coordinates(3), lambda(4), field(3), roots(36)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code, reference
    integer, parameter :: endpoints(2, 6) = reshape( &
        [1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4], [2, 6])
    integer :: edge, component, first, second, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    coordinates = [sym(arena, "x"), sym(arena, "y"), sym(arena, "z")]
    lambda(1) = 1 - coordinates(1) - coordinates(2) - coordinates(3)
    lambda(2:4) = coordinates
    do edge = 1, 6
        first = endpoints(1, edge)
        second = endpoints(2, edge)
        do component = 1, 3
            field(component) = lambda(first)*diff(lambda(second), &
                coordinates(component)) - lambda(second)* &
                diff(lambda(first), coordinates(component))
            roots(3*(edge - 1) + component) = field(component)
        end do
        roots(18 + 3*(edge - 1) + 1) = &
            diff(field(3), coordinates(2)) - diff(field(2), coordinates(3))
        roots(18 + 3*(edge - 1) + 2) = &
            diff(field(1), coordinates(3)) - diff(field(3), coordinates(1))
        roots(18 + 3*(edge - 1) + 3) = &
            diff(field(2), coordinates(1)) - diff(field(1), coordinates(2))
    end do
    do root = 1, size(roots)
        result = engine%simplify(roots(root))
        if (.not. result%ok) error stop "Whitney simplification failed"
        roots(root) = result%value
    end do

    spec%name = str("generated_tetra_whitney")
    spec%module_name = str("fortfem_generated_tetra_whitney")
    spec%mode = KERNEL_SUBROUTINE
    spec%generator = str("gen_tetra_whitney_products")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    spec%pure_procedure = .true.
    allocate(spec%args(3), spec%outputs(2), spec%output_shapes(2))
    allocate(spec%output_references(36))
    spec%args = [str("x"), str("y"), str("z")]
    spec%outputs = [str("values"), str("curls")]
    spec%output_shapes = [str("(3,6)"), str("(3,6)")]
    do edge = 1, 6
        do component = 1, 3
            reference = "("//integer_text(component)//","//integer_text(edge)//")"
            spec%output_references(3*(edge - 1) + component) = &
                str("values"//reference)
            spec%output_references(18 + 3*(edge - 1) + component) = &
                str("curls"//reference)
        end do
    end do
    code = chars(emit_kernel(roots, spec))
    open (newunit=unit, &
        file=generated_path("fortfem_tetra_whitney_products.f90"), &
        status="replace", action="write")
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
end program gen_tetra_whitney_products
