program gen_tetra_rt_candidates
    use fortsym_arena, only: arena_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, operator(*), operator(+), &
        operator(**), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: &
        fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: x, y, z
    integer :: degree

    call arena%init()
    engine = make_native_engine(arena)
    x = sym(arena, "x")
    y = sym(arena, "y")
    z = sym(arena, "z")
    do degree = 0, 4
        call generate_degree(degree)
    end do

contains

    subroutine generate_degree(degree)
        integer, intent(in) :: degree

        type(expr_t), allocatable :: divergences(:), roots(:), values(:, :)
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code, filename, module_name
        character(:), allocatable :: procedure_name
        integer :: candidate, component, dof_count, ios, root, unit

        dof_count = 3*(degree + 1)*(degree + 2)*(degree + 3)/6 + &
            (degree + 1)*(degree + 2)/2
        allocate(values(3, dof_count), divergences(dof_count))
        call build_candidates(degree, values)
        do candidate = 1, dof_count
            divergences(candidate) = diff(values(1, candidate), x) + &
                diff(values(2, candidate), y) + &
                diff(values(3, candidate), z)
        end do
        allocate(roots(4*dof_count))
        root = 0
        do candidate = 1, dof_count
            do component = 1, 3
                root = root + 1
                roots(root) = values(component, candidate)
            end do
        end do
        do candidate = 1, dof_count
            root = root + 1
            roots(root) = divergences(candidate)
        end do
        call simplify_all(roots)

        filename = "fortfem_tetra_rt_candidates_degree_"// &
            integer_text(degree)//".f90"
        module_name = "fortfem_generated_tetra_rt_candidates_degree_"// &
            integer_text(degree)
        procedure_name = "evaluate_rt_candidates_degree_"// &
            integer_text(degree)
        call configure_spec( &
            dof_count, module_name, procedure_name, spec)
        code = chars(emit_kernel(roots, spec))
        open( &
            newunit=unit, file=generated_path(filename), status="replace", &
            action="write", iostat=ios)
        if (ios /= 0) error stop "cannot write generated RT candidate kernel"
        write(unit, "(a)") code(:len(code) - 1)
        close(unit)
    end subroutine generate_degree

    subroutine build_candidates(degree, values)
        integer, intent(in) :: degree
        type(expr_t), intent(out) :: values(:, :)

        type(expr_t) :: monomial, zero
        integer :: candidate, component, total_degree
        integer :: x_degree, y_degree, z_degree

        zero = num(arena, 0)
        values = zero
        candidate = 0
        do component = 1, 3
            do total_degree = 0, degree
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        candidate = candidate + 1
                        values(component, candidate) = &
                            monomial_value(x_degree, y_degree, z_degree)
                    end do
                end do
            end do
        end do
        total_degree = degree
        do x_degree = 0, total_degree
            do y_degree = 0, total_degree - x_degree
                z_degree = total_degree - x_degree - y_degree
                candidate = candidate + 1
                monomial = &
                    monomial_value(x_degree, y_degree, z_degree)
                values(1, candidate) = x*monomial
                values(2, candidate) = y*monomial
                values(3, candidate) = z*monomial
            end do
        end do
        if (candidate /= size(values, 2)) then
            error stop "RT candidate dimension mismatch"
        end if
    end subroutine build_candidates

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        type(engine_result_t) :: result
        integer :: expression

        do expression = 1, size(expressions)
            result = engine%simplify(expressions(expression))
            if (result%ok) expressions(expression) = result%value
        end do
    end subroutine simplify_all

    subroutine configure_spec( &
            dof_count, module_name, procedure_name, spec)
        integer, intent(in) :: dof_count
        character(*), intent(in) :: module_name, procedure_name
        type(kernel_spec_t), intent(out) :: spec

        integer :: candidate, component, root

        spec%name = str(procedure_name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str("t")
        spec%generator = str("gen_tetra_rt_candidates")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        allocate(spec%args(3), spec%outputs(2), spec%output_shapes(2))
        allocate(spec%output_references(4*dof_count))
        spec%args = [str("x"), str("y"), str("z")]
        spec%outputs = [str("values"), str("divergences")]
        spec%output_shapes = [ &
            str("(3,"//integer_text(dof_count)//")"), &
            str("("//integer_text(dof_count)//")")]
        root = 0
        do candidate = 1, dof_count
            do component = 1, 3
                root = root + 1
                spec%output_references(root) = str( &
                    "values("//integer_text(component)//","// &
                    integer_text(candidate)//")")
            end do
        end do
        do candidate = 1, dof_count
            root = root + 1
            spec%output_references(root) = str( &
                "divergences("//integer_text(candidate)//")")
        end do
    end subroutine configure_spec

    function monomial_value(x_degree, y_degree, z_degree) result(value)
        integer, intent(in) :: x_degree, y_degree, z_degree
        type(expr_t) :: value

        value = x**x_degree*y**y_degree*z**z_degree
    end function monomial_value

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text

        character(32) :: buffer

        write(buffer, "(i0)") value
        text = trim(buffer)
    end function integer_text
end program gen_tetra_rt_candidates
