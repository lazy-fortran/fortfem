program gen_tetra_nedelec_candidates
    use fortsym_arena, only: arena_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, operator(*), operator(-), &
        operator(**), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: &
        fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: x, y, z
    integer :: order

    call arena%init()
    engine = make_native_engine(arena)
    x = sym(arena, "x")
    y = sym(arena, "y")
    z = sym(arena, "z")
    do order = 1, 4
        call generate_order(order)
    end do

contains

    subroutine generate_order(order)
        integer, intent(in) :: order

        type(expr_t), allocatable :: curls(:, :), roots(:), values(:, :)
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code, filename, module_name
        character(:), allocatable :: procedure_name
        integer :: dof_count, ios, root, unit

        dof_count = order * (order + 2) * (order + 3) / 2
        allocate(values(3, dof_count), curls(3, dof_count))
        call build_candidates(order, values)
        call differentiate_candidates(values, curls)
        allocate(roots(6 * dof_count))
        call flatten_outputs(values, curls, roots)
        call simplify_all(roots)

        filename = "fortfem_tetra_nedelec_candidates_order_"// &
            integer_text(order)//".f90"
        module_name = &
            "fortfem_generated_tetra_nedelec_candidates_order_"// &
            integer_text(order)
        procedure_name = "evaluate_candidates_order_"//integer_text(order)
        call configure_spec( &
            order, dof_count, module_name, procedure_name, spec)
        code = chars(emit_kernel(roots, spec))
        open( &
            newunit=unit, file=generated_path(filename), status="replace", &
            action="write", iostat=ios)
        if (ios /= 0) error stop "cannot write generated candidate kernel"
        write(unit, "(a)") code(:len(code) - 1)
        close(unit)
    end subroutine generate_order

    subroutine build_candidates(order, values)
        integer, intent(in) :: order
        type(expr_t), intent(out) :: values(:, :)

        type(expr_t) :: monomial, zero
        integer :: candidate, component, total_degree
        integer :: x_degree, y_degree, z_degree

        zero = num(arena, 0)
        values = zero
        candidate = 0
        do component = 1, 3
            do total_degree = 0, order - 1
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

        total_degree = order - 1
        do component = 4, 5
            do x_degree = 0, total_degree
                do y_degree = 0, total_degree - x_degree
                    z_degree = total_degree - x_degree - y_degree
                    candidate = candidate + 1
                    monomial = &
                        monomial_value(x_degree, y_degree, z_degree)
                    if (component == 4) then
                        values(1, candidate) = -y * monomial
                        values(2, candidate) = x * monomial
                    else
                        values(1, candidate) = -z * monomial
                        values(3, candidate) = x * monomial
                    end if
                end do
            end do
        end do
        do y_degree = 0, total_degree
            z_degree = total_degree - y_degree
            candidate = candidate + 1
            monomial = monomial_value(0, y_degree, z_degree)
            values(2, candidate) = -z * monomial
            values(3, candidate) = y * monomial
        end do
        if (candidate /= size(values, 2)) then
            error stop "candidate dimension mismatch"
        end if
    end subroutine build_candidates

    subroutine differentiate_candidates(values, curls)
        type(expr_t), intent(in) :: values(:, :)
        type(expr_t), intent(out) :: curls(:, :)

        integer :: candidate

        do candidate = 1, size(values, 2)
            curls(1, candidate) = &
                diff(values(3, candidate), y) - &
                diff(values(2, candidate), z)
            curls(2, candidate) = &
                diff(values(1, candidate), z) - &
                diff(values(3, candidate), x)
            curls(3, candidate) = &
                diff(values(2, candidate), x) - &
                diff(values(1, candidate), y)
        end do
    end subroutine differentiate_candidates

    subroutine flatten_outputs(values, curls, roots)
        type(expr_t), intent(in) :: values(:, :), curls(:, :)
        type(expr_t), intent(out) :: roots(:)

        integer :: candidate, component, root

        root = 0
        do candidate = 1, size(values, 2)
            do component = 1, 3
                root = root + 1
                roots(root) = values(component, candidate)
            end do
        end do
        do candidate = 1, size(curls, 2)
            do component = 1, 3
                root = root + 1
                roots(root) = curls(component, candidate)
            end do
        end do
    end subroutine flatten_outputs

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
            order, dof_count, module_name, procedure_name, spec)
        integer, intent(in) :: order, dof_count
        character(*), intent(in) :: module_name, procedure_name
        type(kernel_spec_t), intent(out) :: spec

        integer :: candidate, component, root

        spec%name = str(procedure_name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str("t")
        spec%generator = str("gen_tetra_nedelec_candidates")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        allocate(spec%args(3), spec%outputs(2))
        allocate(spec%output_shapes(2))
        allocate(spec%output_references(6 * dof_count))
        spec%args = [str("x"), str("y"), str("z")]
        spec%outputs = [str("values"), str("curls")]
        spec%output_shapes = [ &
            str("(3,"//integer_text(dof_count)//")"), &
            str("(3,"//integer_text(dof_count)//")")]

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
            do component = 1, 3
                root = root + 1
                spec%output_references(root) = str( &
                    "curls("//integer_text(component)//","// &
                    integer_text(candidate)//")")
            end do
        end do
    end subroutine configure_spec

    function monomial_value(x_degree, y_degree, z_degree) result(value)
        integer, intent(in) :: x_degree, y_degree, z_degree
        type(expr_t) :: value

        value = x**x_degree * y**y_degree * z**z_degree
    end function monomial_value

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text

        character(32) :: buffer

        write(buffer, "(i0)") value
        text = trim(buffer)
    end function integer_text

end program gen_tetra_nedelec_candidates
