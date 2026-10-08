program gen_affine_map_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, operator(*), operator(+), operator(-), sym
    use fortsym_matrix, only: from_matrix, to_matrix, matrix_inverse
    use fortsym_products, only: jvp, vjp
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE, &
        KERNEL_SNIPPET, CSE_NONE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t), allocatable :: vertices(:, :), point(:), jacobian(:, :)
    type(expr_t), allocatable :: inverse(:, :), mapped(:), inputs(:), dots(:), bars(:)
    type(expr_t), allocatable :: geometry(:), geometry_bars(:), relative(:)
    type(expr_t) :: matrix, inverse_matrix
    type(str_t) :: why
    character(:), allocatable :: family
    integer :: dimension, row, column, i, unit
    logical :: ok

    call arena%init()
    engine = make_native_engine(arena)
    open(newunit=unit, file=generated_path("fortfem_affine_map_products.f90"), &
        status="replace", action="write")
    do dimension = 2, 3
        family = "triangle"
        if (dimension == 3) family = "tetra"
        allocate(vertices(dimension, dimension + 1), point(dimension))
        allocate(jacobian(dimension, dimension), mapped(dimension))
        allocate(inputs(dimension*(dimension + 2)))
        allocate(dots(size(inputs)))
        allocate(bars(dimension))
        allocate(geometry(dimension*(dimension + 1)))
        allocate(geometry_bars(size(geometry)), relative(dimension))
        i = 0
        do column = 1, dimension + 1
            do row = 1, dimension
                i = i + 1
                vertices(row, column) = sym(arena, indexed("vertices", row, column))
                inputs(i) = vertices(row, column)
                dots(i) = sym(arena, indexed("vertices_dot", row, column))
            end do
        end do
        do row = 1, dimension
            i = i + 1
            point(row) = sym(arena, indexed("point", row))
            inputs(i) = point(row)
            dots(i) = sym(arena, indexed("point_dot", row))
            bars(row) = sym(arena, indexed("reference_bar", row))
        end do
        do column = 1, dimension
            do row = 1, dimension
                jacobian(row, column) = vertices(row, column + 1) - vertices(row, 1)
            end do
        end do
        i = 0
        do column = 1, dimension
            do row = 1, dimension
                i = i + 1
                geometry(i) = jacobian(row, column)
                geometry_bars(i) = sym(arena, indexed("jacobian_bar", row, column))
            end do
        end do
        do row = 1, dimension
            i = i + 1
            geometry(i) = point(row) - vertices(row, 1)
            geometry_bars(i) = sym(arena, indexed("relative_bar", row))
        end do
        call emit(family, "_geometry", geometry, dimension)
        call emit(family, "_geometry_jvp", jvp(geometry, inputs, dots), dimension)
        call emit(family, "_geometry_vjp", &
            vjp(geometry, inputs, geometry_bars), dimension)

        deallocate(inputs, dots)
        allocate(inputs(size(geometry)), dots(size(geometry)))
        i = 0
        do column = 1, dimension
            do row = 1, dimension
                i = i + 1
                jacobian(row, column) = sym(arena, indexed("jacobian", row, column))
                inputs(i) = jacobian(row, column)
                dots(i) = sym(arena, indexed("jacobian_dot", row, column))
            end do
        end do
        do row = 1, dimension
            i = i + 1
            relative(row) = sym(arena, indexed("relative", row))
            inputs(i) = relative(row)
            dots(i) = sym(arena, indexed("relative_dot", row))
        end do
        matrix = from_matrix(arena, jacobian)
        inverse_matrix = matrix_inverse(arena, matrix, ok, why)
        if (.not. ok) error stop "symbolic affine inverse failed"
        call to_matrix(inverse_matrix, inverse, ok)
        if (.not. ok) error stop "symbolic affine inverse shape failed"
        do row = 1, dimension
            mapped(row) = num(arena, 0)
            do column = 1, dimension
                mapped(row) = mapped(row) + inverse(row, column)*relative(column)
            end do
        end do
        call emit(family, "", mapped, dimension)
        call emit(family, "_jvp", jvp(mapped, inputs, dots), dimension)
        call emit(family, "_vjp", vjp(mapped, inputs, bars), dimension)
        deallocate(vertices, point, jacobian, inverse, mapped, inputs, dots, bars)
        deallocate(geometry, geometry_bars, relative)
    end do
    close(unit)
contains
    function indexed(name, row, column) result(text)
        character(*), intent(in) :: name
        integer, intent(in) :: row
        integer, intent(in), optional :: column
        character(:), allocatable :: text
        text = name//"("//integer_text(row)
        if (present(column)) text = text//","//integer_text(column)
        text = text//")"
    end function indexed
    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text
        character(32) :: buffer
        write(buffer, "(i0)") value
        text = trim(buffer)
    end function integer_text
    subroutine emit(family, suffix, expressions, dimension)
        character(*), intent(in) :: family, suffix
        integer, intent(in) :: dimension
        type(expr_t), intent(in) :: expressions(:)
        type(expr_t) :: roots(size(expressions))
        type(engine_result_t) :: simplified
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code, vector_shape, vertices_shape, output
        integer :: k, vertex_count, destination_unit, geometry_unit
        logical :: geometry_mode
        do k = 1, size(expressions)
            simplified = engine%simplify(expressions(k))
            if (.not. simplified%ok) error stop "affine product simplification failed"
            roots(k) = simplified%value
        end do
        vector_shape = "("//integer_text(dimension)//")"
        vertices_shape = "("//integer_text(dimension)//","// &
            integer_text(dimension + 1)//")"
        spec%name = str("generated_"//family//"_affine"//suffix)
        spec%module_name = str("fortfem_generated_"//family//"_affine"//suffix)
        geometry_mode = index(suffix, "_geometry") == 1
        spec%mode = KERNEL_SUBROUTINE
        if (geometry_mode) then
            spec%mode = KERNEL_SNIPPET
            spec%cse_level = CSE_NONE
        end if
        spec%generator = str("gen_affine_map_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        select case(suffix)
        case("")
            spec%args = [str("jacobian"), str("relative")]
            spec%arg_shapes = [str("("//integer_text(dimension)//","// &
                integer_text(dimension)//")"), str(vector_shape)]
            spec%outputs = [str("reference")]
            spec%output_shapes = [str(vector_shape)]
        case("_jvp")
            spec%args = [str("jacobian"), str("relative"), &
                str("jacobian_dot"), str("relative_dot")]
            spec%arg_shapes = [str("("//integer_text(dimension)//","// &
                integer_text(dimension)//")"), str(vector_shape), &
                str("("//integer_text(dimension)//","// &
                integer_text(dimension)//")"), str(vector_shape)]
            spec%outputs = [str("reference_dot")]
            spec%output_shapes = [str(vector_shape)]
        case("_vjp")
            spec%args = [str("jacobian"), str("relative"), str("reference_bar")]
            spec%arg_shapes = [str("("//integer_text(dimension)//","// &
                integer_text(dimension)//")"), str(vector_shape), str(vector_shape)]
            spec%outputs = [str("jacobian_bar"), str("relative_bar")]
            spec%output_shapes = [str("("//integer_text(dimension)//","// &
                integer_text(dimension)//")"), str(vector_shape)]
        case("_geometry")
            spec%args = [str("vertices"), str("point")]
            spec%outputs = [str("jacobian"), str("relative")]
        case("_geometry_jvp")
            spec%args = [str("vertices_dot"), str("point_dot")]
            spec%outputs = [str("jacobian_dot"), str("relative_dot")]
        case("_geometry_vjp")
            spec%args = [str("jacobian_bar"), str("relative_bar")]
            spec%outputs = [str("vertices_bar"), str("point_bar")]
        end select
        allocate(spec%output_references(size(roots)))
        vertex_count = dimension*dimension
        if (suffix == "_geometry_vjp") vertex_count = dimension*(dimension + 1)
        do k = 1, size(roots)
            if (suffix == "_vjp" .or. suffix == "_geometry_vjp" .or. &
                suffix == "_geometry" .or. suffix == "_geometry_jvp") then
                if (k <= vertex_count) then
                    output = indexed(chars(spec%outputs(1)), &
                        mod(k - 1, dimension) + 1, (k - 1)/dimension + 1)
                else
                    output = indexed(chars(spec%outputs(2)), k - vertex_count)
                end if
            else
                output = indexed(chars(spec%outputs(1)), k)
            end if
            spec%output_references(k) = str(output)
        end do
        code = chars(emit_kernel(roots, spec))
        destination_unit = unit
        if (geometry_mode) then
            open(newunit=geometry_unit, &
                file=generated_path("fortfem_"//family//suffix//".inc"), &
                status="replace", action="write")
            destination_unit = geometry_unit
        end if
        write(destination_unit, "(a)") code(:len(code) - 1)
        if (geometry_mode) close(geometry_unit)
    end subroutine emit
end program gen_affine_map_products
