program gen_tetra_piola_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, sym, operator(+), operator(*), operator(/)
    use fortsym_products, only: jvp, vjp
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SNIPPET, CSE_NONE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: j(3, 3), inverse(3, 3), determinant, values(3), curls(3), divergence
    type(expr_t), allocatable :: inputs(:), dots(:), bars(:), mapped(:)
    type(str_t), allocatable :: output_references(:)
    character(:), allocatable :: family, scalar_name
    integer :: family_id, row, column, index

    call arena%init()
    engine = make_native_engine(arena)
    do family_id = 1, 2
        family = "tetra_covariant"
        scalar_name = "curls"
        if (family_id == 2) then
            family = "tetra_contravariant"
            scalar_name = "divergences"
        end if
        if (family_id == 1) then
            allocate(inputs(25), dots(25), output_references(25))
            allocate(mapped(6), bars(6))
        else
            allocate(inputs(14), dots(14), output_references(14))
            allocate(mapped(4), bars(4))
        end if
        index = 0
        if (family_id == 1) then
            do column = 1, 3
                do row = 1, 3
                    inverse(row, column) = &
                        register(indexed("inverse_jacobian", row, column), &
                        indexed("inverse_jacobian_dot", row, column), &
                        indexed("inverse_jacobian_bar_local", row, column))
                end do
            end do
        end if
        do column = 1, 3
            do row = 1, 3
                j(row, column) = register(indexed("jacobian", row, column), &
                    indexed("jacobian_dot", row, column), &
                    indexed("direct_jacobian_bar_local", row, column))
            end do
        end do
        determinant = register("determinant", "determinant_dot", &
            "determinant_bar_local")
        do row = 1, 3
            values(row) = register(indexed("reference_values", row, basis=.true.), &
                indexed("reference_values_dot", row, basis=.true.), &
                indexed("reference_values_bar", row, basis=.true.))
            bars(row) = sym(arena, indexed("physical_values_bar", row, basis=.true.))
        end do
        if (family_id == 1) then
            do row = 1, 3
                curls(row) = register(indexed("reference_curls", row, basis=.true.), &
                    indexed("reference_curls_dot", row, basis=.true.), &
                    indexed("reference_curls_bar", row, basis=.true.))
                bars(3 + row) = sym(arena, &
                    indexed("physical_curls_bar", row, basis=.true.))
            end do
        else
            divergence = register("reference_divergences(basis)", &
                "reference_divergences_dot(basis)", "reference_divergences_bar(basis)")
            bars(4) = sym(arena, "physical_divergences_bar(basis)")
        end if
        if (index /= size(inputs)) error stop "Piola input registration mismatch"
        do row = 1, 3
            mapped(row) = num(arena, 0)
            if (family_id == 1) then
                mapped(3 + row) = num(arena, 0)
                do column = 1, 3
                    mapped(row) = mapped(row) + inverse(column, row)*values(column)
                    mapped(3 + row) = mapped(3 + row) + j(row, column)*curls(column)
                end do
                mapped(3 + row) = mapped(3 + row)/determinant
            else
                do column = 1, 3
                    mapped(row) = mapped(row) + j(row, column)*values(column)
                end do
                mapped(row) = mapped(row)/determinant
            end if
        end do
        if (family_id == 2) mapped(4) = divergence/determinant
        call emit(family, "", mapped, family_id)
        call emit(family, "_jvp", jvp(mapped, inputs, dots), family_id)
        call emit(family, "_vjp", vjp(mapped, inputs, bars), family_id)
        deallocate(inputs, dots, bars, mapped, output_references)
    end do
contains
    function register(name, dot_name, bar_name) result(value)
        character(*), intent(in) :: name, dot_name, bar_name
        type(expr_t) :: value
        index = index + 1
        value = sym(arena, name)
        inputs(index) = value
        dots(index) = sym(arena, dot_name)
        output_references(index) = str(bar_name)
    end function register
    function indexed(name, row, column, basis) result(text)
        character(*), intent(in) :: name
        integer, intent(in) :: row
        integer, intent(in), optional :: column
        logical, intent(in), optional :: basis
        character(:), allocatable :: text
        character(32) :: buffer
        write(buffer, "(i0)") row
        text = name//"("//trim(buffer)
        if (present(column)) then
            write(buffer, "(i0)") column
            text = text//","//trim(buffer)
        end if
        if (present(basis)) then
            if (basis) text = text//",basis"
        end if
        text = text//")"
    end function indexed
    subroutine emit(family, suffix, expressions, family_id)
        character(*), intent(in) :: family, suffix
        integer, intent(in) :: family_id
        type(expr_t), intent(in) :: expressions(:)
        type(expr_t) :: roots(size(expressions))
        type(engine_result_t) :: simplified
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code, output_suffix
        integer :: k, unit

        do k = 1, size(roots)
            simplified = engine%simplify(expressions(k))
            if (.not. simplified%ok) error stop "tetra Piola simplification failed"
            roots(k) = simplified%value
        end do
        spec%mode = KERNEL_SNIPPET
        spec%cse_level = CSE_NONE
        spec%name = str("generated_"//family//suffix)
        spec%generator = str("gen_tetra_piola_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        if (family_id == 1) then
            spec%args = [str("inverse_jacobian"), str("jacobian"), str("determinant"), &
                str("reference_values"), str("reference_curls"), &
                str("inverse_jacobian_dot"), str("jacobian_dot"), str("determinant_dot"), &
                str("reference_values_dot"), str("reference_curls_dot"), &
                str("physical_values_bar"), str("physical_curls_bar")]
        else
            spec%args = [str("jacobian"), str("determinant"), str("reference_values"), &
                str("reference_divergences"), str("jacobian_dot"), &
                str("determinant_dot"), str("reference_values_dot"), &
                str("reference_divergences_dot"), str("physical_values_bar"), &
                str("physical_divergences_bar")]
        end if
        if (suffix == "_vjp") then
            spec%output_references = output_references
            if (family_id == 1) then
                spec%outputs = [str("inverse_jacobian_bar_local"), &
                    str("direct_jacobian_bar_local"), str("determinant_bar_local"), &
                    str("reference_values_bar"), str("reference_curls_bar")]
            else
                spec%outputs = [str("direct_jacobian_bar_local"), &
                    str("determinant_bar_local"), str("reference_values_bar"), &
                    str("reference_divergences_bar")]
            end if
        else
            output_suffix = ""
            if (suffix == "_jvp") output_suffix = "_dot"
            spec%outputs = [str("physical_values"//output_suffix), &
                str("physical_"//scalar_name//output_suffix)]
            allocate(spec%output_references(size(roots)))
            do k = 1, 3
                spec%output_references(k) = &
                    str(indexed("physical_values"//output_suffix, k, basis=.true.))
            end do
            if (family_id == 1) then
                do k = 1, 3
                    spec%output_references(3 + k) = &
                        str(indexed("physical_curls"//output_suffix, k, basis=.true.))
                end do
            else
                spec%output_references(4) = &
                    str("physical_divergences"//output_suffix//"(basis)")
            end if
        end if
        code = chars(emit_kernel(roots, spec))
        open(newunit=unit, file=generated_path("fortfem_"//family//suffix//".inc"), &
            status="replace", action="write")
        write(unit, "(a)") code(:len(code) - 1)
        close(unit)
    end subroutine emit
end program gen_tetra_piola_products
