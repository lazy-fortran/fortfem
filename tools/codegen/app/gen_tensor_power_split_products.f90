program gen_tensor_power_split_products
    !! Emit fixed 3-D tensor-power split products and exact derivatives.
    !!
    !! The runtime wrapper owns shape/finiteness checks; FortSym owns the
    !! symmetric, skew, and total quadratic contractions and their JVP/VJP.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        real_expr, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    integer, parameter :: dimension = 3
    integer, parameter :: tensor_size = dimension*dimension
    integer, parameter :: variable_count = tensor_size + dimension
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: variables(variable_count), variables_dot(variable_count)
    type(expr_t) :: powers(3), powers_dot(3), powers_vjp(variable_count)
    type(expr_t) :: power_bar(3), half
    type(kernel_spec_t) :: spec
    character(:), allocatable :: value_code, jvp_code, vjp_code
    character(len=32) :: name
    integer :: row, column, index, ios, unit

    call arena%init()
    engine = make_native_engine(arena)
    half = real_expr(arena, 0.5d0)
    do row = 1, dimension
        do column = 1, dimension
            index = (row - 1)*dimension + column
            write (name, '("tensor_",i1,i1)') row, column
            variables(index) = sym(arena, trim(name))
            write (name, '("tensor_",i1,i1,"_dot")') row, column
            variables_dot(index) = sym(arena, trim(name))
        end do
    end do
    do row = 1, dimension
        index = tensor_size + row
        write (name, '("vector_",i1)') row
        variables(index) = sym(arena, trim(name))
        write (name, '("vector_",i1,"_dot")') row
        variables_dot(index) = sym(arena, trim(name))
    end do

    powers = [real_expr(arena, 0.0d0), real_expr(arena, 0.0d0), &
        real_expr(arena, 0.0d0)]
    do row = 1, dimension
        do column = 1, dimension
            index = (row - 1)*dimension + column
            powers(1) = powers(1) + half*(variables(index) + &
                variables((column - 1)*dimension + row))* &
                variables(tensor_size + row)*variables(tensor_size + column)
            powers(2) = powers(2) + half*(variables(index) - &
                variables((column - 1)*dimension + row))* &
                variables(tensor_size + row)*variables(tensor_size + column)
            powers(3) = powers(3) + variables(index)* &
                variables(tensor_size + row)*variables(tensor_size + column)
        end do
    end do
    power_bar = [sym(arena, "symmetric_power_bar"), &
        sym(arena, "skew_power_bar"), sym(arena, "total_power_bar")]
    powers_dot = jvp(powers, variables, variables_dot)
    powers_vjp = vjp(powers, variables, power_bar)
    call simplify_all(powers)
    call simplify_all(powers_dot)
    call simplify_all(powers_vjp)

    call initialize_spec(spec, "generated_tensor_power_split", &
        "fortfem_generated_tensor_power_split", variable_count, 3)
    spec%args = primal_arguments()
    spec%outputs = [str("symmetric_power"), str("skew_power"), &
        str("total_power")]
    value_code = chars(emit_kernel(powers, spec))

    call initialize_spec(spec, "generated_tensor_power_split_jvp", &
        "fortfem_generated_tensor_power_split_jvp", 2*variable_count, 3)
    spec%args = jvp_arguments()
    spec%outputs = [str("symmetric_power_dot"), str("skew_power_dot"), &
        str("total_power_dot")]
    jvp_code = chars(emit_kernel(powers_dot, spec))

    call initialize_spec(spec, "generated_tensor_power_split_vjp", &
        "fortfem_generated_tensor_power_split_vjp", variable_count + 3, &
        variable_count)
    spec%args = vjp_arguments()
    spec%outputs = vjp_outputs()
    vjp_code = chars(emit_kernel(powers_vjp, spec))

    open (newunit=unit, file=generated_path( &
        "fortfem_tensor_power_split_products.f90"), status="replace", &
        action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write tensor power split products"
    write (unit, "(a)") value_code(:len(value_code) - 1)
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    subroutine initialize_spec(kernel_spec, name_value, module_name, &
            argument_count, output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name_value, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name_value)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_tensor_power_split_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str( &
            "cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate(kernel_spec%args(argument_count), kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        integer :: expression_index

        do expression_index = 1, size(expressions)
            result = engine%simplify(expressions(expression_index))
            if (result%ok) expressions(expression_index) = result%value
        end do
    end subroutine simplify_all

    function primal_arguments() result(arguments)
        type(str_t) :: arguments(variable_count)
        integer :: argument_index, tensor_row, tensor_column

        do tensor_row = 1, dimension
            do tensor_column = 1, dimension
                argument_index = (tensor_row - 1)*dimension + tensor_column
                write (name, '("tensor_",i1,i1)') tensor_row, tensor_column
                arguments(argument_index) = str(trim(name))
            end do
        end do
        do tensor_row = 1, dimension
            write (name, '("vector_",i1)') tensor_row
            arguments(tensor_size + tensor_row) = str(trim(name))
        end do
    end function primal_arguments

    function jvp_arguments() result(arguments)
        type(str_t) :: arguments(2*variable_count)
        integer :: argument_index, tensor_row, tensor_column

        arguments(1:variable_count) = primal_arguments()
        do tensor_row = 1, dimension
            do tensor_column = 1, dimension
                argument_index = variable_count + &
                    (tensor_row - 1)*dimension + tensor_column
                write (name, '("tensor_",i1,i1,"_dot")') tensor_row, tensor_column
                arguments(argument_index) = str(trim(name))
            end do
        end do
        do tensor_row = 1, dimension
            write (name, '("vector_",i1,"_dot")') tensor_row
            arguments(variable_count + tensor_size + tensor_row) = str(trim(name))
        end do
    end function jvp_arguments

    function vjp_arguments() result(arguments)
        type(str_t) :: arguments(variable_count + 3)

        arguments(1:variable_count) = primal_arguments()
        arguments(variable_count + 1) = str("symmetric_power_bar")
        arguments(variable_count + 2) = str("skew_power_bar")
        arguments(variable_count + 3) = str("total_power_bar")
    end function vjp_arguments

    function vjp_outputs() result(outputs)
        type(str_t) :: outputs(variable_count)
        integer :: output_index, tensor_row, tensor_column

        do tensor_row = 1, dimension
            do tensor_column = 1, dimension
                output_index = (tensor_row - 1)*dimension + tensor_column
                write (name, '("tensor_",i1,i1,"_bar")') tensor_row, tensor_column
                outputs(output_index) = str(trim(name))
            end do
        end do
        do tensor_row = 1, dimension
            write (name, '("vector_",i1,"_bar")') tensor_row
            outputs(tensor_size + tensor_row) = str(trim(name))
        end do
    end function vjp_outputs

end program gen_tensor_power_split_products
