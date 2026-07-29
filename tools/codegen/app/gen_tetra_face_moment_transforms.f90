program gen_tetra_face_moment_transforms
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_arena, only: arena_t, NK_INT, NK_RAT
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, operator(*), operator(+), rat
    use fortsym_linalg, only: &
        exact_linear_system_result_t, solve_exact_linear_system
    use fortfem_codegen_provenance, only: &
        fortsym_revision, generated_path
    implicit none

    integer, parameter :: permutations(3, 6) = reshape([ &
        1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    integer :: order, output_unit

    call arena%init()
    engine = make_native_engine(arena)
    call write_module_header()
    do order = 2, 4
        call generate_order(order)
    end do
    call write_module_footer()

contains

    subroutine generate_order(order)
        integer, intent(in) :: order

        type(exact_linear_system_result_t) :: solution
        type(expr_t), allocatable :: canonical(:, :), local(:, :)
        type(expr_t), allocatable :: left(:, :), right(:, :)
        integer :: dof_count, permutation

        dof_count = order * (order - 1)
        allocate( &
            canonical(dof_count, dof_count), local(dof_count, dof_count), &
            left(dof_count, dof_count), right(dof_count, dof_count))
        call build_moment_matrix( &
            order, permutations(:, 1), canonical)
        call write_order_header( &
            "transform_order_"//integer_text(order), dof_count)
        do permutation = 1, size(permutations, 2)
            call build_moment_matrix( &
                order, permutations(:, permutation), local)
            left = transpose(local)
            right = transpose(canonical)
            solution = solve_exact_linear_system(engine, left, right)
            if (.not. solution%ok) then
                error stop "exact face moment transform solve failed"
            end if
            call write_transform_case( &
                permutation, transpose(solution%values))
        end do
        call write_order_footer( &
            "transform_order_"//integer_text(order))

        call write_order_header( &
            "basis_to_local_order_"//integer_text(order), dof_count)
        do permutation = 1, size(permutations, 2)
            call build_moment_matrix( &
                order, permutations(:, permutation), local)
            left = transpose(canonical)
            right = transpose(local)
            solution = solve_exact_linear_system(engine, left, right)
            if (.not. solution%ok) then
                error stop "exact face basis transform solve failed"
            end if
            call write_transform_case( &
                permutation, transpose(solution%values))
        end do
        call write_order_footer( &
            "basis_to_local_order_"//integer_text(order))
    end subroutine generate_order

    subroutine build_moment_matrix(order, permutation, matrix)
        integer, intent(in) :: order, permutation(3)
        type(expr_t), intent(out) :: matrix(:, :)

        integer :: affine(2, 2), offset(2)
        integer :: basis_component, basis_x, basis_y, column
        integer :: moment_component, moment_x, moment_y, row
        integer :: basis_degree, moment_degree

        call permutation_map(permutation, offset, affine)
        matrix = num(arena, 0)
        column = 0
        do basis_component = 1, 2
            do basis_degree = 0, order - 2
                do basis_x = 0, basis_degree
                    basis_y = basis_degree - basis_x
                    column = column + 1
                    row = 0
                    do moment_component = 1, 2
                        do moment_degree = 0, order - 2
                            do moment_x = 0, moment_degree
                                moment_y = moment_degree - moment_x
                                row = row + 1
                                matrix(row, column) = num( &
                                    arena, affine( &
                                    basis_component, moment_component)) * &
                                    transformed_monomial_integral( &
                                    basis_x, basis_y, moment_x, moment_y, &
                                    offset, affine)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine build_moment_matrix

    function transformed_monomial_integral( &
            basis_x, basis_y, moment_x, moment_y, offset, affine) &
            result(value)
        integer, intent(in) :: basis_x, basis_y, moment_x, moment_y
        integer, intent(in) :: offset(2), affine(2, 2)
        type(expr_t) :: value

        integer(int64) :: coefficients(0:2, 0:2)
        integer :: x_degree, y_degree

        coefficients = 0_int64
        coefficients(0, 0) = 1_int64
        call multiply_affine_power( &
            coefficients, basis_x, offset(1), affine(1, 1), affine(1, 2))
        call multiply_affine_power( &
            coefficients, basis_y, offset(2), affine(2, 1), affine(2, 2))
        value = num(arena, 0)
        do x_degree = 0, 2
            do y_degree = 0, 2 - x_degree
                if (coefficients(x_degree, y_degree) == 0_int64) cycle
                value = value + num( &
                    arena, coefficients(x_degree, y_degree)) * &
                    triangle_monomial_integral( &
                    x_degree + moment_x, y_degree + moment_y)
            end do
        end do
    end function transformed_monomial_integral

    subroutine multiply_affine_power( &
            coefficients, power, constant, x_coefficient, y_coefficient)
        integer(int64), intent(inout) :: coefficients(0:2, 0:2)
        integer, intent(in) :: power, constant, x_coefficient, y_coefficient

        integer(int64) :: next(0:2, 0:2)
        integer :: factor, x_degree, y_degree

        do factor = 1, power
            next = 0_int64
            do x_degree = 0, 2
                do y_degree = 0, 2 - x_degree
                    next(x_degree, y_degree) = &
                        next(x_degree, y_degree) + &
                        int(constant, int64) * &
                        coefficients(x_degree, y_degree)
                    if (x_degree < 2) then
                        next(x_degree + 1, y_degree) = &
                            next(x_degree + 1, y_degree) + &
                            int(x_coefficient, int64) * &
                            coefficients(x_degree, y_degree)
                    end if
                    if (y_degree < 2) then
                        next(x_degree, y_degree + 1) = &
                            next(x_degree, y_degree + 1) + &
                            int(y_coefficient, int64) * &
                            coefficients(x_degree, y_degree)
                    end if
                end do
            end do
            coefficients = next
        end do
    end subroutine multiply_affine_power

    function triangle_monomial_integral(x_degree, y_degree) result(value)
        integer, intent(in) :: x_degree, y_degree
        type(expr_t) :: value

        value = rat( &
            arena, factorial(x_degree) * factorial(y_degree), &
            factorial(x_degree + y_degree + 2))
    end function triangle_monomial_integral

    pure function factorial(argument) result(value)
        integer, intent(in) :: argument
        integer(int64) :: value

        integer :: factor

        value = 1_int64
        do factor = 2, argument
            value = value * int(factor, int64)
        end do
    end function factorial

    pure subroutine permutation_map(permutation, offset, affine)
        integer, intent(in) :: permutation(3)
        integer, intent(out) :: offset(2), affine(2, 2)

        integer, parameter :: vertices(2, 3) = reshape([ &
            0, 0, 1, 0, 0, 1], [2, 3])

        offset = vertices(:, permutation(1))
        affine(:, 1) = vertices(:, permutation(2)) - offset
        affine(:, 2) = vertices(:, permutation(3)) - offset
    end subroutine permutation_map

    subroutine write_module_header()
        integer :: ios

        open( &
            newunit=output_unit, file=generated_path( &
            "fortfem_tetra_face_moment_transforms.f90"), &
            status="replace", action="write", iostat=ios)
        if (ios /= 0) error stop "cannot write generated face transforms"
        write(output_unit, "(a)") "! Generated by fortsym. Do not edit."
        write(output_unit, "(a)") &
            "! Generator: gen_tetra_face_moment_transforms"
        write(output_unit, "(a)") &
            "! Generator revision: "//fortsym_revision()
        write(output_unit, "(a)") &
            "! Regenerate with: cd tools/codegen && ./generate.sh"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") &
            "module fortfem_generated_tetra_face_moment_transforms"
        write(output_unit, "(a)") &
            "    use, intrinsic :: iso_fortran_env, only: real64"
        write(output_unit, "(a)") "    implicit none"
        write(output_unit, "(a)") "    private"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") &
            "    public :: transform_tetra_face_moments"
        write(output_unit, "(a)") &
            "    public :: map_tetra_face_basis_to_local"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") "contains"
        write(output_unit, "(a)") ""
        call write_dispatch_procedure()
        call write_basis_dispatch_procedure()
    end subroutine write_module_header

    subroutine write_dispatch_procedure()
        write(output_unit, "(a)") &
            "    pure subroutine transform_tetra_face_moments( &"
        write(output_unit, "(a)") &
            "            order, permutation, local, canonical, status)"
        write(output_unit, "(a)") &
            "        integer, intent(in) :: order, permutation(3)"
        write(output_unit, "(a)") &
            "        real(real64), intent(in) :: local(:)"
        write(output_unit, "(a)") &
            "        real(real64), intent(out) :: canonical(:)"
        write(output_unit, "(a)") "        integer, intent(out) :: status"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") "        integer :: permutation_index"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") "        canonical = 0.0_real64"
        write(output_unit, "(a)") "        status = 1"
        write(output_unit, "(a)") &
            "        permutation_index = face_permutation_index(permutation)"
        write(output_unit, "(a)") &
            "        if (permutation_index == 0) return"
        write(output_unit, "(a)") "        select case (order)"
        write(output_unit, "(a)") "        case (2)"
        write(output_unit, "(a)") &
            "            call transform_order_2( &"
        write(output_unit, "(a)") &
            "                permutation_index, local, canonical, status)"
        write(output_unit, "(a)") "        case (3)"
        write(output_unit, "(a)") &
            "            call transform_order_3( &"
        write(output_unit, "(a)") &
            "                permutation_index, local, canonical, status)"
        write(output_unit, "(a)") "        case (4)"
        write(output_unit, "(a)") &
            "            call transform_order_4( &"
        write(output_unit, "(a)") &
            "                permutation_index, local, canonical, status)"
        write(output_unit, "(a)") "        end select"
        write(output_unit, "(a)") &
            "    end subroutine transform_tetra_face_moments"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") &
            "    pure function face_permutation_index(permutation) result(index_)"
        write(output_unit, "(a)") &
            "        integer, intent(in) :: permutation(3)"
        write(output_unit, "(a)") "        integer :: index_"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") "        index_ = 0"
        write(output_unit, "(a)") &
            "        if (any(permutation < 1) .or. any(permutation > 3)) return"
        write(output_unit, "(a)") &
            "        if (permutation(1) == permutation(2) .or. &"
        write(output_unit, "(a)") &
            "            permutation(1) == permutation(3) .or. &"
        write(output_unit, "(a)") &
            "            permutation(2) == permutation(3)) return"
        write(output_unit, "(a)") &
            "        select case (100 * permutation(1) + &"
        write(output_unit, "(a)") &
            "            10 * permutation(2) + permutation(3))"
        write(output_unit, "(a)") "        case (123); index_ = 1"
        write(output_unit, "(a)") "        case (132); index_ = 2"
        write(output_unit, "(a)") "        case (213); index_ = 3"
        write(output_unit, "(a)") "        case (231); index_ = 4"
        write(output_unit, "(a)") "        case (312); index_ = 5"
        write(output_unit, "(a)") "        case (321); index_ = 6"
        write(output_unit, "(a)") "        end select"
        write(output_unit, "(a)") &
            "    end function face_permutation_index"
        write(output_unit, "(a)") ""
    end subroutine write_dispatch_procedure

    subroutine write_basis_dispatch_procedure()
        write(output_unit, "(a)") &
            "    pure subroutine map_tetra_face_basis_to_local( &"
        write(output_unit, "(a)") &
            "            order, permutation, canonical, local, status)"
        write(output_unit, "(a)") &
            "        integer, intent(in) :: order, permutation(3)"
        write(output_unit, "(a)") &
            "        real(real64), intent(in) :: canonical(:)"
        write(output_unit, "(a)") &
            "        real(real64), intent(out) :: local(:)"
        write(output_unit, "(a)") "        integer, intent(out) :: status"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") "        integer :: permutation_index"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") "        local = 0.0_real64"
        write(output_unit, "(a)") "        status = 1"
        write(output_unit, "(a)") &
            "        permutation_index = face_permutation_index(permutation)"
        write(output_unit, "(a)") &
            "        if (permutation_index == 0) return"
        write(output_unit, "(a)") "        select case (order)"
        write(output_unit, "(a)") "        case (2)"
        write(output_unit, "(a)") &
            "            call basis_to_local_order_2( &"
        write(output_unit, "(a)") &
            "                permutation_index, canonical, local, status)"
        write(output_unit, "(a)") "        case (3)"
        write(output_unit, "(a)") &
            "            call basis_to_local_order_3( &"
        write(output_unit, "(a)") &
            "                permutation_index, canonical, local, status)"
        write(output_unit, "(a)") "        case (4)"
        write(output_unit, "(a)") &
            "            call basis_to_local_order_4( &"
        write(output_unit, "(a)") &
            "                permutation_index, canonical, local, status)"
        write(output_unit, "(a)") "        end select"
        write(output_unit, "(a)") &
            "    end subroutine map_tetra_face_basis_to_local"
        write(output_unit, "(a)") ""
    end subroutine write_basis_dispatch_procedure

    subroutine write_order_header(procedure_name, dof_count)
        character(*), intent(in) :: procedure_name
        integer, intent(in) :: dof_count

        write(output_unit, "(a)") "    pure subroutine "// &
            procedure_name//"(permutation, input, output, status)"
        write(output_unit, "(a)") &
            "        integer, intent(in) :: permutation"
        write(output_unit, "(a)") "        real(real64), intent(in) :: input(:)"
        write(output_unit, "(a)") &
            "        real(real64), intent(out) :: output(:)"
        write(output_unit, "(a)") "        integer, intent(out) :: status"
        write(output_unit, "(a)") "        real(real64) :: transform("// &
            integer_text(dof_count)//","//integer_text(dof_count)//")"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") "        output = 0.0_real64"
        write(output_unit, "(a)") "        status = 1"
        write(output_unit, "(a)") "        if (size(input) /= "// &
            integer_text(dof_count)//") return"
        write(output_unit, "(a)") "        if (size(output) /= "// &
            integer_text(dof_count)//") return"
        write(output_unit, "(a)") "        select case (permutation)"
    end subroutine write_order_header

    subroutine write_transform_case(permutation, transform)
        integer, intent(in) :: permutation
        type(expr_t), intent(in) :: transform(:, :)

        integer :: column, entry, row

        write(output_unit, "(a)") "        case ("// &
            integer_text(permutation)//")"
        write(output_unit, "(a)") "            transform = reshape([ &"
        entry = 0
        do column = 1, size(transform, 2)
            do row = 1, size(transform, 1)
                entry = entry + 1
                if (entry < size(transform)) then
                    write(output_unit, "(a)") "                "// &
                        rational_literal(transform(row, column))//", &"
                else
                    write(output_unit, "(a)") "                "// &
                        rational_literal(transform(row, column))//" &"
                end if
            end do
        end do
        write(output_unit, "(a)") "            ], shape(transform))"
    end subroutine write_transform_case

    subroutine write_order_footer(procedure_name)
        character(*), intent(in) :: procedure_name

        write(output_unit, "(a)") "        end select"
        write(output_unit, "(a)") &
            "        output = matmul(transform, input)"
        write(output_unit, "(a)") "        status = 0"
        write(output_unit, "(a)") "    end subroutine "//procedure_name
        write(output_unit, "(a)") ""
    end subroutine write_order_footer

    subroutine write_module_footer()
        write(output_unit, "(a)") &
            "end module fortfem_generated_tetra_face_moment_transforms"
        close(output_unit)
    end subroutine write_module_footer

    function rational_literal(expression) result(literal)
        type(expr_t), intent(in) :: expression
        character(:), allocatable :: literal

        select case (expression%kind())
        case (NK_INT)
            literal = integer_text_64(expression%int_value())//".0_real64"
        case (NK_RAT)
            literal = integer_text_64(expression%int_value())// &
                ".0_real64 / "// &
                integer_text_64(expression%den_value())//".0_real64"
        case default
            error stop "face transform is not exact rational"
        end select
    end function rational_literal

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text

        character(32) :: buffer

        write(buffer, "(i0)") value
        text = trim(buffer)
    end function integer_text

    function integer_text_64(value) result(text)
        integer(int64), intent(in) :: value
        character(:), allocatable :: text

        character(32) :: buffer

        write(buffer, "(i0)") value
        text = trim(buffer)
    end function integer_text_64

end program gen_tetra_face_moment_transforms
