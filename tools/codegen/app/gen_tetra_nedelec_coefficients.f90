program gen_tetra_nedelec_coefficients
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_arena, only: arena_t, NK_INT, NK_RAT
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, operator(*), operator(+), rat
    use fortsym_linalg, only: &
        exact_linear_system_result_t, solve_exact_linear_system
    use fortfem_codegen_provenance, only: &
        fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    integer :: order, output_unit

    call arena%init()
    engine = make_native_engine(arena)
    call write_module_header()
    do order = 1, 4
        call generate_order(order)
    end do
    call write_module_footer()

contains

    subroutine generate_order(order)
        integer, intent(in) :: order

        type(exact_linear_system_result_t) :: solution
        type(expr_t), allocatable :: identity(:, :), moment_matrix(:, :)
        integer :: diagonal, dof_count

        dof_count = order * (order + 2) * (order + 3) / 2
        allocate(moment_matrix(dof_count, dof_count))
        allocate(identity(dof_count, dof_count))
        call build_exact_moment_matrix(order, moment_matrix)
        identity = num(arena, 0)
        do diagonal = 1, dof_count
            identity(diagonal, diagonal) = num(arena, 1)
        end do
        solution = solve_exact_linear_system( &
            engine, moment_matrix, identity)
        if (.not. solution%ok) then
            error stop "exact tetrahedral Nedelec coefficient solve failed"
        end if
        call write_order_case(order, solution%values)
    end subroutine generate_order

    subroutine build_exact_moment_matrix(order, matrix)
        integer, intent(in) :: order
        type(expr_t), intent(out) :: matrix(:, :)

        integer, allocatable :: candidate_coefficients(:, :)
        integer, allocatable :: candidate_powers(:, :, :)
        integer :: candidate, component, edge, face, moment
        integer :: total_degree, x_degree, y_degree, z_degree

        allocate(candidate_coefficients(3, size(matrix, 2)))
        allocate(candidate_powers(3, 3, size(matrix, 2)))
        call build_candidates( &
            order, candidate_coefficients, candidate_powers)
        matrix = num(arena, 0)
        moment = 0
        do edge = 1, 6
            do total_degree = 0, order - 1
                moment = moment + 1
                do candidate = 1, size(matrix, 2)
                    do component = 1, 3
                        if (candidate_coefficients( &
                            component, candidate) == 0) cycle
                        matrix(moment, candidate) = &
                            matrix(moment, candidate) + num( &
                            arena, candidate_coefficients( &
                            component, candidate)) * &
                            edge_monomial_integral( &
                            candidate_powers(:, component, candidate), &
                            edge, component, total_degree)
                    end do
                end do
            end do
        end do

        do face = 1, 4
            do component = 1, 2
                do total_degree = 0, order - 2
                    do x_degree = 0, total_degree
                        y_degree = total_degree - x_degree
                        moment = moment + 1
                        do candidate = 1, size(matrix, 2)
                            matrix(moment, candidate) = &
                                face_candidate_moment( &
                                candidate_coefficients(:, candidate), &
                                candidate_powers(:, :, candidate), face, &
                                component, x_degree, y_degree)
                        end do
                    end do
                end do
            end do
        end do

        do component = 1, 3
            do total_degree = 0, order - 3
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        moment = moment + 1
                        do candidate = 1, size(matrix, 2)
                            if (candidate_coefficients( &
                                component, candidate) == 0) cycle
                            matrix(moment, candidate) = num( &
                                arena, candidate_coefficients( &
                                component, candidate)) * &
                                tetra_monomial_integral( &
                                candidate_powers(:, component, candidate) + &
                                [x_degree, y_degree, z_degree])
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= size(matrix, 1)) then
            error stop "exact moment matrix dimension mismatch"
        end if
    end subroutine build_exact_moment_matrix

    subroutine build_candidates(order, coefficients, powers)
        integer, intent(in) :: order
        integer, intent(out) :: coefficients(:, :), powers(:, :, :)

        integer :: candidate, component, total_degree
        integer :: x_degree, y_degree, z_degree

        coefficients = 0
        powers = 0
        candidate = 0
        do component = 1, 3
            do total_degree = 0, order - 1
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        candidate = candidate + 1
                        coefficients(component, candidate) = 1
                        powers(:, component, candidate) = &
                            [x_degree, y_degree, z_degree]
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
                    if (component == 4) then
                        coefficients(1, candidate) = -1
                        powers(:, 1, candidate) = &
                            [x_degree, y_degree + 1, z_degree]
                        coefficients(2, candidate) = 1
                        powers(:, 2, candidate) = &
                            [x_degree + 1, y_degree, z_degree]
                    else
                        coefficients(1, candidate) = -1
                        powers(:, 1, candidate) = &
                            [x_degree, y_degree, z_degree + 1]
                        coefficients(3, candidate) = 1
                        powers(:, 3, candidate) = &
                            [x_degree + 1, y_degree, z_degree]
                    end if
                end do
            end do
        end do
        do y_degree = 0, total_degree
            z_degree = total_degree - y_degree
            candidate = candidate + 1
            coefficients(2, candidate) = -1
            powers(:, 2, candidate) = [0, y_degree, z_degree + 1]
            coefficients(3, candidate) = 1
            powers(:, 3, candidate) = [0, y_degree + 1, z_degree]
        end do
        if (candidate /= size(coefficients, 2)) then
            error stop "candidate dimension mismatch"
        end if
    end subroutine build_candidates

    function edge_monomial_integral( &
            powers, edge, component, legendre_degree) result(value)
        integer, intent(in) :: powers(3), edge, component, legendre_degree
        type(expr_t) :: value

        integer(int64) :: monomial(0:8), product(0:12)
        integer(int64) :: legendre(0:3)
        integer :: degree, legendre_power, monomial_power
        integer :: tangent(3), vertices(3, 4), edge_vertices(2, 6)

        call reference_topology(vertices, edge_vertices)
        tangent = vertices(:, edge_vertices(2, edge)) - &
            vertices(:, edge_vertices(1, edge))
        if (tangent(component) == 0) then
            value = num(arena, 0)
            return
        end if
        monomial = 0_int64
        monomial(0) = 1_int64
        do degree = 1, 3
            call multiply_affine_1d( &
                monomial, powers(degree), &
                vertices(degree, edge_vertices(1, edge)), tangent(degree))
        end do
        call shifted_legendre_coefficients(legendre_degree, legendre)
        product = 0_int64
        do monomial_power = 0, 8
            do legendre_power = 0, legendre_degree
                product(monomial_power + legendre_power) = &
                    product(monomial_power + legendre_power) + &
                    monomial(monomial_power) * legendre(legendre_power)
            end do
        end do
        value = num(arena, 0)
        do degree = 0, 12
            if (product(degree) == 0_int64) cycle
            value = value + num(arena, product(degree)) * &
                rat(arena, 1_int64, int(degree + 1, int64))
        end do
        value = num(arena, tangent(component)) * value
    end function edge_monomial_integral

    function face_candidate_moment( &
            coefficients, powers, face, tangent_component, moment_x, &
            moment_y) result(value)
        integer, intent(in) :: coefficients(3), powers(3, 3)
        integer, intent(in) :: face, tangent_component, moment_x, moment_y
        type(expr_t) :: value

        integer :: component, face_vertices(3, 4), tangents(3, 2)
        integer :: vertices(3, 4)

        call reference_topology(vertices)
        call reference_faces(face_vertices)
        tangents(:, 1) = &
            vertices(:, face_vertices(2, face)) - &
            vertices(:, face_vertices(1, face))
        tangents(:, 2) = &
            vertices(:, face_vertices(3, face)) - &
            vertices(:, face_vertices(1, face))
        value = num(arena, 0)
        do component = 1, 3
            if (coefficients(component) == 0) cycle
            if (tangents(component, tangent_component) == 0) cycle
            value = value + num( &
                arena, coefficients(component) * &
                tangents(component, tangent_component)) * &
                face_monomial_integral( &
                powers(:, component), face, moment_x, moment_y)
        end do
    end function face_candidate_moment

    function face_monomial_integral( &
            powers, face, moment_x, moment_y) result(value)
        integer, intent(in) :: powers(3), face, moment_x, moment_y
        type(expr_t) :: value

        integer(int64) :: polynomial(0:4, 0:4)
        integer :: affine(3, 2), component, face_vertices(3, 4)
        integer :: offset(3), vertices(3, 4), x_degree, y_degree

        call reference_topology(vertices)
        call reference_faces(face_vertices)
        offset = vertices(:, face_vertices(1, face))
        affine(:, 1) = vertices(:, face_vertices(2, face)) - offset
        affine(:, 2) = vertices(:, face_vertices(3, face)) - offset
        polynomial = 0_int64
        polynomial(0, 0) = 1_int64
        do component = 1, 3
            call multiply_affine_2d( &
                polynomial, powers(component), offset(component), &
                affine(component, 1), affine(component, 2))
        end do
        value = num(arena, 0)
        do x_degree = 0, 4
            do y_degree = 0, 4 - x_degree
                if (polynomial(x_degree, y_degree) == 0_int64) cycle
                value = value + num( &
                    arena, polynomial(x_degree, y_degree)) * &
                    triangle_monomial_integral( &
                    x_degree + moment_x, y_degree + moment_y)
            end do
        end do
    end function face_monomial_integral

    subroutine multiply_affine_1d( &
            coefficients, power, constant, linear)
        integer(int64), intent(inout) :: coefficients(0:8)
        integer, intent(in) :: power, constant, linear

        integer(int64) :: next(0:8)
        integer :: degree, factor

        do factor = 1, power
            next = 0_int64
            do degree = 0, 8
                next(degree) = next(degree) + &
                    int(constant, int64) * coefficients(degree)
                if (degree < 8) then
                    next(degree + 1) = next(degree + 1) + &
                        int(linear, int64) * coefficients(degree)
                end if
            end do
            coefficients = next
        end do
    end subroutine multiply_affine_1d

    subroutine multiply_affine_2d( &
            coefficients, power, constant, x_coefficient, y_coefficient)
        integer(int64), intent(inout) :: coefficients(0:4, 0:4)
        integer, intent(in) :: power, constant, x_coefficient, y_coefficient

        integer(int64) :: next(0:4, 0:4)
        integer :: factor, x_degree, y_degree

        do factor = 1, power
            next = 0_int64
            do x_degree = 0, 4
                do y_degree = 0, 4 - x_degree
                    next(x_degree, y_degree) = &
                        next(x_degree, y_degree) + &
                        int(constant, int64) * &
                        coefficients(x_degree, y_degree)
                    if (x_degree + y_degree < 4) then
                        next(x_degree + 1, y_degree) = &
                            next(x_degree + 1, y_degree) + &
                            int(x_coefficient, int64) * &
                            coefficients(x_degree, y_degree)
                        next(x_degree, y_degree + 1) = &
                            next(x_degree, y_degree + 1) + &
                            int(y_coefficient, int64) * &
                            coefficients(x_degree, y_degree)
                    end if
                end do
            end do
            coefficients = next
        end do
    end subroutine multiply_affine_2d

    pure subroutine shifted_legendre_coefficients(degree, coefficients)
        integer, intent(in) :: degree
        integer(int64), intent(out) :: coefficients(0:3)

        coefficients = 0_int64
        select case (degree)
        case (0)
            coefficients(0) = 1_int64
        case (1)
            coefficients(0:1) = [-1_int64, 2_int64]
        case (2)
            coefficients(0:2) = [1_int64, -6_int64, 6_int64]
        case (3)
            coefficients = [-1_int64, 12_int64, -30_int64, 20_int64]
        case default
            error stop "unsupported shifted Legendre degree"
        end select
    end subroutine shifted_legendre_coefficients

    function triangle_monomial_integral( &
            x_degree, y_degree) result(value)
        integer, intent(in) :: x_degree, y_degree
        type(expr_t) :: value

        value = rat( &
            arena, factorial(x_degree) * factorial(y_degree), &
            factorial(x_degree + y_degree + 2))
    end function triangle_monomial_integral

    function tetra_monomial_integral(powers) result(value)
        integer, intent(in) :: powers(3)
        type(expr_t) :: value

        value = rat( &
            arena, factorial(powers(1)) * factorial(powers(2)) * &
            factorial(powers(3)), factorial(sum(powers) + 3))
    end function tetra_monomial_integral

    pure function factorial(argument) result(value)
        integer, intent(in) :: argument
        integer(int64) :: value

        integer :: factor

        value = 1_int64
        do factor = 2, argument
            value = value * int(factor, int64)
        end do
    end function factorial

    pure subroutine reference_topology(vertices, edge_vertices)
        integer, intent(out) :: vertices(3, 4)
        integer, intent(out), optional :: edge_vertices(2, 6)

        vertices(:, 1) = [0, 0, 0]
        vertices(:, 2) = [1, 0, 0]
        vertices(:, 3) = [0, 1, 0]
        vertices(:, 4) = [0, 0, 1]
        if (present(edge_vertices)) then
            edge_vertices(:, 1) = [1, 2]
            edge_vertices(:, 2) = [1, 3]
            edge_vertices(:, 3) = [1, 4]
            edge_vertices(:, 4) = [2, 3]
            edge_vertices(:, 5) = [2, 4]
            edge_vertices(:, 6) = [3, 4]
        end if
    end subroutine reference_topology

    pure subroutine reference_faces(face_vertices)
        integer, intent(out) :: face_vertices(3, 4)

        face_vertices(:, 1) = [1, 2, 3]
        face_vertices(:, 2) = [1, 2, 4]
        face_vertices(:, 3) = [1, 3, 4]
        face_vertices(:, 4) = [2, 3, 4]
    end subroutine reference_faces

    subroutine write_module_header()
        integer :: ios

        open( &
            newunit=output_unit, file=generated_path( &
            "fortfem_tetra_nedelec_coefficients.f90"), &
            status="replace", action="write", iostat=ios)
        if (ios /= 0) error stop "cannot write generated coefficients"
        write(output_unit, "(a)") "! Generated by fortsym. Do not edit."
        write(output_unit, "(a)") &
            "! Generator: gen_tetra_nedelec_coefficients"
        write(output_unit, "(a)") &
            "! Generator revision: "//fortsym_revision()
        write(output_unit, "(a)") &
            "! Regenerate with: cd tools/codegen && ./generate.sh"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") &
            "module fortfem_generated_tetra_nedelec_coefficients"
        write(output_unit, "(a)") &
            "    use, intrinsic :: iso_fortran_env, only: real64"
        write(output_unit, "(a)") "    implicit none"
        write(output_unit, "(a)") "    private"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") &
            "    public :: load_tetra_nedelec_coefficients"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") "contains"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") &
            "    subroutine load_tetra_nedelec_coefficients( &"
        write(output_unit, "(a)") &
            "            order, coefficients, status)"
        write(output_unit, "(a)") "        integer, intent(in) :: order"
        write(output_unit, "(a)") &
            "        real(real64), allocatable, intent(out) :: coefficients(:, :)"
        write(output_unit, "(a)") "        integer, intent(out) :: status"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") "        status = 1"
        write(output_unit, "(a)") "        select case (order)"
    end subroutine write_module_header

    subroutine write_order_case(order, coefficients)
        integer, intent(in) :: order
        type(expr_t), intent(in) :: coefficients(:, :)

        integer :: column, row

        write(output_unit, "(a)") "        case ("//integer_text(order)//")"
        write(output_unit, "(a)") "            allocate(coefficients("// &
            integer_text(size(coefficients, 1))//","// &
            integer_text(size(coefficients, 2))//"))"
        write(output_unit, "(a)") &
            "            coefficients = 0.0_real64"
        do column = 1, size(coefficients, 2)
            do row = 1, size(coefficients, 1)
                if (coefficients(row, column)%int_value() == 0_int64) cycle
                write(output_unit, "(a)") "            coefficients("// &
                    integer_text(row)//","//integer_text(column)//") = "// &
                    rational_literal(coefficients(row, column))
            end do
        end do
        write(output_unit, "(a)") "            status = 0"
    end subroutine write_order_case

    subroutine write_module_footer()
        write(output_unit, "(a)") "        end select"
        write(output_unit, "(a)") &
            "    end subroutine load_tetra_nedelec_coefficients"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") &
            "end module fortfem_generated_tetra_nedelec_coefficients"
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
            error stop "basis coefficient is not exact rational"
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

end program gen_tetra_nedelec_coefficients
