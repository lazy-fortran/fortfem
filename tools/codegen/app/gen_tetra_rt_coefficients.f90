program gen_tetra_rt_coefficients
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
    integer :: degree, output_unit

    call arena%init()
    engine = make_native_engine(arena)
    call write_module_header()
    do degree = 0, 4
        call generate_degree(degree)
    end do
    call write_module_footer()

contains

    subroutine generate_degree(degree)
        integer, intent(in) :: degree

        type(exact_linear_system_result_t) :: solution
        type(expr_t), allocatable :: identity(:, :), moment_matrix(:, :)
        integer :: diagonal, dof_count

        dof_count = 3*(degree + 1)*(degree + 2)*(degree + 3)/6 + &
            (degree + 1)*(degree + 2)/2
        allocate(moment_matrix(dof_count, dof_count))
        allocate(identity(dof_count, dof_count))
        call build_exact_moment_matrix(degree, moment_matrix)
        identity = num(arena, 0)
        do diagonal = 1, dof_count
            identity(diagonal, diagonal) = num(arena, 1)
        end do
        solution = solve_exact_linear_system( &
            engine, moment_matrix, identity)
        if (.not. solution%ok) then
            error stop "exact tetrahedral RT coefficient solve failed"
        end if
        call write_degree_case(degree, solution%values)
    end subroutine generate_degree

    subroutine build_exact_moment_matrix(degree, matrix)
        integer, intent(in) :: degree
        type(expr_t), intent(out) :: matrix(:, :)

        integer, allocatable :: coefficients(:, :)
        integer, allocatable :: powers(:, :, :)
        integer :: candidate, component, face, moment, total
        integer :: x_degree, y_degree, z_degree

        allocate(coefficients(3, size(matrix, 2)))
        allocate(powers(3, 3, size(matrix, 2)))
        call build_candidates(degree, coefficients, powers)
        matrix = num(arena, 0)
        moment = 0
        do face = 1, 4
            do total = 0, degree
                do x_degree = 0, total
                    y_degree = total - x_degree
                    moment = moment + 1
                    do candidate = 1, size(matrix, 2)
                        matrix(moment, candidate) = face_candidate_moment( &
                            coefficients(:, candidate), &
                            powers(:, :, candidate), face, &
                            x_degree, y_degree)
                    end do
                end do
            end do
        end do
        do component = 1, 3
            do total = 0, degree - 1
                do x_degree = 0, total
                    do y_degree = 0, total - x_degree
                        z_degree = total - x_degree - y_degree
                        moment = moment + 1
                        do candidate = 1, size(matrix, 2)
                            if (coefficients(component, candidate) == 0) cycle
                            matrix(moment, candidate) = num( &
                                arena, coefficients(component, candidate))* &
                                tetra_monomial_integral( &
                                powers(:, component, candidate) + &
                                [x_degree, y_degree, z_degree])
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= size(matrix, 1)) then
            error stop "exact RT moment matrix dimension mismatch"
        end if
    end subroutine build_exact_moment_matrix

    subroutine build_candidates(degree, coefficients, powers)
        integer, intent(in) :: degree
        integer, intent(out) :: coefficients(:, :), powers(:, :, :)

        integer :: candidate, component, total
        integer :: x_degree, y_degree, z_degree

        coefficients = 0
        powers = 0
        candidate = 0
        do component = 1, 3
            do total = 0, degree
                do x_degree = 0, total
                    do y_degree = 0, total - x_degree
                        z_degree = total - x_degree - y_degree
                        candidate = candidate + 1
                        coefficients(component, candidate) = 1
                        powers(:, component, candidate) = &
                            [x_degree, y_degree, z_degree]
                    end do
                end do
            end do
        end do
        total = degree
        do x_degree = 0, total
            do y_degree = 0, total - x_degree
                z_degree = total - x_degree - y_degree
                candidate = candidate + 1
                coefficients(:, candidate) = 1
                powers(:, 1, candidate) = &
                    [x_degree + 1, y_degree, z_degree]
                powers(:, 2, candidate) = &
                    [x_degree, y_degree + 1, z_degree]
                powers(:, 3, candidate) = &
                    [x_degree, y_degree, z_degree + 1]
            end do
        end do
        if (candidate /= size(coefficients, 2)) then
            error stop "RT candidate dimension mismatch"
        end if
    end subroutine build_candidates

    function face_candidate_moment( &
            coefficients, powers, face, moment_x, moment_y) result(value)
        integer, intent(in) :: coefficients(3), powers(3, 3)
        integer, intent(in) :: face, moment_x, moment_y
        type(expr_t) :: value

        integer :: component, area_normal(3)

        call face_geometry(face, area_normal=area_normal)
        value = num(arena, 0)
        do component = 1, 3
            if (coefficients(component) == 0 .or. &
                area_normal(component) == 0) cycle
            value = value + num( &
                arena, coefficients(component)*area_normal(component))* &
                face_monomial_integral( &
                powers(:, component), face, moment_x, moment_y)
        end do
    end function face_candidate_moment

    function face_monomial_integral( &
            powers, face, moment_x, moment_y) result(value)
        integer, intent(in) :: powers(3), face, moment_x, moment_y
        type(expr_t) :: value

        integer(int64) :: polynomial(0:8, 0:8)
        integer :: affine(3, 2), component, offset(3)
        integer :: x_degree, y_degree

        call face_geometry(face, offset, affine)
        polynomial = 0_int64
        polynomial(0, 0) = 1_int64
        do component = 1, 3
            call multiply_affine_2d( &
                polynomial, powers(component), offset(component), &
                affine(component, 1), affine(component, 2))
        end do
        value = num(arena, 0)
        do x_degree = 0, 8
            do y_degree = 0, 8 - x_degree
                if (polynomial(x_degree, y_degree) == 0_int64) cycle
                value = value + num( &
                    arena, polynomial(x_degree, y_degree))* &
                    triangle_monomial_integral( &
                    x_degree + moment_x, y_degree + moment_y)
            end do
        end do
    end function face_monomial_integral

    subroutine multiply_affine_2d( &
            coefficients, power, constant, x_coefficient, y_coefficient)
        integer(int64), intent(inout) :: coefficients(0:8, 0:8)
        integer, intent(in) :: power, constant, x_coefficient, y_coefficient

        integer(int64) :: next(0:8, 0:8)
        integer :: factor, x_degree, y_degree

        do factor = 1, power
            next = 0_int64
            do x_degree = 0, 8
                do y_degree = 0, 8 - x_degree
                    next(x_degree, y_degree) = &
                        next(x_degree, y_degree) + &
                        int(constant, int64)* &
                        coefficients(x_degree, y_degree)
                    if (x_degree + y_degree < 8) then
                        next(x_degree + 1, y_degree) = &
                            next(x_degree + 1, y_degree) + &
                            int(x_coefficient, int64)* &
                            coefficients(x_degree, y_degree)
                        next(x_degree, y_degree + 1) = &
                            next(x_degree, y_degree + 1) + &
                            int(y_coefficient, int64)* &
                            coefficients(x_degree, y_degree)
                    end if
                end do
            end do
            coefficients = next
        end do
    end subroutine multiply_affine_2d

    pure subroutine face_geometry( &
            face, offset, affine, area_normal)
        integer, intent(in) :: face
        integer, intent(out), optional :: offset(3), affine(3, 2)
        integer, intent(out), optional :: area_normal(3)

        integer :: local_offset(3), local_affine(3, 2), local_normal(3)

        select case (face)
        case (1)
            local_offset = [0, 0, 0]
            local_affine(:, 1) = [0, 1, 0]
            local_affine(:, 2) = [0, 0, 1]
            local_normal = [-1, 0, 0]
        case (2)
            local_offset = [0, 0, 0]
            local_affine(:, 1) = [1, 0, 0]
            local_affine(:, 2) = [0, 0, 1]
            local_normal = [0, -1, 0]
        case (3)
            local_offset = [0, 0, 0]
            local_affine(:, 1) = [0, 1, 0]
            local_affine(:, 2) = [1, 0, 0]
            local_normal = [0, 0, -1]
        case (4)
            local_offset = [1, 0, 0]
            local_affine(:, 1) = [-1, 1, 0]
            local_affine(:, 2) = [-1, 0, 1]
            local_normal = [1, 1, 1]
        case default
            error stop "invalid reference tetrahedron face"
        end select
        if (present(offset)) offset = local_offset
        if (present(affine)) affine = local_affine
        if (present(area_normal)) area_normal = local_normal
    end subroutine face_geometry

    function triangle_monomial_integral( &
            x_degree, y_degree) result(value)
        integer, intent(in) :: x_degree, y_degree
        type(expr_t) :: value

        value = rat( &
            arena, factorial(x_degree)*factorial(y_degree), &
            factorial(x_degree + y_degree + 2))
    end function triangle_monomial_integral

    function tetra_monomial_integral(powers) result(value)
        integer, intent(in) :: powers(3)
        type(expr_t) :: value

        value = rat( &
            arena, factorial(powers(1))*factorial(powers(2))* &
            factorial(powers(3)), factorial(sum(powers) + 3))
    end function tetra_monomial_integral

    pure integer(int64) function factorial(argument) result(value)
        integer, intent(in) :: argument
        integer :: factor

        value = 1_int64
        do factor = 2, argument
            value = value*int(factor, int64)
        end do
    end function factorial

    subroutine write_module_header()
        integer :: ios

        open( &
            newunit=output_unit, file=generated_path( &
            "fortfem_tetra_rt_coefficients.f90"), &
            status="replace", action="write", iostat=ios)
        if (ios /= 0) error stop "cannot write generated RT coefficients"
        write(output_unit, "(a)") "! Generated by fortsym. Do not edit."
        write(output_unit, "(a)") "! Generator: gen_tetra_rt_coefficients"
        write(output_unit, "(a)") &
            "! Generator revision: "//fortsym_revision()
        write(output_unit, "(a)") &
            "! Regenerate with: cd tools/codegen && ./generate.sh"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") &
            "module fortfem_generated_tetra_rt_coefficients"
        write(output_unit, "(a)") &
            "    use, intrinsic :: iso_fortran_env, only: real64"
        write(output_unit, "(a)") "    implicit none"
        write(output_unit, "(a)") "    private"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") &
            "    public :: load_tetra_rt_coefficients"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") "contains"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") &
            "    subroutine load_tetra_rt_coefficients( &"
        write(output_unit, "(a)") "            degree, coefficients, status)"
        write(output_unit, "(a)") "        integer, intent(in) :: degree"
        write(output_unit, "(a)") &
            "        real(real64), allocatable, intent(out) :: coefficients(:, :)"
        write(output_unit, "(a)") "        integer, intent(out) :: status"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") "        status = 1"
        write(output_unit, "(a)") "        select case (degree)"
    end subroutine write_module_header

    subroutine write_degree_case(degree, coefficients)
        integer, intent(in) :: degree
        type(expr_t), intent(in) :: coefficients(:, :)

        integer :: column, row

        write(output_unit, "(a)") "        case ("//integer_text(degree)//")"
        write(output_unit, "(a)") "            allocate(coefficients("// &
            integer_text(size(coefficients, 1))//","// &
            integer_text(size(coefficients, 2))//"))"
        write(output_unit, "(a)") "            coefficients = 0.0_real64"
        do column = 1, size(coefficients, 2)
            do row = 1, size(coefficients, 1)
                if (coefficients(row, column)%int_value() == 0_int64) cycle
                write(output_unit, "(a)") "            coefficients("// &
                    integer_text(row)//","//integer_text(column)//") = "// &
                    rational_literal(coefficients(row, column))
            end do
        end do
        write(output_unit, "(a)") "            status = 0"
    end subroutine write_degree_case

    subroutine write_module_footer()
        write(output_unit, "(a)") "        end select"
        write(output_unit, "(a)") &
            "    end subroutine load_tetra_rt_coefficients"
        write(output_unit, "(a)") ""
        write(output_unit, "(a)") &
            "end module fortfem_generated_tetra_rt_coefficients"
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
            error stop "RT basis coefficient is not exact rational"
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
end program gen_tetra_rt_coefficients
