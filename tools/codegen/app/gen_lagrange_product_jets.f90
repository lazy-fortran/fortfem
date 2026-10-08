program gen_lagrange_product_jets
    !! The application supplies only affine barycentric coordinates and a
    !! Taylor polynomial. FortSym derives every multiplication/derivative jet.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/)
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    integer :: unit

    call arena%init()
    engine = make_native_engine(arena)
    open(newunit=unit, file=generated_path('fortfem_lagrange_product_jets.f90'), &
        status='replace', action='write')
    call product_jet(1, 1)
    call cardinal_product_jet(1)
    call cardinal_product_jet(2)
    call factor_product_jet(2, 1)
    call factor_product_jet(3, 1)
    call factor_product_jet(3, 2)
    call barycentric_jet(2)
    call barycentric_jet(3)
    close(unit)
contains
    subroutine initialize(spec, name, arguments, outputs)
        type(kernel_spec_t), intent(out) :: spec
        character(*), intent(in) :: name
        type(str_t), intent(in) :: arguments(:), outputs(:)
        spec%name = str(name)
        spec%module_name = str('fortfem_'//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str('gen_lagrange_product_jets')
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str('cd tools/codegen && ./generate.sh')
        spec%pure_procedure = .true.
        spec%args = arguments
        spec%outputs = outputs
    end subroutine initialize

    subroutine emit(expressions, spec)
        type(expr_t), intent(inout) :: expressions(:)
        type(kernel_spec_t), intent(in) :: spec
        type(engine_result_t) :: result
        character(:), allocatable :: code
        integer :: i
        do i = 1, size(expressions)
            result = engine%simplify(expressions(i))
            if (.not. result%ok) error stop 'Lagrange jet simplification failed'
            expressions(i) = result%value
        end do
        code = chars(emit_kernel(expressions, spec))
        write(unit, '(a)') code(:len(code) - 1)
    end subroutine emit

    function indexed(prefix, i, j) result(name)
        character(*), intent(in) :: prefix
        integer, intent(in) :: i
        integer, intent(in), optional :: j
        character(:), allocatable :: name
        character(4) :: a, b
        write(a, '(i0)') i
        name = prefix//trim(a)
        if (present(j)) then
            write(b, '(i0)') j
            name = name//trim(b)
        end if
    end function indexed

    subroutine product_jet(dimension, order)
        integer, intent(in) :: dimension, order
        type(kernel_spec_t) :: spec
        type(expr_t) :: t(dimension), zero(dimension), polynomial, affine
        type(expr_t), allocatable :: roots(:)
        type(str_t), allocatable :: args(:), names(:)
        character(:), allocatable :: label, coefficient
        integer :: i, j, count, index, nargs
        count = 1 + dimension
        if (order == 2) count = count + dimension*(dimension + 1)/2
        nargs = count + 1 + dimension
        allocate(args(nargs), names(count), roots(count))
        args(1) = str('a_value'); names(1) = str('next_value')
        polynomial = sym(arena, 'a_value')
        affine = sym(arena, 'b_value')
        args(count + 1) = str('b_value')
        do i = 1, dimension
            t(i) = sym(arena, indexed('t_', i))
            zero(i) = num(arena, 0)
            coefficient = indexed('a_g', i)
            args(1 + i) = str(coefficient)
            names(1 + i) = str(indexed('next_g', i))
            polynomial = polynomial + sym(arena, coefficient)*t(i)
            coefficient = indexed('b_g', i)
            args(count + 1 + i) = str(coefficient)
            affine = affine + sym(arena, coefficient)*t(i)
        end do
        index = 1 + dimension
        if (order == 2) then
            do j = 1, dimension
                do i = 1, j
                    index = index + 1
                    coefficient = indexed('a_h', i, j)
                    args(index) = str(coefficient)
                    names(index) = str(indexed('next_h', i, j))
                    if (i == j) then
                        polynomial = polynomial + &
                            sym(arena, coefficient)*t(i)*t(j)/2
                    else
                        polynomial = polynomial + sym(arena, coefficient)*t(i)*t(j)
                    end if
                end do
            end do
        end if
        polynomial = polynomial*affine
        roots(1) = subs_many(polynomial, t, zero)
        do i = 1, dimension
            roots(1 + i) = subs_many(diff(polynomial, t(i)), t, zero)
        end do
        index = 1 + dimension
        if (order == 2) then
            do j = 1, dimension
                do i = 1, j
                    index = index + 1
                    roots(index) = subs_many(diff(diff(polynomial, t(i)), t(j)), &
                        t, zero)
                end do
            end do
        end if
        label = indexed('generated_product_jet', dimension)// &
            indexed('_order', order)
        call initialize(spec, label, args, names)
        call emit(roots, spec)
    end subroutine product_jet

    subroutine barycentric_jet(dimension)
        integer, intent(in) :: dimension
        type(kernel_spec_t) :: spec
        type(expr_t) :: x(dimension), lambda(dimension + 1)
        type(expr_t) :: roots((dimension + 1)*(dimension + 1))
        type(str_t) :: args(dimension), outputs(size(roots))
        integer :: i, j, index
        lambda(1) = num(arena, 1)
        do i = 1, dimension
            x(i) = sym(arena, indexed('point_', i))
            args(i) = str(indexed('point_', i))
            lambda(1) = lambda(1) - x(i)
            lambda(i + 1) = x(i)
        end do
        index = 0
        do j = 1, dimension + 1
            index = index + 1
            roots(index) = lambda(j)
            outputs(index) = str(indexed('lambda_', j))
            do i = 1, dimension
                index = index + 1
                roots(index) = diff(lambda(j), x(i))
                outputs(index) = str(indexed('gradient_', i, j))
            end do
        end do
        call initialize(spec, indexed('generated_barycentric_jet', dimension), &
            args, outputs)
        call emit(roots, spec)
    end subroutine barycentric_jet

    subroutine cardinal_product_jet(order)
        integer, intent(in) :: order
        type(kernel_spec_t) :: spec
        type(expr_t) :: t(1), zero(1), polynomial, factor
        type(expr_t), allocatable :: roots(:)
        type(str_t), allocatable :: args(:), outputs(:)
        integer :: i
        allocate(roots(order + 1), outputs(order + 1), args(order + 4))
        args(:2) = [str('value'), str('derivative')]
        outputs(:2) = [str('next_value'), str('next_derivative')]
        t(1) = sym(arena, 'increment'); zero(1) = num(arena, 0)
        polynomial = sym(arena, 'value') + sym(arena, 'derivative')*t(1)
        if (order == 2) then
            args(3) = str('second_derivative')
            outputs(3) = str('next_second_derivative')
            polynomial = polynomial + sym(arena, 'second_derivative')*t(1)*t(1)/2
        end if
        args(order + 2:) = [str('degree'), str('factor_index'), str('lambda')]
        factor = sym(arena, 'degree')*(sym(arena, 'lambda') + t(1)) - &
            sym(arena, 'factor_index')
        polynomial = polynomial*factor
        roots(1) = subs_many(polynomial, t, zero)
        do i = 1, order
            polynomial = diff(polynomial, t(1))
            roots(i + 1) = subs_many(polynomial, t, zero)
        end do
        call initialize(spec, indexed('generated_cardinal_product_jet1_order', &
            order), args, outputs)
        call emit(roots, spec)
    end subroutine cardinal_product_jet

    subroutine factor_product_jet(dimension, order)
        integer, intent(in) :: dimension, order
        type(kernel_spec_t) :: spec
        type(expr_t) :: t(dimension), zero(dimension), increments(dimension + 1)
        type(expr_t) :: polynomial, factor, directional
        type(expr_t) :: roots(dimension + 1)
        type(str_t), allocatable :: args(:)
        type(str_t) :: outputs(dimension + 1)
        integer :: i, j, count, index
        count = (dimension + 1)*(order + 1)
        if (order == 2) count = count + dimension
        allocate(args(count))
        increments(1) = num(arena, 0)
        do i = 1, dimension
            t(i) = sym(arena, indexed('factor_increment_', i))
            zero(i) = num(arena, 0)
            increments(1) = increments(1) - t(i)
            increments(i + 1) = t(i)
        end do
        polynomial = num(arena, 1)
        index = 0
        do j = 1, dimension + 1
            index = index + 1
            args(index) = str(indexed('factor_value_', j))
            factor = sym(arena, chars(args(index)))
            index = index + 1
            args(index) = str(indexed('factor_derivative_', j))
            factor = factor + sym(arena, chars(args(index)))*increments(j)
            if (order == 2) then
                index = index + 1
                args(index) = str(indexed('factor_second_derivative_', j))
                factor = factor + &
                    sym(arena, chars(args(index)))*increments(j)*increments(j)/2
            end if
            polynomial = polynomial*factor
        end do
        if (order == 2) then
            directional = num(arena, 0)
            do i = 1, dimension
                index = index + 1
                args(index) = str(indexed('direction_', i))
                directional = directional + diff(polynomial, t(i))* &
                    sym(arena, chars(args(index)))
            end do
            polynomial = directional
            outputs(1) = str('value_dot')
            do i = 1, dimension
                outputs(i + 1) = str(indexed('gradient_dot_', i))
            end do
        else
            outputs(1) = str('value')
            do i = 1, dimension
                outputs(i + 1) = str(indexed('gradient_', i))
            end do
        end if
        roots(1) = subs_many(polynomial, t, zero)
        do i = 1, dimension
            roots(i + 1) = subs_many(diff(polynomial, t(i)), t, zero)
        end do
        call initialize(spec, indexed('generated_factor_product_jet', dimension)// &
            indexed('_order', order), args, outputs)
        call emit(roots, spec)
    end subroutine factor_product_jet
end program gen_lagrange_product_jets
