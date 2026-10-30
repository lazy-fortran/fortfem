program gen_p1_elasticity_products
    !! Vector P1 time-harmonic elasticity: the 6x6 triangle matrix K - inertia*M
    !! and the 4x4 impedance (absorbing) edge matrix, with JVP and VJP.
    !!
    !! Triangle geometry uses the affine-map entries of gen_element_precision_products,
    !! a = x2 - x1, b = y2 - y1, c = x3 - x1, d = y3 - y1, and the local dof order
    !! (u_x, u_y) per vertex. The area is orientation*det/2: orientation = sign(det)
    !! is a constant input instead of abs(det), so the kernels stay polynomial in the
    !! geometry and d|det| = orientation*d(det) holds exactly for a nondegenerate
    !! triangle.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, i_expr, num, sym, operator(+), &
        operator(-), operator(*), operator(/)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, &
        KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str, str_t
    use fortsym_subs, only: subs_many
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    use fortfem_codegen_triangle_basis, only: triangle_scalar_bases
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: cell(36), cell_jvp(36), cell_vjp(7)
    type(expr_t) :: cell_inputs(7), cell_dots(7), cell_bars(36)
    type(expr_t) :: edge(16), edge_jvp(16), edge_vjp(6)
    type(expr_t) :: edge_inputs(6), edge_dots(6), edge_bars(16)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code
    character(*), parameter :: cell_names(7) = [character(11) :: "a", "b", "c", &
        "d", "lambda", "mu", "inertia"]
    character(*), parameter :: edge_names(6) = [character(17) :: "length", &
        "normal_1", "normal_2", "p_impedance", "s_impedance", "angular_frequency"]
    integer :: k, unit

    call arena%init()
    engine = make_native_engine(arena)
    call build_cell()
    call build_edge()
    do k = 1, 36
        cell_bars(k) = sym(arena, matrix_entry("cell_bar", 6, k))
    end do
    do k = 1, 16
        edge_bars(k) = sym(arena, matrix_entry("edge_bar", 4, k))
    end do
    do k = 1, 7
        cell_dots(k) = sym(arena, trim(cell_names(k))//"_dot")
    end do
    do k = 1, 6
        edge_dots(k) = sym(arena, trim(edge_names(k))//"_dot")
    end do
    cell_jvp = jvp(cell, cell_inputs, cell_dots)
    cell_vjp = vjp(cell, cell_inputs, cell_bars)
    edge_jvp = jvp(edge, edge_inputs, edge_dots)
    edge_vjp = vjp(edge, edge_inputs, edge_bars)
    call simplify(cell)
    call simplify(cell_jvp)
    call simplify(cell_vjp)
    call simplify(edge)
    call simplify(edge_jvp)
    call simplify(edge_vjp)

    open (newunit=unit, file=generated_path("fortfem_p1_elasticity_products.f90"), &
        status="replace", action="write")
    call initialize("p1_elasticity_cell")
    call set_cell_args(.false.)
    spec%outputs = [str("cell")]
    call set_matrix_output(6)
    call write_kernel(cell)

    call initialize("p1_elasticity_cell_jvp")
    call set_cell_args(.true.)
    spec%outputs = [str("cell_dot")]
    call set_matrix_output(6)
    call write_kernel(cell_jvp)

    call initialize("p1_elasticity_cell_vjp")
    spec%args = [named(cell_names(1:4)), str("orientation"), &
        named(cell_names(5:7)), str("cell_bar")]
    spec%arg_types = [real_types(5), complex_types(3), complex_types(1)]
    spec%arg_shapes = [no_shapes(8), str("(6,6)")]
    spec%outputs = [(str(trim(cell_names(k))//"_bar"), k=1, 7)]
    call write_kernel(cell_vjp)

    call initialize("p1_elastic_impedance_edge")
    call set_edge_args(.false.)
    spec%outputs = [str("edge")]
    call set_matrix_output(4)
    call write_kernel(edge)

    call initialize("p1_elastic_impedance_edge_jvp")
    call set_edge_args(.true.)
    spec%outputs = [str("edge_dot")]
    call set_matrix_output(4)
    call write_kernel(edge_jvp)

    call initialize("p1_elastic_impedance_edge_vjp")
    spec%args = [named(edge_names), str("edge_bar")]
    spec%arg_types = [real_types(3), complex_types(2), real_types(1), &
        complex_types(1)]
    spec%arg_shapes = [no_shapes(6), str("(4,4)")]
    spec%outputs = [(str(trim(edge_names(k))//"_bar"), k=1, 6)]
    call write_kernel(edge_vjp)
    close (unit)

contains

    subroutine build_cell()
        type(expr_t) :: a, b, c, d, orientation, lambda, mu, inertia, det, area
        type(expr_t) :: xi, eta, l(3), p2(6), grad(3, 2), qx(3), qy(3)
        type(expr_t) :: mass, stiffness
        integer :: i, j, p, q, point

        a = sym(arena, "a"); b = sym(arena, "b")
        c = sym(arena, "c"); d = sym(arena, "d")
        orientation = sym(arena, "orientation")
        lambda = sym(arena, "lambda"); mu = sym(arena, "mu")
        inertia = sym(arena, "inertia")
        cell_inputs = [a, b, c, d, lambda, mu, inertia]
        xi = sym(arena, "xi"); eta = sym(arena, "eta")
        det = a*d - b*c
        area = orientation*det/2
        call triangle_scalar_bases(xi, eta, l, p2)
        do i = 1, 3
            grad(i, 1) = (d*diff(l(i), xi) - b*diff(l(i), eta))/det
            grad(i, 2) = (-c*diff(l(i), xi) + a*diff(l(i), eta))/det
        end do
        ! Degree-2 rule on the reference triangle (weights 1/6): exact for l_i l_j.
        qx = [num(arena, 1)/6, num(arena, 2)/3, num(arena, 1)/6]
        qy = [num(arena, 1)/6, num(arena, 1)/6, num(arena, 2)/3]
        do j = 1, 3
            do i = 1, 3
                mass = num(arena, 0)
                do point = 1, 3
                    mass = mass + subs_many(l(i)*l(j), [xi, eta], &
                        [qx(point), qy(point)])/6
                end do
                mass = orientation*det*mass
                do q = 1, 2
                    do p = 1, 2
                        ! Row (i, p), column (j, q) of eps(v):D:eps(u), D isotropic.
                        stiffness = lambda*grad(i, p)*grad(j, q) + mu*( &
                            grad(i, q)*grad(j, p))
                        if (p == q) stiffness = stiffness + mu*( &
                            grad(i, 1)*grad(j, 1) + grad(i, 2)*grad(j, 2))
                        stiffness = area*stiffness
                        if (p == q) stiffness = stiffness - inertia*mass
                        cell(6*(2*j - 2 + q - 1) + 2*i - 2 + p) = stiffness
                    end do
                end do
            end do
        end do
    end subroutine build_cell

    subroutine build_edge()
        type(expr_t) :: length, normal(2), p_impedance, s_impedance, frequency
        type(expr_t) :: coefficient, mass
        integer :: i, j, p, q

        length = sym(arena, "length")
        normal = [sym(arena, "normal_1"), sym(arena, "normal_2")]
        p_impedance = sym(arena, "p_impedance")
        s_impedance = sym(arena, "s_impedance")
        frequency = sym(arena, "angular_frequency")
        edge_inputs = [length, normal, p_impedance, s_impedance, frequency]
        do j = 1, 2
            do i = 1, 2
                if (i == j) then
                    mass = length/3
                else
                    mass = length/6
                end if
                do q = 1, 2
                    do p = 1, 2
                        ! i*omega*M_ij*(Z_s I + (Z_p - Z_s) n n^T)_pq
                        coefficient = (p_impedance - s_impedance)* &
                            normal(p)*normal(q)
                        if (p == q) coefficient = coefficient + s_impedance
                        edge(4*(2*j - 2 + q - 1) + 2*i - 2 + p) = &
                            i_expr(arena)*frequency*mass*coefficient
                    end do
                end do
            end do
        end do
    end subroutine build_edge

    subroutine set_cell_args(tangent)
        logical, intent(in) :: tangent

        if (tangent) then
            spec%args = [named(cell_names(1:4)), str("orientation"), &
                named(cell_names(5:7)), &
                (str(trim(cell_names(k))//"_dot"), k=1, 7)]
            spec%arg_types = [real_types(5), complex_types(3), real_types(4), &
                complex_types(3)]
        else
            spec%args = [named(cell_names(1:4)), str("orientation"), &
                named(cell_names(5:7))]
            spec%arg_types = [real_types(5), complex_types(3)]
        end if
    end subroutine set_cell_args

    subroutine set_edge_args(tangent)
        logical, intent(in) :: tangent

        if (tangent) then
            spec%args = [named(edge_names), &
                (str(trim(edge_names(k))//"_dot"), k=1, 6)]
            spec%arg_types = [real_types(3), complex_types(2), real_types(1), &
                real_types(3), complex_types(2), real_types(1)]
        else
            spec%args = named(edge_names)
            spec%arg_types = [real_types(3), complex_types(2), real_types(1)]
        end if
    end subroutine set_edge_args

    subroutine set_matrix_output(order)
        integer, intent(in) :: order
        character(8) :: shape
        integer :: entry

        write (shape, "('(',i0,',',i0,')')") order, order
        spec%output_shapes = [str(trim(shape))]
        allocate (spec%output_references(order*order))
        do entry = 1, order*order
            spec%output_references(entry) = str(matrix_entry( &
                chars(spec%outputs(1)), order, entry))
        end do
    end subroutine set_matrix_output

    function matrix_entry(name, order, entry) result(text)
        character(*), intent(in) :: name
        integer, intent(in) :: order, entry
        character(:), allocatable :: text
        character(16) :: indices

        write (indices, "('(',i0,',',i0,')')") mod(entry - 1, order) + 1, &
            (entry - 1)/order + 1
        text = name//trim(indices)
    end function matrix_entry

    function named(names) result(strings)
        character(*), intent(in) :: names(:)
        type(str_t) :: strings(size(names))
        integer :: entry

        do entry = 1, size(names)
            strings(entry) = str(trim(names(entry)))
        end do
    end function named

    function real_types(count) result(types)
        integer, intent(in) :: count
        type(str_t) :: types(count)

        types = str("real(dp)")
    end function real_types

    function complex_types(count) result(types)
        integer, intent(in) :: count
        type(str_t) :: types(count)

        types = str("complex(dp)")
    end function complex_types

    function no_shapes(count) result(shapes)
        integer, intent(in) :: count
        type(str_t) :: shapes(count)

        shapes = str("")
    end function no_shapes

    subroutine initialize(family)
        character(*), intent(in) :: family

        spec = kernel_spec_t()
        spec%name = str("generated_"//family)
        spec%module_name = str("fortfem_generated_"//family)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_p1_elasticity_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%scalar_type = str("complex(dp)")
        spec%pure_procedure = .true.
    end subroutine initialize

    subroutine write_kernel(roots)
        type(expr_t), intent(in) :: roots(:)
        character(:), allocatable :: message
        logical :: ok

        code = chars(emit_kernel(roots, spec, ok=ok, message=message))
        if (.not. ok) then
            write (*, "(a)") chars(spec%name)//": "//message
            error stop "elasticity kernel emission failed"
        end if
        write (unit, "(a)") code(:len(code) - 1)
    end subroutine write_kernel

    subroutine simplify(roots)
        type(expr_t), intent(inout) :: roots(:)
        type(engine_result_t) :: result
        integer :: root

        do root = 1, size(roots)
            result = engine%simplify(roots(root))
            if (.not. result%ok) error stop "elasticity simplification failed"
            roots(root) = result%value
        end do
    end subroutine simplify

end program gen_p1_elasticity_products
