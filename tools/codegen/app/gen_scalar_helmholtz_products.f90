program gen_scalar_helmholtz_products
    !! P1 scalar Helmholtz element matrices with diagonal complex coefficients,
    !! sum_d g_d (d_d v, d_d u) - k**2 m (v, u), on intervals, triangles and
    !! tetrahedra, with JVP and VJP. g_d and m are the coefficients of a diagonal
    !! complex-coordinate stretch (fortfem_cartesian_helmholtz_pml); g = m = 1 is
    !! the plain Helmholtz operator. The wave number k is real.
    !!
    !! Interval geometry is the length h. Triangle geometry uses the affine-map
    !! entries a = x2 - x1, b = y2 - y1, c = x3 - x1, d = y3 - y1 as in
    !! gen_p1_elasticity_products; tetrahedron geometry the Jacobian columns
    !! j(:, k) = x_(k+1) - x_1, passed as j11, j21, j31, j12, ... The measure is
    !! orientation*det/(dim!) with orientation = sign(det) a constant input, so the
    !! kernels are rational in the geometry and exact for a nondegenerate cell.
    !! Gradients are adj(J)^T/det, written as numerators over det.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, sym, operator(+), operator(-), &
        operator(*), operator(/)
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code
    character(*), parameter :: interval_geometry(1) = [character(6) :: "length"]
    character(*), parameter :: triangle_geometry(4) = [character(1) :: "a", "b", &
        "c", "d"]
    character(*), parameter :: tetra_geometry(9) = [character(3) :: "j11", &
        "j21", "j31", "j12", "j22", "j32", "j13", "j23", "j33"]
    integer :: unit

    call arena%init()
    engine = make_native_engine(arena)
    open (newunit=unit, file=generated_path("fortfem_scalar_helmholtz_products.f90"), &
        status="replace", action="write")
    call emit_family("p1_helmholtz_interval", interval_geometry, .false., 1, 2)
    call emit_family("p1_helmholtz_triangle", triangle_geometry, .true., 2, 3)
    call emit_family("p1_helmholtz_tetra", tetra_geometry, .true., 3, 4)
    close (unit)

contains

    subroutine emit_family(family, geometry_names, oriented, dimension, order)
        !! Primal, JVP and VJP of one element family. Inputs: geometry,
        !! [orientation], gradient_1..dimension, mass, wave_number.
        character(*), intent(in) :: family, geometry_names(:)
        logical, intent(in) :: oriented
        integer, intent(in) :: dimension, order

        type(expr_t), allocatable :: geometry(:), inputs(:), dots(:), bars(:)
        type(expr_t), allocatable :: cell(:), cell_jvp(:), cell_vjp(:)
        type(expr_t) :: orientation, gradient(3), mass, wave_number
        type(str_t), allocatable :: names(:), base_args(:), base_types(:)
        integer :: k, n_geometry, n_inputs

        n_geometry = size(geometry_names)
        n_inputs = n_geometry + dimension + 2
        allocate (geometry(n_geometry), names(n_inputs))
        do k = 1, n_geometry
            names(k) = str(trim(geometry_names(k)))
            geometry(k) = sym(arena, trim(geometry_names(k)))
        end do
        do k = 1, dimension
            names(n_geometry + k) = str(gradient_name(k))
            gradient(k) = sym(arena, gradient_name(k))
        end do
        names(n_inputs - 1) = str("mass")
        names(n_inputs) = str("wave_number")
        mass = sym(arena, "mass")
        wave_number = sym(arena, "wave_number")
        orientation = sym(arena, "orientation")
        inputs = [geometry, gradient(1:dimension), mass, wave_number]

        select case (dimension)
        case (1)
            call interval_cell(geometry(1), gradient(1), mass, wave_number, cell)
        case (2)
            call triangle_cell(geometry, orientation, gradient(1:2), mass, &
                wave_number, cell)
        case default
            call tetra_cell(geometry, orientation, gradient(1:3), mass, &
                wave_number, cell)
        end select

        allocate (dots(n_inputs), bars(order*order))
        do k = 1, n_inputs
            dots(k) = sym(arena, chars(names(k))//"_dot")
        end do
        do k = 1, order*order
            bars(k) = sym(arena, matrix_entry("cell_bar", order, k))
        end do
        cell_jvp = jvp(cell, inputs, dots)
        cell_vjp = vjp(cell, inputs, bars)
        call simplify(cell)
        call simplify(cell_jvp)
        call simplify(cell_vjp)

        ! Argument lists: geometry (real), [orientation (real)], coefficients
        ! (complex), wave number (real).
        if (oriented) then
            base_args = [names(1:n_geometry), str("orientation"), &
                names(n_geometry + 1:n_inputs)]
            base_types = [real_types(n_geometry + 1), complex_types(dimension + 1), &
                real_types(1)]
        else
            base_args = names
            base_types = [real_types(n_geometry), complex_types(dimension + 1), &
                real_types(1)]
        end if

        call initialize(family)
        spec%args = base_args
        spec%arg_types = base_types
        spec%outputs = [str("cell")]
        call set_matrix_output(order)
        call write_kernel(cell)

        call initialize(family//"_jvp")
        spec%args = [base_args, (str(chars(names(k))//"_dot"), k=1, n_inputs)]
        spec%arg_types = [base_types, real_types(n_geometry), &
            complex_types(dimension + 1), real_types(1)]
        spec%outputs = [str("cell_dot")]
        call set_matrix_output(order)
        call write_kernel(cell_jvp)

        call initialize(family//"_vjp")
        spec%args = [base_args, str("cell_bar")]
        spec%arg_types = [base_types, complex_types(1)]
        spec%arg_shapes = [no_shapes(size(base_args)), str(matrix_shape(order))]
        spec%outputs = [(str(chars(names(k))//"_bar"), k=1, n_inputs)]
        call write_kernel(cell_vjp)
    end subroutine emit_family

    subroutine interval_cell(length, gradient, mass, wave_number, cell)
        type(expr_t), intent(in) :: length, gradient, mass, wave_number
        type(expr_t), allocatable, intent(out) :: cell(:)
        integer :: i, j

        allocate (cell(4))
        do j = 1, 2
            do i = 1, 2
                if (i == j) then
                    cell(2*(j - 1) + i) = gradient/length - &
                        wave_number*wave_number*mass*length/3
                else
                    cell(2*(j - 1) + i) = -gradient/length - &
                        wave_number*wave_number*mass*length/6
                end if
            end do
        end do
    end subroutine interval_cell

    subroutine triangle_cell(geometry, orientation, gradient, mass, wave_number, &
            cell)
        type(expr_t), intent(in) :: geometry(4), orientation, gradient(2)
        type(expr_t), intent(in) :: mass, wave_number
        type(expr_t), allocatable, intent(out) :: cell(:)

        type(expr_t) :: a, b, c, d, det, numerator(2, 3), stiffness
        integer :: i, j

        a = geometry(1); b = geometry(2); c = geometry(3); d = geometry(4)
        det = a*d - b*c
        ! grad(lambda_i) = numerator(:, i)/det
        numerator(:, 1) = [b - d, c - a]
        numerator(:, 2) = [d, num(arena, 0) - c]
        numerator(:, 3) = [num(arena, 0) - b, a]
        allocate (cell(9))
        do j = 1, 3
            do i = 1, 3
                ! area*grad_i.G.grad_j with area = orientation*det/2
                stiffness = orientation*(gradient(1)*numerator(1, i)* &
                    numerator(1, j) + gradient(2)*numerator(2, i)* &
                    numerator(2, j))/(2*det)
                cell(3*(j - 1) + i) = stiffness - wave_number*wave_number*mass* &
                    orientation*det*mass_weight(i, j, 24)
            end do
        end do
    end subroutine triangle_cell

    subroutine tetra_cell(geometry, orientation, gradient, mass, wave_number, cell)
        type(expr_t), intent(in) :: geometry(9), orientation, gradient(3)
        type(expr_t), intent(in) :: mass, wave_number
        type(expr_t), allocatable, intent(out) :: cell(:)

        type(expr_t) :: e(3, 3), det, numerator(3, 4), stiffness
        integer :: i, j, p

        e = reshape(geometry, [3, 3])
        ! Rows of J^-1 are the cross products of the Jacobian columns over det.
        numerator(:, 2) = cross(e(:, 2), e(:, 3))
        numerator(:, 3) = cross(e(:, 3), e(:, 1))
        numerator(:, 4) = cross(e(:, 1), e(:, 2))
        do p = 1, 3
            numerator(p, 1) = num(arena, 0) - numerator(p, 2) - numerator(p, 3) - &
                numerator(p, 4)
        end do
        det = e(1, 1)*numerator(1, 2) + e(2, 1)*numerator(2, 2) + &
            e(3, 1)*numerator(3, 2)
        allocate (cell(16))
        do j = 1, 4
            do i = 1, 4
                ! volume*grad_i.G.grad_j with volume = orientation*det/6
                stiffness = gradient(1)*numerator(1, i)*numerator(1, j)
                do p = 2, 3
                    stiffness = stiffness + gradient(p)*numerator(p, i)* &
                        numerator(p, j)
                end do
                stiffness = orientation*stiffness/(6*det)
                cell(4*(j - 1) + i) = stiffness - wave_number*wave_number*mass* &
                    orientation*det*mass_weight(i, j, 120)
            end do
        end do
    end subroutine tetra_cell

    function cross(u, v) result(w)
        type(expr_t), intent(in) :: u(3), v(3)
        type(expr_t) :: w(3)

        w(1) = u(2)*v(3) - u(3)*v(2)
        w(2) = u(3)*v(1) - u(1)*v(3)
        w(3) = u(1)*v(2) - u(2)*v(1)
    end function cross

    function mass_weight(i, j, denominator) result(weight)
        !! P1 mass |det|*(1 + delta_ij)/denominator (24 triangle, 120 tetra).
        integer, intent(in) :: i, j, denominator
        type(expr_t) :: weight

        if (i == j) then
            weight = num(arena, 2)/denominator
        else
            weight = num(arena, 1)/denominator
        end if
    end function mass_weight

    function gradient_name(k) result(name)
        integer, intent(in) :: k
        character(:), allocatable :: name
        character(16) :: buffer

        write (buffer, "('gradient_',i0)") k
        name = trim(buffer)
    end function gradient_name

    function matrix_shape(order) result(text)
        integer, intent(in) :: order
        character(:), allocatable :: text
        character(8) :: buffer

        write (buffer, "('(',i0,',',i0,')')") order, order
        text = trim(buffer)
    end function matrix_shape

    subroutine set_matrix_output(order)
        integer, intent(in) :: order
        integer :: entry

        spec%output_shapes = [str(matrix_shape(order))]
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
        spec%generator = str("gen_scalar_helmholtz_products")
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
            error stop "scalar Helmholtz kernel emission failed"
        end if
        write (unit, "(a)") code(:len(code) - 1)
    end subroutine write_kernel

    subroutine simplify(roots)
        type(expr_t), intent(inout) :: roots(:)
        type(engine_result_t) :: result
        integer :: root

        do root = 1, size(roots)
            result = engine%simplify(roots(root))
            if (.not. result%ok) error stop "scalar Helmholtz simplification failed"
            roots(root) = result%value
        end do
    end subroutine simplify

end program gen_scalar_helmholtz_products
