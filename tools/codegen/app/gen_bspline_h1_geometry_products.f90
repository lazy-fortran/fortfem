program gen_bspline_h1_geometry_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(+), operator(-), operator(*), &
        operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: geometry(4), geometry_dot(4), output(1)
    type(expr_t) :: output_bar(1), jvp_roots(1), vjp_roots(4)
    type(expr_t) :: row_reference(2), column_reference(2)
    type(expr_t) :: row_physical(2), column_physical(2)
    type(expr_t) :: determinant, basis_row, basis_column
    type(expr_t) :: stiffness, mass, quadrature_weight
    type(kernel_spec_t) :: jvp_spec, vjp_spec
    character(:), allocatable :: jvp_code, vjp_code
    integer :: ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    geometry = [ &
        sym(arena, "jacobian_11"), sym(arena, "jacobian_21"), &
        sym(arena, "jacobian_12"), sym(arena, "jacobian_22")]
    geometry_dot = [ &
        sym(arena, "jacobian_11_dot"), sym(arena, "jacobian_21_dot"), &
        sym(arena, "jacobian_12_dot"), sym(arena, "jacobian_22_dot")]
    row_reference = [ &
        sym(arena, "gradient_row_1"), sym(arena, "gradient_row_2")]
    column_reference = [ &
        sym(arena, "gradient_column_1"), &
        sym(arena, "gradient_column_2")]
    basis_row = sym(arena, "basis_row")
    basis_column = sym(arena, "basis_column")
    stiffness = sym(arena, "stiffness_coefficient")
    mass = sym(arena, "mass_coefficient")
    quadrature_weight = sym(arena, "quadrature_weight")
    output_bar(1) = sym(arena, "contribution_bar")

    determinant = geometry(1)*geometry(4) - geometry(3)*geometry(2)
    row_physical(1) = ( &
        geometry(4)*row_reference(1) - &
        geometry(2)*row_reference(2))/determinant
    row_physical(2) = ( &
        -geometry(3)*row_reference(1) + &
        geometry(1)*row_reference(2))/determinant
    column_physical(1) = ( &
        geometry(4)*column_reference(1) - &
        geometry(2)*column_reference(2))/determinant
    column_physical(2) = ( &
        -geometry(3)*column_reference(1) + &
        geometry(1)*column_reference(2))/determinant
    output(1) = quadrature_weight*determinant*( &
        stiffness*(row_physical(1)*column_physical(1) + &
        row_physical(2)*column_physical(2)) + &
        mass*basis_row*basis_column)
    jvp_roots = jvp(output, geometry, geometry_dot)
    vjp_roots = vjp(output, geometry, output_bar)
    call simplify_all(jvp_roots)
    call simplify_all(vjp_roots)

    jvp_spec%name = str("generated_bspline_h1_geometry_jvp")
    jvp_spec%module_name = str("fortfem_generated_bspline_h1_geometry_jvp")
    jvp_spec%mode = KERNEL_SUBROUTINE
    jvp_spec%generator = str("gen_bspline_h1_geometry_products")
    jvp_spec%generator_revision = str(fortsym_revision())
    jvp_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    jvp_spec%pure_procedure = .true.
    allocate (jvp_spec%args(17), jvp_spec%outputs(1))
    jvp_spec%args = [ &
        str("jacobian_11"), str("jacobian_21"), str("jacobian_12"), &
        str("jacobian_22"), str("gradient_row_1"), &
        str("gradient_row_2"), str("gradient_column_1"), &
        str("gradient_column_2"), str("basis_row"), str("basis_column"), &
        str("stiffness_coefficient"), str("mass_coefficient"), &
        str("quadrature_weight"), str("jacobian_11_dot"), &
        str("jacobian_21_dot"), str("jacobian_12_dot"), &
        str("jacobian_22_dot")]
    jvp_spec%outputs = [str("contribution_dot")]

    vjp_spec%name = str("generated_bspline_h1_geometry_vjp")
    vjp_spec%module_name = str("fortfem_generated_bspline_h1_geometry_vjp")
    vjp_spec%mode = KERNEL_SUBROUTINE
    vjp_spec%generator = str("gen_bspline_h1_geometry_products")
    vjp_spec%generator_revision = str(fortsym_revision())
    vjp_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    vjp_spec%pure_procedure = .true.
    allocate (vjp_spec%args(14), vjp_spec%outputs(4))
    vjp_spec%args = [ &
        str("jacobian_11"), str("jacobian_21"), str("jacobian_12"), &
        str("jacobian_22"), str("gradient_row_1"), &
        str("gradient_row_2"), str("gradient_column_1"), &
        str("gradient_column_2"), str("basis_row"), str("basis_column"), &
        str("stiffness_coefficient"), str("mass_coefficient"), &
        str("quadrature_weight"), str("contribution_bar")]
    vjp_spec%outputs = [ &
        str("jacobian_11_bar"), str("jacobian_21_bar"), &
        str("jacobian_12_bar"), str("jacobian_22_bar")]

    jvp_code = chars(emit_kernel(jvp_roots, jvp_spec))
    vjp_code = chars(emit_kernel(vjp_roots, vjp_spec))
    open (newunit=unit, &
        file=generated_path("fortfem_bspline_h1_geometry_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write B-spline H1 geometry products"
    write (unit, "(a)") jvp_code(:len(jvp_code) - 1)
    write (unit, "(a)") vjp_code(:len(vjp_code) - 1)
    close (unit)

contains

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (result%ok) expressions(root) = result%value
        end do
    end subroutine simplify_all

end program gen_bspline_h1_geometry_products
