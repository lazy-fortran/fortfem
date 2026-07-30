program gen_magnetic_curvilinear_coefficients
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: &
        expr_t, operator(*), operator(-), operator(/), operator(**), sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: determinant, g11, g12, g22, g33
    type(expr_t) :: jacobian, reluctivity, roots(5)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code
    integer :: ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    g11 = sym(arena, "metric_11")
    g12 = sym(arena, "metric_12")
    g22 = sym(arena, "metric_22")
    g33 = sym(arena, "metric_33")
    reluctivity = sym(arena, "reluctivity")

    determinant = g11 * g22 - g12**2
    jacobian = sqrt(determinant * g33)
    roots(1) = reluctivity * g33 / jacobian
    roots(2) = reluctivity * jacobian * g22 / (g33 * determinant)
    roots(3) = -reluctivity * jacobian * g12 / (g33 * determinant)
    roots(4) = roots(3)
    roots(5) = reluctivity * jacobian * g11 / (g33 * determinant)
    do root = 1, size(roots)
        result = engine%simplify(roots(root))
        if (result%ok) roots(root) = result%value
    end do

    spec%name = str("generated_magnetic_curvilinear_coefficients_2d")
    spec%module_name = &
        str("fortfem_generated_magnetic_curvilinear_coefficients_2d")
    spec%mode = KERNEL_SUBROUTINE
    spec%generator = str("gen_magnetic_curvilinear_coefficients")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    spec%pure_procedure = .true.
    allocate(spec%args(5), spec%outputs(2), spec%output_shapes(2))
    allocate(spec%output_references(5))
    spec%args = [ &
        str("metric_11"), str("metric_12"), str("metric_22"), &
        str("metric_33"), str("reluctivity")]
    spec%outputs = [str("curl_weight"), str("mass_tensor")]
    spec%output_shapes = [str(""), str("(2,2)")]
    spec%output_references = [ &
        str("curl_weight"), str("mass_tensor(1,1)"), &
        str("mass_tensor(1,2)"), str("mass_tensor(2,1)"), &
        str("mass_tensor(2,2)")]

    code = chars(emit_kernel(roots, spec))
    open (newunit=unit, &
        file=generated_path( &
        "fortfem_magnetic_curvilinear_coefficients_2d.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write magnetic curvilinear kernel"
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)
end program gen_magnetic_curvilinear_coefficients
