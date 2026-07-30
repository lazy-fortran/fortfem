program gen_planar_helmholtz_dtn_products
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(-), operator(**), sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: result
    type(expr_t) :: variables(2), variables_dot(2), beta_bar(1)
    type(expr_t) :: propagating(1), evanescent(1)
    type(expr_t) :: propagating_jvp(1), propagating_vjp(2)
    type(expr_t) :: evanescent_jvp(1), evanescent_vjp(2)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: propagating_jvp_code, propagating_vjp_code
    character(:), allocatable :: evanescent_jvp_code, evanescent_vjp_code
    integer :: ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    variables = [sym(arena, "wavenumber"), sym(arena, "tangential_wavenumber")]
    variables_dot = [ &
        sym(arena, "wavenumber_dot"), sym(arena, "tangential_wavenumber_dot")]
    beta_bar(1) = sym(arena, "beta_bar")
    propagating(1) = sqrt(variables(1)**2 - variables(2)**2)
    evanescent(1) = sqrt(variables(2)**2 - variables(1)**2)
    propagating_jvp = jvp(propagating, variables, variables_dot)
    propagating_vjp = vjp(propagating, variables, beta_bar)
    evanescent_jvp = jvp(evanescent, variables, variables_dot)
    evanescent_vjp = vjp(evanescent, variables, beta_bar)
    call simplify_all(propagating_jvp)
    call simplify_all(propagating_vjp)
    call simplify_all(evanescent_jvp)
    call simplify_all(evanescent_vjp)

    call initialize_spec( &
        spec, "generated_planar_helmholtz_propagating_jvp", &
        "fortfem_generated_planar_helmholtz_propagating_jvp", 4, 1)
    spec%args = [ &
        str("wavenumber"), str("tangential_wavenumber"), &
        str("wavenumber_dot"), str("tangential_wavenumber_dot")]
    spec%outputs = [str("beta_dot")]
    propagating_jvp_code = chars(emit_kernel(propagating_jvp, spec))

    call initialize_spec( &
        spec, "generated_planar_helmholtz_propagating_vjp", &
        "fortfem_generated_planar_helmholtz_propagating_vjp", 3, 2)
    spec%args = [ &
        str("wavenumber"), str("tangential_wavenumber"), str("beta_bar")]
    spec%outputs = [str("wavenumber_bar"), str("tangential_wavenumber_bar")]
    propagating_vjp_code = chars(emit_kernel(propagating_vjp, spec))

    call initialize_spec( &
        spec, "generated_planar_helmholtz_evanescent_jvp", &
        "fortfem_generated_planar_helmholtz_evanescent_jvp", 4, 1)
    spec%args = [ &
        str("wavenumber"), str("tangential_wavenumber"), &
        str("wavenumber_dot"), str("tangential_wavenumber_dot")]
    spec%outputs = [str("beta_dot")]
    evanescent_jvp_code = chars(emit_kernel(evanescent_jvp, spec))

    call initialize_spec( &
        spec, "generated_planar_helmholtz_evanescent_vjp", &
        "fortfem_generated_planar_helmholtz_evanescent_vjp", 3, 2)
    spec%args = [ &
        str("wavenumber"), str("tangential_wavenumber"), str("beta_bar")]
    spec%outputs = [str("wavenumber_bar"), str("tangential_wavenumber_bar")]
    evanescent_vjp_code = chars(emit_kernel(evanescent_vjp, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_planar_helmholtz_dtn_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write planar Helmholtz DtN products"
    write (unit, "(a)") &
        propagating_jvp_code(:len(propagating_jvp_code) - 1)
    write (unit, "(a)") &
        propagating_vjp_code(:len(propagating_vjp_code) - 1)
    write (unit, "(a)") evanescent_jvp_code(:len(evanescent_jvp_code) - 1)
    write (unit, "(a)") evanescent_vjp_code(:len(evanescent_vjp_code) - 1)
    close (unit)

contains

    subroutine initialize_spec( &
            kernel_spec, name, module_name, argument_count, output_count)
        type(kernel_spec_t), intent(out) :: kernel_spec
        character(*), intent(in) :: name, module_name
        integer, intent(in) :: argument_count, output_count

        kernel_spec%name = str(name)
        kernel_spec%module_name = str(module_name)
        kernel_spec%mode = KERNEL_SUBROUTINE
        kernel_spec%generator = str("gen_planar_helmholtz_dtn_products")
        kernel_spec%generator_revision = str(fortsym_revision())
        kernel_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        kernel_spec%pure_procedure = .true.
        allocate ( &
            kernel_spec%args(argument_count), &
            kernel_spec%outputs(output_count))
    end subroutine initialize_spec

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        do root = 1, size(expressions)
            result = engine%simplify(expressions(root))
            if (result%ok) expressions(root) = result%value
        end do
    end subroutine simplify_all

end program gen_planar_helmholtz_dtn_products
