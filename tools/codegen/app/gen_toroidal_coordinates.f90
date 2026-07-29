program gen_toroidal_coordinates
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    use fortsym_toroidal, only: make_toroidal_chart
    use fortsym_chart, only: chart_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(chart_t) :: chart
    type(expr_t) :: eta, phi, scale, theta
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code
    integer :: ios, unit

    call arena%init()
    eta = sym(arena, "eta")
    theta = sym(arena, "theta")
    phi = sym(arena, "phi")
    scale = sym(arena, "scale")
    chart = make_toroidal_chart(arena, eta, theta, phi, scale)

    spec%name = str("generated_toroidal_point_to_cartesian")
    spec%module_name = str("fortfem_generated_toroidal_coordinates")
    spec%mode = KERNEL_SUBROUTINE
    spec%generator = str("gen_toroidal_coordinates")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    spec%pure_procedure = .true.
    allocate(spec%args(4), spec%outputs(1), spec%output_shapes(1))
    allocate(spec%output_references(3))
    spec%args = [str("scale"), str("eta"), str("theta"), str("phi")]
    spec%outputs = [str("point")]
    spec%output_shapes = [str("(3)")]
    spec%output_references = [str("point(1)"), str("point(2)"), str("point(3)")]

    code = chars(emit_kernel(chart%x, spec))
    open (newunit=unit, &
        file=generated_path("fortfem_toroidal_coordinates.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write toroidal coordinate kernel"
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)
end program gen_toroidal_coordinates
