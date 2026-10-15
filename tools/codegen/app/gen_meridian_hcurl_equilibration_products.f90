program gen_meridian_hcurl_equilibration_products
    !! Pointwise residual of the patchwise equilibrated H(curl_n) bound
    !! (src/estimators/equilibrated_hcurl_meridian_n.f90).  With the
    !! equilibrated field H = (H_pol, s), s continuous, and the Ampere law
    !! curl_n^* H = J in the tested (gauge a = 0) sense,
    !!
    !!   H_pol = (1/(i n)) (grad(R s) + R rot90 J),   rot90 v = (-v_Z, v_R),
    !!
    !! the patch residual against the discrete field B_h = (Bp, sh) is
    !!
    !!   r_pol = (1/(i n)) (grad(R s) + psi R rot90 J - R sh grad psi) - psi Bp
    !!   r_phi = s - psi sh.
    !!
    !! psi = 1 gives the global residual H - B_h; psi a vertex hat function
    !! gives the patch residual whose sum over patches is the global one
    !! (sum psi = 1, sum grad psi = 0).  Complex quantities are split into
    !! real and imaginary parts.  R is the metric radius with gradient
    !! (grr1, grr2).
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, operator(*), operator(+), operator(-), &
        operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    integer, parameter :: nin = 23, nout = 6
    character(len=*), parameter :: in_names(nin) = [character(len=5) :: &
        "n", "rr", "grr1", "grr2", "psi", "gpsi1", "gpsi2", "sre", "sim", &
        "gsre1", "gsim1", "gsre2", "gsim2", "jre1", "jim1", "jre2", "jim2", &
        "shre", "shim", "bre1", "bim1", "bre2", "bim2"]
    character(len=*), parameter :: out_names(nout) = [character(len=5) :: &
        "r1re", "r1im", "r2re", "r2im", "rphre", "rphim"]

    type(arena_t), target :: arena
    type(engine_result_t) :: result
    type(native_engine_t) :: engine
    type(expr_t) :: v(nin), outputs(nout)
    type(expr_t) :: x1re, x1im, x2re, x2im
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code
    integer :: ios, root, unit, i

    call arena%init()
    engine = make_native_engine(arena)
    do i = 1, nin
        v(i) = sym(arena, trim(in_names(i)))
    end do

    ! X = grad(R s) + psi R rot90 J - R sh grad psi (complex 2-vector)
    ! rot90 J = (-J_Z, J_R); s = (sre, sim), J = (jre1 + i jim1, jre2 + i jim2)
    x1re = v(3)*v(8) + v(2)*v(10) - v(5)*v(2)*v(16) &
        - v(2)*(v(18)*v(6))
    x1im = v(3)*v(9) + v(2)*v(11) - v(5)*v(2)*v(17) &
        - v(2)*(v(19)*v(6))
    x2re = v(4)*v(8) + v(2)*v(12) + v(5)*v(2)*v(14) &
        - v(2)*(v(18)*v(7))
    x2im = v(4)*v(9) + v(2)*v(13) + v(5)*v(2)*v(15) &
        - v(2)*(v(19)*v(7))
    ! r_pol = (-i/n) X - psi Bp: (-i)(a + i b) = b - i a
    outputs(1) = x1im/v(1) - v(5)*v(20)
    outputs(2) = -x1re/v(1) - v(5)*v(21)
    outputs(3) = x2im/v(1) - v(5)*v(22)
    outputs(4) = -x2re/v(1) - v(5)*v(23)
    outputs(5) = v(8) - v(5)*v(18)
    outputs(6) = v(9) - v(5)*v(19)
    do root = 1, nout
        result = engine%simplify(outputs(root))
        if (result%ok) outputs(root) = result%value
    end do

    spec%name = str("generated_meridian_hcurl_residual")
    spec%module_name = str("fortfem_generated_meridian_hcurl_equilibration")
    spec%mode = KERNEL_SUBROUTINE
    spec%generator = str("gen_meridian_hcurl_equilibration_products")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    spec%pure_procedure = .true.
    allocate (spec%args(nin), spec%outputs(nout))
    do i = 1, nin
        spec%args(i) = str(trim(in_names(i)))
    end do
    do i = 1, nout
        spec%outputs(i) = str(trim(out_names(i)))
    end do
    code = chars(emit_kernel(outputs, spec))

    open (newunit=unit, &
        file=generated_path("fortfem_meridian_hcurl_equilibration_products.f90"), &
        status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write meridian H(curl) products"
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)
end program gen_meridian_hcurl_equilibration_products
