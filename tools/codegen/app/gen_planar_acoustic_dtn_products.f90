program gen_planar_acoustic_dtn_products
    !! Mode multiplier of the planar acoustic displacement-to-pressure map,
    !! -i rho omega**2/k_n with k_n = sqrt((omega/c)**2 - (2 pi m/period)**2)
    !! (principal root, time convention exp(-i omega t)), with JVP and VJP in
    !! (omega, c, rho, period). The radicand is a separate primal leaf so the
    !! caller can reject a mode at cutoff before dividing by k_n.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, i_expr, pi_expr, sqrt, sym, operator(+), &
        operator(-), operator(*), operator(/), operator(**)
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(kernel_spec_t) :: spec
    type(expr_t) :: omega, c, rho, period, mode, radicand, multiplier(1)
    type(expr_t) :: inputs(4), dots(4), bar(1), tangent(1), cotangent(4)
    type(expr_t) :: roots(1)
    character(*), parameter :: names(4) = [character(17) :: "angular_frequency", &
        "sound_speed", "fluid_density", "period"]
    integer :: k, unit

    call arena%init()
    engine = make_native_engine(arena)
    omega = sym(arena, "angular_frequency")
    c = sym(arena, "sound_speed")
    rho = sym(arena, "fluid_density")
    period = sym(arena, "period")
    mode = sym(arena, "mode")
    radicand = (omega/c)**2 - (2*pi_expr(arena)*mode/period)**2
    multiplier(1) = -i_expr(arena)*rho*omega**2/sqrt(radicand)
    inputs = [omega, c, rho, period]
    do k = 1, 4
        dots(k) = sym(arena, trim(names(k))//"_dot")
    end do
    bar(1) = sym(arena, "multiplier_bar")
    tangent = jvp(multiplier, inputs, dots)
    cotangent = vjp(multiplier, inputs, bar)

    open (newunit=unit, file=generated_path( &
        "fortfem_planar_acoustic_dtn_products.f90"), &
        status="replace", action="write")
    call initialize("planar_acoustic_dtn_radicand")
    spec%args = [str("angular_frequency"), str("sound_speed"), str("period"), &
        str("mode")]
    spec%arg_types = [str("real(dp)"), str("complex(dp)"), str("real(dp)"), &
        str("real(dp)")]
    spec%outputs = [str("radicand")]
    roots(1) = radicand
    call write_kernel(roots)

    call initialize("planar_acoustic_dtn_multiplier")
    call set_args(.false.)
    spec%outputs = [str("multiplier")]
    call write_kernel(multiplier)

    call initialize("planar_acoustic_dtn_multiplier_jvp")
    call set_args(.true.)
    spec%outputs = [str("multiplier_dot")]
    call write_kernel(tangent)

    call initialize("planar_acoustic_dtn_multiplier_vjp")
    call set_args(.false.)
    spec%args = [spec%args, str("multiplier_bar")]
    spec%arg_types = [spec%arg_types, str("complex(dp)")]
    spec%outputs = [(str(trim(names(k))//"_bar"), k=1, 4)]
    call write_kernel(cotangent)
    close (unit)

contains

    subroutine set_args(tangent_args)
        logical, intent(in) :: tangent_args
        type(str_t) :: types(5)

        types = [str("real(dp)"), str("complex(dp)"), str("real(dp)"), &
            str("real(dp)"), str("real(dp)")]
        spec%args = [(str(trim(names(k))), k=1, 4), str("mode")]
        spec%arg_types = types
        if (tangent_args) then
            spec%args = [spec%args, (str(trim(names(k))//"_dot"), k=1, 4)]
            spec%arg_types = [types, types(1:4)]
        end if
    end subroutine set_args

    subroutine initialize(family)
        character(*), intent(in) :: family

        spec = kernel_spec_t()
        spec%name = str("generated_"//family)
        spec%module_name = str("fortfem_generated_"//family)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_planar_acoustic_dtn_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%scalar_type = str("complex(dp)")
        spec%pure_procedure = .true.
    end subroutine initialize

    subroutine write_kernel(expressions)
        type(expr_t), intent(in) :: expressions(:)
        type(expr_t), allocatable :: work(:)
        type(engine_result_t) :: result
        character(:), allocatable :: code, message
        logical :: ok
        integer :: root

        work = expressions
        do root = 1, size(work)
            result = engine%simplify(work(root))
            if (.not. result%ok) error stop "planar DtN simplification failed"
            work(root) = result%value
        end do
        code = chars(emit_kernel(work, spec, ok=ok, message=message))
        if (.not. ok) then
            write (*, "(a)") chars(spec%name)//": "//message
            error stop "planar DtN kernel emission failed"
        end if
        write (unit, "(a)") code(:len(code) - 1)
    end subroutine write_kernel

end program gen_planar_acoustic_dtn_products
