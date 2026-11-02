program gen_elastic_acoustic_coupling_products
    !! Pointwise leaves of the coupled elasticity-acoustic solvers:
    !! - isotropic material: Lame parameters, inertia rho*omega**2 and the P and S
    !!   impedances sqrt((lambda + 2 mu) rho), sqrt(mu rho) from (E, nu, rho, omega);
    !! - segment length of an edge from its coordinate differences;
    !! - one row of the periodic P1 mass on equispaced samples, h/6 (1, 4, 1), and
    !!   the incident-pressure load -n h/6 (1, 4, 1) p it gives at an interface node;
    !! - the interface block n_a n_b^T F_ab that contracts a scalar acoustic form
    !!   entry with the normals of its two nodes, and the same contraction
    !!   (n_a·e_p)(n_b·e_q) F_ab of one normal-component pair, which lets the
    !!    solver store the coupling on the active components only;
    !! - the panel coupling of the curved NtD form: length, normal trace n/2,
    !!   endpoint mass length/2 and the P1 traction mass n_c length (2 or 1)/6.
    !! JVP and VJP are emitted for the leaves that the solver derivatives use.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, sqrt, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_products, only: jvp, vjp
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(kernel_spec_t) :: spec
    integer :: unit

    call arena%init()
    engine = make_native_engine(arena)
    open (newunit=unit, file=generated_path( &
        "fortfem_elastic_acoustic_coupling_products.f90"), status="replace", &
        action="write")
    call material()
    call segment_length()
    call periodic_mass()
    call incident_load()
    call normal_block()
    call normal_component_block()
    call panel_coupling()
    close (unit)

contains

    subroutine material()
        type(expr_t) :: e, nu, rho, omega, lambda, mu, roots(5)

        e = sym(arena, "young_modulus")
        nu = sym(arena, "poisson_ratio")
        rho = sym(arena, "solid_density")
        omega = sym(arena, "angular_frequency")
        mu = e/(2*(1 + nu))
        lambda = e*nu/((1 + nu)*(1 - 2*nu))
        roots = [lambda, mu, rho*omega**2, sqrt((lambda + 2*mu)*rho), &
            sqrt(mu*rho)]
        call emit_family("isotropic_elastic_material", roots, [e, nu, rho, omega], &
            [str("young_modulus"), str("poisson_ratio"), str("solid_density"), &
            str("angular_frequency")], [str("complex(dp)"), real_types(3)], &
            [.true., .true., .true., .true.], &
            [str("lambda"), str("mu"), str("inertia"), str("p_impedance"), &
            str("s_impedance")], no_shapes(5), [1, 2, 3, 4, 5], no_shapes(5), &
            "complex(dp)", .true.)
    end subroutine material

    subroutine segment_length()
        type(expr_t) :: dx, dy

        dx = sym(arena, "dx")
        dy = sym(arena, "dy")
        call emit_family("segment_length", [sqrt(dx**2 + dy**2)], [dx, dy], &
            [str("dx"), str("dy")], real_types(2), [.true., .true.], &
            [str("length")], no_shapes(1), [1], no_shapes(1), "real(dp)", .true.)
    end subroutine segment_length

    subroutine periodic_mass()
        type(expr_t) :: period, count, previous, centre, next

        period = sym(arena, "period")
        count = sym(arena, "sample_count")
        previous = sym(arena, "previous")
        centre = sym(arena, "centre")
        next = sym(arena, "next")
        call emit_family("periodic_p1_mass", &
            [period/count/6*(previous + 4*centre + next)], &
            [period, count, previous, centre, next], &
            [str("period"), str("sample_count"), str("previous"), str("centre"), &
            str("next")], [real_types(2), complex_types(3)], &
            [.true., .false., .true., .true., .true.], [str("value")], &
            no_shapes(1), [1], no_shapes(1), "complex(dp)", .true.)
    end subroutine periodic_mass

    subroutine incident_load()
        type(expr_t) :: period, count, previous, centre, next, normal(2), weak

        period = sym(arena, "period")
        count = sym(arena, "sample_count")
        previous = sym(arena, "previous")
        centre = sym(arena, "centre")
        next = sym(arena, "next")
        normal = [sym(arena, "normal_1"), sym(arena, "normal_2")]
        weak = period/count/6*(previous + 4*centre + next)
        ! The incident pressure pushes against the outward normal.
        call emit_family("planar_incident_load", [-normal(1)*weak, &
            -normal(2)*weak], [period, count, previous, centre, next, normal], &
            [str("period"), str("sample_count"), str("previous"), str("centre"), &
            str("next"), str("normal_1"), str("normal_2")], &
            [real_types(2), complex_types(3), real_types(2)], &
            [.true., .false., .true., .true., .true., .false., .false.], &
            [str("load")], [str("(2)")], [1, 1], [str("(1)"), str("(2)")], &
            "complex(dp)", .true.)
    end subroutine incident_load

    subroutine normal_block()
        type(expr_t) :: a(2), b(2), value, roots(4)
        type(str_t) :: refs(4)
        integer :: p, q

        a = [sym(arena, "normal_a_1"), sym(arena, "normal_a_2")]
        b = [sym(arena, "normal_b_1"), sym(arena, "normal_b_2")]
        value = sym(arena, "value")
        do q = 1, 2
            do p = 1, 2
                roots(p + 2*(q - 1)) = a(p)*value*b(q)
                refs(p + 2*(q - 1)) = str(index_text([p, q]))
            end do
        end do
        call emit_family("interface_normal_block", roots, [a, b, value], &
            [str("normal_a_1"), str("normal_a_2"), str("normal_b_1"), &
            str("normal_b_2"), str("value")], [real_types(4), complex_types(1)], &
            [.false., .false., .false., .false., .true.], [str("block")], &
            [str("(2,2)")], [1, 1, 1, 1], refs, "complex(dp)", .true.)
    end subroutine normal_block

    subroutine normal_component_block()
        !! One component (p, q) of the interface block: the solver assembles
        !! the coupling on the normal-component dofs only, so it needs the
        !! scalar (n_a·e_p)(n_b·e_q) F_ab as its own kernel.
        type(expr_t) :: a, b, value

        a = sym(arena, "normal_a")
        b = sym(arena, "normal_b")
        value = sym(arena, "value")
        call emit_family("interface_normal_component_block", [a*value*b], &
            [a, b, value], [str("normal_a"), str("normal_b"), str("value")], &
            [real_types(2), complex_types(1)], [.false., .false., .true.], &
            [str("block_value")], no_shapes(1), [1], no_shapes(1), &
            "complex(dp)", .true.)
    end subroutine normal_component_block

    subroutine panel_coupling()
        type(expr_t) :: dx, dy, length, normal(2), roots(12)
        type(str_t) :: refs(12)
        integer :: c, first, second, k, owners(12)

        dx = sym(arena, "dx")
        dy = sym(arena, "dy")
        length = sqrt(dx**2 + dy**2)
        ! Right normal of a counter-clockwise panel: out of the bounded domain.
        normal = [dy/length, -dx/length]
        roots(1) = length
        roots(2) = length/2
        roots(3) = normal(1)/2
        roots(4) = normal(2)/2
        owners(1:4) = [1, 2, 3, 3]
        refs(1:2) = str("")
        refs(3:4) = [str("(1)"), str("(2)")]
        k = 4
        do second = 1, 2
            do first = 1, 2
                do c = 1, 2
                    k = k + 1
                    if (first == second) then
                        roots(k) = normal(c)*length/3
                    else
                        roots(k) = normal(c)*length/6
                    end if
                    owners(k) = 4
                    refs(k) = str(index_text([c, first, second]))
                end do
            end do
        end do
        ! Primal only: the BEM operators have no geometric derivative yet.
        call emit_family("ntd_panel_coupling", roots, [dx, dy], &
            [str("dx"), str("dy")], real_types(2), [.true., .true.], &
            [str("length"), str("endpoint_mass"), str("trace"), str("traction")], &
            [str(""), str(""), str("(2)"), str("(2,2,2)")], owners, refs, &
            "real(dp)", .false.)
    end subroutine panel_coupling

    subroutine emit_family(family, roots, inputs, names, types, active, &
            outputs, output_shapes, owners, suffixes, scalar_type, derivatives)
        !! Primal kernel and, if requested, the JVP and VJP with respect to the
        !! active inputs. Root k is written to outputs(owners(k))//suffixes(k).
        !! The VJP is the plain transpose (no conjugation), as in the other
        !! complex fortfem kernels.
        character(*), intent(in) :: family, scalar_type
        type(expr_t), intent(in) :: roots(:), inputs(:)
        type(str_t), intent(in) :: names(:), types(:), outputs(:)
        type(str_t), intent(in) :: output_shapes(:), suffixes(:)
        logical, intent(in) :: active(:), derivatives
        integer, intent(in) :: owners(:)

        type(expr_t), allocatable :: dots(:), bars(:), tangent(:), cotangent(:)
        type(expr_t), allocatable :: variables(:), work(:)
        type(str_t), allocatable :: active_names(:), active_types(:)
        type(str_t), allocatable :: dot_names(:), bar_names(:)
        integer :: k

        variables = pack(inputs, active)
        active_names = pack(names, active)
        active_types = pack(types, active)
        call initialize(family, scalar_type)
        spec%args = names
        spec%arg_types = types
        call set_outputs(outputs, output_shapes, owners, suffixes, "")
        work = roots
        call simplify(work)
        call write_kernel(work)
        if (.not. derivatives) return

        allocate(dots(size(variables)), dot_names(size(variables)))
        do k = 1, size(variables)
            dot_names(k) = str(chars(active_names(k))//"_dot")
            dots(k) = sym(arena, chars(dot_names(k)))
        end do
        tangent = jvp(roots, variables, dots)
        call simplify(tangent)
        call initialize(family//"_jvp", scalar_type)
        spec%args = [names, dot_names]
        spec%arg_types = [types, active_types]
        call set_outputs(outputs, output_shapes, owners, suffixes, "_dot")
        call write_kernel(tangent)

        allocate(bars(size(roots)))
        do k = 1, size(roots)
            bars(k) = sym(arena, chars(outputs(owners(k)))//"_bar"// &
                chars(suffixes(k)))
        end do
        cotangent = vjp(roots, variables, bars)
        call simplify(cotangent)
        allocate(bar_names(size(outputs)))
        do k = 1, size(outputs)
            bar_names(k) = str(chars(outputs(k))//"_bar")
        end do
        call initialize(family//"_vjp", scalar_type)
        spec%args = [names, bar_names]
        spec%arg_types = [types, (str(scalar_type), k=1, size(outputs))]
        spec%arg_shapes = [no_shapes(size(names)), output_shapes]
        spec%outputs = [(str(chars(active_names(k))//"_bar"), &
            k=1, size(active_names))]
        call write_kernel(cotangent)
    end subroutine emit_family

    subroutine set_outputs(outputs, output_shapes, owners, suffixes, tag)
        type(str_t), intent(in) :: outputs(:), output_shapes(:), suffixes(:)
        integer, intent(in) :: owners(:)
        character(*), intent(in) :: tag
        integer :: k

        spec%outputs = [(str(chars(outputs(k))//tag), k=1, size(outputs))]
        spec%output_shapes = output_shapes
        allocate(spec%output_references(size(owners)))
        do k = 1, size(owners)
            spec%output_references(k) = str(chars(outputs(owners(k)))//tag// &
                chars(suffixes(k)))
        end do
    end subroutine set_outputs

    function index_text(indices) result(text)
        integer, intent(in) :: indices(:)
        character(:), allocatable :: text
        character(8) :: part
        integer :: k

        text = "("
        do k = 1, size(indices)
            write (part, "(i0)") indices(k)
            text = text//trim(part)
            if (k < size(indices)) text = text//","
        end do
        text = text//")"
    end function index_text

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

    subroutine initialize(family, scalar_type)
        character(*), intent(in) :: family, scalar_type

        spec = kernel_spec_t()
        spec%name = str("generated_"//family)
        spec%module_name = str("fortfem_generated_"//family)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_elastic_acoustic_coupling_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%scalar_type = str(scalar_type)
        spec%pure_procedure = .true.
    end subroutine initialize

    subroutine write_kernel(roots)
        type(expr_t), intent(in) :: roots(:)
        character(:), allocatable :: code, message
        logical :: ok

        code = chars(emit_kernel(roots, spec, ok=ok, message=message))
        if (.not. ok) then
            write (*, "(a)") chars(spec%name)//": "//message
            error stop "elastic-acoustic kernel emission failed"
        end if
        write (unit, "(a)") code(:len(code) - 1)
    end subroutine write_kernel

    subroutine simplify(roots)
        type(expr_t), intent(inout) :: roots(:)
        type(engine_result_t) :: result
        integer :: root

        do root = 1, size(roots)
            result = engine%simplify(roots(root))
            if (.not. result%ok) error stop "elastic-acoustic simplification failed"
            roots(root) = result%value
        end do
    end subroutine simplify

end program gen_elastic_acoustic_coupling_products
