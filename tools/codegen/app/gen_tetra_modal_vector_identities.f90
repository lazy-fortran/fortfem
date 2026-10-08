program gen_tetra_modal_vector_identities
    use fortsym_arena, only: arena_t
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_products, only: jvp
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, operator(*), operator(+), &
        operator(-), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: component_curls(3, 3), cross_curls(3, 3)
    type(expr_t) :: cross_curls_dot(3, 3), cross_values(3, 3)
    type(expr_t) :: cross_values_dot(3, 3), roots(27), roots_dot(18)
    type(expr_t) :: dx, dy, dz, phi, x, y, z, zero
    type(expr_t) :: dx_dot, dy_dot, dz_dot, phi_dot
    type(expr_t) :: x_dot, y_dot, z_dot
    type(expr_t) :: shift(3), shift_zero(3), local_phi, local_field(3), local_curl(3)
    type(expr_t) :: component_roots_dot(9), radial_roots(4), radial_dot(4)
    type(expr_t) :: local_divergence
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code, filename
    integer :: column, component, ios, root, unit

    call arena%init()
    engine = make_native_engine(arena)
    x = sym(arena, "x")
    y = sym(arena, "y")
    z = sym(arena, "z")
    phi = sym(arena, "phi")
    dx = sym(arena, "dx")
    dy = sym(arena, "dy")
    dz = sym(arena, "dz")
    x_dot = sym(arena, "x_dot")
    y_dot = sym(arena, "y_dot")
    z_dot = sym(arena, "z_dot")
    phi_dot = sym(arena, "phi_dot")
    dx_dot = sym(arena, "dx_dot")
    dy_dot = sym(arena, "dy_dot")
    dz_dot = sym(arena, "dz_dot")
    zero = num(arena, 0)

    shift=[sym(arena,'shift_x'),sym(arena,'shift_y'),sym(arena,'shift_z')]
    shift_zero=zero
    local_phi=phi+dx*shift(1)+dy*shift(2)+dz*shift(3)
    do column=1,3
        local_field=zero
        local_field(column)=local_phi
        local_curl=spatial_curl(local_field)
        do component=1,3
            component_curls(component,column)=subs_many(local_curl(component),shift,shift_zero)
        end do
        select case(column)
        case(1)
            local_field=[-(y+shift(2))*local_phi,(x+shift(1))*local_phi,zero]
        case(2)
            local_field=[-(z+shift(3))*local_phi,zero,(x+shift(1))*local_phi]
        case(3)
            local_field=[zero,-(z+shift(3))*local_phi,(y+shift(2))*local_phi]
        end select
        local_curl=spatial_curl(local_field)
        do component=1,3
            cross_values(component,column)=subs_many(local_field(component),shift,shift_zero)
            cross_curls(component,column)=subs_many(local_curl(component),shift,shift_zero)
        end do
    end do

    root = 0
    call append_matrix(component_curls)
    call append_matrix(cross_values)
    call append_matrix(cross_curls)
    call simplify_all(roots)

    spec%name = str("evaluate_tetra_modal_vector_identities")
    spec%module_name = str("fortfem_generated_tetra_modal_vector_identities")
    spec%mode = KERNEL_SUBROUTINE
    spec%temp_prefix = str("t")
    spec%generator = str("gen_tetra_modal_vector_identities")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    spec%pure_procedure = .true.
    allocate(spec%args(7), spec%outputs(3), spec%output_shapes(3))
    allocate(spec%output_references(27))
    spec%args = [ &
        str("x"), str("y"), str("z"), str("phi"), &
        str("dx"), str("dy"), str("dz")]
    spec%outputs = [ &
        str("component_curls"), str("cross_values"), str("cross_curls")]
    spec%output_shapes = [str("(3,3)"), str("(3,3)"), str("(3,3)")]
    root = 0
    call append_references("component_curls")
    call append_references("cross_values")
    call append_references("cross_curls")

    filename = generated_path( &
        "fortfem_tetra_modal_vector_identities.f90")
    open( &
        newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write modal vector identity kernel"
    code = chars(emit_kernel(roots, spec))
    write(unit, "(a)") code(:len(code) - 1)
    close(unit)

    do column = 1, 3
        do component = 1, 3
            cross_values_dot(component, column) = &
                directional_derivative(cross_values(component, column))
            cross_curls_dot(component, column) = &
                directional_derivative(cross_curls(component, column))
        end do
    end do
    root = 0
    do column = 1, 3
        do component = 1, 3
            root = root + 1
            roots_dot(root) = cross_values_dot(component, column)
        end do
    end do
    do column = 1, 3
        do component = 1, 3
            root = root + 1
            roots_dot(root) = cross_curls_dot(component, column)
        end do
    end do
    call simplify_all(roots_dot)
    spec = kernel_spec_t()
    spec%name = str("evaluate_tetra_modal_vector_identities_jvp")
    spec%module_name = str( &
        "fortfem_generated_tetra_modal_vector_identities_jvp")
    spec%mode = KERNEL_SUBROUTINE
    spec%temp_prefix = str("t")
    spec%generator = str("gen_tetra_modal_vector_identities")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
    spec%pure_procedure = .true.
    allocate(spec%args(14), spec%outputs(2), spec%output_shapes(2))
    allocate(spec%output_references(18))
    spec%args = [ &
        str("x"), str("y"), str("z"), str("phi"), &
        str("dx"), str("dy"), str("dz"), str("x_dot"), str("y_dot"), &
        str("z_dot"), str("phi_dot"), str("dx_dot"), str("dy_dot"), &
        str("dz_dot")]
    spec%outputs = [str("cross_values_dot"), str("cross_curls_dot")]
    spec%output_shapes = [str("(3,3)"), str("(3,3)")]
    root = 0
    call append_references("cross_values_dot")
    call append_references("cross_curls_dot")
    filename = generated_path( &
        "fortfem_tetra_modal_vector_identities_jvp.f90")
    open( &
        newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write modal vector identity JVP kernel"
    code = chars(emit_kernel(roots_dot, spec))
    write(unit, "(a)") code(:len(code) - 1)
    close(unit)

    root = 0
    do column = 1, 3
        do component = 1, 3
            root = root + 1
            component_roots_dot(root) = directional_derivative( &
                component_curls(component, column))
        end do
    end do
    call write_product(component_roots_dot, &
        "tetra_modal_component_curls_jvp", &
        [str("dx_dot"), str("dy_dot"), str("dz_dot")], &
        [str("component_curls_dot")], [str("(3,3)")], &
        [str("component_curls_dot(1,1)"), str("component_curls_dot(2,1)"), &
         str("component_curls_dot(3,1)"), str("component_curls_dot(1,2)"), &
         str("component_curls_dot(2,2)"), str("component_curls_dot(3,2)"), &
         str("component_curls_dot(1,3)"), str("component_curls_dot(2,3)"), &
         str("component_curls_dot(3,3)")])

    ! Derive the radial field and its divergence from the same spatial jet.
    local_field = [(x + shift(1))*local_phi, (y + shift(2))*local_phi, &
        (z + shift(3))*local_phi]
    local_divergence = zero
    do component = 1, 3
        radial_roots(component) = subs_many(local_field(component), shift, shift_zero)
        local_divergence = local_divergence + &
            diff(local_field(component), shift(component))
    end do
    radial_roots(4) = subs_many(local_divergence, shift, shift_zero)
    call write_product(radial_roots, "tetra_modal_radial_products", &
        [str("x"), str("y"), str("z"), str("phi"), &
         str("dx"), str("dy"), str("dz")], &
        [str("radial_values"), str("radial_divergence")], &
        [str("(3)"), str("")], &
        [str("radial_values(1)"), str("radial_values(2)"), &
         str("radial_values(3)"), str("radial_divergence")])
    do component = 1, 4
        radial_dot(component) = directional_derivative(radial_roots(component))
    end do
    call write_product(radial_dot, "tetra_modal_radial_products_jvp", &
        [str("x"), str("y"), str("z"), str("phi"), &
         str("dx"), str("dy"), str("dz"), str("x_dot"), str("y_dot"), &
         str("z_dot"), str("phi_dot"), str("dx_dot"), str("dy_dot"), str("dz_dot")], &
        [str("radial_values_dot"), str("radial_divergence_dot")], &
        [str("(3)"), str("")], &
        [str("radial_values_dot(1)"), str("radial_values_dot(2)"), &
         str("radial_values_dot(3)"), str("radial_divergence_dot")])

contains
    subroutine write_product(expressions, name, args, outputs, shapes, references)
        use fortsym_string, only: str_t
        type(expr_t), intent(inout) :: expressions(:)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:), shapes(:), references(:)
        type(kernel_spec_t) :: product_spec
        character(:), allocatable :: product_code, product_filename, message
        logical :: ok
        integer :: product_unit, product_ios

        call simplify_all(expressions)
        product_spec%name = str("evaluate_"//name)
        product_spec%module_name = str("fortfem_generated_"//name)
        product_spec%mode = KERNEL_SUBROUTINE
        product_spec%temp_prefix = str("t")
        product_spec%generator = str("gen_tetra_modal_vector_identities")
        product_spec%generator_revision = str(fortsym_revision())
        product_spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        product_spec%pure_procedure = .true.
        product_spec%args = args
        product_spec%outputs = outputs
        product_spec%output_shapes = shapes
        product_spec%output_references = references
        product_code = chars(emit_kernel( &
            expressions, product_spec, ok=ok, message=message))
        if (.not. ok) then
            print *, message
            error stop "native modal differential product emission failed"
        end if
        product_filename = generated_path("fortfem_"//name//".f90")
        open(newunit=product_unit, file=product_filename, status="replace", &
            action="write", iostat=product_ios)
        if (product_ios /= 0) error stop "cannot write modal differential product"
        write(product_unit, "(a)") product_code(:len(product_code) - 1)
        close(product_unit)
    end subroutine write_product

    function spatial_curl(field) result(curl)
        type(expr_t),intent(in)::field(3)
        type(expr_t)::curl(3)
        curl=[diff(field(3),shift(2))-diff(field(2),shift(3)), &
            diff(field(1),shift(3))-diff(field(3),shift(1)), &
            diff(field(2),shift(1))-diff(field(1),shift(2))]
    end function spatial_curl

    subroutine append_matrix(matrix)
        type(expr_t), intent(in) :: matrix(3, 3)

        do column = 1, 3
            do component = 1, 3
                root = root + 1
                roots(root) = matrix(component, column)
            end do
        end do
    end subroutine append_matrix

    subroutine append_references(name)
        character(*), intent(in) :: name

        do column = 1, 3
            do component = 1, 3
                root = root + 1
                spec%output_references(root) = str( &
                    name//"("//integer_text(component)//","// &
                    integer_text(column)//")")
            end do
        end do
    end subroutine append_references

    function directional_derivative(expression) result(value)
        type(expr_t), intent(in) :: expression
        type(expr_t) :: value

        type(expr_t)::roots(1)
        roots=jvp([expression],[x,y,z,phi,dx,dy,dz], &
            [x_dot,y_dot,z_dot,phi_dot,dx_dot,dy_dot,dz_dot])
        value=roots(1)
    end function directional_derivative

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)

        type(engine_result_t) :: result
        integer :: expression

        do expression = 1, size(expressions)
            result = engine%simplify(expressions(expression))
            if (.not. result%ok) error stop "native modal simplification failed"
            expressions(expression) = result%value
        end do
    end subroutine simplify_all

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text
        character(16) :: buffer

        write(buffer, "(i0)") value
        text = trim(buffer)
    end function integer_text

end program gen_tetra_modal_vector_identities
