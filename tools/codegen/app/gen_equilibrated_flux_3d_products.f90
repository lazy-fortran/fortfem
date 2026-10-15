program gen_equilibrated_flux_3d_products
    !! Pointwise leaves of the 3D patchwise equilibrated-flux estimator:
    !! tetrahedral geometry (Jacobian, determinant, hat gradients), monomial
    !! bases P_k for k = 0..4, squared edge lengths, and the patch load, Gram,
    !! flux right-hand side, multiplier, projection and defect integrands.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, operator(*), operator(+), operator(-), &
        operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE, &
        CSE_NONE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    integer :: unit

    call arena%init()
    engine = make_native_engine(arena)
    open (newunit=unit, file=generated_path("fortfem_equilibrated_flux_3d_products.f90"), &
        status="replace", action="write")
    call emit_geometry()
    call emit_monomials()
    call emit_edge_lengths()
    call emit_patch_load()
    call emit_patch_gram()
    call emit_patch_flux_rhs()
    call emit_patch_multiplier()
    call emit_projection_terms()
    call emit_flux_defect()
    call emit_divergence_defect()
    call emit_oscillation_defect()
    close (unit)
contains

    function sym_(name) result(value)
        character(*), intent(in) :: name
        type(expr_t) :: value

        value = sym(arena, name)
    end function sym_

    function at(name, i) result(text)
        character(*), intent(in) :: name
        integer, intent(in) :: i
        character(:), allocatable :: text

        text = name//"("//itext(i)//")"
    end function at

    function at2(name, i, j) result(text)
        character(*), intent(in) :: name
        integer, intent(in) :: i, j
        character(:), allocatable :: text

        text = name//"("//itext(i)//","//itext(j)//")"
    end function at2

    function itext(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text
        character(32) :: buffer

        write (buffer, "(i0)") value
        text = trim(buffer)
    end function itext

    subroutine cross(a, b, c)
        type(expr_t), intent(in) :: a(3), b(3)
        type(expr_t), intent(out) :: c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end subroutine cross

    subroutine write_kernel(roots, spec)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t), intent(inout) :: spec
        type(expr_t), allocatable :: simplified(:)
        type(engine_result_t) :: result
        character(:), allocatable :: code
        integer :: k

        allocate (simplified(size(roots)))
        do k = 1, size(roots)
            result = engine%simplify(roots(k))
            if (.not. result%ok) error stop "equilibrated leaf simplification failed"
            simplified(k) = result%value
        end do
        spec%generator = str("gen_equilibrated_flux_3d_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%mode = KERNEL_SUBROUTINE
        spec%cse_level = CSE_NONE
        spec%pure_procedure = .true.
        code = chars(emit_kernel(simplified, spec))
        write (unit, "(a)") code(:len(code) - 1)
    end subroutine write_kernel

    subroutine output_spec(spec, name, module_name, args, shapes, output_name, &
            output_shape, count)
        type(kernel_spec_t), intent(inout) :: spec
        character(*), intent(in) :: name, module_name, output_name, output_shape
        character(*), intent(in) :: args(:), shapes(:)
        integer, intent(in) :: count
        integer :: k

        if (allocated(spec%args)) deallocate (spec%args)
        if (allocated(spec%arg_shapes)) deallocate (spec%arg_shapes)
        if (allocated(spec%outputs)) deallocate (spec%outputs)
        if (allocated(spec%output_shapes)) deallocate (spec%output_shapes)
        if (allocated(spec%output_references)) deallocate (spec%output_references)
        spec%name = str(name)
        spec%module_name = str(module_name)
        allocate (spec%args(size(args)), spec%arg_shapes(size(args)))
        do k = 1, size(args)
            spec%args(k) = str(trim(args(k)))
            spec%arg_shapes(k) = str(trim(shapes(k)))
        end do
        allocate (spec%outputs(1), spec%output_shapes(1), spec%output_references(count))
        spec%outputs(1) = str(output_name)
        spec%output_shapes(1) = str(trim(output_shape))
        do k = 1, count
            spec%output_references(k) = str(output_name//"("//itext(k)//")")
        end do
    end subroutine output_spec

    subroutine emit_geometry()
        !! Outputs product(1:9) = Jacobian (column-major), product(10) =
        !! determinant, product(11:22) = hat gradients (3 x 4, column-major).
        type(kernel_spec_t) :: spec
        type(expr_t) :: vertices(3, 4), jac(3, 3), c1(3), c2(3), c3(3)
        type(expr_t) :: det, roots(22), grad(3, 4), inv_rows(3, 3), n1(3), n2(3), n3(3)
        integer :: r, c, s, k

        do c = 1, 4
            do r = 1, 3
                vertices(r, c) = sym_(at2("vertices", r, c))
            end do
        end do
        do c = 1, 3
            do r = 1, 3
                jac(r, c) = vertices(r, c + 1) - vertices(r, 1)
            end do
        end do
        c1 = jac(:, 1)
        c2 = jac(:, 2)
        c3 = jac(:, 3)
        call cross(c2, c3, n1)
        call cross(c3, c1, n2)
        call cross(c1, c2, n3)
        det = c1(1)*n1(1) + c1(2)*n1(2) + c1(3)*n1(3)
        do r = 1, 3
            inv_rows(r, 1) = n1(r)/det
            inv_rows(r, 2) = n2(r)/det
            inv_rows(r, 3) = n3(r)/det
        end do
        do s = 1, 4
            if (s == 1) then
                do r = 1, 3
                    grad(r, s) = num(arena, 0) - inv_rows(r, 1) &
                        - inv_rows(r, 2) - inv_rows(r, 3)
                end do
            else
                do r = 1, 3
                    grad(r, s) = inv_rows(r, s - 1)
                end do
            end if
        end do
        k = 0
        do c = 1, 3
            do r = 1, 3
                k = k + 1
                roots(k) = jac(r, c)
            end do
        end do
        roots(10) = det
        k = 10
        do s = 1, 4
            do r = 1, 3
                k = k + 1
                roots(k) = grad(r, s)
            end do
        end do
        call output_spec(spec, "generated_tetra_simplex_geometry", &
            "fortfem_generated_tetra_simplex_geometry", ["vertices"], ["(3,4)"], &
            "product", "(22)", 22)
        call write_kernel(roots, spec)
    end subroutine emit_geometry

    subroutine emit_monomials()
        !! Monomials x^i y^j z^l, i + j + l <= degree, ordered with i outermost.
        integer :: degree, i, j, l, n, nk, r
        type(kernel_spec_t) :: spec
        type(expr_t), allocatable :: roots(:)
        type(expr_t) :: point(3), term
        character(:), allocatable :: name

        do degree = 0, 4
            nk = (degree + 1)*(degree + 2)*(degree + 3)/6
            allocate (roots(nk))
            do r = 1, 3
                point(r) = sym_(at("point", r))
            end do
            n = 0
            do i = 0, degree
                do j = 0, degree - i
                    do l = 0, degree - i - j
                        n = n + 1
                        term = num(arena, 1)
                        term = term*power(point(1), i)*power(point(2), j) &
                            *power(point(3), l)
                        roots(n) = term
                    end do
                end do
            end do
            name = "tetra_monomials_degree"//itext(degree)
            call output_spec(spec, "generated_"//name, &
                "fortfem_generated_"//name, ["point"], ["(3)"], "product", &
                "("//itext(nk)//")", nk)
            call write_kernel(roots, spec)
            deallocate (roots)
        end do
    end subroutine emit_monomials

    recursive function power(base, exponent) result(value)
        type(expr_t), intent(in) :: base
        integer, intent(in) :: exponent
        type(expr_t) :: value

        if (exponent == 0) then
            value = num(arena, 1)
        else
            value = base*power(base, exponent - 1)
        end if
    end function power

    subroutine emit_edge_lengths()
        type(kernel_spec_t) :: spec
        type(expr_t) :: vertices(3, 4), roots(6), d
        integer :: r, c, a, b, k, pair(2, 6)

        pair(:, 1) = [1, 2]
        pair(:, 2) = [1, 3]
        pair(:, 3) = [1, 4]
        pair(:, 4) = [2, 3]
        pair(:, 5) = [2, 4]
        pair(:, 6) = [3, 4]
        do c = 1, 4
            do r = 1, 3
                vertices(r, c) = sym_(at2("vertices", r, c))
            end do
        end do
        do k = 1, 6
            a = pair(1, k)
            b = pair(2, k)
            d = num(arena, 0)
            do r = 1, 3
                d = d + (vertices(r, a) - vertices(r, b))*(vertices(r, a) - vertices(r, b))
            end do
            roots(k) = d
        end do
        call output_spec(spec, "generated_tetra_edge_lengths_squared", &
            "fortfem_generated_tetra_edge_lengths_squared", ["vertices"], &
            ["(3,4)"], "product", "(6)", 6)
        call write_kernel(roots, spec)
    end subroutine emit_edge_lengths

    subroutine emit_patch_load()
        !! product(s) = source psi_s - A grad u_h . grad psi_s, s = 1..4,
        !! product(4 + s) = psi_s at the reference point (barycentric).
        type(kernel_spec_t) :: spec
        type(expr_t) :: point(3), source, flux(3), grad(3, 4), psi(4), roots(8)
        type(expr_t) :: dot
        integer :: r, s

        do r = 1, 3
            point(r) = sym_(at("point", r))
            flux(r) = sym_(at("flux", r))
        end do
        source = sym_("source")
        do s = 1, 4
            do r = 1, 3
                grad(r, s) = sym_(at2("hat_grad", r, s))
            end do
        end do
        psi(1) = num(arena, 1) - point(1) - point(2) - point(3)
        psi(2) = point(1)
        psi(3) = point(2)
        psi(4) = point(3)
        do s = 1, 4
            dot = flux(1)*grad(1, s) + flux(2)*grad(2, s) + flux(3)*grad(3, s)
            roots(s) = source*psi(s) - dot
            roots(4 + s) = psi(s)
        end do
        call output_spec(spec, "generated_tetra_patch_load", &
            "fortfem_generated_tetra_patch_load", &
            ["point      ", "source     ", "flux       ", "hat_grad   "], &
            ["(3)  ", "     ", "(3)  ", "(3,4)"], "product", "(8)", 8)
        call write_kernel(roots, spec)
    end subroutine emit_patch_load

    subroutine emit_patch_gram()
        !! product(1) = wdet ainv phi_j . phi_k.
        type(kernel_spec_t) :: spec
        type(expr_t) :: phi_j(3), phi_k(3), ainv, wdet, roots(1)
        integer :: r

        do r = 1, 3
            phi_j(r) = sym_(at("phi_j", r))
            phi_k(r) = sym_(at("phi_k", r))
        end do
        ainv = sym_("ainv")
        wdet = sym_("wdet")
        roots(1) = wdet*ainv*(phi_j(1)*phi_k(1) + phi_j(2)*phi_k(2) &
            + phi_j(3)*phi_k(3))
        call output_spec(spec, "generated_tetra_patch_gram", &
            "fortfem_generated_tetra_patch_gram", &
            ["phi_j", "phi_k", "ainv ", "wdet "], ["(3)", "(3)", "   ", "   "], &
            "product", "(1)", 1)
        call write_kernel(roots, spec)
    end subroutine emit_patch_gram

    subroutine emit_patch_flux_rhs()
        !! product(1) = wdet psi ainv (A grad u_h) . phi.
        type(kernel_spec_t) :: spec
        type(expr_t) :: phi(3), flux(3), psi, ainv, wdet, roots(1)
        integer :: r

        do r = 1, 3
            phi(r) = sym_(at("phi", r))
            flux(r) = sym_(at("flux", r))
        end do
        psi = sym_("psi")
        ainv = sym_("ainv")
        wdet = sym_("wdet")
        roots(1) = wdet*psi*ainv*(flux(1)*phi(1) + flux(2)*phi(2) + flux(3)*phi(3))
        call output_spec(spec, "generated_tetra_patch_flux_rhs", &
            "fortfem_generated_tetra_patch_flux_rhs", &
            ["phi  ", "flux ", "psi  ", "ainv ", "wdet "], &
            ["(3)", "(3)", "   ", "   ", "   "], "product", "(1)", 1)
        call write_kernel(roots, spec)
    end subroutine emit_patch_flux_rhs

    subroutine emit_patch_multiplier()
        !! product(1) = wdet div_k poly_l (constraint coupling),
        !! product(2) = wdet load poly_l (constraint right-hand side).
        type(kernel_spec_t) :: spec
        type(expr_t) :: div, load, poly, wdet, roots(2)

        div = sym_("div")
        load = sym_("load")
        poly = sym_("poly")
        wdet = sym_("wdet")
        roots(1) = wdet*div*poly
        roots(2) = wdet*load*poly
        call output_spec(spec, "generated_tetra_patch_multiplier", &
            "fortfem_generated_tetra_patch_multiplier", &
            ["div  ", "load ", "poly ", "wdet "], ["   ", "   ", "   ", "   "], &
            "product", "(2)", 2)
        call write_kernel(roots, spec)
    end subroutine emit_patch_multiplier

    subroutine emit_projection_terms()
        !! product(1) = wdet poly_i poly_j (mass), product(2) = wdet f poly_i.
        type(kernel_spec_t) :: spec
        type(expr_t) :: poly_i, poly_j, source, wdet, roots(2)

        poly_i = sym_("poly_i")
        poly_j = sym_("poly_j")
        source = sym_("source")
        wdet = sym_("wdet")
        roots(1) = wdet*poly_i*poly_j
        roots(2) = wdet*source*poly_i
        call output_spec(spec, "generated_tetra_projection_terms", &
            "fortfem_generated_tetra_projection_terms", &
            ["poly_i", "poly_j", "source", "wdet  "], &
            ["     ", "     ", "     ", "     "], "product", "(2)", 2)
        call write_kernel(roots, spec)
    end subroutine emit_projection_terms

    subroutine emit_flux_defect()
        !! product(1) = wdet |sigma - A grad u_h|^2 / A.
        type(kernel_spec_t) :: spec
        type(expr_t) :: sigma(3), flux(3), coefficient, wdet, roots(1), d
        integer :: r

        do r = 1, 3
            sigma(r) = sym_(at("sigma", r))
            flux(r) = sym_(at("flux", r))
        end do
        coefficient = sym_("coefficient")
        wdet = sym_("wdet")
        d = num(arena, 0)
        do r = 1, 3
            d = d + (sigma(r) - flux(r))*(sigma(r) - flux(r))
        end do
        roots(1) = wdet*d/coefficient
        call output_spec(spec, "generated_tetra_flux_defect", &
            "fortfem_generated_tetra_flux_defect", &
            ["sigma      ", "flux       ", "coefficient", "wdet       "], &
            ["(3)", "(3)", "   ", "   "], "product", "(1)", 1)
        call write_kernel(roots, spec)
    end subroutine emit_flux_defect

    subroutine emit_divergence_defect()
        !! product(1) = wdet (div sigma + projection)^2.
        type(kernel_spec_t) :: spec
        type(expr_t) :: div_sigma, projection, wdet, roots(1)

        div_sigma = sym_("div_sigma")
        projection = sym_("projection")
        wdet = sym_("wdet")
        roots(1) = wdet*(div_sigma + projection)*(div_sigma + projection)
        call output_spec(spec, "generated_tetra_divergence_defect", &
            "fortfem_generated_tetra_divergence_defect", &
            ["div_sigma ", "projection", "wdet      "], ["   ", "   ", "   "], &
            "product", "(1)", 1)
        call write_kernel(roots, spec)
    end subroutine emit_divergence_defect

    subroutine emit_oscillation_defect()
        !! product(1) = wdet (source - projection)^2.
        type(kernel_spec_t) :: spec
        type(expr_t) :: source, projection, wdet, roots(1)

        source = sym_("source")
        projection = sym_("projection")
        wdet = sym_("wdet")
        roots(1) = wdet*(source - projection)*(source - projection)
        call output_spec(spec, "generated_tetra_oscillation_defect", &
            "fortfem_generated_tetra_oscillation_defect", &
            ["source    ", "projection", "wdet      "], ["   ", "   ", "   "], &
            "product", "(1)", 1)
        call write_kernel(roots, spec)
    end subroutine emit_oscillation_defect

end program gen_equilibrated_flux_3d_products
