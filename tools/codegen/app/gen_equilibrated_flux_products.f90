program gen_equilibrated_flux_products
    !! Pointwise leaves of the patchwise equilibrated-flux estimator: cell
    !! geometry (Lagrange P1/P2/P3 Jacobians, determinant, inverse transpose,
    !! hat gradients), reference monomials, the quadrature integrands
    !! (densities per unit reference area) and the Neumann edge leaves
    !! (shifted Legendre moments, edge speed, normal trace, edge integrands).
    !! Contractions over runtime sizes, patch topology and the dense solves
    !! stay handwritten in the estimator.
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_expr, only: expr_t, num, sqrt, operator(*), operator(+), &
        operator(-), operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    integer :: unit, degree, l
    character(:), allocatable :: text
    logical :: mapped_output

    mapped_output = .true.
    call arena%init()
    engine = make_native_engine(arena)
    open(newunit=unit, file=generated_path("fortfem_mapped_oscillation_products.f90"), &
        status="replace", action="write")

    call emit("equilibrated_bernstein_midpoint", &
        [str("midpoint"), str("left"), str("right")], &
        [str("(2,2)"), str("(2,2)"), str("(2,2)")], &
        [str("control")], [str("(2,2)")], &
        [str("control(1,1)"), str("control(2,1)"), &
        str("control(1,2)"), str("control(2,2)")], midpoint_roots())
    call emit("equilibrated_relative_jacobian_norm", &
        [str("jacobian"), str("inverse_transpose")], &
        [str("(2,2)"), str("(2,2)")], &
        [str("value")], [str("")], [str("value")], relative_norm_roots())
    call emit("equilibrated_mapped_poincare", &
        [str("diameter"), str("norm_bound"), str("det_min"), &
        str("lower"), str("pi")], &
        [str(""), str(""), str(""), str(""), str("")], &
        [str("value")], [str("")], [str("value")], poincare_roots())

    call determinant_control_kernel()

    call emit("equilibrated_reference_oscillation", &
        [str("source"), str("projection"), str("determinant")], &
        [str(""), str(""), str("")], &
        [str("value")], [str("")], [str("value")], reference_oscillation_roots())

    close(unit)
    mapped_output = .false.
    open(newunit=unit, file=generated_path("fortfem_equilibrated_flux_products.f90"), &
        status="replace", action="write")

    ! Lagrange geometry Jacobian for P1 (affine), P2 and P3 isoparametric
    ! cells (vertices first).
    call jacobian_kernel(3, "p1")
    call jacobian_kernel(6, "p2")
    call jacobian_kernel(10, "p3")

    call emit("equilibrated_determinant", &
        [str("jacobian")], [str("(2,2)")], &
        [str("determinant")], [str("")], [str("determinant")], determinant_roots())

    call emit("equilibrated_inverse_transpose", &
        [str("jacobian"), str("determinant")], [str("(2,2)"), str("")], &
        [str("inverse_transpose")], [str("(2,2)")], &
        [str("inverse_transpose(1,1)"), str("inverse_transpose(2,1)"), &
        str("inverse_transpose(1,2)"), str("inverse_transpose(2,2)")], &
        inverse_transpose_roots())

    call emit("equilibrated_hat_gradients", &
        [str("inverse_transpose")], [str("(2,2)")], &
        [str("hat_gradient")], [str("(2,3)")], &
        [str("hat_gradient(1,1)"), str("hat_gradient(2,1)"), &
        str("hat_gradient(1,2)"), str("hat_gradient(2,2)"), &
        str("hat_gradient(1,3)"), str("hat_gradient(2,3)")], &
        hat_gradient_roots())

    call emit("equilibrated_hat_values", &
        [str("xi"), str("eta")], [str(""), str("")], &
        [str("hat")], [str("(3)")], &
        [str("hat(1)"), str("hat(2)"), str("hat(3)")], &
        hat_value_roots())

    do degree = 0, 6
        call monomial_kernel(degree)
    end do

    call emit("equilibrated_load_integrand", &
        [str("source"), str("hat"), str("flux"), str("hat_gradient")], &
        [str(""), str(""), str("(2)"), str("(2)")], &
        [str("value")], [str("")], [str("value")], load_roots())

    call emit("equilibrated_flux_defect_integrand", &
        [str("residual"), str("coefficient"), str("determinant")], &
        [str("(2)"), str(""), str("")], &
        [str("value")], [str("")], [str("value")], flux_defect_roots())

    ! Symmetric positive definite tensor coefficient A (anisotropic case):
    ! A^(-1), its action, the A^(-1)-weighted defect density and the smaller
    ! eigenvalue (pointwise check of the caller's lower bound).
    call emit("equilibrated_tensor_inverse", &
        [str("coefficient")], [str("(2,2)")], &
        [str("inverse")], [str("(2,2)")], &
        [str("inverse(1,1)"), str("inverse(2,1)"), str("inverse(1,2)"), &
        str("inverse(2,2)")], tensor_inverse_roots())

    call emit("equilibrated_tensor_apply", &
        [str("matrix"), str("vector")], [str("(2,2)"), str("(2)")], &
        [str("value")], [str("(2)")], [str("value(1)"), str("value(2)")], &
        tensor_apply_roots())

    call emit("equilibrated_tensor_flux_defect_integrand", &
        [str("residual"), str("inverse"), str("determinant")], &
        [str("(2)"), str("(2,2)"), str("")], &
        [str("value")], [str("")], [str("value")], tensor_flux_defect_roots())

    call emit("equilibrated_tensor_min_eigenvalue", &
        [str("coefficient")], [str("(2,2)")], &
        [str("value")], [str("")], [str("value")], tensor_min_eigenvalue_roots())

    call emit("equilibrated_oscillation_integrand", &
        [str("source"), str("projection"), str("determinant")], &
        [str(""), str(""), str("")], &
        [str("value")], [str("")], [str("value")], oscillation_roots())

    call emit("equilibrated_divergence_integrand", &
        [str("divergence"), str("projection"), str("determinant")], &
        [str(""), str(""), str("")], &
        [str("value")], [str("")], [str("value")], divergence_roots())

    call emit("equilibrated_edge_legendre", &
        [str("s")], [str("")], [str("values")], [str("(7)")], &
        [str("values(1)"), str("values(2)"), str("values(3)"), str("values(4)"), &
        str("values(5)"), str("values(6)"), str("values(7)")], legendre_roots())

    call emit("equilibrated_edge_speed", &
        [str("jacobian"), str("reference_tangent")], [str("(2,2)"), str("(2)")], &
        [str("speed")], [str("")], [str("speed")], edge_speed_roots())

    call emit("equilibrated_reference_normal", &
        [str("reference_tangent")], [str("(2)")], &
        [str("normal")], [str("(2)")], &
        [str("normal(1)"), str("normal(2)")], reference_normal_roots())

    call emit("equilibrated_flux_normal_component", &
        [str("reference_value"), str("normal")], [str("(2)"), str("(2)")], &
        [str("value")], [str("")], [str("value")], normal_component_roots())

    call emit("equilibrated_neumann_target", &
        [str("psi"), str("data"), str("speed")], [str(""), str(""), str("")], &
        [str("value")], [str("")], [str("value")], neumann_target_roots())

    call emit("equilibrated_neumann_defect_integrand", &
        [str("data"), str("reference_flux"), str("speed")], &
        [str(""), str(""), str("")], [str("value")], [str("")], &
        [str("value")], neumann_defect_roots())

    close(unit)

contains

    function ex(name) result(e)
        character(*), intent(in) :: name
        type(expr_t) :: e
        e = sym(arena, name)
    end function ex

    function ix(name, i, j) result(text)
        character(*), intent(in) :: name
        integer, intent(in) :: i
        integer, intent(in), optional :: j
        character(:), allocatable :: text
    logical :: mapped_output
        text = name//"("//itxt(i)
        if (present(j)) text = text//","//itxt(j)
        text = text//")"
    end function ix

    function itxt(v) result(text)
        integer, intent(in) :: v
        character(:), allocatable :: text
    logical :: mapped_output
        character(16) :: buffer
        write(buffer, "(i0)") v
        text = trim(buffer)
    end function itxt

    subroutine jacobian_kernel(count, label)
        !! jacobian(i,c) = sum_j nodes(i,j) grad(j,c) for a P_k cell.
        integer, intent(in) :: count
        character(*), intent(in) :: label
        type(expr_t) :: roots(4), acc
        type(str_t) :: refs(4)
        integer :: r, c, k

        do c = 1, 2
            do r = 1, 2
                acc = num(arena, 0)
                do k = 1, count
                    acc = acc + ex(ix("nodes", r, k))*ex(ix("grad", k, c))
                end do
                roots(2*(c - 1) + r) = acc
                refs(2*(c - 1) + r) = str(ix("jacobian", r, c))
            end do
        end do
        call emit_roots("equilibrated_jacobian_"//label, &
            [str("nodes"), str("grad")], &
            [str("(2,"//itxt(count)//")"), str("("//itxt(count)//",2)")], &
            [str("jacobian")], [str("(2,2)")], refs, roots)
    end subroutine jacobian_kernel

    function midpoint_roots() result(roots)
        ! Quadratic Bernstein edge control from midpoint and endpoint values.
        type(expr_t) :: roots(4)
        integer :: i, j, n
        n = 0
        do j = 1, 2
            do i = 1, 2
                n = n + 1
                roots(n) = num(arena, 2)*ex(ix("midpoint", i, j)) &
                    - (ex(ix("left", i, j)) + ex(ix("right", i, j)))/num(arena, 2)
            end do
        end do
    end function midpoint_roots

    function relative_norm_roots() result(roots)
        ! Spectral norm of J J0^(-1), from the larger eigenvalue of M^T M.
        type(expr_t) :: roots(1), entry(2, 2), aa, bb, cc
        integer :: i, j, m
        do j = 1, 2
            do i = 1, 2
                entry(i, j) = num(arena, 0)
                do m = 1, 2
                    entry(i, j) = entry(i, j) + ex(ix("jacobian", i, m)) &
                        *ex(ix("inverse_transpose", j, m))
                end do
            end do
        end do
        aa = entry(1,1)*entry(1,1) + entry(2,1)*entry(2,1)
        bb = entry(1,1)*entry(1,2) + entry(2,1)*entry(2,2)
        cc = entry(1,2)*entry(1,2) + entry(2,2)*entry(2,2)
        roots(1) = sqrt((aa + cc + sqrt((aa - cc)*(aa - cc) &
            + num(arena, 4)*bb*bb))/num(arena, 2))
    end function relative_norm_roots

    function poincare_roots() result(roots)
        ! Pullback of Payne--Weinberger on the straight vertex triangle:
        ! h/pi * sup ||DG|| / sqrt(det_min A_min), for the reference residual.
        type(expr_t) :: roots(1)
        roots(1) = ex("diameter")*ex("norm_bound") &
            /(ex("pi")*sqrt(ex("det_min")*ex("lower")))
    end function poincare_roots

    subroutine determinant_control_kernel()
        ! Exact Bernstein product in barycentric coordinates: a quadratic
        ! Jacobian has a quartic determinant. Controls follow vertices then
        ! edge midpoints; determinant controls use lexicographic powers.
        type(expr_t) :: roots(15), term
        type(str_t) :: refs(15)
        integer :: powers(3, 6), weights(6), a, b, i, j, n, weight4
        powers = reshape([2,0,0, 0,2,0, 0,0,2, 1,1,0, 0,1,1, 1,0,1], [3,6])
        weights = [1, 1, 1, 2, 2, 2]
        n = 0
        do a = 0, 4
            do b = 0, 4 - a
                n = n + 1
                refs(n) = str("controls("//itxt(n)//")")
                roots(n) = num(arena, 0)
                weight4 = factorial(4)/(factorial(a)*factorial(b)*factorial(4-a-b))
                do i = 1, 6
                    do j = 1, 6
                        if (any(powers(:, i) + powers(:, j) /= [a, b, 4-a-b])) cycle
                        term = ex("jacobians(1,1,"//itxt(i)//")") &
                            *ex("jacobians(2,2,"//itxt(j)//")") &
                            - ex("jacobians(1,2,"//itxt(i)//")") &
                            *ex("jacobians(2,1,"//itxt(j)//")")
                        roots(n) = roots(n) + num(arena, weights(i)*weights(j)) &
                            *term/num(arena, weight4)
                    end do
                end do
            end do
        end do
        call emit("equilibrated_determinant_controls", &
            [str("jacobians")], [str("(2,2,6)")], &
            [str("controls")], [str("(15)")], refs, roots)
    end subroutine determinant_control_kernel

    integer function factorial(n) result(value)
        integer, intent(in) :: n
        integer :: i
        value = 1
        do i = 2, n
            value = value*i
        end do
    end function factorial

    function determinant_roots() result(roots)
        type(expr_t) :: roots(1)
        roots(1) = ex("jacobian(1,1)")*ex("jacobian(2,2)") &
            - ex("jacobian(1,2)")*ex("jacobian(2,1)")
    end function determinant_roots

    function inverse_transpose_roots() result(roots)
        !! J^(-T) = [[j22, -j21], [-j12, j11]] / det; roots follow the
        !! output references (1,1), (2,1), (1,2), (2,2).
        type(expr_t) :: roots(4)
        roots(1) = ex("jacobian(2,2)")/ex("determinant")
        roots(2) = -ex("jacobian(1,2)")/ex("determinant")
        roots(3) = -ex("jacobian(2,1)")/ex("determinant")
        roots(4) = ex("jacobian(1,1)")/ex("determinant")
    end function inverse_transpose_roots

    function hat_gradient_roots() result(roots)
        !! grad(hat_1) = -(e1+e2), grad(hat_2) = e1, grad(hat_3) = e2 applied
        !! through inverse_transpose.
        type(expr_t) :: roots(6)
        integer :: r
        do r = 1, 2
            roots(r) = -(ex(ix("inverse_transpose", r, 1)) &
                + ex(ix("inverse_transpose", r, 2)))
            roots(2 + r) = ex(ix("inverse_transpose", r, 1))
            roots(4 + r) = ex(ix("inverse_transpose", r, 2))
        end do
    end function hat_gradient_roots

    function hat_value_roots() result(roots)
        type(expr_t) :: roots(3)
        roots(1) = num(arena, 1) - ex("xi") - ex("eta")
        roots(2) = ex("xi")
        roots(3) = ex("eta")
    end function hat_value_roots

    subroutine monomial_kernel(deg)
        !! values(l) = xi**i eta**j in the estimator's monomial order.
        integer, intent(in) :: deg
        type(expr_t), allocatable :: roots(:)
        type(str_t), allocatable :: refs(:)
        type(expr_t) :: term
        integer :: ii, jj, p, count

        count = (deg + 1)*(deg + 2)/2
        allocate (roots(count), refs(count))
        l = 0
        do ii = 0, deg
            do jj = 0, deg - ii
                l = l + 1
                term = num(arena, 1)
                do p = 1, ii
                    term = term*ex("xi")
                end do
                do p = 1, jj
                    term = term*ex("eta")
                end do
                roots(l) = term
                refs(l) = str(ix("values", l))
            end do
        end do
        call emit_roots("equilibrated_monomials_degree"//itxt(deg), &
            [str("xi"), str("eta")], [str(""), str("")], &
            [str("values")], [str("("//itxt(count)//")")], refs, roots)
    end subroutine monomial_kernel

    function load_roots() result(roots)
        type(expr_t) :: roots(1)
        roots(1) = ex("source")*ex("hat") - ex("flux(1)")*ex("hat_gradient(1)") &
            - ex("flux(2)")*ex("hat_gradient(2)")
    end function load_roots

    function flux_defect_roots() result(roots)
        !! A^(-1) |sigma - A grad u_h|^2 det: energy density per reference area.
        type(expr_t) :: roots(1)
        roots(1) = (ex("residual(1)")*ex("residual(1)") + &
            ex("residual(2)")*ex("residual(2)"))*ex("determinant")/ex("coefficient")
    end function flux_defect_roots

    function tensor_inverse_roots() result(roots)
        !! Inverse of the symmetrised 2x2 coefficient, column-major order.
        type(expr_t) :: roots(4), b, det
        b = (ex("coefficient(1,2)") + ex("coefficient(2,1)"))/num(arena, 2)
        det = ex("coefficient(1,1)")*ex("coefficient(2,2)") - b*b
        roots(1) = ex("coefficient(2,2)")/det
        roots(2) = -b/det
        roots(3) = -b/det
        roots(4) = ex("coefficient(1,1)")/det
    end function tensor_inverse_roots

    function tensor_apply_roots() result(roots)
        type(expr_t) :: roots(2)
        roots(1) = ex("matrix(1,1)")*ex("vector(1)") + ex("matrix(1,2)")*ex("vector(2)")
        roots(2) = ex("matrix(2,1)")*ex("vector(1)") + ex("matrix(2,2)")*ex("vector(2)")
    end function tensor_apply_roots

    function tensor_flux_defect_roots() result(roots)
        !! r . A^(-1) r det, r = sigma - A grad u_h: energy density per
        !! reference area of the A^(-1)-weighted flux defect.
        type(expr_t) :: roots(1), r1, r2
        r1 = ex("residual(1)")
        r2 = ex("residual(2)")
        roots(1) = (ex("inverse(1,1)")*r1*r1 + (ex("inverse(1,2)") &
            + ex("inverse(2,1)"))*r1*r2 + ex("inverse(2,2)")*r2*r2) &
            *ex("determinant")
    end function tensor_flux_defect_roots

    function tensor_min_eigenvalue_roots() result(roots)
        !! Smaller eigenvalue (a + c)/2 - sqrt(((a - c)/2)^2 + b^2) of the
        !! symmetrised coefficient [[a, b], [b, c]].
        type(expr_t) :: roots(1), b, h
        b = (ex("coefficient(1,2)") + ex("coefficient(2,1)"))/num(arena, 2)
        h = (ex("coefficient(1,1)") - ex("coefficient(2,2)"))/num(arena, 2)
        roots(1) = (ex("coefficient(1,1)") + ex("coefficient(2,2)"))/num(arena, 2) &
            - sqrt(h*h + b*b)
    end function tensor_min_eigenvalue_roots

    function oscillation_roots() result(roots)
        !! (f - Pi f)^2 det with Pi f = projection/det, the projection being the
        !! reference-cell L2 projection of det f onto P_k.
        type(expr_t) :: roots(1), d
        d = ex("determinant")*ex("source") - ex("projection")
        roots(1) = d*d/ex("determinant")
    end function oscillation_roots

    function reference_oscillation_roots() result(roots)
        type(expr_t) :: roots(1), d
        d = ex("determinant")*ex("source") - ex("projection")
        roots(1) = d*d
    end function reference_oscillation_roots

    function divergence_roots() result(roots)
        !! (div sigma + Pi f)^2 det, div sigma the physical divergence.
        type(expr_t) :: roots(1), d
        d = ex("determinant")*ex("divergence") + ex("projection")
        roots(1) = d*d/ex("determinant")
    end function divergence_roots

    function legendre_roots() result(roots)
        !! Shifted Legendre polynomials P_0..P_6 on [0,1] by Bonnet's recurrence
        !! in x = 2 s - 1, the edge-moment basis of the Nedelec/RT dofs.
        type(expr_t) :: roots(7), x
        integer :: j

        x = num(arena, 2)*ex("s") - num(arena, 1)
        roots(1) = num(arena, 1)
        roots(2) = x
        do j = 1, 5
            roots(j + 2) = (num(arena, 2*j + 1)*x*roots(j + 1) &
                - num(arena, j)*roots(j))/num(arena, j + 1)
        end do
    end function legendre_roots

    function edge_speed_roots() result(roots)
        !! |J t| for the reference edge tangent t: edge length per parameter.
        type(expr_t) :: roots(1), t1, t2
        t1 = ex("jacobian(1,1)")*ex("reference_tangent(1)") + &
            ex("jacobian(1,2)")*ex("reference_tangent(2)")
        t2 = ex("jacobian(2,1)")*ex("reference_tangent(1)") + &
            ex("jacobian(2,2)")*ex("reference_tangent(2)")
        roots(1) = sqrt(t1*t1 + t2*t2)
    end function edge_speed_roots

    function reference_normal_roots() result(roots)
        !! Outward reference normal (t_y, -t_x) of a counter-clockwise edge with
        !! tangent t: the direction of the RT normal-trace moments.
        type(expr_t) :: roots(2)
        roots(1) = ex("reference_tangent(2)")
        roots(2) = -ex("reference_tangent(1)")
    end function reference_normal_roots

    function normal_component_roots() result(roots)
        type(expr_t) :: roots(1)
        roots(1) = ex("reference_value(1)")*ex("normal(1)") + &
            ex("reference_value(2)")*ex("normal(2)")
    end function normal_component_roots

    function neumann_target_roots() result(roots)
        !! psi g |x'(s)|: the flux per unit edge parameter to be projected.
        type(expr_t) :: roots(1)
        roots(1) = ex("psi")*ex("data")*ex("speed")
    end function neumann_target_roots

    function neumann_defect_roots() result(roots)
        !! (g - q/speed)^2 speed: squared L2(edge) density per unit parameter,
        !! q the reconstructed flux per unit parameter.
        type(expr_t) :: roots(1), d
        d = ex("data") - ex("reference_flux")/ex("speed")
        roots(1) = d*d*ex("speed")
    end function neumann_defect_roots

    subroutine emit(name, args, arg_shapes, outs, out_shapes, refs, roots)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), arg_shapes(:), outs(:), out_shapes(:)
        type(str_t), intent(in) :: refs(:)
        type(expr_t), intent(in) :: roots(:)

        call emit_roots(name, args, arg_shapes, outs, out_shapes, refs, roots)
    end subroutine emit

    subroutine emit_roots(name, args, arg_shapes, outs, out_shapes, refs, roots)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), arg_shapes(:), outs(:), out_shapes(:)
        type(str_t), intent(in) :: refs(:)
        type(expr_t), intent(in) :: roots(:)
        type(expr_t) :: simplified(size(roots))
        type(engine_result_t) :: result
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code
        integer :: k

        do k = 1, size(roots)
            result = engine%simplify(roots(k))
            if (.not. result%ok) error stop "equilibrated simplification failed"
            simplified(k) = result%value
        end do
        spec%name = str("generated_"//name)
        spec%module_name = str("fortfem_generated_"//name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_equilibrated_flux_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%pure_procedure = .true.
        spec%args = args
        spec%arg_shapes = arg_shapes
        spec%outputs = outs
        spec%output_shapes = out_shapes
        spec%output_references = refs
        code = chars(emit_kernel(simplified, spec))
        if (mapped_output) then
            call write_wrapped(code(:len(code) - 1))
        else
            write(unit, "(a)") code(:len(code) - 1)
        end if
    end subroutine emit_roots

    subroutine write_wrapped(code)
        ! The pinned emitter uses 100 columns. Wrap the new products at token
        ! boundaries, leaving margin for the standard formatter's whitespace.
        character(*), intent(in) :: code
        character(:), allocatable :: line
        integer :: first, last, cut
        first = 1
        do while (first <= len(code))
            last = index(code(first:), achar(10))
            if (last == 0) then
                last = len(code)
            else
                last = first + last - 2
            end if
            line = code(first:last)
            do while (len_trim(line) > 80)
                cut = scan(line(:78), ' ', back=.true.)
                if (cut <= len(line) - len_trim(adjustl(line))) &
                    error stop "no generated token boundary"
                write(unit, '(a)') line(:cut - 1)//' &'
                line = '            '//trim(adjustl(line(cut + 1:)))
            end do
            write(unit, '(a)') line
            first = last + 2
        end do
    end subroutine write_wrapped

end program gen_equilibrated_flux_products
