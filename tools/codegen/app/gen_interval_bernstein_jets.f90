program gen_interval_bernstein_jets
    !! Bernstein polynomials on [0, 1] with first and second derivatives,
    !! B_i^p(t) = binom(p, i) t^i (1 - t)^(p - i), i = 0..p, one kernel per
    !! degree. Spline cells evaluate every smooth spline space through these
    !! jets and a cell extraction matrix (Borden et al. 2011, IJNME 87:15).
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_engine, only: engine_result_t
    use fortsym_expr, only: expr_t, sym, num, operator(-), operator(*), &
        operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t, itoa
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    integer, parameter :: max_degree = 8
    type(native_engine_t) :: engine
    integer :: unit, degree

    engine = make_native_engine()
    open(newunit=unit, file=generated_path("fortfem_interval_bernstein_jets.f90"), &
        status="replace", action="write")
    do degree = 0, max_degree
        call write_jet(unit, degree)
    end do
    close(unit)

contains

    subroutine write_jet(unit, degree)
        integer, intent(in) :: unit, degree

        type(kernel_spec_t) :: spec
        type(engine_result_t) :: result
        type(expr_t) :: t, basis, first
        type(expr_t), allocatable :: roots(:)
        type(str_t), allocatable :: refs(:)
        character(:), allocatable :: code
        integer :: i, n, k

        n = degree + 1
        allocate(roots(3*n), refs(3*n))
        t = sym("t")
        do i = 0, degree
            basis = num(binomial(degree, i))*t**i*(num(1) - t)**(degree - i)
            first = diff(basis, t)
            roots(i + 1) = basis
            roots(n + i + 1) = first
            roots(2*n + i + 1) = diff(first, t)
            refs(i + 1) = str("value("//itoa(i + 1)//")")
            refs(n + i + 1) = str("first("//itoa(i + 1)//")")
            refs(2*n + i + 1) = str("second("//itoa(i + 1)//")")
        end do
        do k = 1, size(roots)
            result = engine%simplify(roots(k))
            if (.not. result%ok) error stop "Bernstein jet simplification failed"
            roots(k) = result%value
        end do
        spec%name = str("generated_interval_bernstein_jet_"//itoa(degree))
        spec%module_name = str("fortfem_generated_interval_bernstein_jet_"// &
            itoa(degree))
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%generator = str("gen_interval_bernstein_jets")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        spec%args = [str("t")]
        spec%outputs = [str("value"), str("first"), str("second")]
        spec%output_shapes = [str("("//itoa(n)//")"), str("("//itoa(n)//")"), &
            str("("//itoa(n)//")")]
        spec%output_references = refs
        code = chars(emit_kernel(roots, spec))
        write(unit, "(a)") code(:len(code) - 1)
    end subroutine write_jet

    pure integer function binomial(n, k) result(value)
        integer, intent(in) :: n, k
        integer :: j

        value = 1
        do j = 1, k
            value = value*(n - k + j)/j
        end do
    end function binomial

end program gen_interval_bernstein_jets
