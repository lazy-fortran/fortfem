program gen_lab_grad_shafranov_source_products
    ! Consumer-side pointwise kernel for the NEQ1 lab benchmark (not a library
    ! kernel): the axisymmetric Grad-Shafranov source in weak form,
    !   div((1/R) grad psi) = -g,  g = mu0 R dp/dpsi + F dF/dpsi / R,
    ! with normalised-flux profiles p = pa (1 - psin)^ap and
    ! F^2 = F0^2 + cf (1 - psin)^af, psin = (psi - psia)/(psib - psia).
    ! Outputs g and its derivatives with respect to the local flux psi and the
    ! free-axis value psia, which carries the nonlocal Newton term.
    use fortsym_arena, only: arena_t
    use fortsym_diff, only: diff
    use fortsym_expr, only: expr_t, sym, real_expr, &
                            operator(+), operator(-), operator(*), operator(/), &
                            operator(**)
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: str, chars
    use fortfem_codegen_provenance, only: fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: r, psi, psia, psib, mu0, pa, cf, psin, pressure, fsq, g, &
                    roots(3), one, half
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code, directory, all_code
    character(4096) :: buffer
    character(32) :: tag
    integer :: unit, length, status, variant
    integer, parameter :: exponents(2, 2) = reshape([1, 1, 2, 3], [2, 2])

    call arena%init()
    one = real_expr(arena, 1.0d0)
    half = real_expr(arena, 0.5d0)
    r = sym(arena, 'r')
    psi = sym(arena, 'psi')
    psia = sym(arena, 'psia')
    psib = sym(arena, 'psib')
    mu0 = sym(arena, 'mu0')
    pa = sym(arena, 'pa')
    cf = sym(arena, 'cf')
    psin = (psi - psia)/(psib - psia)
    all_code = ''
    do variant = 1, size(exponents, 2)
        pressure = pa*(one - psin)**exponents(1, variant)
        fsq = cf*(one - psin)**exponents(2, variant)
        g = mu0*r*diff(pressure, psi) + half*diff(fsq, psi)/r
        roots(1) = g
        roots(2) = diff(g, psi)
        roots(3) = diff(g, psia)
        write (tag, '(a,i0,a,i0)') 'p', exponents(1, variant), 'f', &
            exponents(2, variant)
        if (allocated(spec%args)) deallocate (spec%args)
        if (allocated(spec%outputs)) deallocate (spec%outputs)
        allocate (spec%args(7), spec%outputs(3))
        spec%args = [str('r'), str('psi'), str('psia'), str('psib'), str('mu0'), &
                     str('pa'), str('cf')]
        spec%outputs = [str('g'), str('g_psi'), str('g_psia')]
        spec%name = str('lab_gs_source_'//trim(tag))
        spec%module_name = str('fortfem_lab_gs_source_'//trim(tag))
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str('gen_lab_grad_shafranov_source_products')
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str('tools/codegen/generate.sh')
        spec%pure_procedure = .true.
        code = chars(emit_kernel(roots, spec))
        all_code = all_code//code
    end do
    call get_environment_variable('FORTFEM_LAB_CODEGEN_OUTPUT_DIR', buffer, &
                                  length=length, status=status)
    if (status /= 0 .or. length == 0) then
        directory = '../../benchmark/lab/neq1/generated/'
    else
        directory = buffer(:length)//'/'
    end if
    open (newunit=unit, file=directory//'fortfem_lab_gs_source_products.f90', &
          status='replace', action='write')
    write (unit, '(a)') all_code(:len(all_code) - 1)
    close (unit)
end program gen_lab_grad_shafranov_source_products
