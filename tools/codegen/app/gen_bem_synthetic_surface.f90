program gen_bem_synthetic_surface
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, cos, sin, operator(+), operator(-), &
        operator(*)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_string, only: str, chars
    use fortfem_codegen_provenance, only: generated_path, fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: theta, phi, major, minor, ripple, nfp, radial, roots(9)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: code
    integer :: i, unit
    call arena%init()
    theta=sym(arena,'theta'); phi=sym(arena,'phi')
    major=sym(arena,'major'); minor=sym(arena,'minor')
    ripple=sym(arena,'ripple'); nfp=sym(arena,'nfp')
    radial=major+minor*cos(theta)+ripple*cos(theta-nfp*phi)
    roots(1:3)=[radial*cos(phi),radial*sin(phi), &
        minor*sin(theta)+ripple*sin(theta-nfp*phi)]
    do i=1,3
        roots(i+3)=diff(roots(i),theta)
        roots(i+6)=diff(roots(i),phi)
    end do
    spec%name=str('bem_synthetic_surface')
    spec%module_name=str('fortfem_generated_bem_synthetic_surface')
    spec%generator=str('gen_bem_synthetic_surface')
    spec%generator_revision=str(fortsym_revision())
    spec%regenerate_command=str('cd tools/codegen && ./generate.sh')
    spec%mode=KERNEL_SUBROUTINE; spec%pure_procedure=.true.
    spec%args=[str('theta'),str('phi'),str('major'),str('minor'), &
        str('ripple'),str('nfp')]
    spec%outputs=[str('x'),str('y'),str('z'),str('dt1'),str('dt2'),str('dt3'), &
        str('dp1'),str('dp2'),str('dp3')]
    code=chars(emit_kernel(roots,spec))
    open(newunit=unit,file=generated_path('fortfem_bem_synthetic_surface.f90'), &
        status='replace',action='write')
    write(unit,'(a)') code(:len(code)-1)
    close(unit)
end program
