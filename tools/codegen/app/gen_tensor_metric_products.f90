program gen_tensor_metric_products
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, &
        operator(+), operator(-), operator(*), operator(/)
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: str, chars
    use fortfem_codegen_provenance, only: fortsym_revision, generated_path
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: j(3,3), k(3,3), adj(3,3), det, w, c, roots(11), zero
    type(kernel_spec_t) :: spec
    character(8) :: name
    character(:), allocatable :: code
    integer :: a, b, r, s, t, unit
    call arena%init()
    zero = real_expr(arena, 0.0d0)
    allocate(spec%args(20), spec%outputs(11))
    t = 0
    do b = 1, 3
        do a = 1, 3
            t = t + 1
            write(name, '(a,i1,i1)') 'j', a, b
            j(a,b) = sym(arena, trim(name))
            spec%args(t) = str(trim(name))
        end do
    end do
    do b = 1, 3
        do a = 1, 3
            t = t + 1
            write(name, '(a,i1,i1)') 'k', a, b
            k(a,b) = sym(arena, trim(name))
            spec%args(t) = str(trim(name))
        end do
    end do
    w = sym(arena, 'w')
    c = sym(arena, 'c')
    spec%args(19:20) = [str('w'), str('c')]
    ! adj(a,b) is the cofactor transpose, so J maps reference to physical.
    do a = 1, 3
        do b = 1, 3
            r = modulo(b,3) + 1
            s = modulo(a,3) + 1
            adj(a,b) = j(r,s)*j(modulo(r,3)+1,modulo(s,3)+1) - &
                j(r,modulo(s,3)+1)*j(modulo(r,3)+1,s)
        end do
    end do
    det = j(1,1)*adj(1,1) + j(1,2)*adj(2,1) + j(1,3)*adj(3,1)
    t = 0
    do b = 1, 3
        do a = 1, 3
            t = t + 1
            roots(t) = zero
            do s = 1, 3
                do r = 1, 3
                    roots(t) = roots(t) + w*adj(a,r)*k(r,s)*adj(b,s)/det
                end do
            end do
            write(name, '(a,i1,i1)') 'g', a, b
            spec%outputs(t) = str(trim(name))
        end do
    end do
    roots(10) = c*w*det
    roots(11) = det
    spec%outputs(10:11) = [str('mass'), str('det')]
    spec%name = str('generated_tensor_metric')
    spec%module_name = str('fortfem_generated_tensor_metric')
    spec%mode = KERNEL_SUBROUTINE
    spec%generator = str('gen_tensor_metric_products')
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str('tools/codegen/generate.sh')
    spec%pure_procedure = .true.
    code = chars(emit_kernel(roots, spec))
    open(newunit=unit, file=generated_path('fortfem_tensor_metric_products.f90'), &
        status='replace', action='write')
    write(unit,'(a)') code(:len(code)-1)
    close(unit)
end program gen_tensor_metric_products
