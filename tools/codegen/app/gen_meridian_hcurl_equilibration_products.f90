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
    use fortsym_products, only: jvp, vjp
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

    ! Scaled potential s = sigma/R (sheet and trace loads prescribe R s):
    ! s and grad s from sigma, grad sigma, R and grad R.
    block
        type(expr_t) :: sv(6), so(3)
        type(kernel_spec_t) :: sspec
        character(len=*), parameter :: sin_names(6) = [character(len=5) :: &
            "rr", "grr1", "grr2", "sig", "gsig1", "gsig2"]
        character(len=*), parameter :: sout_names(3) = [character(len=3) :: &
            "s", "gs1", "gs2"]
        do i = 1, 6
            sv(i) = sym(arena, trim(sin_names(i)))
        end do
        so(1) = sv(4)/sv(1)
        so(2) = sv(5)/sv(1) - sv(4)*sv(2)/(sv(1)*sv(1))
        so(3) = sv(6)/sv(1) - sv(4)*sv(3)/(sv(1)*sv(1))
        do root = 1, 3
            result = engine%simplify(so(root))
            if (result%ok) so(root) = result%value
        end do
        sspec%name = str("generated_meridian_hcurl_scaled_potential")
        sspec%module_name = str("fortfem_generated_meridian_hcurl_scaled_potential")
        sspec%mode = KERNEL_SUBROUTINE
        sspec%generator = str("gen_meridian_hcurl_equilibration_products")
        sspec%generator_revision = str(fortsym_revision())
        sspec%regenerate_command = str("cd tools/codegen && ./generate.sh")
        sspec%pure_procedure = .true.
        allocate (sspec%args(6), sspec%outputs(3))
        do i = 1, 6
            sspec%args(i) = str(trim(sin_names(i)))
        end do
        do i = 1, 3
            sspec%outputs(i) = str(trim(sout_names(i)))
        end do
        code = chars(emit_kernel(so, sspec))
        write (unit, "(a)") code(:len(code) - 1)
    end block
    block
        type(expr_t) :: lv(8), lo(3), tau
        character(len=4), parameter :: names(8) = &
            ["la  ", "lb  ", "ga1 ", "ga2 ", "gb1 ", "gb2 ", "r   ", "dr  "]
        do i = 1, 8
            lv(i) = sym(arena, trim(names(i)))
        end do
        tau = lv(2)/(lv(1) + lv(2))
        call emit_family('edge_coordinate', lv(:2), names(:2), [tau])
        ! sigma = -(lambda_a + lambda_b)^4 r(tau), tau = lambda_b/(la+lb).
        ! r vanishes at the endpoints: the lift is zero on the other edges.
        ! Fourth-power attenuation keeps the oscillatory interpolation
        ! residual near its edge instead of extending it through the cell.
        lo(1) = -(lv(1) + lv(2))*(lv(1) + lv(2)) &
            *(lv(1) + lv(2))*(lv(1) + lv(2))*lv(7)
        lo(2) = -(lv(1) + lv(2))*(lv(1) + lv(2))*(lv(1) + lv(2)) &
            *(4*(lv(3) + lv(5))*lv(7) &
            + ((1 - tau)*lv(5) - tau*lv(3))*lv(8))
        lo(3) = -(lv(1) + lv(2))*(lv(1) + lv(2))*(lv(1) + lv(2)) &
            *(4*(lv(4) + lv(6))*lv(7) &
            + ((1 - tau)*lv(6) - tau*lv(4))*lv(8))
        call emit_family('edge_lift', lv, names, lo)
    end block
    close (unit)
contains
    subroutine emit_family(name, inputs, names, roots)
        character(len=*), intent(in) :: name, names(:)
        type(expr_t), intent(in) :: inputs(:), roots(:)
        type(expr_t) :: dots(size(inputs)), bars(size(roots))
        type(str_t) :: args(2*size(inputs)), adjargs(size(inputs) + size(roots))
        integer :: j
        do j = 1, size(inputs)
            args(j) = str(trim(names(j)))
            args(size(inputs) + j) = str(trim(names(j))//'_dot')
            dots(j) = sym(arena, chars(args(size(inputs) + j)))
        end do
        adjargs(:size(inputs)) = args(:size(inputs))
        do j = 1, size(roots)
            bars(j) = sym(arena, 'bar'//achar(iachar('0') + j))
            adjargs(size(inputs) + j) = str('bar'//achar(iachar('0') + j))
        end do
        call emit_one(name, roots, args(:size(inputs)))
        call emit_one(name//'_jvp', jvp(roots, inputs, dots), args)
        call emit_one(name//'_vjp', vjp(roots, inputs, bars), adjargs)
    end subroutine emit_family

    subroutine emit_one(name, roots, args)
        character(len=*), intent(in) :: name
        type(expr_t), intent(in) :: roots(:)
        type(str_t), intent(in) :: args(:)
        type(kernel_spec_t) :: local_spec
        integer :: j
        local_spec%name = str('generated_meridian_'//name)
        local_spec%module_name = str('fortfem_generated_meridian_'//name)
        local_spec%mode = KERNEL_SUBROUTINE
        local_spec%generator = str('gen_meridian_hcurl_equilibration_products')
        local_spec%generator_revision = str(fortsym_revision())
        local_spec%regenerate_command = str('cd tools/codegen && ./generate.sh')
        local_spec%pure_procedure = .true.
        local_spec%args = args
        allocate(local_spec%outputs(size(roots)))
        do j = 1, size(roots)
            local_spec%outputs(j) = str('out'//achar(iachar('0') + j))
        end do
        code = chars(emit_kernel(roots, local_spec))
        call write_wrapped(code(:len(code) - 1))
    end subroutine emit_one

    subroutine write_wrapped(source)
        character(len=*), intent(in) :: source
        character(len=:), allocatable :: line
        integer :: first, next, cut, k
        first = 1
        do while (first <= len(source))
            next = index(source(first:), new_line('a'))
            if (next == 0) next = len(source) - first + 2
            line = source(first:first + next - 2)
            do while (len_trim(line) > 84)
                if (index(adjustl(line), '!') == 1) exit
                cut = 0
                do k = 1, 80
                    if (line(k:k) == ',' .or. line(k:k) == ' ') cut = k
                end do
                if (cut <= 8) error stop 'cannot wrap generated token'
                write(unit, '(a)') trim(line(:cut))//' &'
                line = '        '//trim(adjustl(line(cut + 1:)))
            end do
            write(unit, '(a)') line
            first = first + next
        end do
    end subroutine write_wrapped
end program gen_meridian_hcurl_equilibration_products
