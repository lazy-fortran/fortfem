module meridian_fourier_curl_fields
    !! Phi-dependent tensor and polynomial-trigonometric fields with known
    !! analytic components, for the mode-coupled curl-complex forms.
    use fortfem_kinds, only: dp
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t
    implicit none
    private
    public :: identity_t, coupled_t, field_t, coupled_tensor
    public :: series, field_a, field_b, zc, ratio_k3

    real(dp), parameter :: zc = 0.5_dp
    real(dp), parameter :: ratio_k3 = 1.0_dp

    type, extends(meridian_fourier_tensor_t) :: identity_t
    contains
        procedure :: tensor => identity_tensor
    end type

    !> Symmetric positive definite tensor, Fourier degree <= 2 in phi.
    type, extends(meridian_fourier_tensor_t) :: coupled_t
    contains
        procedure :: tensor => coupled_eval
    end type

    !> A = grad_pol(R^2 Z) T1(phi) + grad(R^2 (1 + Z) T3(phi))
    !!     + e_phi R (1 + zc Z) S(phi); weights per Fourier block.
    type :: field_t
        real(dp), allocatable :: w1(:), w3(:), ws(:)
    end type

contains

    subroutine identity_tensor(self, t, x, k, status)
        class(identity_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: status
        integer :: p
        k = 0.0_dp
        do p = 1, size(x, 2)
            k(1, 1, p) = 1.0_dp; k(2, 2, p) = 1.0_dp; k(3, 3, p) = 1.0_dp
        end do
        status = 0
    end subroutine identity_tensor

    pure subroutine coupled_tensor(x, k)
        !! x = (R, Z, phi).
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: k(3, 3)
        k = 0.0_dp
        k(1, 1) = 2.0_dp + 0.2_dp*x(1) + 0.5_dp*cos(x(3))
        k(2, 2) = 1.5_dp + 0.1_dp*x(2) + 0.4_dp*cos(2*x(3))
        k(3, 3) = 1.0_dp + 0.1_dp*x(1) + 0.3_dp*sin(x(3))
        k(1, 2) = 0.3_dp*sin(x(3)); k(2, 1) = k(1, 2)
        k(1, 3) = 0.1_dp + 0.05_dp*x(1); k(3, 1) = k(1, 3)
        k(2, 3) = 0.2_dp*cos(x(3)); k(3, 2) = k(2, 3)
    end subroutine coupled_tensor

    subroutine coupled_eval(self, t, x, k, status)
        class(coupled_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: status
        integer :: p
        do p = 1, size(x, 2)
            call coupled_tensor(x(:, p), k(:, :, p))
        end do
        status = 0
    end subroutine coupled_eval

    pure real(dp) function series(w, phi, d) result(s)
        !! sum_b w(b) trig_b(phi) (d = 0) or its phi derivative (d = 1);
        !! block 1 is the constant, 2n cos(n phi), 2n + 1 sin(n phi).
        real(dp), intent(in) :: w(:), phi
        integer, intent(in) :: d
        integer :: n

        s = 0.0_dp
        if (d == 0) s = w(1)
        do n = 1, (size(w) - 1)/2
            if (d == 0) then
                s = s + w(2*n)*cos(n*phi) + w(2*n + 1)*sin(n*phi)
            else
                s = s + n*(-w(2*n)*sin(n*phi) + w(2*n + 1)*cos(n*phi))
            end if
        end do
    end function series

    pure function field_a(f, x, zcoef) result(a)
        !! Physical components (R, phi, Z) of A at x = (R, Z, phi).
        type(field_t), intent(in) :: f
        real(dp), intent(in) :: x(3), zcoef
        real(dp) :: a(3), t1, t3, dt3, s
        t1 = series(f%w1, x(3), 0)
        t3 = series(f%w3, x(3), 0)
        dt3 = series(f%w3, x(3), 1)
        s = series(f%ws, x(3), 0)
        a(1) = 2*x(1)*x(2)*t1 + 2*x(1)*(1 + x(2))*t3
        a(2) = x(1)*(1 + zcoef*x(2))*s + x(1)*(1 + x(2))*dt3
        a(3) = x(1)**2*t1 + x(1)**2*t3
    end function field_a

    pure function field_b(f, x, zcoef) result(b)
        !! Physical components (R, phi, Z) of curl A.
        type(field_t), intent(in) :: f
        real(dp), intent(in) :: x(3), zcoef
        real(dp) :: b(3), dt1, s
        dt1 = series(f%w1, x(3), 1)
        s = series(f%ws, x(3), 0)
        b(1) = x(1)*(dt1 - zcoef*s)
        b(2) = 0.0_dp
        b(3) = 2*(1 + zcoef*x(2))*s - 2*x(2)*dt1
    end function field_b

end module meridian_fourier_curl_fields

program test_meridian_fourier_curl
    !! Mode-coupled real-Fourier forms on the curl complex:
    !!  1. a phi-independent coefficient reproduces the per-n forms exactly
    !!     (block diagonal; half weight for cos/sin);
    !!  2. a phi-dependent, spatially varying tensor gives the exact bilinear
    !!     values on polynomial-trigonometric fields (mass on V1 and V2,
    !!     curl-curl, V1 x V0 gradient coupling) and couples only modes
    !!     within the Fourier content of the coefficient;
    !!  3. symmetry and positive definiteness (semi-definiteness of the
    !!     curl-curl form);
    !!  4. curl grad = 0 and curl-curl grad = 0 through the coupled assembly.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, &
        meridian_mesh_from_map
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_lagrange_nodes
    use fortfem_meridian_fourier_curl, only: meridian_fourier_mass, &
        meridian_fourier_curl_curl, meridian_fourier_gradient_coupling, &
        meridian_fourier_gradient, meridian_fourier_curl, &
        meridian_fourier_identity_forms
    use fortsparse, only: csc_t, csc_matvec
    use meridian_fourier_curl_fields, only: identity_t, coupled_t, field_t, &
        coupled_tensor, series, field_a, field_b, zc
    use meridian_test_meshes, only: make_warped_rectangle_mesh
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    logical :: all_passed

    all_passed = .true.
    call test_reduction_symmetry_derham()
    call test_exact_forms()
    call test_identity_forms()
    call check_summary("Meridian real Fourier curl-complex forms")
    if (.not. all_passed) error stop 1

contains

    subroutine test_identity_forms()
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(identity_t) :: iso
        type(csc_t) :: c, s, mg, reference, g
        real(dp), allocatable :: exact(:, :), d(:, :), v(:), cv(:)
        integer :: st, n, k, j
        call make_warped_rectangle_mesh(2, 2, mm)
        do k = 2, 3
            call meridian_feec_build(mm, k, fs, st)
            call record(st == 0, "identity forms complex builds")
            do n = 0, 2
                call meridian_fourier_identity_forms(fs, n, c, s, mg, st)
                call record(st == 0, "identity Fourier forms assemble")
                call meridian_fourier_curl_curl(fs, n, 4*n + 1, iso, reference, st)
                call record(st == 0, "quadrature identity curl-curl assembles")
                exact = to_dense(reference); d = to_dense(s)
                call record(maxval(abs(d - exact)) <= &
                    1.0e-13_dp*maxval(abs(exact)), "analytic identity curl-curl")
                call meridian_fourier_gradient_coupling(fs, n, 4*n + 1, iso, &
                    reference, st)
                call record(st == 0, "quadrature identity gradient assembles")
                exact = to_dense(reference); d = to_dense(mg)
                call record(maxval(abs(d - exact)) <= &
                    1.0e-13_dp*maxval(abs(exact)), "analytic identity gradient")
                call meridian_fourier_gradient(fs, n, g, st)
                v = [(sin(0.37_dp*j), j=1, g%ncol)]
                cv = csc_matvec(c, csc_matvec(g, v))
                call record(norm2(cv) < 1.0e-10_dp*norm2(v), &
                    "identity forms retain curl grad nullspace")
            end do
        end do
    end subroutine

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    function to_dense(a) result(d)
        type(csc_t), intent(in) :: a
        real(dp) :: d(a%nrow, a%ncol)
        integer :: j, p
        d = 0.0_dp
        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                d(a%row_idx(p), j) = d(a%row_idx(p), j) + a%val(p)
            end do
        end do
    end function to_dense

    function block_diag(m, modes) result(d)
        !! Phi-mean mass of K = I: M, then M/2 for every cos and sin block.
        real(dp), intent(in) :: m(:, :)
        integer, intent(in) :: modes
        real(dp), allocatable :: d(:, :)
        real(dp) :: f
        integer :: b, nm

        nm = size(m, 1)
        allocate (d(nm*(2*modes + 1), nm*(2*modes + 1)))
        d = 0.0_dp
        do b = 0, 2*modes
            f = merge(1.0_dp, 0.5_dp, b == 0)
            d(b*nm + 1:(b + 1)*nm, b*nm + 1:(b + 1)*nm) = f*m
        end do
    end function block_diag

    function block_op(a, b, modes) result(d)
        !! Dense real-block operator of a + i n b.
        real(dp), intent(in) :: a(:, :), b(:, :)
        integer, intent(in) :: modes
        real(dp), allocatable :: d(:, :)
        integer :: n, nr, nc

        nr = size(a, 1); nc = size(a, 2)
        allocate (d(nr*(2*modes + 1), nc*(2*modes + 1)))
        d = 0.0_dp
        d(1:nr, 1:nc) = a
        do n = 1, modes
            d((2*n - 1)*nr + 1:2*n*nr, (2*n - 1)*nc + 1:2*n*nc) = a
            d(2*n*nr + 1:(2*n + 1)*nr, 2*n*nc + 1:(2*n + 1)*nc) = a
            d((2*n - 1)*nr + 1:2*n*nr, 2*n*nc + 1:(2*n + 1)*nc) = n*b
            d(2*n*nr + 1:(2*n + 1)*nr, (2*n - 1)*nc + 1:2*n*nc) = -n*b
        end do
    end function block_op

    logical function cholesky_ok(m) result(ok)
        real(dp), intent(in) :: m(:, :)
        real(dp), allocatable :: l(:, :)
        integer :: i, j, n
        real(dp) :: s

        n = size(m, 1)
        allocate (l(n, n))
        l = 0.0_dp
        ok = .false.
        do j = 1, n
            s = m(j, j) - dot_product(l(j, :j - 1), l(j, :j - 1))
            if (s <= 0.0_dp) return
            l(j, j) = sqrt(s)
            do i = j + 1, n
                l(i, j) = (m(i, j) - dot_product(l(i, :j - 1), l(j, :j - 1)))/l(j, j)
            end do
        end do
        ok = .true.
    end function cholesky_ok

    subroutine test_reduction_symmetry_derham()
        integer, parameter :: modes = 3, nphi = 12
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(identity_t) :: iso
        type(coupled_t) :: cpd
        type(csc_t) :: m1, m2, s, cpl, g, c
        real(dp), allocatable :: dm1(:, :), dm2(:, :), ds(:, :), dcpl(:, :)
        real(dp), allocatable :: dg(:, :), dc(:, :), da(:, :), db(:, :), ex(:, :)
        real(dp), allocatable :: lam(:), gl(:), x(:), sx(:)
        real(dp) :: err, scale
        integer :: st, n0, n1, n2, nb, i, trial, b, b2, nm
        logical :: pos

        call make_warped_rectangle_mesh(3, 2, mm)
        call meridian_feec_build(mm, 2, fs, st)
        call record(st == 0, "build P2 complex on the curved mesh")
        n0 = fs%n0
        n1 = fs%sp(1)%ndof
        n2 = fs%sp(2)%ndof
        nb = 2*modes + 1

        ! 1. phi-independent coefficient: per-n forms, block diagonal.
        call meridian_fourier_mass(fs, 1, modes, nphi, iso, m1, st)
        call record(st == 0, "V1 mass assembly")
        call meridian_fourier_mass(fs, 2, modes, nphi, iso, m2, st)
        call record(st == 0, "V2 mass assembly")
        dm1 = to_dense(m1)
        dm2 = to_dense(m2)
        ex = block_diag(to_dense(fs%mass(1)), modes)
        scale = maxval(abs(ex))
        err = maxval(abs(dm1 - ex))/scale
        print '(a, es10.3)', "  K = I V1 mass defect: ", err
        call record(err <= 1.0e-14_dp, "V1 mass equals block-diagonal per-n mass")
        ex = block_diag(to_dense(fs%mass(2)), modes)
        scale = maxval(abs(ex))
        err = maxval(abs(dm2 - ex))/scale
        print '(a, es10.3)', "  K = I V2 mass defect: ", err
        call record(err <= 1.0e-14_dp, "V2 mass equals block-diagonal per-n mass")

        da = to_dense(fs%curl%a)
        db = to_dense(fs%curl%b)
        dc = block_op(da, db, modes)
        call meridian_fourier_curl(fs, modes, c, st)
        err = maxval(abs(to_dense(c) - dc))
        call record(st == 0 .and. err == 0.0_dp, "block curl operator")
        ex = matmul(transpose(dc), matmul(block_diag(to_dense(fs%mass(2)), &
            modes), dc))
        call meridian_fourier_curl_curl(fs, modes, nphi, iso, s, st)
        call record(st == 0, "curl-curl assembly (K = I)")
        ds = to_dense(s)
        scale = maxval(abs(ex))
        err = maxval(abs(ds - ex))/scale
        print '(a, es10.3)', "  K = I curl-curl defect: ", err
        ! C^T M C sums large cancelling entries: a few 1e-14 of roundoff.
        call record(err <= 1.0e-13_dp, "curl-curl equals per-n forms")

        dg = block_op(to_dense(fs%grad%a), to_dense(fs%grad%b), modes)
        ex = matmul(block_diag(to_dense(fs%mass(1)), modes), dg)
        call meridian_fourier_gradient_coupling(fs, modes, nphi, iso, cpl, st)
        call record(st == 0, "gradient coupling assembly (K = I)")
        dcpl = to_dense(cpl)
        scale = maxval(abs(ex))
        err = maxval(abs(dcpl - ex))/scale
        print '(a, es10.3)', "  K = I gradient coupling defect: ", err
        call record(err <= 1.0e-14_dp, "gradient coupling equals per-n forms")

        ! 3. symmetry and positivity, phi-dependent tensor.
        call meridian_fourier_mass(fs, 1, modes, nphi, cpd, m1, st)
        call meridian_fourier_mass(fs, 2, modes, nphi, cpd, m2, st)
        call meridian_fourier_curl_curl(fs, modes, nphi, cpd, s, st)
        call record(st == 0, "coupled assemblies")
        dm1 = to_dense(m1); dm2 = to_dense(m2); ds = to_dense(s)
        err = maxval(abs(dm1 - transpose(dm1)))/maxval(abs(dm1))
        err = max(err, maxval(abs(dm2 - transpose(dm2)))/maxval(abs(dm2)))
        err = max(err, maxval(abs(ds - transpose(ds)))/maxval(abs(ds)))
        print '(a, es10.3)', "  asymmetry: ", err
        call record(err <= 1.0e-14_dp, "mass and curl-curl blocks are symmetric")
        call record(cholesky_ok(dm1), "V1 mass is positive definite")
        call record(cholesky_ok(dm2), "V2 mass is positive definite")
        allocate (x(n1*nb))
        pos = .true.
        do trial = 1, 20
            do i = 1, n1*nb
                x(i) = sin(0.37_dp*i*trial) + 0.3_dp*cos(1.1_dp*i + trial)
            end do
            pos = pos .and. dot_product(x, matmul(ds, x)) > 0.0_dp
        end do
        call record(pos, "curl-curl form is positive on generic vectors")

        ! Coupling only through the Fourier content (degree 2) of K.
        nm = n1
        err = 0.0_dp
        do b2 = 0, nb - 1
            do b = 0, nb - 1
                if (abs((b + 1)/2 - (b2 + 1)/2) > 2 .and. &
                    (b + 1)/2 + (b2 + 1)/2 > 2) then
                    err = max(err, maxval(abs(dm1(b*nm + 1:(b + 1)*nm, &
                        b2*nm + 1:(b2 + 1)*nm))))
                end if
            end do
        end do
        print '(a, es10.3)', "  out-of-band mode coupling: ", err
        call record(err <= 1.0e-14_dp*maxval(abs(dm1)), &
            "modes couple only within the coefficient's Fourier content")
        call record(maxval(abs(dm1(nm + 1:2*nm, 1:nm))) > 1.0e-3_dp, &
            "phi dependence couples neighbouring modes")

        ! 4. de Rham consistency through the coupled assembly.
        call meridian_fourier_gradient(fs, modes, g, st)
        call record(st == 0, "block gradient operator")
        dg = to_dense(g)
        err = maxval(abs(matmul(dc, dg)))/(maxval(abs(dc))*maxval(abs(dg)))
        print '(a, es10.3)', "  curl grad: ", err
        call record(err <= 1.0e-14_dp, "curl_F grad_F = 0")
        allocate (lam(n0*nb), gl(n1*nb), sx(n1*nb))
        do i = 1, n0*nb
            lam(i) = sin(0.53_dp*i) + 0.2_dp*cos(0.91_dp*i*i)
        end do
        gl = csc_matvec(g, lam)
        sx = csc_matvec(s, gl)
        do i = 1, n1*nb
            x(i) = sin(0.37_dp*i) + 0.3_dp*cos(1.1_dp*i + 1)
        end do
        err = norm2(sx)/(norm2(gl)*maxval(abs(ds)))
        print '(a, es10.3)', "  curl-curl grad: ", err
        call record(err <= 1.0e-13_dp .and. norm2(gl) > 0.0_dp, &
            "curl-curl annihilates gradients (phi-dependent K)")
        call meridian_fourier_gradient_coupling(fs, modes, nphi, cpd, cpl, st)
        sx = csc_matvec(cpl, lam)
        call record(norm2(sx) > 1.0e-3_dp*norm2(lam), &
            "gradient coupling is nonzero (phi-dependent K)")
    end subroutine test_reduction_symmetry_derham

    subroutine test_exact_forms()
        !! Exact bilinear values on the rectangle [2, 3] x [-1/2, 1/2].
        integer, parameter :: modes = 3, nphi = 16, nb = 2*modes + 1
        type(mesh_2d_t) :: pmesh, topo
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(coupled_t) :: cpd
        type(csc_t) :: m1, m2, s, cpl, g, c
        type(field_t) :: fu, fv
        real(dp), allocatable :: pv(:, :, :), xy(:, :), u(:), v(:), bu(:), bv(:)
        real(dp), allocatable :: lam(:), wl(:)
        real(dp) :: ref, val, scale
        integer :: st, t, i, n0, n1

        call pmesh%create_rectangular(5, 4, 2.0_dp, 3.0_dp, -0.5_dp, 0.5_dp)
        topo = pmesh
        allocate (pv(2, 3, topo%n_triangles))
        do t = 1, topo%n_triangles
            pv(:, :, t) = pmesh%vertices(:, pmesh%triangles(:, t))
        end do
        call meridian_mesh_from_map(topo, pv, 1, mm, st, ident_map)
        call record(st == 0, "build rectangle mesh")
        call meridian_feec_build(mm, 3, fs, st)
        call record(st == 0, "build P3 complex")
        n0 = fs%n0
        n1 = fs%sp(1)%ndof
        allocate (fu%w1(nb), fu%w3(nb), fu%ws(nb), fv%w1(nb), fv%w3(nb), &
            fv%ws(nb))
        do i = 1, nb
            fu%w1(i) = 0.8_dp*sin(1.3_dp*i) + 0.2_dp
            fu%w3(i) = cos(0.7_dp*i*i) - 0.1_dp
            fu%ws(i) = 0.5_dp*sin(2.1_dp*i + 0.4_dp)
            fv%w1(i) = cos(0.9_dp*i) + 0.3_dp
            fv%w3(i) = 0.6_dp*sin(0.5_dp*i*i + 1.0_dp)
            fv%ws(i) = cos(1.7_dp*i) - 0.2_dp
        end do
        u = field_dofs(fs, fu)
        v = field_dofs(fs, fv)

        call meridian_fourier_mass(fs, 1, modes, nphi, cpd, m1, st)
        call record(st == 0, "coupled V1 mass assembly")
        val = dot_product(v, csc_matvec(m1, u))
        ref = oracle(fv, fu, 1, zc)
        scale = abs(ref)
        print '(a, 2es22.14)', "  V1 mass value, exact: ", val, ref
        call record(abs(val - ref) <= 1.0e-12_dp*scale, &
            "V1 mass exact on polynomial-trigonometric fields")

        call meridian_fourier_curl(fs, modes, c, st)
        bu = csc_matvec(c, u)
        bv = csc_matvec(c, v)
        call meridian_fourier_mass(fs, 2, modes, nphi, cpd, m2, st)
        call record(st == 0, "coupled V2 mass assembly")
        val = dot_product(bv, csc_matvec(m2, bu))
        ref = oracle(fv, fu, 2, zc)
        print '(a, 2es22.14)', "  V2 mass value, exact: ", val, ref
        call record(abs(val - ref) <= 1.0e-12_dp*abs(ref), &
            "V2 mass exact on curl of polynomial-trigonometric fields")

        call meridian_fourier_curl_curl(fs, modes, nphi, cpd, s, st)
        call record(st == 0, "coupled curl-curl assembly")
        val = dot_product(v, csc_matvec(s, u))
        print '(a, 2es22.14)', "  curl-curl value, exact: ", val, ref
        call record(abs(val - ref) <= 1.0e-11_dp*abs(ref), &
            "curl-curl exact on polynomial-trigonometric fields")

        call meridian_fourier_gradient_coupling(fs, modes, nphi, cpd, cpl, st)
        call record(st == 0, "coupled gradient coupling assembly")
        call meridian_lagrange_nodes(fs, xy)
        allocate (lam(n0*nb), wl(nb))
        do i = 1, nb
            wl(i) = sin(0.8_dp*i + 0.3_dp) + 0.1_dp*i
        end do
        do i = 0, nb - 1
            lam(i*n0 + 1:(i + 1)*n0) = wl(i + 1)*xy(1, :)*xy(2, :)
        end do
        val = dot_product(v, csc_matvec(cpl, lam))
        ref = oracle_grad(fv, wl)
        print '(a, 2es22.14)', "  gradient coupling value, exact: ", val, ref
        call record(abs(val - ref) <= 1.0e-12_dp*abs(ref), &
            "gradient coupling exact on polynomial-trigonometric fields")
    end subroutine test_exact_forms

    subroutine ident_map(p, x)
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: x(2)
        x = p
    end subroutine ident_map

    function field_dofs(fs, f) result(x)
        !! Fourier-V1 dofs of the field of meridian_fourier_curl_fields.
        type(meridian_feec_t), intent(in) :: fs
        type(field_t), intent(in) :: f
        real(dp), allocatable :: x(:)
        type(csc_t) :: g
        real(dp), allocatable :: xy(:, :), l1(:), lam(:), gl(:)
        integer :: n0, n1, nn, b, st, nb, i

        nb = size(f%w1)
        n0 = fs%n0
        n1 = fs%sp(1)%ndof
        nn = n1 - n0
        call meridian_lagrange_nodes(fs, xy)
        allocate (x(n1*nb), lam(n0*nb))
        x = 0.0_dp
        ! grad_pol(R^2 Z) T1: Nedelec part only (the a rows of grad%a vanish).
        l1 = xy(1, :)**2*xy(2, :)
        gl = csc_matvec(fs%grad%a, l1)
        do b = 0, nb - 1
            x(b*n1 + 1:b*n1 + nn) = f%w1(b + 1)*gl(1:nn)
            ! a = R A_phi = R^2 (1 + zc Z) S
            do i = 1, n0
                x(b*n1 + nn + i) = f%ws(b + 1)*xy(1, i)**2*(1 + zc*xy(2, i))
            end do
            lam(b*n0 + 1:(b + 1)*n0) = f%w3(b + 1)*xy(1, :)**2*(1 + xy(2, :))
        end do
        call meridian_fourier_gradient(fs, (nb - 1)/2, g, st)
        x = x + csc_matvec(g, lam)
    end function field_dofs

    real(dp) function oracle(fv, fu, space, zcoef) result(r)
        !! int <v . K u> over [2,3]x[-1/2,1/2] with A (space 1) or curl A
        !! (space 2): Gauss-Legendre in (R, Z), uniform phi mean.
        type(field_t), intent(in) :: fv, fu
        integer, intent(in) :: space
        real(dp), intent(in) :: zcoef
        integer, parameter :: ng = 8, nph = 64
        real(dp) :: xg(ng), wg(ng), x(3), k(3, 3), a(3), b(3)
        integer :: i, j, l

        call gauss_legendre(ng, xg, wg)
        r = 0.0_dp
        do i = 1, ng
            do j = 1, ng
                do l = 0, nph - 1
                    x = [2.5_dp + 0.5_dp*xg(i), 0.5_dp*xg(j), 2*pi*l/nph]
                    call coupled_tensor(x, k)
                    if (space == 1) then
                        a = field_a(fu, x, zcoef)
                        b = field_a(fv, x, zcoef)
                    else
                        a = field_b(fu, x, zcoef)
                        b = field_b(fv, x, zcoef)
                    end if
                    r = r + 0.25_dp*wg(i)*wg(j)*x(1)/nph*dot_product(b, &
                        matmul(k, a))
                end do
            end do
        end do
    end function oracle

    real(dp) function oracle_grad(fv, wl) result(r)
        !! int <A_v . K grad lambda>, lambda = R Z W(phi).
        type(field_t), intent(in) :: fv
        real(dp), intent(in) :: wl(:)
        integer, parameter :: ng = 8, nph = 64
        real(dp) :: xg(ng), wg(ng), x(3), k(3, 3), gr(3)
        integer :: i, j, l

        call gauss_legendre(ng, xg, wg)
        r = 0.0_dp
        do i = 1, ng
            do j = 1, ng
                do l = 0, nph - 1
                    x = [2.5_dp + 0.5_dp*xg(i), 0.5_dp*xg(j), 2*pi*l/nph]
                    call coupled_tensor(x, k)
                    gr = [x(2)*series(wl, x(3), 0), x(2)*series(wl, x(3), 1), &
                        x(1)*series(wl, x(3), 0)]
                    r = r + 0.25_dp*wg(i)*wg(j)*x(1)/nph*dot_product( &
                        field_a(fv, x, zc), matmul(k, gr))
                end do
            end do
        end do
    end function oracle_grad

    subroutine gauss_legendre(n, x, w)
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        real(dp) :: z, p1, p2, p3, pp
        integer :: i, j, it

        do i = 1, n
            z = cos(pi*(i - 0.25_dp)/(n + 0.5_dp))
            do it = 1, 100
                p1 = 1.0_dp
                p2 = 0.0_dp
                do j = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = ((2*j - 1)*z*p2 - (j - 1)*p3)/j
                end do
                pp = n*(z*p1 - p2)/(z*z - 1.0_dp)
                z = z - p1/pp
            end do
            x(i) = z
            w(i) = 2.0_dp/((1.0_dp - z*z)*pp*pp)
        end do
    end subroutine gauss_legendre

end program test_meridian_fourier_curl
