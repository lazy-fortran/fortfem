program test_dg_curved_2d
    !! DG upwind advection, SIPG diffusion and their anisotropic complex
    !! combination on curved isoparametric triangle meshes.  Oracles: a
    !! manufactured solution (L2 convergence orders), the divergence theorem
    !! and the telescoping of interior fluxes (conservation), and the
    !! symmetry of the SIPG form.
    use fortfem_kinds, only: dp
    use fortfem_dg_curved_2d, only: DG_FORM_ADVECTIVE, DG_FORM_CONSERVATIVE, &
        dg_curved_assemble, dg_curved_qpoint, dg_curved_rhs, &
        dg_curved_space_init, dg_curved_space_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use dg_curved_test_problems, only: anisotropic_tensor, l2_error, &
        mms_coef_t, mms_data_t
    use meridian_test_meshes, only: make_annulus_mesh, make_warped_rectangle_mesh
    use fortfem_sparse_direct, only: sparse_direct_solve_once
    use fortsparse, only: csc_is_valid, csc_matvec, csc_z_t, fortsparse_status_t, &
        status_ok
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    logical :: all_passed
    type(mms_coef_t) :: coef
    complex(dp) :: eye(2, 2)
    integer :: k

    all_passed = .true.
    eye = reshape([one, (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), one], [2, 2])

    print '(a)', "Upwind advection-reaction, conservative form, warped mesh"
    coef%kk = (0.0_dp, 0.0_dp)
    coef%c = (1.0_dp, 0.5_dp)
    do k = 1, 3
        call convergence(coef, DG_FORM_CONSERVATIVE, .false., k, &
            real(k, dp) + 0.5_dp)
    end do
    print '(a)', "Upwind advection-reaction, advective form, warped mesh"
    call convergence(coef, DG_FORM_ADVECTIVE, .false., 2, 2.5_dp)

    print '(a)', "SIPG diffusion, complex isotropic K, curved annulus"
    coef%kk = (1.0_dp, 0.3_dp)*eye
    coef%c = (0.0_dp, 0.0_dp)
    coef%bscale = 0.0_dp
    do k = 1, 3
        call convergence(coef, DG_FORM_CONSERVATIVE, .true., k, real(k, dp) + 0.8_dp)
    end do

    print '(a)', "Advection-diffusion, anisotropic complex K (eps = 1e-4)"
    coef%kk = anisotropic_tensor((1.0_dp, 0.2_dp), 0.5_dp, 1.0e-4_dp)
    coef%c = (0.5_dp, 1.0_dp)
    coef%bscale = 1.0_dp
    do k = 1, 3
        call convergence(coef, DG_FORM_CONSERVATIVE, .true., k, &
            real(k, dp) + 0.5_dp)
    end do

    call check_conservation()
    call check_sipg_symmetry()

    if (.not. all_passed) error stop 1
    print '(a)', "All curved DG tests passed"

contains

    subroutine build_mesh(annulus, level, k, mm)
        logical, intent(in) :: annulus
        integer, intent(in) :: level, k
        type(meridian_mesh_t), intent(out) :: mm

        if (annulus) then
            call make_annulus_mesh(2*level, 8*level, 0.5_dp, 1.0_dp, 3.0_dp, &
                max(k, 2), mm)
        else
            call make_warped_rectangle_mesh(4*level, max(k, 2), mm)
        end if
    end subroutine build_mesh

    subroutine solve(space, c, form, u, ok)
        type(dg_curved_space_t), intent(in) :: space
        type(mms_coef_t), intent(in) :: c
        integer, intent(in) :: form
        complex(dp), allocatable, intent(out) :: u(:)
        logical, intent(out) :: ok

        type(csc_z_t) :: a
        type(mms_data_t) :: data
        type(fortsparse_status_t) :: fst
        complex(dp), allocatable :: rhs(:)
        integer :: st

        data%coef = c
        data%form = form
        allocate (rhs(space%ndof), u(space%ndof))
        call dg_curved_assemble(space, c, form, a, st)
        ok = st == 0
        if (.not. ok) return
        ok = csc_is_valid(a)
        if (.not. ok) return
        call dg_curved_rhs(space, c, data, form, rhs, st)
        ok = st == 0
        if (.not. ok) return
        call sparse_direct_solve_once(a, rhs, u, fst)
        ok = status_ok(fst)
    end subroutine solve

    subroutine convergence(c, form, annulus, k, min_rate)
        type(mms_coef_t), intent(in) :: c
        integer, intent(in) :: form, k
        logical, intent(in) :: annulus
        real(dp), intent(in) :: min_rate

        type(meridian_mesh_t) :: mm
        type(dg_curved_space_t) :: space
        complex(dp), allocatable :: u(:)
        real(dp) :: err(3), rate
        integer :: level, st
        logical :: ok

        do level = 1, 3
            call build_mesh(annulus, 2**(level - 1), k, mm)
            call dg_curved_space_init(mm, k, space, st)
            if (st /= 0) error stop "space init failed"
            call solve(space, c, form, u, ok)
            if (.not. ok) error stop "assembly or solve failed"
            err(level) = l2_error(space, u)
        end do
        rate = log(err(2)/err(3))/log(2.0_dp)
        print '(a, i0, a, 3es11.3, a, f6.2, a, f5.2, a)', "  P", k, " L2 errors", &
            err, "  rate ", rate, " (min ", min_rate, ")"
        if (rate < min_rate) then
            print '(a)', "  FAILED: convergence rate below the bound"
            all_passed = .false.
        end if
    end subroutine convergence

    subroutine check_conservation()
        !! Conservative form: (i) interior-element columns of A sum to zero
        !! for c = 0 (interior fluxes telescope), while those of the
        !! advective form are -int (div b) phi_j; (ii) for u = 1 and
        !! g = 1, sum(A 1 - F) = int_dOmega b.n = 0.1 |Omega_h| (divergence
        !! theorem), and |Omega_h| is the annulus area to geometry accuracy.
        type(meridian_mesh_t) :: mm
        type(dg_curved_space_t) :: space
        type(mms_coef_t) :: c
        type(mms_data_t) :: data
        type(csc_z_t) :: a
        complex(dp), allocatable :: ones(:), au(:), rhs(:)
        logical, allocatable :: interior(:)
        real(dp) :: defect(2), scale, balance, area, exact_area, target
        real(dp), allocatable :: mass(:)
        integer :: form, st, f, t, j

        print '(a)', "Conservation on the curved annulus (P2, geometry P3)"
        call make_annulus_mesh(3, 24, 0.5_dp, 1.0_dp, 3.0_dp, 3, mm)
        call dg_curved_space_init(mm, 2, space, st)
        if (st /= 0) error stop "space init failed"
        allocate (interior(space%nt))
        interior = .true.
        do f = space%nint + 1, space%nface
            interior(space%face_t(1, f)) = .false.
        end do
        c%kk = anisotropic_tensor((1.0_dp, 0.4_dp), 0.3_dp, 1.0e-3_dp)
        c%c = (0.0_dp, 0.0_dp)
        call basis_integrals(space, mass)
        do form = DG_FORM_CONSERVATIVE, DG_FORM_ADVECTIVE
            call dg_curved_assemble(space, c, form, a, st)
            if (st /= 0) error stop "assembly failed"
            ! Column sums: 0 (conservative); -int (div b) phi_j (advective).
            target = 0.0_dp
            if (form == DG_FORM_ADVECTIVE) target = -0.1_dp
            defect(form) = 0.0_dp
            scale = maxval(abs(a%val))
            do j = 1, space%ndof
                t = (j - 1)/space%nl + 1
                if (.not. interior(t)) cycle
                defect(form) = max(defect(form), abs(sum(a%val(a%col_ptr(j): &
                    a%col_ptr(j + 1) - 1)) - target*mass(j))/scale)
            end do
        end do
        print '(a, 2es11.3)', "  interior column-sum defects (cons., adv.):", &
            defect
        if (maxval(defect) > 1.0e-12_dp) then
            print '(a)', "  FAILED: interior fluxes do not telescope"
            all_passed = .false.
        end if

        c%kk = (0.0_dp, 0.0_dp)
        data%coef = c
        data%unit = .true.
        call dg_curved_assemble(space, c, DG_FORM_CONSERVATIVE, a, st)
        if (st /= 0) error stop "assembly failed"
        allocate (ones(space%ndof), au(space%ndof), rhs(space%ndof))
        ones = (1.0_dp, 0.0_dp)
        call dg_curved_rhs(space, c, data, DG_FORM_CONSERVATIVE, rhs, st)
        if (st /= 0) error stop "rhs failed"
        au = csc_matvec(a, ones)
        area = sum(space%area)
        balance = abs(sum(au - rhs) - 0.1_dp*area)
        exact_area = pi*(1.0_dp**2 - 0.5_dp**2)
        print '(a, es11.3, a, es11.3)', "  |sum(A 1 - F) - 0.1 |Omega_h||:", &
            balance, "  area defect:", abs(area - exact_area)/exact_area
        if (balance > 1.0e-12_dp .or. abs(area - exact_area) > 1.0e-5_dp) then
            print '(a)', "  FAILED: global balance"
            all_passed = .false.
        end if
    end subroutine check_conservation

    subroutine basis_integrals(space, mass)
        !! mass(j) = int phi_j dx.
        type(dg_curved_space_t), intent(in) :: space
        real(dp), allocatable, intent(out) :: mass(:)

        real(dp) :: x(2), w, val(space%nl), grad(2, space%nl)
        integer :: t, q, st, base

        allocate (mass(space%ndof))
        mass = 0.0_dp
        do t = 1, space%nt
            base = (t - 1)*space%nl
            do q = 1, space%nq
                call dg_curved_qpoint(space, t, q, x, w, val, grad, st)
                if (st /= 0) error stop "bad geometry"
                mass(base + 1:base + space%nl) = mass(base + 1:base + space%nl) &
                    + w*val
            end do
        end do
    end subroutine basis_integrals

    subroutine check_sipg_symmetry()
        !! For a real symmetric K and b = 0 the SIPG matrix is symmetric.
        type(meridian_mesh_t) :: mm
        type(dg_curved_space_t) :: space
        type(mms_coef_t) :: c
        type(csc_z_t) :: a
        complex(dp), allocatable :: x(:), y(:), ax(:), ay(:)
        real(dp) :: asym
        integer :: st, i

        print '(a)', "SIPG symmetry (P3 on the curved annulus)"
        call make_annulus_mesh(2, 16, 0.5_dp, 1.0_dp, 3.0_dp, 3, mm)
        call dg_curved_space_init(mm, 3, space, st)
        if (st /= 0) error stop "space init failed"
        c%kk = anisotropic_tensor(one, 0.7_dp, 0.05_dp)
        c%bscale = 0.0_dp
        call dg_curved_assemble(space, c, DG_FORM_CONSERVATIVE, a, st)
        if (st /= 0) error stop "assembly failed"
        allocate (x(space%ndof), y(space%ndof), ax(space%ndof), ay(space%ndof))
        do i = 1, space%ndof
            x(i) = cmplx(sin(0.37_dp*i), cos(1.3_dp*i), dp)
            y(i) = cmplx(cos(0.11_dp*i), sin(2.1_dp*i), dp)
        end do
        ax = csc_matvec(a, x)
        ay = csc_matvec(a, y)
        asym = abs(sum(y*ax) - sum(x*ay))/abs(sum(y*ax))
        print '(a, es11.3)', "  |y.Ax - x.Ay| / |y.Ax|:", asym
        if (asym > 1.0e-12_dp) then
            print '(a)', "  FAILED: SIPG form is not symmetric"
            all_passed = .false.
        end if
    end subroutine check_sipg_symmetry

end program test_dg_curved_2d
