module fortfem_meridian_feec_n
    !! Per-toroidal-mode FEEC complex on a curved meridional (R, Z) mesh.
    !!
    !! A field is u(R, Z) exp(i n phi) with complex amplitudes.  All spaces
    !! carry the cylindrical measure R dR dZ (the factor 2 pi of the phi
    !! integral is dropped).  With physical components (R, phi, Z):
    !!
    !!   V0  f                         Lagrange P_k
    !!   V1  A = (A_pol, a)            A_pol Nedelec N_k (covariant Piola),
    !!                                 a = R A_phi Lagrange P_k
    !!   V2  B = (q, s)                q = R B_pol Raviart-Thomas RT_k
    !!                                 (contravariant Piola, the rotated N_k),
    !!                                 s = B_phi in DG_{k-1} (density Piola)
    !!   V3  g = R div B               DG_{k-1} (density Piola)
    !!
    !! Operators (rot90 v = (-v_Z, v_R), rot A = d_R A_Z - d_Z A_R):
    !!
    !!   grad_n f = (grad f, i n f)                       [A_pol, a]
    !!   curl_n A = (rot90(grad a - i n A_pol), -rot A_pol)   [q, s]
    !!   div_n B  = div q + i n s                            [g]
    !!
    !! These are the physical operators: (grad f)_phi = i n f / R,
    !! B_pol = q / R = (1/R) rot90(grad a - i n A_pol) and B_phi = d_Z A_R -
    !! d_R A_Z = -rot A_pol, div B = g / R.  In these variables the operators
    !! are metric free: each is a + i n b with real sparse a, b, so
    !! curl_n grad_n = 0 and div_n curl_n = 0 hold as exact matrix identities
    !! for every n (they reduce to the 2D identities rot90 grad = grad^perp
    !! with rot grad = 0 and div rot90 = -rot).  The metric (1/R and R
    !! factors, curved Jacobians) lives only in the mass matrices
    !! M_j = (phi_i, phi_j)_{R dR dZ} of the physical fields, which are real
    !! and independent of n.
    !!
    !! Correspondence with Albert, Lainer and Biro, Comput. Phys. Commun. 277
    !! (2022) 108401, doi:10.1016/j.cpc.2022.108401 (the paper), for the
    !! right-handed coordinates (x1, x2, x3) = (Z, R, phi) and Fourier
    !! convention f = sum_n f_n exp(i n x3) (paper Eq. (31)):
    !!   grad_n  = Eq. (34): grad_t U + i n U e^3 (covariant a = A_3 = R A_phi);
    !!   curl_n  = Eq. (36) with Eqs. (38)-(39): B^3 = curl_t A, b = -i n
    !!             eps_t A + curl_t A_3 (q = R B_pol is the contravariant
    !!             density, s = B_phi = B^3 the scalar density);
    !!   div_n   = Eq. (35): div_t b + i n B^3.
    !! With nu = 1 the curl-curl form K = curl_n^H M_2 curl_n equals the
    !! weak form Eq. (69) in the gauge A_3 = 0 (Eq. (37)): nu_33 = nu R,
    !! bar nu_t = nu / R (Eqs. (28), (44)), i.e. int R |B_phi|^2 + (n^2 / R)
    !! |A_pol|^2 for the unweighted measure dR dZ, and the same form as Eq.
    !! (11) of Lainer et al., J. Phys.: Conf. Ser. 2397 (2022) 012014,
    !! doi:10.1088/1742-6596/2397/1/012014 (MEPHIT).  The n = 0 Poisson form
    !! Eq. (47), -div_t(bar nu grad_t A_3) = J^3, is meridian_gs_stiffness
    !! (bar nu = nu / R).  Spaces: A in H(curl_t) by Nedelec elements, b and
    !! J in H(div_t) by Raviart-Thomas elements (paper Sec. 2.5).  The
    !! mode-by-mode reduction is also presented in Albert et al., 8th FreeFem++
    !! Workshop, 2016 (Eqs. (17)-(20) of the slides).
    !!
    !! Periodic-cylinder chart.  A straight cylinder with axial coordinate z
    !! periodic with period 2 pi R0 is the same complex with the constant
    !! metric radius R0 in place of R: phi = z / R0 is the ignorable angle,
    !! the mode exp(i n phi) has axial wavenumber k = n / R0, and the meridian
    !! plane (x1, x2) = (R, Z) slots hold the Cartesian cross-section (x, y),
    !! which may contain the axis x = y = 0 and negative x1.  The metric is
    !! ds^2 = dx^2 + dy^2 + R0^2 dphi^2, so (grad f)_z = i n f / R0, the
    !! measure is R0 dx dy, and in the variables above
    !!
    !!   a = R0 A_z,  q = R0 B_pol,  s = B_z,  g = R0 div B
    !!
    !! the operators grad_n, curl_n and div_n are unchanged (they are metric
    !! free); only the masses, loads, error norms and the Piola 1/R factors
    !! use R0.  Select it with the optional cylinder_r0 argument of
    !! meridian_feec_build (cylinder_r0 > 0); the default is the axisymmetric
    !! chart of the (R, Z) meridian plane.
    !!
    !! Element-parallel (OpenMP) assembly into per-element triplet slices.
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_feec_operators, only: build_triangle_discrete_gradient
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        evaluate_triangle_lagrange_basis, initialize_triangle_lagrange_basis, &
        triangle_lagrange_basis_t, triangle_lagrange_nodes
    use fortfem_triangle_nedelec_arbitrary_order, only: &
        evaluate_triangle_nedelec_first_kind, &
        initialize_triangle_nedelec_first_kind, &
        triangle_nedelec_first_kind_t
    use fortfem_triangle_rt_arbitrary_order, only: &
        evaluate_triangle_raviart_thomas, initialize_triangle_raviart_thomas, &
        triangle_rt_basis_t
    use fortsparse, only: csc_conjugate_transpose, csc_from_triplet, &
        csc_matmul, csc_t, csc_z_t, fortsparse_status_t, status_ok
    implicit none
    private

    !> Global numbering of one space of the complex.
    type, public :: meridian_space_t
        integer :: ndof = 0
        integer :: nloc = 0
        integer, allocatable :: dof(:, :)      !! (nloc, ntriangles)
        integer, allocatable :: sgn(:, :)      !! orientation signs
        logical, allocatable :: bnd(:)         !! dofs of the boundary trace
    end type meridian_space_t

    !> Operator a + i n b with real sparse a and b.
    type, public :: meridian_op_t
        type(csc_t) :: a
        type(csc_t) :: b
    end type meridian_op_t

    type, public :: meridian_feec_t
        type(meridian_mesh_t) :: mm
        integer :: k = 0                       !! polynomial order
        integer :: g = 1                       !! geometry order
        integer :: nq = 0
        real(dp) :: cyl_r0 = 0.0_dp          !! > 0: periodic-cylinder chart, R = R0
        integer :: n0 = 0, nn = 0, nd = 0      !! Lagrange, Nedelec/RT, DG dofs
        integer :: nl0 = 0, nln = 0, nld = 0   !! local counts
        type(meridian_space_t) :: sp(0:3)
        real(dp), allocatable :: qw(:)
        real(dp), allocatable :: tl(:, :), tlg(:, :, :)
        real(dp), allocatable :: tn(:, :, :), tr(:, :, :), td(:, :)
        real(dp), allocatable :: tgeo(:, :), tgeog(:, :, :)
        type(meridian_op_t) :: grad          !! V0 -> V1
        type(meridian_op_t) :: curl          !! V1 -> V2
        type(meridian_op_t) :: div           !! V2 -> V3
        type(csc_t) :: mass(0:3)             !! physical L2(R dR dZ) masses
    end type meridian_feec_t

    abstract interface
        subroutine meridian_field_i(x, v)
            !! Physical components v = (v_R, v_phi, v_Z) at x = (R, Z).
            import :: dp
            real(dp), intent(in) :: x(2)
            complex(dp), intent(out) :: v(3)
        end subroutine meridian_field_i
        function meridian_weight_i(x) result(w)
            import :: dp
            real(dp), intent(in) :: x(2)
            real(dp) :: w
        end function meridian_weight_i
    end interface

    public :: meridian_field_i, meridian_weight_i
    public :: meridian_feec_build
    public :: meridian_mass_weighted
    public :: meridian_gs_stiffness
    public :: meridian_op_complex
    public :: meridian_stiffness
    public :: meridian_load
    public :: meridian_error_norm
    public :: meridian_zcsc_matmul
    public :: meridian_real_to_complex
    public :: meridian_lagrange_nodes

    real(dp), parameter :: tiny_entry = 1.0e-12_dp
    !> Triangle count below which assembly loops stay serial.
    integer, parameter :: par_min = 128

contains

    ! ------------------------------------------------------------------
    ! Construction
    ! ------------------------------------------------------------------

    subroutine meridian_feec_build(mm, k, fs, status, quad_degree, &
            cylinder_r0)
        !! Build the complex of order k on the meridional mesh mm.  With
        !! cylinder_r0 > 0 the mesh is a cross-section (x, y) of a straight
        !! cylinder of axial period 2 pi cylinder_r0 (see the module header);
        !! otherwise it is an (R, Z) half-plane mesh with R > 0.
        type(meridian_mesh_t), intent(in) :: mm
        integer, intent(in) :: k
        type(meridian_feec_t), intent(out) :: fs
        integer, intent(out) :: status
        integer, intent(in), optional :: quad_degree
        real(dp), intent(in), optional :: cylinder_r0

        integer :: nt, qdeg, i

        status = 1
        if (k < 1 .or. k > 6) return
        if (.not. allocated(mm%nodes)) return
        if (present(cylinder_r0)) then
            if (.not. (cylinder_r0 > 0.0_dp)) return
            fs%cyl_r0 = cylinder_r0
        end if
        fs%mm = mm
        fs%k = k
        fs%g = mm%geo_order
        nt = fs%mm%topology%n_triangles
        fs%nl0 = (k + 1)*(k + 2)/2
        fs%nln = k*(k + 2)
        fs%nld = k*(k + 1)/2
        fs%nd = nt*fs%nld

        call build_numbering(fs, status)
        if (status /= 0) return
        qdeg = 2*k + 2*fs%g + 2
        if (present(quad_degree)) qdeg = quad_degree
        call build_tables(fs, qdeg, status)
        if (status /= 0) return
        call build_operators(fs, status)
        if (status /= 0) return
        do i = 0, 3
            call meridian_mass_weighted(fs, i, fs%mass(i), status)
            if (status /= 0) return
        end do
    end subroutine meridian_feec_build

    subroutine build_numbering(fs, status)
        type(meridian_feec_t), intent(inout) :: fs
        integer, intent(out) :: status

        integer, allocatable :: ned_dof(:, :), ned_sgn(:, :), lag_dof(:, :)
        integer, allocatable :: ecount(:)
        logical, allocatable :: lag_bnd(:)
        integer :: nt, k, nn, ne, t, i, l, e0, nl, m

        status = 0
        k = fs%k
        nt = fs%mm%topology%n_triangles
        call build_edge_numbering(fs%mm%topology%triangles, &
            fs%mm%topology%n_vertices, k, ned_dof, ned_sgn, ne, ecount)
        nn = ne*k + nt*k*(k - 1)
        fs%nn = nn
        call build_lagrange_map(fs, ned_dof, ned_sgn, ecount, ne, lag_dof, &
            lag_bnd)

        ! V0
        nl = fs%nl0
        fs%sp(0)%ndof = fs%n0
        fs%sp(0)%nloc = nl
        fs%sp(0)%dof = lag_dof
        allocate (fs%sp(0)%sgn(nl, nt))
        fs%sp(0)%sgn = 1
        fs%sp(0)%bnd = lag_bnd
        ! V1
        nl = fs%nln + fs%nl0
        fs%sp(1)%ndof = nn + fs%n0
        fs%sp(1)%nloc = nl
        allocate (fs%sp(1)%dof(nl, nt), fs%sp(1)%sgn(nl, nt))
        fs%sp(1)%dof(1:fs%nln, :) = ned_dof
        fs%sp(1)%sgn(1:fs%nln, :) = ned_sgn
        fs%sp(1)%dof(fs%nln + 1:nl, :) = nn + lag_dof
        fs%sp(1)%sgn(fs%nln + 1:nl, :) = 1
        allocate (fs%sp(1)%bnd(nn + fs%n0))
        fs%sp(1)%bnd = .false.
        call mark_edge_bnd(fs, ecount, fs%sp(1)%bnd(1:nn))
        fs%sp(1)%bnd(nn + 1:) = lag_bnd
        ! V2
        nl = fs%nln + fs%nld
        fs%sp(2)%ndof = nn + fs%nd
        fs%sp(2)%nloc = nl
        allocate (fs%sp(2)%dof(nl, nt), fs%sp(2)%sgn(nl, nt))
        fs%sp(2)%dof(1:fs%nln, :) = ned_dof
        fs%sp(2)%sgn(1:fs%nln, :) = ned_sgn
        do t = 1, nt
            do l = 1, fs%nld
                fs%sp(2)%dof(fs%nln + l, t) = nn + (t - 1)*fs%nld + l
            end do
        end do
        fs%sp(2)%sgn(fs%nln + 1:nl, :) = 1
        allocate (fs%sp(2)%bnd(nn + fs%nd))
        fs%sp(2)%bnd = .false.
        call mark_edge_bnd(fs, ecount, fs%sp(2)%bnd(1:nn))
        ! V3
        nl = fs%nld
        fs%sp(3)%ndof = fs%nd
        fs%sp(3)%nloc = nl
        allocate (fs%sp(3)%dof(nl, nt), fs%sp(3)%sgn(nl, nt))
        do t = 1, nt
            do m = 1, nl
                fs%sp(3)%dof(m, t) = (t - 1)*nl + m
            end do
        end do
        fs%sp(3)%sgn = 1
        allocate (fs%sp(3)%bnd(fs%nd))
        fs%sp(3)%bnd = .false.
        status = 0
    end subroutine build_numbering

    subroutine build_edge_numbering(tri, nv, k, ned_dof, ned_sgn, ne, ecount)
        !! Nedelec/RT dof numbering with edge moments (e0 k + m) and element
        !! interior moments, edges found by a linear-time sort of the
        !! half-edges.  An odd edge moment changes sign when the local edge
        !! direction opposes the global one (lower to higher vertex).
        integer, intent(in) :: tri(:, :), nv, k
        integer, allocatable, intent(out) :: ned_dof(:, :), ned_sgn(:, :)
        integer, intent(out) :: ne
        integer, allocatable, intent(out) :: ecount(:)

        integer, allocatable :: lo(:), hi(:), perm(:), tmp(:), eid(:)
        integer :: nt, nh, h, i, j, t, e0, m, nloc, nint, c, a, b, prev_lo, prev_hi

        nt = size(tri, 2)
        nh = 3*nt
        allocate (lo(nh), hi(nh), perm(nh), tmp(nh), eid(nh))
        do t = 1, nt
            do i = 1, 3
                j = mod(i, 3) + 1
                h = 3*(t - 1) + i
                lo(h) = min(tri(i, t), tri(j, t))
                hi(h) = max(tri(i, t), tri(j, t))
            end do
        end do
        ! stable counting sort by hi, then by lo
        call counting_sort(hi, nv, [(h, h=1, nh)], tmp)
        call counting_sort(lo, nv, tmp, perm)
        ne = 0
        prev_lo = 0
        prev_hi = 0
        do h = 1, nh
            a = lo(perm(h))
            b = hi(perm(h))
            if (a /= prev_lo .or. b /= prev_hi) then
                ne = ne + 1
                prev_lo = a
                prev_hi = b
            end if
            eid(perm(h)) = ne - 1
        end do
        allocate (ecount(0:ne - 1))
        ecount = 0
        do h = 1, nh
            ecount(eid(h)) = ecount(eid(h)) + 1
        end do
        nloc = k*(k + 2)
        nint = k*(k - 1)
        allocate (ned_dof(nloc, nt), ned_sgn(nloc, nt))
        ned_sgn = 1
        do t = 1, nt
            do i = 1, 3
                h = 3*(t - 1) + i
                e0 = eid(h)
                do m = 1, k
                    ned_dof((i - 1)*k + m, t) = e0*k + m
                    if (tri(i, t) /= lo(h) .and. mod(m, 2) == 1) &
                        ned_sgn((i - 1)*k + m, t) = -1
                end do
            end do
            do c = 1, nint
                ned_dof(3*k + c, t) = ne*k + (t - 1)*nint + c
            end do
        end do
    end subroutine build_edge_numbering

    subroutine counting_sort(key, nkey, order_in, order_out)
        !! order_out = order_in stably sorted by key(order_in), keys 1..nkey.
        integer, intent(in) :: key(:), nkey, order_in(:)
        integer, intent(out) :: order_out(:)
        integer, allocatable :: start(:)
        integer :: h, n

        n = size(order_in)
        allocate (start(nkey + 1))
        start = 0
        do h = 1, n
            start(key(order_in(h)) + 1) = start(key(order_in(h)) + 1) + 1
        end do
        start(1) = 1
        do h = 2, nkey + 1
            start(h) = start(h) + start(h - 1)
        end do
        do h = 1, n
            order_out(start(key(order_in(h)))) = order_in(h)
            start(key(order_in(h))) = start(key(order_in(h))) + 1
        end do
    end subroutine counting_sort

    subroutine mark_edge_bnd(fs, ecount, mask)
        !! Nedelec/RT dofs on boundary edges (one adjacent triangle).
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: ecount(0:)
        logical, intent(inout) :: mask(:)
        integer :: e0, j

        do e0 = 0, size(ecount) - 1
            if (ecount(e0) == 1) then
                do j = 1, fs%k
                    mask(e0*fs%k + j) = .true.
                end do
            end if
        end do
    end subroutine mark_edge_bnd

    subroutine lagrange_local_nodes(k, pw)
        !! Barycentric lattice weights (a1, a2, a3), a1 + a2 + a3 = k, of the
        !! local nodes (weight a1 belongs to vertex 1 = (0, 0), a2 to
        !! (1, 0) and a3 to (0, 1)).
        integer, intent(in) :: k
        integer, allocatable, intent(out) :: pw(:, :)
        type(triangle_lagrange_basis_t) :: basis
        real(dp), allocatable :: nodes(:, :)
        integer :: i, status

        call initialize_triangle_lagrange_basis(k, basis, status)
        call triangle_lagrange_nodes(basis, nodes)
        allocate (pw(3, size(nodes, 2)))
        do i = 1, size(nodes, 2)
            if (k == 0) then
                pw(:, i) = 0
            else
                pw(2, i) = nint(nodes(1, i)*k)
                pw(3, i) = nint(nodes(2, i)*k)
                pw(1, i) = k - pw(2, i) - pw(3, i)
            end if
        end do
    end subroutine lagrange_local_nodes

    subroutine build_lagrange_map(fs, ned_dof, ned_sgn, ecount, ne, lag_dof, &
            lag_bnd)
        type(meridian_feec_t), intent(inout) :: fs
        integer, intent(in) :: ned_dof(:, :), ned_sgn(:, :), ecount(0:)
        integer, intent(in) :: ne
        integer, allocatable, intent(out) :: lag_dof(:, :)
        logical, allocatable, intent(out) :: lag_bnd(:)

        integer, allocatable :: pw(:, :), tri(:, :)
        integer :: k, nt, nv, nint, t, m, i, j, third, e0, ori, mpos, cnt

        k = fs%k
        nt = fs%mm%topology%n_triangles
        nv = fs%mm%topology%n_vertices
        nint = (k - 1)*(k - 2)/2
        fs%n0 = nv + ne*(k - 1) + nt*nint
        call lagrange_local_nodes(k, pw)
        tri = fs%mm%topology%triangles
        allocate (lag_dof(fs%nl0, nt), lag_bnd(fs%n0))
        lag_bnd = .false.
        do t = 1, nt
            cnt = 0
            do m = 1, fs%nl0
                if (any(pw(:, m) == k)) then
                    i = maxloc(pw(:, m), 1)
                    lag_dof(m, t) = tri(i, t)
                    j = mod(i + 1, 3) + 1      ! other edge at vertex i
                    if (ecount((ned_dof((i - 1)*k + 1, t) - 1)/k) == 1 .or. &
                        ecount((ned_dof((j - 1)*k + 1, t) - 1)/k) == 1) then
                        lag_bnd(lag_dof(m, t)) = .true.
                    end if
                else if (any(pw(:, m) == 0)) then
                    third = minloc(pw(:, m), 1)
                    i = mod(third, 3) + 1        ! edge i: vertex i -> j
                    j = mod(i, 3) + 1
                    e0 = (ned_dof((i - 1)*k + 1, t) - 1)/k
                    ori = ned_sgn((i - 1)*k + 1, t)
                    ! position measured from the lower global vertex
                    if (ori == 1) then
                        mpos = pw(j, m)
                    else
                        mpos = pw(i, m)
                    end if
                    lag_dof(m, t) = nv + e0*(k - 1) + mpos
                    if (ecount(e0) == 1) lag_bnd(lag_dof(m, t)) = .true.
                else
                    cnt = cnt + 1
                    lag_dof(m, t) = nv + ne*(k - 1) + (t - 1)*nint + cnt
                end if
            end do
        end do
    end subroutine build_lagrange_map

    subroutine build_tables(fs, qdeg, status)
        type(meridian_feec_t), intent(inout) :: fs
        integer, intent(in) :: qdeg
        integer, intent(out) :: status

        type(triangle_lagrange_basis_t) :: lb, db, gb
        type(triangle_nedelec_first_kind_t) :: nb
        type(triangle_rt_basis_t) :: rb
        real(dp), allocatable :: xi(:), eta(:), cur(:), dgr(:, :)
        integer :: q, ngeo

        call triangle_duffy_quadrature(qdeg, xi, eta, fs%qw, status)
        if (status /= 0) return
        fs%nq = size(xi)
        ngeo = size(fs%mm%nodes, 2)
        call initialize_triangle_lagrange_basis(fs%k, lb, status)
        if (status /= 0) return
        call initialize_triangle_lagrange_basis(fs%k - 1, db, status)
        if (status /= 0) return
        call initialize_triangle_lagrange_basis(fs%g, gb, status)
        if (status /= 0) return
        call initialize_triangle_nedelec_first_kind(fs%k, nb, status)
        if (status /= 0) return
        call initialize_triangle_raviart_thomas(fs%k - 1, rb, status)
        if (status /= 0) return
        allocate (fs%tl(fs%nl0, fs%nq), fs%tlg(2, fs%nl0, fs%nq))
        allocate (fs%tn(2, fs%nln, fs%nq), fs%tr(2, fs%nln, fs%nq))
        allocate (fs%td(fs%nld, fs%nq))
        allocate (fs%tgeo(ngeo, fs%nq), fs%tgeog(2, ngeo, fs%nq))
        allocate (cur(fs%nln), dgr(2, fs%nld))
        do q = 1, fs%nq
            call evaluate_triangle_lagrange_basis(lb, xi(q), eta(q), &
                fs%tl(:, q), fs%tlg(:, :, q), status)
            if (status /= 0) return
            call evaluate_triangle_lagrange_basis(db, xi(q), eta(q), &
                fs%td(:, q), dgr, status)
            if (status /= 0) return
            call evaluate_triangle_lagrange_basis(gb, xi(q), eta(q), &
                fs%tgeo(:, q), fs%tgeog(:, :, q), status)
            if (status /= 0) return
            call evaluate_triangle_nedelec_first_kind(nb, xi(q), eta(q), &
                fs%tn(:, :, q), cur, status)
            if (status /= 0) return
            call evaluate_triangle_raviart_thomas(rb, xi(q), eta(q), &
                fs%tr(:, :, q), cur, status)
            if (status /= 0) return
        end do
        status = 0
    end subroutine build_tables

    ! ------------------------------------------------------------------
    ! Geometry and physical basis functions
    ! ------------------------------------------------------------------

    pure function metric_radius(fs, x) result(r)
        !! The metric radius R: the coordinate R, or R0 on the cylinder chart.
        type(meridian_feec_t), intent(in) :: fs
        real(dp), intent(in) :: x(2)
        real(dp) :: r

        if (fs%cyl_r0 > 0.0_dp) then
            r = fs%cyl_r0
        else
            r = x(1)
        end if
    end function metric_radius

    pure subroutine point_geometry(fs, t, q, x, jac, det)
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: t, q
        real(dp), intent(out) :: x(2), jac(2, 2), det
        integer :: i

        x = 0.0_dp
        jac = 0.0_dp
        do i = 1, size(fs%tgeo, 1)
            x = x + fs%tgeo(i, q)*fs%mm%nodes(:, i, t)
            jac(:, 1) = jac(:, 1) + fs%tgeog(1, i, q)*fs%mm%nodes(:, i, t)
            jac(:, 2) = jac(:, 2) + fs%tgeog(2, i, q)*fs%mm%nodes(:, i, t)
        end do
        det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
    end subroutine point_geometry

    pure subroutine point_basis(fs, space, q, x, jac, det, phi)
        !! Physical components (R, phi, Z) of the local basis functions.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: space, q
        real(dp), intent(in) :: x(2), jac(2, 2), det
        real(dp), intent(out) :: phi(:, :)
        real(dp) :: jit(2, 2), v(2), rr
        integer :: i, nln

        rr = metric_radius(fs, x)
        nln = fs%nln
        phi = 0.0_dp
        select case (space)
        case (0)
            phi(1, :) = fs%tl(:, q)
        case (1)
            jit(1, 1) = jac(2, 2)/det
            jit(2, 1) = -jac(1, 2)/det
            jit(1, 2) = -jac(2, 1)/det
            jit(2, 2) = jac(1, 1)/det
            do i = 1, nln
                v = matmul(jit, fs%tn(:, i, q))
                phi(1, i) = v(1)
                phi(3, i) = v(2)
            end do
            do i = 1, fs%nl0
                phi(2, nln + i) = fs%tl(i, q)/rr
            end do
        case (2)
            do i = 1, nln
                v = matmul(jac, fs%tr(:, i, q))/(det*rr)
                phi(1, i) = v(1)
                phi(3, i) = v(2)
            end do
            do i = 1, fs%nld
                phi(2, nln + i) = fs%td(i, q)/det
            end do
        case (3)
            phi(1, :) = fs%td(:, q)/(det*rr)
        end select
    end subroutine point_basis

    ! ------------------------------------------------------------------
    ! Differential operators
    ! ------------------------------------------------------------------

    subroutine build_operators(fs, status)
        type(meridian_feec_t), intent(inout) :: fs
        integer, intent(out) :: status

        type(triangle_nedelec_first_kind_t) :: nb
        type(triangle_rt_basis_t) :: rb
        type(triangle_lagrange_basis_t) :: db
        type(csc_t) :: gpol, pmat, pg
        type(fortsparse_status_t) :: fst
        real(dp), allocatable :: gref(:, :), cref(:, :), dref(:, :)
        real(dp), allocatable :: nodes(:, :), cur(:), vals(:, :), dv(:)
        integer, allocatable :: owner(:), pi_(:), rows(:), cols(:)
        real(dp), allocatable :: sig(:), tv(:)
        integer :: k, nn, nt, n0, nd, nln, nl0, nld, t, j, m, l, idx, half
        integer :: ncap, pos, i, row, e3
        real(dp) :: v

        k = fs%k
        nn = fs%nn
        nt = fs%mm%topology%n_triangles
        n0 = fs%n0
        nd = fs%nd
        nln = fs%nln
        nl0 = fs%nl0
        nld = fs%nld

        ! reference matrices
        call build_triangle_discrete_gradient(k, gref, status)
        if (status /= 0) return
        call initialize_triangle_nedelec_first_kind(k, nb, status)
        call initialize_triangle_raviart_thomas(k - 1, rb, status)
        call initialize_triangle_lagrange_basis(k - 1, db, status)
        call triangle_lagrange_nodes(db, nodes)
        allocate (cref(nld, nln), dref(nld, nln), vals(2, nln), cur(nln))
        allocate (dv(nln))
        do l = 1, nld
            call evaluate_triangle_nedelec_first_kind(nb, nodes(1, l), &
                nodes(2, l), vals, cur, status)
            if (status /= 0) return
            cref(l, :) = cur
            call evaluate_triangle_raviart_thomas(rb, nodes(1, l), &
                nodes(2, l), vals, dv, status)
            if (status /= 0) return
            dref(l, :) = dv
        end do

        ! signed permutation RT local <- Nedelec local (rot90 N = -sig RT)
        allocate (pi_(nln), sig(nln))
        half = k*(k - 1)/2
        e3 = 3*k
        do l = 1, nln
            if (l <= e3) then
                pi_(l) = l
                sig(l) = 1.0_dp
            else if (l <= e3 + half) then
                pi_(l) = l + half
                sig(l) = 1.0_dp
            else
                pi_(l) = l - half
                sig(l) = -1.0_dp
            end if
        end do

        ! owner element of every Nedelec dof
        allocate (owner(nn))
        owner = 0
        do t = 1, nt
            do j = 1, nln
                row = fs%sp(1)%dof(j, t)
                if (owner(row) == 0) owner(row) = t
            end do
        end do

        ! G_pol: Lagrange -> Nedelec (rows owned by one element)
        ncap = nt*nln*nl0
        allocate (rows(ncap), cols(ncap), tv(ncap))
        !$omp parallel do if(nt >= par_min) default(shared) &
        !$omp private(t, j, m, idx, row, v) &
        !$omp schedule(static)
        do t = 1, nt
            do j = 1, nln
                row = fs%sp(1)%dof(j, t)
                do m = 1, nl0
                    idx = ((t - 1)*nln + j - 1)*nl0 + m
                    rows(idx) = row
                    cols(idx) = fs%sp(0)%dof(m, t)
                    v = real(fs%sp(1)%sgn(j, t), dp)*gref(j, m)
                    if (owner(row) /= t .or. abs(v) < tiny_entry) v = 0.0_dp
                    tv(idx) = v
                end do
            end do
        end do
        !$omp end parallel do
        call triplets_to_csc(nn, n0, rows, cols, tv, gpol, status)
        if (status /= 0) return
        deallocate (rows, cols, tv)

        ! P: beta coefficients -> q coefficients, q = rot90 beta
        ncap = nt*nln
        allocate (rows(ncap), cols(ncap), tv(ncap))
        !$omp parallel do if(nt >= par_min) default(shared) &
        !$omp private(t, l, idx, row, v) &
        !$omp schedule(static)
        do t = 1, nt
            do l = 1, nln
                idx = (t - 1)*nln + l
                row = fs%sp(1)%dof(l, t)
                rows(idx) = row
                cols(idx) = fs%sp(1)%dof(pi_(l), t)
                v = -sig(l)*real(fs%sp(1)%sgn(pi_(l), t)* &
                    fs%sp(1)%sgn(l, t), dp)
                if (owner(row) /= t) v = 0.0_dp
                tv(idx) = v
            end do
        end do
        !$omp end parallel do
        call triplets_to_csc(nn, nn, rows, cols, tv, pmat, status)
        if (status /= 0) return
        deallocate (rows, cols, tv)
        call csc_matmul(pmat, gpol, pg, fst)
        if (.not. status_ok(fst)) then
            status = 2
            return
        end if

        ! grad = [G_pol; 0] + i n [0; I]
        call csc_to_triplets_offset(gpol, 0, 0, 1.0_dp, rows, cols, tv, 0)
        call triplets_to_csc(nn + n0, n0, rows, cols, tv, fs%grad%a, status)
        if (status /= 0) return
        deallocate (rows, cols, tv)
        allocate (rows(n0), cols(n0), tv(n0))
        do i = 1, n0
            rows(i) = nn + i
            cols(i) = i
            tv(i) = 1.0_dp
        end do
        call triplets_to_csc(nn + n0, n0, rows, cols, tv, fs%grad%b, status)
        if (status /= 0) return
        deallocate (rows, cols, tv)

        ! curl.a = [[0, P G_pol], [-C_pol, 0]],  curl.b = [[-P, 0], [0, 0]]
        ncap = pg%nnz + nt*nld*nln
        allocate (rows(ncap), cols(ncap), tv(ncap))
        call csc_to_triplets_offset(pg, 0, nn, 1.0_dp, rows, cols, tv, 0)
        pos = pg%nnz
        do t = 1, nt
            do l = 1, nld
                do j = 1, nln
                    pos = pos + 1
                    rows(pos) = nn + (t - 1)*nld + l
                    cols(pos) = fs%sp(1)%dof(j, t)
                    v = -real(fs%sp(1)%sgn(j, t), dp)*cref(l, j)
                    if (abs(v) < tiny_entry) v = 0.0_dp
                    tv(pos) = v
                end do
            end do
        end do
        call triplets_to_csc(nn + nd, nn + n0, rows, cols, tv(1:pos), &
            fs%curl%a, status, pos)
        if (status /= 0) return
        deallocate (rows, cols, tv)
        allocate (rows(pmat%nnz), cols(pmat%nnz), tv(pmat%nnz))
        call csc_to_triplets_offset(pmat, 0, 0, -1.0_dp, rows, cols, tv, 0)
        call triplets_to_csc(nn + nd, nn + n0, rows, cols, tv, fs%curl%b, &
            status)
        if (status /= 0) return
        deallocate (rows, cols, tv)

        ! div.a = [D_pol, 0],  div.b = [0, I]
        ncap = nt*nld*nln
        allocate (rows(ncap), cols(ncap), tv(ncap))
        pos = 0
        do t = 1, nt
            do l = 1, nld
                do j = 1, nln
                    pos = pos + 1
                    rows(pos) = (t - 1)*nld + l
                    cols(pos) = fs%sp(2)%dof(j, t)
                    v = real(fs%sp(2)%sgn(j, t), dp)*dref(l, j)
                    if (abs(v) < tiny_entry) v = 0.0_dp
                    tv(pos) = v
                end do
            end do
        end do
        call triplets_to_csc(nd, nn + nd, rows, cols, tv, fs%div%a, status)
        if (status /= 0) return
        deallocate (rows, cols, tv)
        allocate (rows(nd), cols(nd), tv(nd))
        do i = 1, nd
            rows(i) = i
            cols(i) = nn + i
            tv(i) = 1.0_dp
        end do
        call triplets_to_csc(nd, nn + nd, rows, cols, tv, fs%div%b, status)
    end subroutine build_operators

    subroutine triplets_to_csc(nrow, ncol, rows, cols, vals, a, status, count)
        !! CSC matrix from triplets, dropping exact zeros.
        integer, intent(in) :: nrow, ncol
        integer, intent(in) :: rows(:), cols(:)
        real(dp), intent(in) :: vals(:)
        type(csc_t), intent(out) :: a
        integer, intent(out) :: status
        integer, intent(in), optional :: count

        type(fortsparse_status_t) :: fst
        integer, allocatable :: r(:), c(:)
        real(dp), allocatable :: v(:)
        integer :: n, i, cnt

        n = size(vals)
        if (present(count)) n = count
        cnt = count_nonzero(vals(1:n))
        allocate (r(cnt), c(cnt), v(cnt))
        cnt = 0
        do i = 1, n
            if (vals(i) /= 0.0_dp) then
                cnt = cnt + 1
                r(cnt) = rows(i)
                c(cnt) = cols(i)
                v(cnt) = vals(i)
            end if
        end do
        call csc_from_triplet(nrow, ncol, r, c, v, a, fst)
        status = 0
        if (.not. status_ok(fst)) status = 3
    end subroutine triplets_to_csc

    pure integer function count_nonzero(v) result(c)
        real(dp), intent(in) :: v(:)
        c = count(v /= 0.0_dp)
    end function count_nonzero

    subroutine csc_to_triplets_offset(a, ro, co, scale, rows, cols, vals, pos)
        !! Append the entries of a, shifted by (ro, co), at index pos + 1.
        !! The arrays are (re)allocated when pos = 0.
        type(csc_t), intent(in) :: a
        integer, intent(in) :: ro, co, pos
        real(dp), intent(in) :: scale
        integer, allocatable, intent(inout) :: rows(:), cols(:)
        real(dp), allocatable, intent(inout) :: vals(:)
        integer :: j, p, n

        if (pos == 0 .and. .not. (allocated(rows) .and. &
            size(rows) >= a%nnz)) then
            if (allocated(rows)) deallocate (rows, cols, vals)
            allocate (rows(a%nnz), cols(a%nnz), vals(a%nnz))
        end if
        n = pos
        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                n = n + 1
                rows(n) = a%row_idx(p) + ro
                cols(n) = j + co
                vals(n) = scale*a%val(p)
            end do
        end do
    end subroutine csc_to_triplets_offset

    ! ------------------------------------------------------------------
    ! Weak forms
    ! ------------------------------------------------------------------

    subroutine meridian_mass_weighted(fs, space, m, status, weight)
        !! M_ij = int phi_i . phi_j w(R, Z) R dR dZ of the physical basis
        !! fields of the given space (weight w = 1 by default).
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: space
        type(csc_t), intent(out) :: m
        integer, intent(out) :: status
        procedure(meridian_weight_i), optional :: weight

        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:), phi(:, :), ml(:, :)
        real(dp) :: x(2), jac(2, 2), det, wq, extra
        integer :: nt, nl, t, q, i, j, base, bad

        status = 1
        if (space < 0 .or. space > 3) return
        nt = fs%mm%topology%n_triangles
        nl = fs%sp(space)%nloc
        allocate (rows(nl*nl*nt), cols(nl*nl*nt), vals(nl*nl*nt))
        bad = 0
        !$omp parallel if(nt >= par_min) default(shared) &
        !$omp private(t, q, i, j, base, phi, ml, x, jac, det, wq, extra) &
        !$omp reduction(max:bad)
        allocate (phi(3, nl), ml(nl, nl))
        !$omp do schedule(static)
        do t = 1, nt
            ml = 0.0_dp
            do q = 1, fs%nq
                call point_geometry(fs, t, q, x, jac, det)
                if (det <= 0.0_dp .or. metric_radius(fs, x) <= 0.0_dp) then
                    bad = 1
                    cycle
                end if
                call point_basis(fs, space, q, x, jac, det, phi)
                extra = 1.0_dp
                if (present(weight)) extra = weight(x)
                wq = fs%qw(q)*det*metric_radius(fs, x)*extra
                do j = 1, nl
                    do i = 1, nl
                        ml(i, j) = ml(i, j) + wq*(phi(1, i)*phi(1, j) &
                            + phi(2, i)*phi(2, j) + phi(3, i)*phi(3, j))
                    end do
                end do
            end do
            base = (t - 1)*nl*nl
            do j = 1, nl
                do i = 1, nl
                    base = base + 1
                    rows(base) = fs%sp(space)%dof(i, t)
                    cols(base) = fs%sp(space)%dof(j, t)
                    vals(base) = real(fs%sp(space)%sgn(i, t)* &
                        fs%sp(space)%sgn(j, t), dp)*ml(i, j)
                    if (abs(vals(base)) < 1.0e-300_dp) vals(base) = 0.0_dp
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        if (bad /= 0) return
        call triplets_to_csc(fs%sp(space)%ndof, fs%sp(space)%ndof, rows, &
            cols, vals, m, status)
    end subroutine meridian_mass_weighted

    subroutine meridian_gs_stiffness(fs, k, status)
        !! Direct axisymmetric Grad-Shafranov stiffness
        !! K_ij = int (1/R) grad L_i . grad L_j dR dZ on V0.
        type(meridian_feec_t), intent(in) :: fs
        type(csc_t), intent(out) :: k
        integer, intent(out) :: status

        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:), ml(:, :), gp(:, :)
        real(dp) :: x(2), jac(2, 2), det, wq, jit(2, 2)
        integer :: nt, nl, t, q, i, j, base, bad

        nt = fs%mm%topology%n_triangles
        nl = fs%nl0
        allocate (rows(nl*nl*nt), cols(nl*nl*nt), vals(nl*nl*nt))
        bad = 0
        !$omp parallel if(nt >= par_min) default(shared) &
        !$omp private(t, q, i, j, base, gp, ml, x, jac, det, wq, jit) &
        !$omp reduction(max:bad)
        allocate (gp(2, nl), ml(nl, nl))
        !$omp do schedule(static)
        do t = 1, nt
            ml = 0.0_dp
            do q = 1, fs%nq
                call point_geometry(fs, t, q, x, jac, det)
                if (det <= 0.0_dp .or. metric_radius(fs, x) <= 0.0_dp) then
                    bad = 1
                    cycle
                end if
                jit(1, 1) = jac(2, 2)/det
                jit(2, 1) = -jac(1, 2)/det
                jit(1, 2) = -jac(2, 1)/det
                jit(2, 2) = jac(1, 1)/det
                gp = matmul(jit, fs%tlg(:, :, q))
                wq = fs%qw(q)*det/metric_radius(fs, x)
                do j = 1, nl
                    do i = 1, nl
                        ml(i, j) = ml(i, j) + wq*(gp(1, i)*gp(1, j) &
                            + gp(2, i)*gp(2, j))
                    end do
                end do
            end do
            base = (t - 1)*nl*nl
            do j = 1, nl
                do i = 1, nl
                    base = base + 1
                    rows(base) = fs%sp(0)%dof(i, t)
                    cols(base) = fs%sp(0)%dof(j, t)
                    vals(base) = ml(i, j)
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        status = 1
        if (bad /= 0) return
        call triplets_to_csc(fs%n0, fs%n0, rows, cols, vals, k, status)
    end subroutine meridian_gs_stiffness

    subroutine meridian_lagrange_nodes(fs, xy)
        !! Physical coordinates xy(:, i) of the V0 (Lagrange) nodes.
        type(meridian_feec_t), intent(in) :: fs
        real(dp), allocatable, intent(out) :: xy(:, :)

        type(triangle_lagrange_basis_t) :: lb, gb
        real(dp), allocatable :: nodes(:, :), gv(:), gg(:, :)
        integer :: t, m, i, status

        call initialize_triangle_lagrange_basis(fs%k, lb, status)
        call initialize_triangle_lagrange_basis(fs%g, gb, status)
        call triangle_lagrange_nodes(lb, nodes)
        allocate (xy(2, fs%n0), gv(size(fs%mm%nodes, 2)))
        allocate (gg(2, size(gv)))
        do t = 1, fs%mm%topology%n_triangles
            do m = 1, fs%nl0
                call evaluate_triangle_lagrange_basis(gb, nodes(1, m), &
                    nodes(2, m), gv, gg, status)
                xy(:, fs%sp(0)%dof(m, t)) = 0.0_dp
                do i = 1, size(gv)
                    xy(:, fs%sp(0)%dof(m, t)) = xy(:, fs%sp(0)%dof(m, t)) &
                        + gv(i)*fs%mm%nodes(:, i, t)
                end do
            end do
        end do
    end subroutine meridian_lagrange_nodes

    subroutine meridian_real_to_complex(a, z)
        type(csc_t), intent(in) :: a
        type(csc_z_t), intent(out) :: z

        z%nrow = a%nrow
        z%ncol = a%ncol
        z%nnz = a%nnz
        z%col_ptr = a%col_ptr
        z%row_idx = a%row_idx
        z%val = cmplx(a%val, 0.0_dp, dp)
    end subroutine meridian_real_to_complex

    subroutine meridian_op_complex(op, n, z, status)
        !! z = a + i n b.
        type(meridian_op_t), intent(in) :: op
        integer, intent(in) :: n
        type(csc_z_t), intent(out) :: z
        integer, intent(out) :: status

        type(fortsparse_status_t) :: fst
        integer, allocatable :: rows(:), cols(:)
        complex(dp), allocatable :: vals(:)
        integer :: j, p, c

        c = op%a%nnz + op%b%nnz
        allocate (rows(c), cols(c), vals(c))
        c = 0
        do j = 1, op%a%ncol
            do p = op%a%col_ptr(j), op%a%col_ptr(j + 1) - 1
                c = c + 1
                rows(c) = op%a%row_idx(p)
                cols(c) = j
                vals(c) = cmplx(op%a%val(p), 0.0_dp, dp)
            end do
            do p = op%b%col_ptr(j), op%b%col_ptr(j + 1) - 1
                c = c + 1
                rows(c) = op%b%row_idx(p)
                cols(c) = j
                vals(c) = cmplx(0.0_dp, real(n, dp)*op%b%val(p), dp)
            end do
        end do
        call csc_from_triplet(op%a%nrow, op%a%ncol, rows, cols, vals, z, fst)
        status = 0
        if (.not. status_ok(fst)) status = 1
    end subroutine meridian_op_complex

    subroutine meridian_zcsc_matmul(a, b, c, status)
        type(csc_z_t), intent(in) :: a, b
        type(csc_z_t), intent(out) :: c
        integer, intent(out) :: status

        type(fortsparse_status_t) :: fst
        integer, allocatable :: rows(:), cols(:)
        complex(dp), allocatable :: vals(:)
        integer :: j, p, r, cnt, total

        status = 1
        if (a%ncol /= b%nrow) return
        total = 0
        do j = 1, b%ncol
            do p = b%col_ptr(j), b%col_ptr(j + 1) - 1
                r = b%row_idx(p)
                total = total + a%col_ptr(r + 1) - a%col_ptr(r)
            end do
        end do
        allocate (rows(total), cols(total), vals(total))
        cnt = 0
        do j = 1, b%ncol
            do p = b%col_ptr(j), b%col_ptr(j + 1) - 1
                r = b%row_idx(p)
                rows(cnt + 1:cnt + a%col_ptr(r + 1) - a%col_ptr(r)) = &
                    a%row_idx(a%col_ptr(r):a%col_ptr(r + 1) - 1)
                vals(cnt + 1:cnt + a%col_ptr(r + 1) - a%col_ptr(r)) = &
                    a%val(a%col_ptr(r):a%col_ptr(r + 1) - 1)*b%val(p)
                cols(cnt + 1:cnt + a%col_ptr(r + 1) - a%col_ptr(r)) = j
                cnt = cnt + a%col_ptr(r + 1) - a%col_ptr(r)
            end do
        end do
        call csc_from_triplet(a%nrow, b%ncol, rows, cols, vals, c, fst)
        if (status_ok(fst)) status = 0
    end subroutine meridian_zcsc_matmul

    subroutine meridian_stiffness(fs, space, n, kmat, status)
        !! K = op_n^H M_{space+1} op_n for op = grad_n (space 0),
        !! curl_n (1) or div_n (2): the per-n Dirichlet-form of the
        !! gradient (resistive diffusion of V0), curl-curl (ideal energy and
        !! resistive operator on V1) and grad-div forms.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: space, n
        type(csc_z_t), intent(out) :: kmat
        integer, intent(out) :: status

        type(csc_z_t) :: z, zh, mz, mm
        integer :: st

        status = 1
        select case (space)
        case (0)
            call meridian_op_complex(fs%grad, n, z, st)
        case (1)
            call meridian_op_complex(fs%curl, n, z, st)
        case (2)
            call meridian_op_complex(fs%div, n, z, st)
        case default
            return
        end select
        if (st /= 0) return
        call meridian_real_to_complex(fs%mass(space + 1), mm)
        call meridian_zcsc_matmul(mm, z, mz, st)
        if (st /= 0) return
        call csc_conj_transpose_z(z, zh, st)
        if (st /= 0) return
        call meridian_zcsc_matmul(zh, mz, kmat, status)
    end subroutine meridian_stiffness

    subroutine csc_conj_transpose_z(a, ah, status)
        type(csc_z_t), intent(in) :: a
        type(csc_z_t), intent(out) :: ah
        integer, intent(out) :: status
        type(fortsparse_status_t) :: fst

        call csc_conjugate_transpose(a, ah, fst)
        status = 0
        if (.not. status_ok(fst)) status = 1
    end subroutine csc_conj_transpose_z

    ! ------------------------------------------------------------------
    ! Loads and errors for analytic physical fields
    ! ------------------------------------------------------------------

    subroutine meridian_load(fs, space, fun, b, status)
        !! b_i = int F . phi_i R dR dZ with the physical field F = fun(R, Z).
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: space
        procedure(meridian_field_i) :: fun
        complex(dp), allocatable, intent(out) :: b(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: bl(:, :)
        real(dp), allocatable :: phi(:, :)
        complex(dp) :: f(3)
        real(dp) :: x(2), jac(2, 2), det, wq
        integer :: nt, nl, t, q, i, bad

        status = 1
        nt = fs%mm%topology%n_triangles
        nl = fs%sp(space)%nloc
        allocate (bl(nl, nt), b(fs%sp(space)%ndof))
        bad = 0
        !$omp parallel if(nt >= par_min) default(shared) &
        !$omp private(t, q, i, phi, f, x, jac, det, wq) reduction(max:bad)
        allocate (phi(3, nl))
        !$omp do schedule(static)
        do t = 1, nt
            bl(:, t) = (0.0_dp, 0.0_dp)
            do q = 1, fs%nq
                call point_geometry(fs, t, q, x, jac, det)
                if (det <= 0.0_dp .or. metric_radius(fs, x) <= 0.0_dp) then
                    bad = 1
                    cycle
                end if
                call point_basis(fs, space, q, x, jac, det, phi)
                call fun(x, f)
                wq = fs%qw(q)*det*metric_radius(fs, x)
                do i = 1, nl
                    bl(i, t) = bl(i, t) + wq*(f(1)*phi(1, i) &
                        + f(2)*phi(2, i) + f(3)*phi(3, i))
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        if (bad /= 0) return
        b = (0.0_dp, 0.0_dp)
        do t = 1, nt
            do i = 1, nl
                b(fs%sp(space)%dof(i, t)) = b(fs%sp(space)%dof(i, t)) &
                    + real(fs%sp(space)%sgn(i, t), dp)*bl(i, t)
            end do
        end do
        status = 0
    end subroutine meridian_load

    subroutine meridian_error_norm(fs, space, coeff, fun, err, status)
        !! err = ||u_h - F||_{L2(R dR dZ)} for the field with coefficients
        !! coeff against the analytic physical field F.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: space
        complex(dp), intent(in) :: coeff(:)
        procedure(meridian_field_i) :: fun
        real(dp), intent(out) :: err
        integer, intent(out) :: status

        real(dp), allocatable :: phi(:, :)
        complex(dp), allocatable :: cl(:)
        complex(dp) :: f(3), u(3)
        real(dp) :: x(2), jac(2, 2), det, acc, tot
        integer :: nt, nl, t, q, i, bad

        status = 1
        err = 0.0_dp
        nt = fs%mm%topology%n_triangles
        nl = fs%sp(space)%nloc
        tot = 0.0_dp
        bad = 0
        !$omp parallel if(nt >= par_min) default(shared) &
        !$omp private(t, q, i, phi, cl, f, u, x, jac, det, acc) &
        !$omp reduction(+:tot) reduction(max:bad)
        allocate (phi(3, nl), cl(nl))
        !$omp do schedule(static)
        do t = 1, nt
            do i = 1, nl
                cl(i) = real(fs%sp(space)%sgn(i, t), dp)* &
                    coeff(fs%sp(space)%dof(i, t))
            end do
            do q = 1, fs%nq
                call point_geometry(fs, t, q, x, jac, det)
                if (det <= 0.0_dp .or. metric_radius(fs, x) <= 0.0_dp) then
                    bad = 1
                    cycle
                end if
                call point_basis(fs, space, q, x, jac, det, phi)
                call fun(x, f)
                u = (0.0_dp, 0.0_dp)
                do i = 1, nl
                    u(1) = u(1) + cl(i)*phi(1, i)
                    u(2) = u(2) + cl(i)*phi(2, i)
                    u(3) = u(3) + cl(i)*phi(3, i)
                end do
                acc = abs(u(1) - f(1))**2 + abs(u(2) - f(2))**2 &
                    + abs(u(3) - f(3))**2
                tot = tot + fs%qw(q)*det*metric_radius(fs, x)*acc
            end do
        end do
        !$omp end do
        !$omp end parallel
        if (bad /= 0) return
        err = sqrt(tot)
        status = 0
    end subroutine meridian_error_norm

end module fortfem_meridian_feec_n
