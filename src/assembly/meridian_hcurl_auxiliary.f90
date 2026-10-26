module fortfem_meridian_hcurl_auxiliary
    !! Auxiliary-space data of the meridional Nedelec space N_k (the A_pol
    !! block of V1 in fortfem_meridian_feec_n) for the Hiptmair-Xu
    !! preconditioner of fortfem_hcurl_auxiliary_space:
    !!
    !!   transfer  E : N_1 -> N_k, the exact embedding of the Whitney space.
    !!             On the reference triangle w_i = lambda_a grad lambda_b -
    !!             lambda_b grad lambda_a (edge i from local vertex a = i to
    !!             b = i + 1) lies in N_k, so its coefficients follow from the
    !!             reference Gram system; the covariant Piola map carries them
    !!             to every (curved) element with the edge and dof signs.
    !!   gradient  G : P1 -> N_1, the signed edge-vertex incidence
    !!             (Whitney dofs are tangential edge integrals, lower to higher
    !!             global vertex).
    !!   interp    Pi_R, Pi_Z : P1 -> N_1, the Whitney interpolation of a
    !!             piecewise linear vector field from its vertex values,
    !!             (Pi v)_e = (v_a + v_b)/2 . (x_b - x_a) with the straight
    !!             chord (exact for straight edges).
    !!   patches   vertex stars: the N_k dofs on the edges at a vertex and in
    !!             the interiors of its triangles.
    !!
    !! The maps depend only on the mesh and the order, not on n or on the
    !! coefficients, so one build serves every toroidal mode.
    use fortfem_kinds, only: dp
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_stiffness
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortsparse, only: csc_t, csc_z_t, csc_from_triplet, fortsparse_status_t, &
        status_ok
    implicit none
    private

    type, public :: meridian_hcurl_aux_t
        integer :: ne = 0, nv = 0
        type(csc_t) :: transfer           !! E (nn x ne)
        type(csc_t) :: gradient           !! G (ne x nv)
        type(csc_t) :: interp(2)          !! Pi_R, Pi_Z (ne x nv)
        logical, allocatable :: edge_fixed(:)
        integer, allocatable :: edge_vertex(:, :)   !! (2, ne): lower, higher
        integer, allocatable :: patch_ptr(:), patch_dof(:)
    end type meridian_hcurl_aux_t

    public :: meridian_hcurl_aux_build, meridian_whitney_coefficients
    public :: meridian_hcurl_mode_matrix


contains

    subroutine meridian_whitney_coefficients(k, c, status)
        !! c(l, i): coefficient of the reference N_k basis function l in the
        !! reference Whitney function w_i of local edge i.  The N_k basis of
        !! fortfem_triangle_nedelec_arbitrary_order is dual to its moments
        !! (edge: int_0^1 P_m(s) v.t ds with shifted Legendre P_m and the
        !! reference_edge tangent; interior: int xi^a eta^b v_c over the
        !! reference triangle, c = 1, 2, a + b <= k - 2), so c(:, i) are the
        !! moments of w_i: no Gram solve and no loss of digits with k.  On
        !! its own edge w_i . t = 1 (P_0 moment 1, higher moments 0), and
        !! w_i . t = 0 on the other edges.
        integer, intent(in) :: k
        real(dp), allocatable, intent(out) :: c(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: xi(:), eta(:), w(:)
        real(dp) :: lam(3), glam(2, 3), wv(2), mono
        integer :: nln, q, i, a, b, comp, td, xd, l

        nln = k*(k + 2)
        allocate (c(nln, 3))
        c = 0.0_dp
        do i = 1, 3
            c((i - 1)*k + 1, i) = 1.0_dp
        end do
        status = 0
        if (k < 2) return
        call triangle_duffy_quadrature(2*k, xi, eta, w, status)
        if (status /= 0) return
        glam(:, 1) = [-1.0_dp, -1.0_dp]
        glam(:, 2) = [1.0_dp, 0.0_dp]
        glam(:, 3) = [0.0_dp, 1.0_dp]
        do q = 1, size(w)
            lam = [1.0_dp - xi(q) - eta(q), xi(q), eta(q)]
            do i = 1, 3
                a = i
                b = mod(i, 3) + 1
                wv = lam(a)*glam(:, b) - lam(b)*glam(:, a)
                l = 3*k
                do comp = 1, 2
                    do td = 0, k - 2
                        do xd = 0, td
                            l = l + 1
                            mono = xi(q)**xd*eta(q)**(td - xd)
                            c(l, i) = c(l, i) + w(q)*mono*wv(comp)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine meridian_whitney_coefficients

    subroutine meridian_hcurl_aux_build(fs, aux, status)
        type(meridian_feec_t), intent(in) :: fs
        type(meridian_hcurl_aux_t), intent(out) :: aux
        integer, intent(out) :: status
        real(dp), allocatable :: c(:, :), tv(:), xv(:, :)
        integer, allocatable :: rows(:), cols(:), owner(:), eid(:, :), cnt(:)
        integer, allocatable :: fill(:)
        type(fortsparse_status_t) :: fst
        integer :: k, nt, nv, nn, nln, nint, t, i, j, l, e, a, b, row, pos, d
        integer :: v, ne
        real(dp) :: val, s

        status = 1
        k = fs%k
        nt = fs%mm%topology%n_triangles
        nv = fs%mm%topology%n_vertices
        nn = fs%nn
        nln = fs%nln
        nint = k*(k - 1)
        ! edges
        allocate (eid(3, nt))
        ne = 0
        do t = 1, nt
            do i = 1, 3
                eid(i, t) = (fs%sp(1)%dof((i - 1)*k + 1, t) - 1)/k + 1
                ne = max(ne, eid(i, t))
            end do
        end do
        if (ne*k + nt*nint /= nn) return
        aux%ne = ne
        aux%nv = nv
        allocate (aux%edge_vertex(2, ne), aux%edge_fixed(ne))
        do t = 1, nt
            do i = 1, 3
                a = fs%mm%topology%triangles(i, t)
                b = fs%mm%topology%triangles(mod(i, 3) + 1, t)
                aux%edge_vertex(:, eid(i, t)) = [min(a, b), max(a, b)]
            end do
        end do
        do e = 1, ne
            aux%edge_fixed(e) = fs%sp(1)%bnd((e - 1)*k + 1)
        end do

        ! E = N_1 -> N_k
        call meridian_whitney_coefficients(k, c, status)
        if (status /= 0) return
        status = 1
        allocate (owner(nn))
        owner = 0
        do t = 1, nt
            do l = 1, nln
                row = fs%sp(1)%dof(l, t)
                if (owner(row) == 0) owner(row) = t
            end do
        end do
        allocate (rows(3*nln*nt), cols(3*nln*nt), tv(3*nln*nt))
        pos = 0
        do t = 1, nt
            do i = 1, 3
                s = 1.0_dp
                if (fs%mm%topology%triangles(i, t) /= aux%edge_vertex(1, eid(i, t))) &
                    s = -1.0_dp
                do l = 1, nln
                    row = fs%sp(1)%dof(l, t)
                    if (owner(row) /= t) cycle
                    val = s*c(l, i)*real(fs%sp(1)%sgn(l, t), dp)
                    if (abs(val) < tiny(1.0_dp)) cycle
                    pos = pos + 1
                    rows(pos) = row
                    cols(pos) = eid(i, t)
                    tv(pos) = val
                end do
            end do
        end do
        call csc_from_triplet(nn, ne, rows(:pos), cols(:pos), tv(:pos), &
            aux%transfer, fst)
        if (.not. status_ok(fst)) return
        deallocate (rows, cols, tv)

        ! G and Pi_d
        allocate (rows(2*ne), cols(2*ne), tv(2*ne))
        do e = 1, ne
            rows(2*e - 1:2*e) = e
            cols(2*e - 1) = aux%edge_vertex(1, e)
            cols(2*e) = aux%edge_vertex(2, e)
            tv(2*e - 1) = -1.0_dp
            tv(2*e) = 1.0_dp
        end do
        call csc_from_triplet(ne, nv, rows, cols, tv, aux%gradient, fst)
        if (.not. status_ok(fst)) return
        xv = fs%mm%topology%vertices
        do d = 1, 2
            do e = 1, ne
                a = aux%edge_vertex(1, e)
                b = aux%edge_vertex(2, e)
                tv(2*e - 1:2*e) = 0.5_dp*(xv(d, b) - xv(d, a))
            end do
            call csc_from_triplet(ne, nv, rows, cols, tv, aux%interp(d), fst)
            if (.not. status_ok(fst)) return
        end do

        ! vertex patches
        allocate (cnt(nv + 1))
        cnt = 0
        do e = 1, ne
            do j = 1, 2
                v = aux%edge_vertex(j, e)
                cnt(v + 1) = cnt(v + 1) + k
            end do
        end do
        do t = 1, nt
            do i = 1, 3
                v = fs%mm%topology%triangles(i, t)
                cnt(v + 1) = cnt(v + 1) + nint
            end do
        end do
        cnt(1) = 1
        do v = 2, nv + 1
            cnt(v) = cnt(v) + cnt(v - 1)
        end do
        aux%patch_ptr = cnt
        allocate (aux%patch_dof(cnt(nv + 1) - 1), fill(nv))
        fill = cnt(1:nv)
        do e = 1, ne
            do j = 1, 2
                v = aux%edge_vertex(j, e)
                do l = 1, k
                    aux%patch_dof(fill(v)) = (e - 1)*k + l
                    fill(v) = fill(v) + 1
                end do
            end do
        end do
        do t = 1, nt
            do i = 1, 3
                v = fs%mm%topology%triangles(i, t)
                do l = 1, nint
                    aux%patch_dof(fill(v)) = fs%sp(1)%dof(3*k + l, t)
                    fill(v) = fill(v) + 1
                end do
            end do
        end do
        status = 0
    end subroutine meridian_hcurl_aux_build

    subroutine meridian_hcurl_mode_matrix(fs, n, kmat, status)
        !! Real nn x nn block K_n = rot^T M rot + n^2 M_(1/R) of the per-mode
        !! curl-curl form on A_pol (gauge R A_phi = 0, n /= 0), the system of
        !! the per-mode vacuum solve.  Rejects a non-negligible imaginary part.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: n
        type(csc_t), intent(out) :: kmat
        integer, intent(out) :: status
        type(csc_z_t) :: z
        type(fortsparse_status_t) :: fst
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        real(dp) :: scale
        integer :: j, p, c, nn

        call meridian_stiffness(fs, 1, n, z, status)
        if (status /= 0) return
        status = 1
        nn = fs%nn
        scale = maxval(abs(z%val(:z%nnz)))
        allocate (rows(z%nnz), cols(z%nnz), vals(z%nnz))
        c = 0
        do j = 1, nn
            do p = z%col_ptr(j), z%col_ptr(j + 1) - 1
                if (z%row_idx(p) > nn) cycle
                if (abs(aimag(z%val(p))) > 1.0e-12_dp*scale) return
                c = c + 1
                rows(c) = z%row_idx(p)
                cols(c) = j
                vals(c) = real(z%val(p), dp)
            end do
        end do
        call csc_from_triplet(nn, nn, rows(:c), cols(:c), vals(:c), kmat, fst)
        if (status_ok(fst)) status = 0
    end subroutine meridian_hcurl_mode_matrix

end module fortfem_meridian_hcurl_auxiliary
