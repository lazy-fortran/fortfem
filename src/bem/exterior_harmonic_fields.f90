module fortfem_exterior_harmonic_fields
    !! Decaying exterior cohomology representatives from closed interior wires.
    !! Unit mu0*current; outward normals point from the excluded region.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortfem_generated_current_segment, only: current_segment
    use fortfem_generated_exterior_gram_density, only: exterior_gram_density
    use fortfem_generated_bem_frame, only: bem_frame
    use fortfem_generated_p1_basis, only: generated_p1_jet
    use fortfem_generated_p2_basis, only: generated_p2_jet
    use fortfem_generated_field_inner_product, only: field_inner_product
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t
    use fortfem_parametric_bem_high_order, only: build_p2_trace_connectivity
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private
    public :: evaluate_current_loops, assemble_exterior_harmonics
    public :: prescribe_exterior_periods, tangent_exterior_harmonics
    real(dp), parameter :: scale = 1/(4*acos(-1.0_dp))
contains
    subroutine evaluate_current_loops(vertices, offsets, x, potential, field, status)
        !! Each [offsets(j),offsets(j+1)) polygon closes implicitly. The caller
        !! supplies g independent interior cycles, disjoint from the exterior.
        real(dp), intent(in) :: vertices(:, :), x(:, :)
        integer, intent(in) :: offsets(:)
        real(dp), intent(out) :: potential(:, :, :), field(:, :, :)
        integer, intent(out) :: status
        real(dp) :: r1(3), r2(3), a(3), b(3), l, s, n1, n2
        integer :: g, nv, nq, j, q, e, next, errors
        status = 1; potential = 0; field = 0
        g = size(offsets) - 1; nv = size(vertices, 2); nq = size(x, 2)
        if (g < 1 .or. nq < 1 .or. size(vertices, 1) /= 3 .or. size(x, 1) /= 3) return
        if (any(shape(potential) /= [3, nq, g])) return
        if (any(shape(field) /= [3, nq, g])) return
        if (offsets(1) /= 1 .or. offsets(g + 1) /= nv + 1) return
        if (any(offsets(2:) - offsets(:g) < 3)) return
        if (.not. all(ieee_is_finite(vertices))) return
        if (.not. all(ieee_is_finite(x))) return
        do j = 1, g
            do e = offsets(j), offsets(j + 1) - 1
                next = e + 1
                if (next == offsets(j + 1)) next = offsets(j)
                if (norm2(vertices(:, next) - vertices(:, e)) <= tiny(1.0_dp)) return
            end do
        end do
        errors = 0
        !$omp parallel do collapse(2) default(none) schedule(static) &
        !$omp shared(g,nq,offsets,vertices,x,potential,field) &
        !$omp private(e,next,r1,r2,a,b,l,s,n1,n2) reduction(+:errors)
        do j = 1, g
            do q = 1, nq
                do e = offsets(j), offsets(j + 1) - 1
                    next = e + 1
                    if (next == offsets(j + 1)) next = offsets(j)
                    r1 = x(:, q) - vertices(:, e); r2 = x(:, q) - vertices(:, next)
                    l = norm2(vertices(:, e) - vertices(:, next))
                    n1 = norm2(r1); n2 = norm2(r2); s = n1 + n2
                    if (min(n1, n2) <= tiny(1.0_dp) .or. &
                        s - l <= 64*epsilon(1.0_dp)*s) then
                        errors = errors + 1
                        cycle
                    end if
                    call current_segment(r1(1), r1(2), r1(3), r2(1), r2(2), r2(3), &
                                         a(1), a(2), a(3), b(1), b(2), b(3))
                    potential(:, q, j) = potential(:, q, j) + scale*a
                    field(:, q, j) = field(:, q, j) + scale*b
                end do
            end do
        end do
        !$omp end parallel do
        status = 2
        if (errors /= 0) then
            potential = 0; field = 0
            return
        end if
        if (.not. all(ieee_is_finite(potential))) return
        if (.not. all(ieee_is_finite(field))) return
        status = 0
    end subroutine

    subroutine prescribe_exterior_periods(periods, prescribed, independent, &
                                          tolerance, coefficients, status)
        !! Select g independent cycles; validate ALL specified surface periods.
        !! Status 3: incompatible data (2g boundary periods need not extend).
        real(dp), intent(in) :: periods(:, :), prescribed(:), tolerance
        integer, intent(in) :: independent(:)
        real(dp), intent(out) :: coefficients(:)
        integer, intent(out) :: status
        real(dp), allocatable :: c(:)
        integer :: g, info, i
        status = 1; coefficients = 0; g = size(periods, 2)
        if (g < 1 .or. size(periods, 1) /= size(prescribed)) return
        if (size(coefficients) /= g .or. size(independent) /= g) return
        if (any(independent < 1) .or. any(independent > size(prescribed))) return
        if (.not. ieee_is_finite(tolerance)) return
        if (tolerance <= 0) return
        if (.not. all(ieee_is_finite(periods))) return
        if (.not. all(ieee_is_finite(prescribed))) return
        allocate (c(g))
        call dense_solve(periods(independent, :), prescribed(independent), c, info)
        status = 2
        if (info /= 0) return
        if (.not. all(ieee_is_finite(c))) return
        status = 3
        do i = 1, size(prescribed)
            if (abs(dot_product(periods(i, :), c) - prescribed(i)) > &
                tolerance*max(1.0_dp, maxval(abs(prescribed)))) return
        end do
        coefficients = c; status = 0
    end subroutine

    subroutine assemble_exterior_harmonics(surface, faces, dofs, degree, &
                            vertices, offsets, normal_load, gram, raw_asymmetry, status)
        !! C=-normal_load, G_ij=-int n.(A_i x B_j). No exterior truncation.
        !! P1 or P2 scalar traces on caller-owned arbitrary curved geometry.
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: faces(:, :), dofs(:, :), degree, offsets(:)
        real(dp), intent(in) :: vertices(:, :)
        real(dp), allocatable, intent(out) :: normal_load(:, :), gram(:, :)
        real(dp), intent(out) :: raw_asymmetry
        integer, intent(out) :: status
        real(dp), allocatable :: xi(:), eta(:), qw(:), x(:, :), normal(:, :), weights(:)
        real(dp), allocatable :: potential(:, :, :), field(:, :, :), basis(:, :)
        real(dp) :: a(3), b(3), n(3), jac, value, total, grad(2), hess(2, 2)
        integer :: nf, nq, g, nb, nt, p, t, q, i, j, k, info, errors, nv, ne
        integer, allocatable :: topology_trace(:, :), topology_flux(:, :)
        status = 1; raw_asymmetry = 0
        nf = size(faces, 2); g = size(offsets) - 1; nb = size(dofs, 1)
        if (size(faces, 1) /= 3 .or. nf < 4 .or. g < 1) return
        if (size(dofs, 2) /= nf .or. any(dofs < 1)) return
        if (nb /= 3 .and. nb /= 6) return
        call build_p2_trace_connectivity(faces, topology_trace, topology_flux, status)
        if (status /= 0) return
        status = 1
        nv = maxval(faces); ne = maxval(topology_trace) - nv
        if (g /= (2 - nv + ne - nf)/2) return
        call triangle_duffy_quadrature(degree, xi, eta, qw, status)
        if (status /= 0) return
        nq = size(qw); nt = nf*nq
        allocate (x(3, nt), normal(3, nt), weights(nt), basis(nb, nq), &
                  potential(3, nt, g), field(3, nt, g))
        do t = 1, nq
            do k = 1, nb
                if (nb == 3) then
                    call generated_p1_jet(k, xi(t), eta(t), basis(k, t), grad, hess)
                else
                    call generated_p2_jet(k, xi(t), eta(t), basis(k, t), grad, hess)
                end if
            end do
        end do
        errors = 0
        !$omp parallel do default(none) schedule(static) &
        !$omp shared(surface,nf,nq,xi,eta,qw,x,normal,weights) &
        !$omp private(t,q,a,b,n,jac,info) reduction(+:errors)
        do p = 1, nf
            do t = 1, nq
                q = (p - 1)*nq + t
                call surface%evaluate(p, [xi(t), eta(t)], x(:, q), a, b, info)
                if (info /= 0) then
                    errors = errors + 1
                    cycle
                end if
                if (.not. all(ieee_is_finite(x(:, q))) .or. &
                    .not. all(ieee_is_finite(a)) .or. &
                    .not. all(ieee_is_finite(b))) then
                    errors = errors + 1
                    cycle
                end if
                if (min(norm2(a), norm2(b)) <= tiny(1.0_dp)) then
                    errors = errors + 1
                    cycle
                end if
                if (abs(dot_product(a, b)) >= &
                    (1 - 64*epsilon(1.0_dp))*norm2(a)*norm2(b)) then
                    errors = errors + 1
                    cycle
                end if
               call bem_frame(a(1), a(2), a(3), b(1), b(2), b(3), jac, n(1), n(2), n(3))
                normal(:, q) = n; weights(q) = qw(t)*jac
            end do
        end do
        !$omp end parallel do
        status = 2
        if (errors /= 0) return
        if (.not. all(ieee_is_finite(normal))) return
        if (.not. all(ieee_is_finite(weights))) return
        if (any(weights <= 0)) return
        call evaluate_current_loops(vertices, offsets, x, potential, field, status)
        if (status /= 0) return
        allocate (normal_load(maxval(dofs), g), gram(g, g))
        normal_load = 0; gram = 0
        do j = 1, g
            do i = 1, g
                total = 0
                !$omp parallel do default(none) schedule(static) &
                !$omp shared(nt,normal,potential,field,weights,i,j) &
                !$omp private(value) reduction(+:total)
                do q = 1, nt
                    call exterior_gram_density(potential(1, q, i), potential(2, q, i), &
                   potential(3, q, i), field(1, q, j), field(2, q, j), field(3, q, j), &
                                        normal(1, q), normal(2, q), normal(3, q), value)
                    total = total + weights(q)*value
                end do
                !$omp end parallel do
                gram(i, j) = total
            end do
        end do
        !$omp parallel do default(none) schedule(static) &
        !$omp shared(g,nf,nq,nb,dofs,basis,normal,field,weights,normal_load) &
        !$omp private(p,t,q,k,value)
        do j = 1, g
            do p = 1, nf
                do t = 1, nq
                    q = (p - 1)*nq + t
                  call field_inner_product(weights(q), field(1, q, j), field(2, q, j), &
                        field(3, q, j), normal(1, q), normal(2, q), normal(3, q), value)
                    do k = 1, nb
                        normal_load(dofs(k, p), j) = normal_load(dofs(k, p), j) + &
                                                     basis(k, t)*value
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        raw_asymmetry = maxval(abs(gram - transpose(gram)))/max(tiny(1.0_dp), &
                                                                maxval(abs(gram)))
        gram = 0.5_dp*(gram + transpose(gram))
        status = 2
        if (.not. all(ieee_is_finite(gram))) return
        if (.not. all(ieee_is_finite(normal_load))) return
        status = 0
    end subroutine

    subroutine tangent_exterior_harmonics(steklov, normal_load, gram, &
                                          correction, tangent_gram, status)
        !! H_t=H+grad(w), S w=N. Periods unchanged; G_t=G-N^T S^-1 N.
        real(dp), intent(in) :: steklov(:, :), normal_load(:, :), gram(:, :)
        real(dp), allocatable, intent(out) :: correction(:, :), tangent_gram(:, :)
        integer, intent(out) :: status
        integer :: n, g
        n = size(steklov, 1); g = size(gram, 1); status = 1
        if (n < 1 .or. g < 1) return
        if (any(shape(steklov) /= [n, n])) return
        if (any(shape(normal_load) /= [n, g])) return
        if (any(shape(gram) /= [g, g])) return
        if (.not. all(ieee_is_finite(steklov))) return
        if (.not. all(ieee_is_finite(normal_load))) return
        if (.not. all(ieee_is_finite(gram))) return
        allocate (correction(n, g))
        call dense_solve(steklov, normal_load, correction, status)
        if (status /= 0) return
        tangent_gram = gram - matmul(transpose(normal_load), correction)
        tangent_gram = 0.5_dp*(tangent_gram + transpose(tangent_gram))
        status = 2
        if (.not. all(ieee_is_finite(tangent_gram))) return
        status = 0
    end subroutine
end module
