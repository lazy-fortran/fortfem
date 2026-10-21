program test_verified_mapped_majorant
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortnum_interval, only: interval_t, interval
    use fortfem_verified_bernstein_majorant, only: &
        verified_mapped_bernstein_majorant, verified_bernstein_majorant
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_triangle_rt_arbitrary_order, only: triangle_rt_basis_t, &
                                                   evaluate_triangle_raviart_thomas
    use equilibrated_test_mesh, only: test_mesh_t, build_test_mesh, cell_point, &
                                      random_rt_flux, edge_jump
    implicit none
    real(dp) :: vertices(2, 3), geometry(2, 10, 1), solution(10), flux(24, 1)
    real(dp) :: weight(10), a, b, exact, gx, gy, j21, j22, x, y, scale
    real(dp) :: previous, error, effectivity
    real(dp), allocatable :: xi(:), eta(:), qw(:)
    integer :: triangles(3, 1), cells(10, 1), status, order, level, q, j, ng
    type(interval_t) :: bound, fs, rs, old_bound, old_fs, old_rs
    real(dp), parameter :: nodes(2, 10) = reshape([ &
                                             0._dp, 0._dp, 1._dp, 0._dp, 0._dp, 1._dp, &
                                     1._dp/3, 0._dp, 2._dp/3, 0._dp, 2._dp/3, 1._dp/3, &
                                     1._dp/3, 2._dp/3, 0._dp, 2._dp/3, 0._dp, 1._dp/3, &
                                                  1._dp/3, 1._dp/3], [2, 10])
    real(dp), parameter :: p2nodes(2, 6) = reshape([ &
                                             0._dp, 0._dp, 1._dp, 0._dp, 0._dp, 1._dp, &
                                      .5_dp, 0._dp, .5_dp, .5_dp, 0._dp, .5_dp], [2, 6])

    triangles(:, 1) = [1, 2, 3]
    cells(:, 1) = [(j, j=1, 10)]
    vertices = nodes(:, :3)
    flux = 0; weight = 1
    call triangle_duffy_quadrature(40, xi, eta, qw, status)
    call check_condition(status == 0, 'independent energy quadrature')
    ! Manufactured exact solution u=0, f=0; the conforming approximation is
    ! scale*xi*eta*(1-xi-eta). Its exact error is its energy seminorm.
    ! Curvature stays nonzero through every level. No solve or fitted value
    ! supplies the oracle. The polynomial Jacobians below are independent.
    do order = 2, 3
        ng = (order + 1)*(order + 2)/2
        a = .125_dp; b = 0
        if (order == 3) b = .0625_dp
        do j = 1, ng
            if (order == 2) then
                x = p2nodes(1, j); y = p2nodes(2, j)
            else
                x = nodes(1, j); y = nodes(2, j)
            end if
            geometry(:, j, 1) = [x, y + a*x*y + b*x*y*(1 - x - y)]
        end do
        previous = 0
        do level = 0, 2
            scale = 0.5_dp**(order*level)
            solution = scale*nodes(1, :)*nodes(2, :) &
                       *(1 - nodes(1, :) - nodes(2, :))
            call verified_mapped_bernstein_majorant(vertices, triangles, 3, &
                              flux, cells, solution, geometry(:, :ng, :), weight(:ng), &
                              zero_source, 0, 0, interval(1._dp), bound, fs, rs, status)
            call check_condition(status == 0, 'curved P2/P3 audit accepted')
            exact = 0
            do q = 1, size(qw)
                x = xi(q); y = eta(q)
                j21 = a*y + b*y*(1 - 2*x - y)
                j22 = 1 + a*x + b*x*(1 - x - 2*y)
                gx = scale*y*(1 - 2*x - y)
                gy = scale*x*(1 - x - 2*y)
                exact = exact + qw(q)*j22 &
                        *((gx - j21*gy/j22)**2 + (gy/j22)**2)
            end do
            error = sqrt(exact); effectivity = bound%hi/error
            call check_condition(fs%lo <= exact .and. exact <= fs%hi, &
                             'rational integral enclosure contains exact error squared')
            call check_condition(effectivity >= 1 .and. effectivity < 1.3_dp, &
                                 'curved manufactured effectivity')
            if (level > 0) then
                call check_condition(abs(log(previous/bound%hi)/log(2._dp) &
                           - order) < 1.e-6_dp, 'manufactured observed amplitude order')
            end if
            previous = bound%hi
            print '(a,2i2,3es14.5)', 'mapped P,level,error,bound,effectivity ', &
                order, level, error, bound%hi, effectivity
        end do
    end do
    geometry(:, :, 1) = nodes
    solution = nodes(1, :)*nodes(2, :)*(1 - nodes(1, :) - nodes(2, :))
    call verified_mapped_bernstein_majorant(vertices, triangles, 3, flux, &
                cells, solution, geometry, weight, zero_source, 0, 0, interval(1._dp), &
                                            bound, fs, rs, status)
    call check_condition(status == 0, 'affine mapped limit')
    call verified_bernstein_majorant(vertices, triangles, 3, flux, cells, &
                             solution, weight(:3), zero_source, 0, 0, interval(1._dp), &
                                     old_bound, old_fs, old_rs, status)
    call check_condition(status == 0, 'original affine checker')
    call check_condition(abs(bound%hi - old_bound%hi) < 1.e-10_dp, &
                         'affine-limit agreement')
    ! Positive at all three vertices, but negative near xi=.37.
    do j = 1, 10
        x = nodes(1, j); y = nodes(2, j)
        geometry(:, j, 1) = [x, y*((x - .37_dp)**2 - .0001_dp)]
    end do
    vertices = geometry(:, :3, 1)
    call verified_mapped_bernstein_majorant(vertices, triangles, 3, flux, &
                cells, solution, geometry, weight, zero_source, 0, 0, interval(1._dp), &
                                            bound, fs, rs, status)
    call check_condition(status /= 0, 'between-vertex Jacobian inversion refused')
    call shared_edge_rejection()
    call nonzero_piola_field()
    call check_summary('verified mapped Bernstein majorant')
contains
    subroutine nonzero_piola_field()
        type(test_mesh_t) :: mesh
        type(triangle_rt_basis_t) :: rt
        real(dp), allocatable :: coeff(:, :), w(:), u(:)
        real(dp) :: rv(2, 8), rd(8), shat(2), sigma(2), pos(2), jac(2, 2), det
        real(dp) :: phi(10), gr(2, 10), flux_integral, residual_integral, jump
        integer :: m, stat, cell, q, i, e, other, oe
        do m = 2, 3
            call build_test_mesh(2, m, .04_dp, mesh, stat)
            call check_condition(stat == 0, 'curved shared-edge test mesh')
            call random_rt_flux(mesh, 1, rt, coeff)
            allocate (w(mesh%nb*mesh%nt), u(mesh%ndof))
            w = reshape(1 + mesh%nodes(1, :, :), [size(w)]); u = 0
            flux_integral = 0; residual_integral = 0; jump = 0
            do cell = 1, mesh%nt
                do e = 1, 3
                    do other = cell + 1, mesh%nt
                        do oe = 1, 3
                            if (mesh%triangles(e, cell) &
                                /= mesh%triangles(mod(oe, 3) + 1, other)) cycle
                            if (mesh%triangles(mod(e, 3) + 1, cell) &
                                /= mesh%triangles(oe, other)) cycle
                            jump = max(jump, edge_jump(mesh, rt, coeff, &
                                                       cell, e, other, oe))
                        end do
                    end do
                end do
                do q = 1, size(qw)
                    call cell_point(mesh, cell, xi(q), eta(q), pos, jac, det, &
                                    phi(:mesh%nb), gr(:, :mesh%nb))
                    call evaluate_triangle_raviart_thomas(rt, xi(q), eta(q), &
                                                          rv, rd, stat)
                    do i = 1, 2
                        shat(i) = dot_product(rv(i, :), coeff(:, cell))
                    end do
                    do i = 1, 2
                        sigma(i) = dot_product(jac(i, :), shat)/det
                    end do
                    flux_integral = flux_integral &
                                    + qw(q)*det*(1 + pos(1))*sum(sigma**2)
                    residual_integral = residual_integral &
                                        + qw(q)*dot_product(rd, coeff(:, cell))**2/det
                end do
            end do
            call check_condition(jump < 1.e-12_dp, 'independent Piola normal trace')
            call verified_mapped_bernstein_majorant(mesh%vertices, mesh%triangles, &
                             1, coeff, mesh%gdof, u, mesh%nodes, w, zero_source, 0, 0, &
                                                   interval(1._dp), bound, fs, rs, stat)
            call check_condition(stat == 0, 'nonzero reference field audit')
            call check_condition(fs%lo <= flux_integral .and. &
                                 fs%hi >= flux_integral, 'nonzero rational flux oracle')
            call check_condition(rs%lo <= residual_integral .and. &
                       rs%hi >= residual_integral, 'nonzero rational divergence oracle')
            deallocate (w, u)
        end do
    end subroutine

    subroutine zero_source(cell, coefficients, ok)
        integer, intent(in) :: cell
        type(interval_t), intent(out) :: coefficients(:)
        logical, intent(out) :: ok
        coefficients = interval(0._dp)
        ok = cell >= 1
    end subroutine

    subroutine shared_edge_rejection()
        real(dp) :: xy(2, 4), geo(2, 6, 2), fc(8, 2), u(4), w(6, 2)
        integer :: tri(3, 2), i, e, next
        xy = reshape([0._dp, 0._dp, 1._dp, 0._dp, &
                      1._dp, 1._dp, 0._dp, 1._dp], [2, 4])
        tri = reshape([1, 2, 3, 1, 3, 4], [3, 2])
        do e = 1, 2
            geo(:, :3, e) = xy(:, tri(:, e))
            do i = 1, 3
                next = mod(i, 3) + 1
                geo(:, 3 + i, e) = .5_dp*(geo(:, i, e) + geo(:, next, e))
            end do
        end do
        fc = 0; u = 0; w = 1
        call verified_mapped_bernstein_majorant(xy, tri, 1, fc, tri, u, geo, &
                                 reshape(w, [12]), zero_source, 0, 0, interval(1._dp), &
                                                bound, fs, rs, status)
        call check_condition(status == 0, 'shared mapped edge accepted')
        geo(1, 6, 1) = geo(1, 6, 1) + .001_dp
        call verified_mapped_bernstein_majorant(xy, tri, 1, fc, tri, u, geo, &
                                 reshape(w, [12]), zero_source, 0, 0, interval(1._dp), &
                                                bound, fs, rs, status)
        call check_condition(status /= 0, 'inconsistent shared-edge map refused')
    end subroutine
end program
