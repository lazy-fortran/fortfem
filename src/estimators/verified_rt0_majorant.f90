module fortfem_verified_rt0_majorant
    !! Independent, outward-rounded audit of a recovered RT flux on affine
    !! triangles. Only its zeroth edge moments are retained. One exactly
    !! represented real number is selected per global edge, so the mathematical
    !! RT0 Piola field has continuous normal trace even when the supplied RTp
    !! coefficients differ through local floating-point assembly. No numerical
    !! inverse of the RTp moment matrix enters this checker.
    !!
    !! For conforming u_h in H1_0, A>0 and any such sigma in H(div), integration
    !! by parts gives ||u-u_h||_A <= ||A^-1/2(sigma-A grad u_h)||
    !! + C_F ||f+div sigma||. C_F must enclose ||v||/||v||_A on H1_0.
    !! Each reference triangle is split into four congruent children per level.
    !! A callback encloses A,f,grad u_h throughout the physical coordinate box;
    !! generated interval leaves enclose both squared integrands throughout
    !! that child. Multiplication by its exact dyadic reference area encloses
    !! its integral. Interval summation and sqrt include arithmetic error.
    !! There is no assumed Gaussian exactness, estimated quadrature remainder,
    !! source projection, or local-solve accuracy requirement.
    !!
    !! The caller guarantees a nonoverlapping, conforming affine mesh and an
    !! H1_0 discrete function. Unsupported/nonfinite interval callbacks fail.
    !! This may be wider than an RTp estimate: it is an audit, not a claim that
    !! the original sampled RTp quadrature has become rigorous.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, operator(+), operator(*), sqrt
    use fortfem_equilibrated_flux_topology, only: topology_t, build_topology
    use fortfem_verified_rt0_products, only: verified_rt0_geometry, &
                                             verified_rt0_integrands
    implicit none
    private
    public :: verified_rt0_majorant, verified_rt0_problem
    abstract interface
        subroutine verified_rt0_problem(cell, x, y, coefficient, source, gradient, ok)
            import :: interval_t
            integer, intent(in) :: cell
            type(interval_t), intent(in) :: x, y
            type(interval_t), intent(out) :: coefficient, source, gradient(2)
            logical, intent(out) :: ok
        end subroutine verified_rt0_problem
    end interface
contains
    subroutine verified_rt0_majorant(vertices, triangles, edge_moments, friedrichs, &
                          problem, bound, flux_squared, residual_squared, status, level)
        real(dp), intent(in) :: vertices(:, :), edge_moments(:, :)
        integer, intent(in) :: triangles(:, :), level
        type(interval_t), intent(in) :: friedrichs
        procedure(verified_rt0_problem) :: problem
        type(interval_t), intent(out) :: bound, flux_squared, residual_squared
        integer, intent(out) :: status
        type(topology_t) :: topo
        type(interval_t), allocatable :: cell_flux(:), cell_residual(:)
        real(dp), allocatable :: shared(:), moments(:, :)
        integer, allocatable :: count_edge(:), orientation_sum(:)
        logical, allocatable :: seen(:)
        integer :: nt, nv, e, t, edge, stat, failed, nedges
        real(dp) :: reference(2, 3)

        status = 1
        bound = interval(0.0_dp, huge(1.0_dp))
        flux_squared = bound
        residual_squared = bound
        if (size(vertices, 1) /= 2 .or. size(triangles, 1) /= 3) return
        nv = size(vertices, 2); nt = size(triangles, 2)
        if (nv < 3 .or. nt < 1 .or. level < 0 .or. level > 5) return
        if (any(triangles < 1) .or. any(triangles > nv)) return
        if (size(edge_moments, 1) /= 3 .or. size(edge_moments, 2) /= nt) return
        if (.not. all(ieee_is_finite(vertices))) return
        if (.not. all(ieee_is_finite(edge_moments))) return
        if (.not. valid(friedrichs) .or. friedrichs%lo < 0) return
        call build_topology(nv, triangles, topo, stat)
        if (stat /= 0) return
        nedges = maxval(topo%cell_edge)
        allocate (shared(nedges), seen(nedges), count_edge(nedges), &
                  orientation_sum(nedges), moments(3, nt))
        seen = .false.; count_edge = 0; orientation_sum = 0
        do t = 1, nt
            do e = 1, 3
                edge = topo%cell_edge(e, t)
                count_edge(edge) = count_edge(edge) + 1
                orientation_sum(edge) = orientation_sum(edge) &
                                        + topo%cell_orientation(e, t)
                if (.not. seen(edge)) then
                    shared(edge) = edge_moments(e, t)*topo%cell_orientation(e, t)
                    seen(edge) = .true.
                end if
                moments(e, t) = shared(edge)*topo%cell_orientation(e, t)
            end do
        end do
        if (any(count_edge > 2)) return
        do edge = 1, nedges
            if (count_edge(edge) == 2 .and. orientation_sum(edge) /= 0) return
        end do
        reference = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                             0.0_dp, 1.0_dp], [2, 3])
        allocate (cell_flux(nt), cell_residual(nt))
        failed = 0
        !$omp parallel do schedule(static) reduction(max:failed) private(stat)
        do t = 1, nt
            cell_flux(t) = interval(0.0_dp)
            cell_residual(t) = interval(0.0_dp)
            call integrate(t, reference, level, 0.5_dp/4.0_dp**level, &
                           cell_flux(t), cell_residual(t), stat)
            if (stat /= 0) failed = 1
        end do
        !$omp end parallel do
        if (failed /= 0) return
        flux_squared = interval(0.0_dp); residual_squared = interval(0.0_dp)
        do t = 1, nt
            flux_squared = flux_squared + cell_flux(t)
            residual_squared = residual_squared + cell_residual(t)
        end do
        flux_squared%lo = max(0.0_dp, flux_squared%lo)
        residual_squared%lo = max(0.0_dp, residual_squared%lo)
        bound = sqrt(flux_squared) + friedrichs*sqrt(residual_squared)
        if (.not. valid(bound)) return
        status = 0
    contains
        recursive subroutine integrate(cell, ref, depth, area, flux, residual, stat)
            integer, intent(in) :: cell, depth
            real(dp), intent(in) :: ref(2, 3), area
            type(interval_t), intent(inout) :: flux, residual
            integer, intent(out) :: stat
            type(interval_t) :: xy(2, 3), xi, eta, x, y, sx, sy, div, det
            type(interval_t) :: coefficient, source, gradient(2), fi, ri
            real(dp) :: midpoint(2, 3), child(2, 3)
            integer :: j, next, child_stat
            logical :: ok
            stat = 1
            if (depth > 0) then
                do j = 1, 3
                    next = mod(j, 3) + 1
                    midpoint(:, j) = (ref(:, j) + ref(:, next))*0.5_dp
                end do
                do j = 1, 4
                    select case (j)
                    case (1)
                        child(:, 1) = ref(:, 1); child(:, 2) = midpoint(:, 1)
                        child(:, 3) = midpoint(:, 3)
                    case (2)
                        child(:, 1) = midpoint(:, 1); child(:, 2) = ref(:, 2)
                        child(:, 3) = midpoint(:, 2)
                    case (3)
                        child(:, 1) = midpoint(:, 3); child(:, 2) = midpoint(:, 2)
                        child(:, 3) = ref(:, 3)
                    case (4)
                        child = midpoint
                    end select
                    call integrate(cell, child, depth - 1, area, flux, residual, &
                                   child_stat)
                    if (child_stat /= 0) return
                end do
                stat = 0
                return
            end if
            do j = 1, 3
                xy(:, j) = interval(vertices(:, triangles(j, cell)))
            end do
            xi = interval(minval(ref(1, :)), maxval(ref(1, :)))
            eta = interval(minval(ref(2, :)), maxval(ref(2, :)))
           call verified_rt0_geometry(xi, eta, xy(1, 1), xy(2, 1), xy(1, 2), xy(2, 2), &
                                       xy(1, 3), xy(2, 3), interval(moments(1, cell)), &
                               interval(moments(2, cell)), interval(moments(3, cell)), &
                                       x, y, sx, sy, div, det)
            if (.not. valid(det)) return
            if (det%lo <= 0) return
            call problem(cell, x, y, coefficient, source, gradient, ok)
            if (.not. ok) return
            if (.not. valid(coefficient)) return
            if (coefficient%lo <= 0) return
            if (.not. valid(source)) return
            if (.not. all(valid(gradient))) return
           call verified_rt0_integrands(coefficient, source, gradient(1), gradient(2), &
                                         sx, sy, div, det, fi, ri)
            if (.not. valid(fi) .or. .not. valid(ri)) return
            fi%lo = max(0.0_dp, fi%lo); ri%lo = max(0.0_dp, ri%lo)
            flux = flux + fi*area
            residual = residual + ri*area
            stat = 0
        end subroutine integrate
    end subroutine verified_rt0_majorant

    elemental logical function valid(value)
        type(interval_t), intent(in) :: value
        valid = ieee_is_finite(value%lo) .and. ieee_is_finite(value%hi) &
                .and. value%lo <= value%hi
    end function valid
end module fortfem_verified_rt0_majorant
