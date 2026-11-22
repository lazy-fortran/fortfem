program test_spline_derham_t_basis
    !! spline_derham_t_basis against the defining harmonic sums: the
    !! trigonometric cardinal basis L_k(t) = (1 + 2 sum_m cos(2 pi m x))/n and
    !! the histopolation basis H_k(t) = 1 + 2 sum_m cos(2 pi m (x - h/2))/
    !! sinc(pi m h), x = t - (k - 1) h, h = 1/n, m = 1..(n - 1)/2, with their
    !! t-derivatives; plus the cardinal property L_k(t_j) = delta_kj, the
    !! partitions of unity sum L_k = 1 and h sum H_k = 1, and L_k' =
    !! H_(k-1) - H_k (the commuting derivative of the de Rham complex).
    use fortfem_kinds, only: dp
    use fortfem_spline_derham, only: spline_derham_t_basis
    use check, only: check_condition, check_summary
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: nt, i, k, j
    real(dp) :: t, err_nodal, err_edge, err_card, err_unity, err_commute, scale
    real(dp), allocatable :: nodal(:, :), edge(:, :), rn(:, :), re(:, :)

    err_nodal = 0
    err_edge = 0
    err_card = 0
    err_unity = 0
    err_commute = 0
    do nt = 1, 61, 2
        allocate (nodal(nt, 2), edge(nt, 2), rn(nt, 2), re(nt, 2))
        scale = 1 + 4*pi*nt**2
        do i = 0, 400
            t = -1.0_dp + 3.0_dp*i/400.0_dp
            if (i == 7) t = 1.0_dp - 1e-15_dp
            if (i == 8) t = 1e-15_dp
            call spline_derham_t_basis(nt, t, nodal, edge)
            call reference(nt, t, rn, re)
            err_nodal = max(err_nodal, maxval(abs(nodal - rn))/scale)
            err_edge = max(err_edge, maxval(abs(edge - re))/(scale*nt))
            err_unity = max(err_unity, abs(sum(nodal(:, 1)) - 1), abs(sum(edge(:, 1))/nt - 1), &
                            abs(sum(nodal(:, 2)))/scale)
            do k = 1, nt
                err_commute = max(err_commute, abs(nodal(k, 2) &
                                  - (edge(modulo(k - 2, nt) + 1, 1) - edge(k, 1)))/(scale*nt))
            end do
        end do
        do j = 1, nt
            call spline_derham_t_basis(nt, (j - 1)/real(nt, dp), nodal, edge)
            do k = 1, nt
                err_card = max(err_card, abs(nodal(k, 1) - merge(1.0_dp, 0.0_dp, k == j)))
            end do
        end do
        deallocate (nodal, edge, rn, re)
    end do
    print '(a,5es10.2)', 'nodal, edge, cardinal, unity, commuting errors', err_nodal, err_edge, &
        err_card, err_unity, err_commute
    call check_condition(err_nodal < 1e-13_dp, 'nodal basis and derivative equal the harmonic sum')
    call check_condition(err_edge < 1e-13_dp, 'edge basis and derivative equal the harmonic sum')
    call check_condition(err_card < 1e-13_dp, 'cardinal property at the nodes')
    call check_condition(err_unity < 1e-13_dp, 'partitions of unity')
    call check_condition(err_commute < 1e-13_dp, 'L_k'' = H_(k-1) - H_k')
    call check_summary('spline_derham_t_basis')

contains

    subroutine reference(n, tt, nod, edg)
        integer, intent(in) :: n
        real(dp), intent(in) :: tt
        real(dp), intent(out) :: nod(n, 2), edg(n, 2)
        real(dp) :: h, x, s
        integer :: kk, m

        h = 1.0_dp/n
        do kk = 1, n
            x = tt - (kk - 1)*h
            nod(kk, :) = [1.0_dp, 0.0_dp]
            edg(kk, :) = [1.0_dp, 0.0_dp]
            do m = 1, (n - 1)/2
                s = sin(pi*m*h)/(pi*m*h)
                nod(kk, 1) = nod(kk, 1) + 2*cos(2*pi*m*x)
                nod(kk, 2) = nod(kk, 2) - 4*pi*m*sin(2*pi*m*x)
                edg(kk, 1) = edg(kk, 1) + 2*cos(2*pi*m*(x - 0.5_dp*h))/s
                edg(kk, 2) = edg(kk, 2) - 4*pi*m*sin(2*pi*m*(x - 0.5_dp*h))/s
            end do
            nod(kk, :) = nod(kk, :)/n
        end do
    end subroutine reference

end program test_spline_derham_t_basis
