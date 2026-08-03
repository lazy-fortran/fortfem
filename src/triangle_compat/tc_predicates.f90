module tc_predicates
    !! Adaptive exact geometric predicates after Shewchuk's "Adaptive
    !! Precision Floating-Point Arithmetic and Fast Robust Geometric
    !! Predicates".  The orientation test reproduces the reference
    !! implementation value-exactly (its result is consumed numerically by
    !! the circumcenter routine); the incircle test reproduces the fast
    !! path bit-exactly (including FMA contraction observed in the oracle
    !! binary) and falls back to a fully exact expansion computation whose
    !! sign matches the reference adaptive evaluation.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_fma
    implicit none
    private

    public :: tc_counterclockwise, tc_incircle

    real(dp), parameter :: EPS_TC = 2.0_dp**(-53)
    real(dp), parameter :: SPLITTER = 134217729.0_dp
    real(dp), parameter :: RESULTERRBOUND = (3.0_dp + 8.0_dp*EPS_TC)*EPS_TC
    real(dp), parameter :: CCWERRBOUND_A = (3.0_dp + 16.0_dp*EPS_TC)*EPS_TC
    real(dp), parameter :: CCWERRBOUND_B = (2.0_dp + 12.0_dp*EPS_TC)*EPS_TC
    real(dp), parameter :: CCWERRBOUND_C = (9.0_dp + 64.0_dp*EPS_TC)*EPS_TC*EPS_TC
    real(dp), parameter :: ICCERRBOUND_A = (10.0_dp + 96.0_dp*EPS_TC)*EPS_TC

contains

    function tc_counterclockwise(ax, ay, bx, by, cx, cy) result(det)
        real(dp), intent(in) :: ax, ay, bx, by, cx, cy
        real(dp) :: det
        real(dp), volatile :: detleft, detright
        real(dp) :: detsum, errbound

        detleft = (ax - cx)*(by - cy)
        detright = (ay - cy)*(bx - cx)
        det = detleft - detright

        if (detleft > 0.0_dp) then
            if (detright <= 0.0_dp) then
                return
            else
                detsum = detleft + detright
            end if
        else if (detleft < 0.0_dp) then
            if (detright >= 0.0_dp) then
                return
            else
                detsum = -detleft - detright
            end if
        else
            return
        end if

        errbound = CCWERRBOUND_A*detsum
        if ((det >= errbound) .or. (-det >= errbound)) return

        det = ccw_adapt(ax, ay, bx, by, cx, cy, detsum)
    end function tc_counterclockwise

    function ccw_adapt(ax, ay, bx, by, cx, cy, detsum) result(det)
        real(dp), intent(in) :: ax, ay, bx, by, cx, cy, detsum
        real(dp) :: det
        real(dp) :: acx, acy, bcx, bcy
        real(dp) :: acxtail, acytail, bcxtail, bcytail
        real(dp) :: detleft, detlefttail, detright, detrighttail
        real(dp) :: errbound
        real(dp) :: b4(4), u(4), c1(8), c2(12), d16(16)
        real(dp), volatile :: pos1, pos2, neg1, neg2, possum, negsum
        integer :: c1len, c2len, dlen
        real(dp) :: s1, s0, t1, t0

        acx = ax - cx
        bcx = bx - cx
        acy = ay - cy
        bcy = by - cy

        call two_product(acx, bcy, detleft, detlefttail)
        call two_product(acy, bcx, detright, detrighttail)
        call two_two_diff(detleft, detlefttail, detright, detrighttail, &
                          b4(4), b4(3), b4(2), b4(1))

        det = estimate(4, b4)
        errbound = CCWERRBOUND_B*detsum
        if ((det >= errbound) .or. (-det >= errbound)) return

        call two_diff_tail(ax, cx, acx, acxtail)
        call two_diff_tail(bx, cx, bcx, bcxtail)
        call two_diff_tail(ay, cy, acy, acytail)
        call two_diff_tail(by, cy, bcy, bcytail)

        if ((acxtail == 0.0_dp) .and. (acytail == 0.0_dp) &
            .and. (bcxtail == 0.0_dp) .and. (bcytail == 0.0_dp)) return

        errbound = CCWERRBOUND_C*detsum + RESULTERRBOUND*abs(det)
        pos1 = acx*bcytail
        pos2 = bcy*acxtail
        neg1 = acy*bcxtail
        neg2 = bcx*acytail
        possum = pos1 + pos2
        negsum = neg1 + neg2
        det = det + (possum - negsum)
        if ((det >= errbound) .or. (-det >= errbound)) return

        call two_product(acxtail, bcy, s1, s0)
        call two_product(acytail, bcx, t1, t0)
        call two_two_diff(s1, s0, t1, t0, u(4), u(3), u(2), u(1))
        c1len = fast_expansion_sum_zeroelim(4, b4, 4, u, c1)

        call two_product(acx, bcytail, s1, s0)
        call two_product(acy, bcxtail, t1, t0)
        call two_two_diff(s1, s0, t1, t0, u(4), u(3), u(2), u(1))
        c2len = fast_expansion_sum_zeroelim(c1len, c1, 4, u, c2)

        call two_product(acxtail, bcytail, s1, s0)
        call two_product(acytail, bcxtail, t1, t0)
        call two_two_diff(s1, s0, t1, t0, u(4), u(3), u(2), u(1))
        dlen = fast_expansion_sum_zeroelim(c2len, c2, 4, u, d16)

        det = d16(dlen)
    end function ccw_adapt

    function tc_incircle(ax, ay, bx, by, cx, cy, dx, dy) result(det)
        real(dp), intent(in) :: ax, ay, bx, by, cx, cy, dx, dy
        real(dp) :: det
        real(dp) :: adx, bdx, cdx, ady, bdy, cdy
        real(dp) :: bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady
        real(dp) :: alift, blift, clift
        real(dp) :: permanent, errbound

        adx = ax - dx
        bdx = bx - dx
        cdx = cx - dx
        ady = ay - dy
        bdy = by - dy
        cdy = cy - dy

        bdxcdy = bdx*cdy
        cdxbdy = cdx*bdy
        alift = ieee_fma(adx, adx, ady*ady)

        cdxady = cdx*ady
        adxcdy = adx*cdy
        blift = ieee_fma(bdx, bdx, bdy*bdy)

        adxbdy = adx*bdy
        bdxady = bdx*ady
        clift = ieee_fma(cdx, cdx, cdy*cdy)

        det = ieee_fma(clift, adxbdy - bdxady, &
                       ieee_fma(alift, bdxcdy - cdxbdy, blift*(cdxady - adxcdy)))

        permanent = ieee_fma(abs(adxbdy) + abs(bdxady), clift, &
                             ieee_fma(abs(bdxcdy) + abs(cdxbdy), alift, &
                                      blift*(abs(cdxady) + abs(adxcdy))))
        errbound = permanent*ICCERRBOUND_A
        if ((det > errbound) .or. (-det > errbound)) return

        det = incircle_exact(ax, ay, bx, by, cx, cy, dx, dy)
    end function tc_incircle

    function incircle_exact(ax, ay, bx, by, cx, cy, dx, dy) result(det)
        !! Fully exact incircle determinant via expansion arithmetic.
        !! Only the sign (and exact zero) of the result is meaningful.
        real(dp), intent(in) :: ax, ay, bx, by, cx, cy, dx, dy
        real(dp) :: det
        real(dp) :: adx(2), ady(2), bdx(2), bdy(2), cdx(2), cdy(2)
        real(dp) :: bc(64), ca(64), ab(64)
        real(dp) :: alift(32), blift(32), clift(32)
        real(dp) :: ta(2048), tb(2048), tc_(2048), tsum(4096)
        integer :: nbc, nca, nab, nal, nbl, ncl, na, nb, nc, ns

        call two_diff_e(ax, dx, adx)
        call two_diff_e(ay, dy, ady)
        call two_diff_e(bx, dx, bdx)
        call two_diff_e(by, dy, bdy)
        call two_diff_e(cx, dx, cdx)
        call two_diff_e(cy, dy, cdy)

        call cross_e(bdx, bdy, cdx, cdy, bc, nbc)
        call cross_e(cdx, cdy, adx, ady, ca, nca)
        call cross_e(adx, ady, bdx, bdy, ab, nab)

        call lift_e(adx, ady, alift, nal)
        call lift_e(bdx, bdy, blift, nbl)
        call lift_e(cdx, cdy, clift, ncl)

        call exp_product(nal, alift, nbc, bc, ta, na)
        call exp_product(nbl, blift, nca, ca, tb, nb)
        call exp_product(ncl, clift, nab, ab, tc_, nc)

        ns = fast_expansion_sum_zeroelim(na, ta, nb, tb, tsum)
        na = fast_expansion_sum_zeroelim(ns, tsum, nc, tc_, ta)
        det = ta(na)
    end function incircle_exact

    subroutine cross_e(ux, uy, vx, vy, h, n)
        !! h = ux*vy - uy*vx for two-component operands, exact.
        real(dp), intent(in) :: ux(2), uy(2), vx(2), vy(2)
        real(dp), intent(out) :: h(:)
        integer, intent(out) :: n
        real(dp) :: p1(16), p2(16), p2n(16)
        integer :: n1, n2, i

        call exp_product(2, ux, 2, vy, p1, n1)
        call exp_product(2, uy, 2, vx, p2, n2)
        do i = 1, n2
            p2n(i) = -p2(i)
        end do
        n = fast_expansion_sum_zeroelim(n1, p1, n2, p2n, h)
    end subroutine cross_e

    subroutine lift_e(ux, uy, h, n)
        !! h = ux*ux + uy*uy for two-component operands, exact.
        real(dp), intent(in) :: ux(2), uy(2)
        real(dp), intent(out) :: h(:)
        integer, intent(out) :: n
        real(dp) :: p1(16), p2(16)
        integer :: n1, n2

        call exp_product(2, ux, 2, ux, p1, n1)
        call exp_product(2, uy, 2, uy, p2, n2)
        n = fast_expansion_sum_zeroelim(n1, p1, n2, p2, h)
    end subroutine lift_e

    subroutine exp_product(ne, e, nf, f, h, nh)
        !! Exact product of two expansions: scale e by each component of f
        !! and accumulate.
        real(dp), intent(in) :: e(:), f(:)
        integer, intent(in) :: ne, nf
        real(dp), intent(out) :: h(:)
        integer, intent(out) :: nh
        real(dp) :: scaled(2*size(e)), acc(size(h)), tmp(size(h))
        integer :: i, nsc, nacc

        nacc = 0
        do i = 1, nf
            nsc = scale_expansion_zeroelim(ne, e, f(i), scaled)
            if (nacc == 0) then
                acc(1:nsc) = scaled(1:nsc)
                nacc = nsc
            else
                nacc = fast_expansion_sum_zeroelim(nacc, acc, nsc, scaled, tmp)
                acc(1:nacc) = tmp(1:nacc)
            end if
        end do
        if (nacc == 0) then
            h(1) = 0.0_dp
            nh = 1
        else
            h(1:nacc) = acc(1:nacc)
            nh = nacc
        end if
    end subroutine exp_product

    subroutine two_diff_e(a, b, x)
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: x(2)
        call two_diff(a, b, x(2), x(1))
    end subroutine two_diff_e

    subroutine two_sum(a, b, x, y)
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: x, y
        real(dp) :: bvirt, avirt, bround, around
        x = a + b
        bvirt = x - a
        avirt = x - bvirt
        bround = b - bvirt
        around = a - avirt
        y = around + bround
    end subroutine two_sum

    subroutine fast_two_sum(a, b, x, y)
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: x, y
        real(dp) :: bvirt
        x = a + b
        bvirt = x - a
        y = b - bvirt
    end subroutine fast_two_sum

    subroutine two_diff(a, b, x, y)
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: x, y
        x = a - b
        call two_diff_tail(a, b, x, y)
    end subroutine two_diff

    subroutine two_diff_tail(a, b, x, y)
        real(dp), intent(in) :: a, b, x
        real(dp), intent(out) :: y
        real(dp) :: bvirt, avirt, bround, around
        bvirt = a - x
        avirt = x + bvirt
        bround = bvirt - b
        around = a - avirt
        y = around + bround
    end subroutine two_diff_tail

    subroutine two_product(a, b, x, y)
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: x, y
        real(dp) :: ahi, alo, bhi, blo, err1, err2, err3
        x = a*b
        call split(a, ahi, alo)
        call split(b, bhi, blo)
        err1 = x - ahi*bhi
        err2 = err1 - alo*bhi
        err3 = err2 - ahi*blo
        y = alo*blo - err3
    end subroutine two_product

    subroutine split(a, ahi, alo)
        real(dp), intent(in) :: a
        real(dp), intent(out) :: ahi, alo
        real(dp) :: c, abig
        c = SPLITTER*a
        abig = c - a
        ahi = c - abig
        alo = a - ahi
    end subroutine split

    subroutine two_one_diff(a1, a0, b, x2, x1, x0)
        real(dp), intent(in) :: a1, a0, b
        real(dp), intent(out) :: x2, x1, x0
        real(dp) :: i_
        call two_diff(a0, b, i_, x0)
        call two_sum(a1, i_, x2, x1)
    end subroutine two_one_diff

    subroutine two_two_diff(a1, a0, b1, b0, x3, x2, x1, x0)
        real(dp), intent(in) :: a1, a0, b1, b0
        real(dp), intent(out) :: x3, x2, x1, x0
        real(dp) :: j_, t0
        call two_one_diff(a1, a0, b0, j_, t0, x0)
        call two_one_diff(j_, t0, b1, x3, x2, x1)
    end subroutine two_two_diff

    function estimate(elen, e) result(q)
        integer, intent(in) :: elen
        real(dp), intent(in) :: e(:)
        real(dp) :: q
        integer :: i
        q = e(1)
        do i = 2, elen
            q = q + e(i)
        end do
    end function estimate

    function fast_expansion_sum_zeroelim(elen, e, flen, f, h) result(hlen)
        integer, intent(in) :: elen, flen
        real(dp), intent(in) :: e(:), f(:)
        real(dp), intent(out) :: h(:)
        integer :: hlen
        real(dp) :: q, qnew, hh, enow, fnow
        integer :: eindex, findex, hindex

        enow = e(1)
        fnow = f(1)
        eindex = 1
        findex = 1
        if ((fnow > enow) .eqv. (fnow > -enow)) then
            q = enow
            eindex = eindex + 1
            if (eindex <= elen) enow = e(eindex)
        else
            q = fnow
            findex = findex + 1
            if (findex <= flen) fnow = f(findex)
        end if
        hindex = 0
        if ((eindex <= elen) .and. (findex <= flen)) then
            if ((fnow > enow) .eqv. (fnow > -enow)) then
                call fast_two_sum(enow, q, qnew, hh)
                eindex = eindex + 1
                if (eindex <= elen) enow = e(eindex)
            else
                call fast_two_sum(fnow, q, qnew, hh)
                findex = findex + 1
                if (findex <= flen) fnow = f(findex)
            end if
            q = qnew
            if (hh /= 0.0_dp) then
                hindex = hindex + 1
                h(hindex) = hh
            end if
            do while ((eindex <= elen) .and. (findex <= flen))
                if ((fnow > enow) .eqv. (fnow > -enow)) then
                    call two_sum(q, enow, qnew, hh)
                    eindex = eindex + 1
                    if (eindex <= elen) enow = e(eindex)
                else
                    call two_sum(q, fnow, qnew, hh)
                    findex = findex + 1
                    if (findex <= flen) fnow = f(findex)
                end if
                q = qnew
                if (hh /= 0.0_dp) then
                    hindex = hindex + 1
                    h(hindex) = hh
                end if
            end do
        end if
        do while (eindex <= elen)
            call two_sum(q, enow, qnew, hh)
            eindex = eindex + 1
            if (eindex <= elen) enow = e(eindex)
            q = qnew
            if (hh /= 0.0_dp) then
                hindex = hindex + 1
                h(hindex) = hh
            end if
        end do
        do while (findex <= flen)
            call two_sum(q, fnow, qnew, hh)
            findex = findex + 1
            if (findex <= flen) fnow = f(findex)
            q = qnew
            if (hh /= 0.0_dp) then
                hindex = hindex + 1
                h(hindex) = hh
            end if
        end do
        if ((q /= 0.0_dp) .or. (hindex == 0)) then
            hindex = hindex + 1
            h(hindex) = q
        end if
        hlen = hindex
    end function fast_expansion_sum_zeroelim

    function scale_expansion_zeroelim(elen, e, b, h) result(hlen)
        integer, intent(in) :: elen
        real(dp), intent(in) :: e(:), b
        real(dp), intent(out) :: h(:)
        integer :: hlen
        real(dp) :: q, sum_, hh, product1, product0, enow
        real(dp) :: bhi, blo, ahi, alo, err1, err2, err3
        integer :: eindex, hindex

        call split(b, bhi, blo)
        product1 = e(1)*b
        call split(e(1), ahi, alo)
        err1 = product1 - ahi*bhi
        err2 = err1 - alo*bhi
        err3 = err2 - ahi*blo
        hh = alo*blo - err3
        q = product1
        hindex = 0
        if (hh /= 0.0_dp) then
            hindex = hindex + 1
            h(hindex) = hh
        end if
        do eindex = 2, elen
            enow = e(eindex)
            product1 = enow*b
            call split(enow, ahi, alo)
            err1 = product1 - ahi*bhi
            err2 = err1 - alo*bhi
            err3 = err2 - ahi*blo
            product0 = alo*blo - err3
            call two_sum(q, product0, sum_, hh)
            if (hh /= 0.0_dp) then
                hindex = hindex + 1
                h(hindex) = hh
            end if
            call fast_two_sum(product1, sum_, q, hh)
            if (hh /= 0.0_dp) then
                hindex = hindex + 1
                h(hindex) = hh
            end if
        end do
        if ((q /= 0.0_dp) .or. (hindex == 0)) then
            hindex = hindex + 1
            h(hindex) = q
        end if
        hlen = hindex
    end function scale_expansion_zeroelim

end module tc_predicates
