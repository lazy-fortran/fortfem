module tc_divconq
    !! Divide-and-conquer Delaunay triangulation with alternating cuts,
    !! following the algorithm of Guibas and Stolfi as engineered in
    !! Shewchuk's mesh generator (ghost-triangle bounding fan, identical
    !! creation and removal order so that pool slot numbering matches).
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use tc_state
    use tc_alloc
    use tc_sort, only: vertexsort, alternateaxes
    use tc_geom, only: ccw_v, incircle_v
    implicit none
    private

    public :: divconq_delaunay

contains

    function divconq_delaunay(s) result(hullsize)
        type(tc_state_t), intent(inout) :: s
        integer(int64) :: hullsize
        integer, allocatable :: sortarray(:)
        integer :: i, j, divider, hullleft, hullright

        allocate (sortarray(s%invertices))
        do i = 1, s%invertices
            sortarray(i) = i
        end do
        call vertexsort(s, sortarray, s%invertices)
        i = 1
        do j = 2, s%invertices
            if ((s%vx(sortarray(i)) == s%vx(sortarray(j))) .and. &
                (s%vy(sortarray(i)) == s%vy(sortarray(j)))) then
                s%vtype(sortarray(j)) = V_UNDEAD
                s%undeads = s%undeads + 1
            else
                i = i + 1
                sortarray(i) = sortarray(j)
            end if
        end do
        divider = shiftr(i, 1)
        if (i - divider >= 2) then
            if (divider >= 2) then
                call alternateaxes(s, sortarray(1:divider), divider, 1)
            end if
            call alternateaxes(s, sortarray(divider + 1:i), i - divider, 1)
        end if
        call divconqrecurse(s, sortarray(1:i), i, 0, hullleft, hullright)
        hullsize = removeghosts(s, hullleft)
    end function divconq_delaunay

    recursive subroutine divconqrecurse(s, a, n, axis, farleft, farright)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: a(:)
        integer, intent(in) :: n, axis
        integer, intent(out) :: farleft, farright
        integer :: midtri, tri1, tri2, tri3
        integer :: innerleft, innerright
        real(dp) :: area
        integer :: divider

        if (n == 2) then
            farleft = make_triangle(s)
            call set_org(s, farleft, a(1))
            call set_dest(s, farleft, a(2))
            farright = make_triangle(s)
            call set_org(s, farright, a(2))
            call set_dest(s, farright, a(1))
            call t_bond(s, farleft, farright)
            farleft = t_lprev(farleft)
            farright = t_lnext(farright)
            call t_bond(s, farleft, farright)
            farleft = t_lprev(farleft)
            farright = t_lnext(farright)
            call t_bond(s, farleft, farright)
            farleft = t_lprev(farright)
            return
        else if (n == 3) then
            midtri = make_triangle(s)
            tri1 = make_triangle(s)
            tri2 = make_triangle(s)
            tri3 = make_triangle(s)
            area = ccw_v(s, a(1), a(2), a(3))
            if (area == 0.0_dp) then
                call set_org(s, midtri, a(1))
                call set_dest(s, midtri, a(2))
                call set_org(s, tri1, a(2))
                call set_dest(s, tri1, a(1))
                call set_org(s, tri2, a(3))
                call set_dest(s, tri2, a(2))
                call set_org(s, tri3, a(2))
                call set_dest(s, tri3, a(3))
                call t_bond(s, midtri, tri1)
                call t_bond(s, tri2, tri3)
                midtri = t_lnext(midtri)
                tri1 = t_lprev(tri1)
                tri2 = t_lnext(tri2)
                tri3 = t_lprev(tri3)
                call t_bond(s, midtri, tri3)
                call t_bond(s, tri1, tri2)
                midtri = t_lnext(midtri)
                tri1 = t_lprev(tri1)
                tri2 = t_lnext(tri2)
                tri3 = t_lprev(tri3)
                call t_bond(s, midtri, tri1)
                call t_bond(s, tri2, tri3)
                farleft = tri1
                farright = tri2
            else
                call set_org(s, midtri, a(1))
                call set_dest(s, tri1, a(1))
                call set_org(s, tri3, a(1))
                if (area > 0.0_dp) then
                    call set_dest(s, midtri, a(2))
                    call set_org(s, tri1, a(2))
                    call set_dest(s, tri2, a(2))
                    call set_apex(s, midtri, a(3))
                    call set_org(s, tri2, a(3))
                    call set_dest(s, tri3, a(3))
                else
                    call set_dest(s, midtri, a(3))
                    call set_org(s, tri1, a(3))
                    call set_dest(s, tri2, a(3))
                    call set_apex(s, midtri, a(2))
                    call set_org(s, tri2, a(2))
                    call set_dest(s, tri3, a(2))
                end if
                call t_bond(s, midtri, tri1)
                midtri = t_lnext(midtri)
                call t_bond(s, midtri, tri2)
                midtri = t_lnext(midtri)
                call t_bond(s, midtri, tri3)
                tri1 = t_lprev(tri1)
                tri2 = t_lnext(tri2)
                call t_bond(s, tri1, tri2)
                tri1 = t_lprev(tri1)
                tri3 = t_lprev(tri3)
                call t_bond(s, tri1, tri3)
                tri2 = t_lnext(tri2)
                tri3 = t_lprev(tri3)
                call t_bond(s, tri2, tri3)
                farleft = tri1
                if (area > 0.0_dp) then
                    farright = tri2
                else
                    farright = t_lnext(farleft)
                end if
            end if
            return
        else
            divider = shiftr(n, 1)
            call divconqrecurse(s, a(1:divider), divider, 1 - axis, &
                                farleft, innerleft)
            call divconqrecurse(s, a(divider + 1:n), n - divider, 1 - axis, &
                                innerright, farright)
            call mergehulls(s, farleft, innerleft, innerright, farright, axis)
        end if
    end subroutine divconqrecurse

    function removeghosts(s, startghost) result(hullsize)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: startghost
        integer(int64) :: hullsize
        integer :: searchedge, dissolveedge, deadtriangle

        searchedge = t_sym(s, t_lprev(startghost))
        s%tnbr(0, 0) = searchedge
        dissolveedge = startghost
        hullsize = 0
        do
            hullsize = hullsize + 1
            deadtriangle = t_lnext(dissolveedge)
            dissolveedge = t_sym(s, t_lprev(dissolveedge))
            call t_dissolve(s, dissolveedge)
            dissolveedge = t_sym(s, deadtriangle)
            call triangle_dealloc(s, t_slot(deadtriangle))
            if (dissolveedge == startghost) exit
        end do
    end function removeghosts

    subroutine mergehulls(s, farleft, innerleft, innerright, farright, axis)
        type(tc_state_t), intent(inout) :: s
        integer, intent(inout) :: farleft, innerleft, innerright, farright
        integer, intent(in) :: axis
        integer :: leftcand, rightcand, baseedge, nextedge
        integer :: sidecasing, topcasing, outercasing, checkedge
        integer :: innerleftdest, innerrightorg, innerleftapex, innerrightapex
        integer :: farleftpt, farrightpt, farleftapex, farrightapex
        integer :: lowerleft, lowerright, upperleft, upperright
        integer :: nextapex, checkvertex
        logical :: changemade, badedge, leftfinished, rightfinished

        innerleftdest = t_dest(s, innerleft)
        innerleftapex = t_apex(s, innerleft)
        innerrightorg = t_org(s, innerright)
        innerrightapex = t_apex(s, innerright)
        if (axis == 1) then
            farleftpt = t_org(s, farleft)
            farleftapex = t_apex(s, farleft)
            farrightpt = t_dest(s, farright)
            farrightapex = t_apex(s, farright)
            do while (s%vy(farleftapex) < s%vy(farleftpt))
                farleft = t_sym(s, t_lnext(farleft))
                farleftpt = farleftapex
                farleftapex = t_apex(s, farleft)
            end do
            checkedge = t_sym(s, innerleft)
            checkvertex = t_apex(s, checkedge)
            do while (s%vy(checkvertex) > s%vy(innerleftdest))
                innerleft = t_lnext(checkedge)
                innerleftapex = innerleftdest
                innerleftdest = checkvertex
                checkedge = t_sym(s, innerleft)
                checkvertex = t_apex(s, checkedge)
            end do
            do while (s%vy(innerrightapex) < s%vy(innerrightorg))
                innerright = t_sym(s, t_lnext(innerright))
                innerrightorg = innerrightapex
                innerrightapex = t_apex(s, innerright)
            end do
            checkedge = t_sym(s, farright)
            checkvertex = t_apex(s, checkedge)
            do while (s%vy(checkvertex) > s%vy(farrightpt))
                farright = t_lnext(checkedge)
                farrightapex = farrightpt
                farrightpt = checkvertex
                checkedge = t_sym(s, farright)
                checkvertex = t_apex(s, checkedge)
            end do
        end if
        do
            changemade = .false.
            if (ccw_v(s, innerleftdest, innerleftapex, innerrightorg) &
                > 0.0_dp) then
                innerleft = t_sym(s, t_lprev(innerleft))
                innerleftdest = innerleftapex
                innerleftapex = t_apex(s, innerleft)
                changemade = .true.
            end if
            if (ccw_v(s, innerrightapex, innerrightorg, innerleftdest) &
                > 0.0_dp) then
                innerright = t_sym(s, t_lnext(innerright))
                innerrightorg = innerrightapex
                innerrightapex = t_apex(s, innerright)
                changemade = .true.
            end if
            if (.not. changemade) exit
        end do
        leftcand = t_sym(s, innerleft)
        rightcand = t_sym(s, innerright)
        baseedge = make_triangle(s)
        call t_bond(s, baseedge, innerleft)
        baseedge = t_lnext(baseedge)
        call t_bond(s, baseedge, innerright)
        baseedge = t_lnext(baseedge)
        call set_org(s, baseedge, innerrightorg)
        call set_dest(s, baseedge, innerleftdest)
        farleftpt = t_org(s, farleft)
        if (innerleftdest == farleftpt) farleft = t_lnext(baseedge)
        farrightpt = t_dest(s, farright)
        if (innerrightorg == farrightpt) farright = t_lprev(baseedge)
        lowerleft = innerleftdest
        lowerright = innerrightorg
        upperleft = t_apex(s, leftcand)
        upperright = t_apex(s, rightcand)
        do
            leftfinished = ccw_v(s, upperleft, lowerleft, lowerright) <= 0.0_dp
            rightfinished = ccw_v(s, upperright, lowerleft, lowerright) <= 0.0_dp
            if (leftfinished .and. rightfinished) then
                nextedge = make_triangle(s)
                call set_org(s, nextedge, lowerleft)
                call set_dest(s, nextedge, lowerright)
                call t_bond(s, nextedge, baseedge)
                nextedge = t_lnext(nextedge)
                call t_bond(s, nextedge, rightcand)
                nextedge = t_lnext(nextedge)
                call t_bond(s, nextedge, leftcand)
                if (axis == 1) then
                    farleftpt = t_org(s, farleft)
                    farleftapex = t_apex(s, farleft)
                    farrightpt = t_dest(s, farright)
                    farrightapex = t_apex(s, farright)
                    checkedge = t_sym(s, farleft)
                    checkvertex = t_apex(s, checkedge)
                    do while (s%vx(checkvertex) < s%vx(farleftpt))
                        farleft = t_lprev(checkedge)
                        farleftapex = farleftpt
                        farleftpt = checkvertex
                        checkedge = t_sym(s, farleft)
                        checkvertex = t_apex(s, checkedge)
                    end do
                    do while (s%vx(farrightapex) > s%vx(farrightpt))
                        farright = t_sym(s, t_lprev(farright))
                        farrightpt = farrightapex
                        farrightapex = t_apex(s, farright)
                    end do
                end if
                return
            end if
            if (.not. leftfinished) then
                nextedge = t_sym(s, t_lprev(leftcand))
                nextapex = t_apex(s, nextedge)
                if (nextapex /= 0) then
                    badedge = incircle_v(s, lowerleft, lowerright, upperleft, &
                                         nextapex) > 0.0_dp
                    do while (badedge)
                        nextedge = t_lnext(nextedge)
                        topcasing = t_sym(s, nextedge)
                        nextedge = t_lnext(nextedge)
                        sidecasing = t_sym(s, nextedge)
                        call t_bond(s, nextedge, topcasing)
                        call t_bond(s, leftcand, sidecasing)
                        leftcand = t_lnext(leftcand)
                        outercasing = t_sym(s, leftcand)
                        nextedge = t_lprev(nextedge)
                        call t_bond(s, nextedge, outercasing)
                        call set_org(s, leftcand, lowerleft)
                        call set_dest(s, leftcand, 0)
                        call set_apex(s, leftcand, nextapex)
                        call set_org(s, nextedge, 0)
                        call set_dest(s, nextedge, upperleft)
                        call set_apex(s, nextedge, nextapex)
                        upperleft = nextapex
                        nextedge = sidecasing
                        nextapex = t_apex(s, nextedge)
                        if (nextapex /= 0) then
                            badedge = incircle_v(s, lowerleft, lowerright, &
                                                 upperleft, nextapex) > 0.0_dp
                        else
                            badedge = .false.
                        end if
                    end do
                end if
            end if
            if (.not. rightfinished) then
                nextedge = t_sym(s, t_lnext(rightcand))
                nextapex = t_apex(s, nextedge)
                if (nextapex /= 0) then
                    badedge = incircle_v(s, lowerleft, lowerright, upperright, &
                                         nextapex) > 0.0_dp
                    do while (badedge)
                        nextedge = t_lprev(nextedge)
                        topcasing = t_sym(s, nextedge)
                        nextedge = t_lprev(nextedge)
                        sidecasing = t_sym(s, nextedge)
                        call t_bond(s, nextedge, topcasing)
                        call t_bond(s, rightcand, sidecasing)
                        rightcand = t_lprev(rightcand)
                        outercasing = t_sym(s, rightcand)
                        nextedge = t_lnext(nextedge)
                        call t_bond(s, nextedge, outercasing)
                        call set_org(s, rightcand, 0)
                        call set_dest(s, rightcand, lowerright)
                        call set_apex(s, rightcand, nextapex)
                        call set_org(s, nextedge, upperright)
                        call set_dest(s, nextedge, 0)
                        call set_apex(s, nextedge, nextapex)
                        upperright = nextapex
                        nextedge = sidecasing
                        nextapex = t_apex(s, nextedge)
                        if (nextapex /= 0) then
                            badedge = incircle_v(s, lowerleft, lowerright, &
                                                 upperright, nextapex) > 0.0_dp
                        else
                            badedge = .false.
                        end if
                    end do
                end if
            end if
            if (leftfinished .or. (.not. rightfinished .and. &
                                   (incircle_v(s, upperleft, lowerleft, &
                                               lowerright, upperright) &
                                    > 0.0_dp))) then
                call t_bond(s, baseedge, rightcand)
                baseedge = t_lprev(rightcand)
                call set_dest(s, baseedge, lowerleft)
                lowerright = upperright
                rightcand = t_sym(s, baseedge)
                upperright = t_apex(s, rightcand)
            else
                call t_bond(s, baseedge, leftcand)
                baseedge = t_lnext(leftcand)
                call set_org(s, baseedge, lowerright)
                lowerleft = upperleft
                leftcand = t_sym(s, baseedge)
                upperleft = t_apex(s, leftcand)
            end if
        end do
    end subroutine mergehulls

end module tc_divconq
