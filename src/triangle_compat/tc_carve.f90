module tc_carve
    !! Hole and concavity carving by virus infection, in the exact
    !! traversal and deallocation order of the reference implementation.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tc_state
    use tc_alloc
    use tc_geom, only: ccw_p
    use tc_locate, only: locate, LOC_OUTSIDE
    implicit none
    private

    public :: carveholes

contains

    subroutine virus_push(s, slot)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: slot
        integer, allocatable :: i1(:)
        integer :: oldcap

        oldcap = size(s%vir)
        if (s%vir_n + 1 > oldcap) then
            allocate (i1(2*oldcap))
            i1(1:oldcap) = s%vir
            call move_alloc(i1, s%vir)
        end if
        s%vir_n = s%vir_n + 1
        s%vir(s%vir_n) = slot
    end subroutine virus_push

    subroutine carveholes(s, holes)
        type(tc_state_t), intent(inout) :: s
        real(dp), intent(in) :: holes(:, :)
        integer :: i, searchtri, searchorg, searchdest, intersect

        call infecthull(s)

        do i = 1, size(holes, 2)
            if ((holes(1, i) >= s%xmin) .and. (holes(1, i) <= s%xmax) &
                .and. (holes(2, i) >= s%ymin) .and. (holes(2, i) <= s%ymax)) then
                searchtri = s%tnbr(0, 0)
                searchorg = t_org(s, searchtri)
                searchdest = t_dest(s, searchtri)
                if (ccw_p(s, searchorg, searchdest, holes(1, i), holes(2, i)) &
                    > 0.0_dp) then
                    intersect = locate(s, holes(1, i), holes(2, i), searchtri)
                    if ((intersect /= LOC_OUTSIDE) .and. &
                        (.not. s%tinf(t_slot(searchtri)))) then
                        s%tinf(t_slot(searchtri)) = .true.
                        call virus_push(s, t_slot(searchtri))
                    end if
                end if
            end if
        end do

        if (s%vir_n > 0) call plague(s)
    end subroutine carveholes

    subroutine infecthull(s)
        type(tc_state_t), intent(inout) :: s
        integer :: hulltri, nexttri, starttri, hullsubseg
        integer :: horg, hdest

        hulltri = s%tnbr(0, 0)
        starttri = hulltri
        do
            if (.not. s%tinf(t_slot(hulltri))) then
                hullsubseg = ts_pivot(s, hulltri)
                if (s_slot(hullsubseg) == 0) then
                    if (.not. s%tinf(t_slot(hulltri))) then
                        s%tinf(t_slot(hulltri)) = .true.
                        call virus_push(s, t_slot(hulltri))
                    end if
                else
                    if (s%smark(s_slot(hullsubseg)) == 0) then
                        s%smark(s_slot(hullsubseg)) = 1
                        horg = t_org(s, hulltri)
                        hdest = t_dest(s, hulltri)
                        if (s%vmark(horg) == 0) s%vmark(horg) = 1
                        if (s%vmark(hdest) == 0) s%vmark(hdest) = 1
                    end if
                end if
            end if
            hulltri = t_lnext(hulltri)
            nexttri = t_oprev(s, hulltri)
            do while (t_slot(nexttri) /= 0)
                hulltri = nexttri
                nexttri = t_oprev(s, hulltri)
            end do
            if (hulltri == starttri) exit
        end do
    end subroutine infecthull

    subroutine plague(s)
        type(tc_state_t), intent(inout) :: s
        integer :: i, orient, testtri, neighbor, neighborsubseg
        integer :: testvertex, norg, ndest
        logical :: killorg

        ! Spread the virus.
        i = 1
        do while (i <= s%vir_n)
            testtri = t_make(s%vir(i), 0)
            s%tinf(s%vir(i)) = .false.
            do orient = 0, 2
                testtri = t_make(s%vir(i), orient)
                neighbor = t_sym(s, testtri)
                neighborsubseg = ts_pivot(s, testtri)
                if ((t_slot(neighbor) == 0) .or. &
                    s%tinf(t_slot(neighbor))) then
                    if (s_slot(neighborsubseg) /= 0) then
                        call subseg_dealloc(s, s_slot(neighborsubseg))
                        if (t_slot(neighbor) /= 0) then
                            s%tinf(t_slot(neighbor)) = .false.
                            call ts_dissolve(s, neighbor)
                            s%tinf(t_slot(neighbor)) = .true.
                        end if
                    end if
                else
                    if (s_slot(neighborsubseg) == 0) then
                        s%tinf(t_slot(neighbor)) = .true.
                        call virus_push(s, t_slot(neighbor))
                    else
                        call st_dissolve(s, neighborsubseg)
                        if (s%smark(s_slot(neighborsubseg)) == 0) then
                            s%smark(s_slot(neighborsubseg)) = 1
                        end if
                        norg = t_org(s, neighbor)
                        ndest = t_dest(s, neighbor)
                        if (s%vmark(norg) == 0) s%vmark(norg) = 1
                        if (s%vmark(ndest) == 0) s%vmark(ndest) = 1
                    end if
                end if
            end do
            s%tinf(s%vir(i)) = .true.
            i = i + 1
        end do

        ! Kill the infected triangles and orphaned vertices.
        do i = 1, s%vir_n
            do orient = 0, 2
                testtri = t_make(s%vir(i), orient)
                testvertex = t_org(s, testtri)
                if (testvertex /= 0) then
                    killorg = .true.
                    call set_org(s, testtri, 0)
                    neighbor = t_onext(s, testtri)
                    do while ((t_slot(neighbor) /= 0) .and. &
                              (neighbor /= testtri))
                        if (s%tinf(t_slot(neighbor))) then
                            call set_org(s, neighbor, 0)
                        else
                            killorg = .false.
                        end if
                        neighbor = t_onext(s, neighbor)
                    end do
                    if (t_slot(neighbor) == 0) then
                        neighbor = t_oprev(s, testtri)
                        do while (t_slot(neighbor) /= 0)
                            if (s%tinf(t_slot(neighbor))) then
                                call set_org(s, neighbor, 0)
                            else
                                killorg = .false.
                            end if
                            neighbor = t_oprev(s, neighbor)
                        end do
                    end if
                    if (killorg) then
                        s%vtype(testvertex) = V_UNDEAD
                        s%undeads = s%undeads + 1
                    end if
                end if
            end do
            do orient = 0, 2
                testtri = t_make(s%vir(i), orient)
                neighbor = t_sym(s, testtri)
                if (t_slot(neighbor) == 0) then
                    s%hullsize = s%hullsize - 1
                else
                    call t_dissolve(s, neighbor)
                    s%hullsize = s%hullsize + 1
                end if
            end do
            call triangle_dealloc(s, s%vir(i))
        end do
        s%vir_n = 0
    end subroutine plague

end module tc_carve
