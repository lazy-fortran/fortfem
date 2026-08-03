module tc_flip
    !! Edge flip and its exact inverse on the oriented-triangle structure.
    use tc_state
    implicit none
    private

    public :: flip, unflip

contains

    subroutine flip(s, flipedge)
        type(tc_state_t), intent(inout) :: s
        integer, intent(inout) :: flipedge
        integer :: botleft, botright, topleft, topright, top
        integer :: botlcasing, botrcasing, toplcasing, toprcasing
        integer :: botlsubseg, botrsubseg, toplsubseg, toprsubseg
        integer :: leftvertex, rightvertex, botvertex, farvertex

        rightvertex = t_org(s, flipedge)
        leftvertex = t_dest(s, flipedge)
        botvertex = t_apex(s, flipedge)
        top = t_sym(s, flipedge)
        farvertex = t_apex(s, top)

        topleft = t_lprev(top)
        toplcasing = t_sym(s, topleft)
        topright = t_lnext(top)
        toprcasing = t_sym(s, topright)
        botleft = t_lnext(flipedge)
        botlcasing = t_sym(s, botleft)
        botright = t_lprev(flipedge)
        botrcasing = t_sym(s, botright)
        call t_bond(s, topleft, botlcasing)
        call t_bond(s, botleft, botrcasing)
        call t_bond(s, botright, toprcasing)
        call t_bond(s, topright, toplcasing)

        if (s%checksegments) then
            toplsubseg = ts_pivot(s, topleft)
            botlsubseg = ts_pivot(s, botleft)
            botrsubseg = ts_pivot(s, botright)
            toprsubseg = ts_pivot(s, topright)
            if (s_slot(toplsubseg) == 0) then
                call ts_dissolve(s, topright)
            else
                call ts_bond(s, topright, toplsubseg)
            end if
            if (s_slot(botlsubseg) == 0) then
                call ts_dissolve(s, topleft)
            else
                call ts_bond(s, topleft, botlsubseg)
            end if
            if (s_slot(botrsubseg) == 0) then
                call ts_dissolve(s, botleft)
            else
                call ts_bond(s, botleft, botrsubseg)
            end if
            if (s_slot(toprsubseg) == 0) then
                call ts_dissolve(s, botright)
            else
                call ts_bond(s, botright, toprsubseg)
            end if
        end if

        call set_org(s, flipedge, farvertex)
        call set_dest(s, flipedge, botvertex)
        call set_apex(s, flipedge, rightvertex)
        call set_org(s, top, botvertex)
        call set_dest(s, top, farvertex)
        call set_apex(s, top, leftvertex)
    end subroutine flip

    subroutine unflip(s, flipedge)
        type(tc_state_t), intent(inout) :: s
        integer, intent(inout) :: flipedge
        integer :: botleft, botright, topleft, topright, top
        integer :: botlcasing, botrcasing, toplcasing, toprcasing
        integer :: botlsubseg, botrsubseg, toplsubseg, toprsubseg
        integer :: leftvertex, rightvertex, botvertex, farvertex

        rightvertex = t_org(s, flipedge)
        leftvertex = t_dest(s, flipedge)
        botvertex = t_apex(s, flipedge)
        top = t_sym(s, flipedge)
        farvertex = t_apex(s, top)

        topleft = t_lprev(top)
        toplcasing = t_sym(s, topleft)
        topright = t_lnext(top)
        toprcasing = t_sym(s, topright)
        botleft = t_lnext(flipedge)
        botlcasing = t_sym(s, botleft)
        botright = t_lprev(flipedge)
        botrcasing = t_sym(s, botright)
        call t_bond(s, topleft, toprcasing)
        call t_bond(s, botleft, toplcasing)
        call t_bond(s, botright, botlcasing)
        call t_bond(s, topright, botrcasing)

        if (s%checksegments) then
            toplsubseg = ts_pivot(s, topleft)
            botlsubseg = ts_pivot(s, botleft)
            botrsubseg = ts_pivot(s, botright)
            toprsubseg = ts_pivot(s, topright)
            if (s_slot(toplsubseg) == 0) then
                call ts_dissolve(s, botleft)
            else
                call ts_bond(s, botleft, toplsubseg)
            end if
            if (s_slot(botlsubseg) == 0) then
                call ts_dissolve(s, botright)
            else
                call ts_bond(s, botright, botlsubseg)
            end if
            if (s_slot(botrsubseg) == 0) then
                call ts_dissolve(s, topright)
            else
                call ts_bond(s, topright, botrsubseg)
            end if
            if (s_slot(toprsubseg) == 0) then
                call ts_dissolve(s, topleft)
            else
                call ts_bond(s, topleft, toprsubseg)
            end if
        end if

        call set_org(s, flipedge, botvertex)
        call set_dest(s, flipedge, farvertex)
        call set_apex(s, flipedge, leftvertex)
        call set_org(s, top, farvertex)
        call set_dest(s, top, botvertex)
        call set_apex(s, top, rightvertex)
    end subroutine unflip

end module tc_flip
