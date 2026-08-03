module tc_geom
    !! Vertex-id wrappers around the exact predicates.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tc_state, only: tc_state_t
    use tc_predicates, only: tc_counterclockwise, tc_incircle
    implicit none
    private

    public :: ccw_v, incircle_v, ccw_p

    interface ccw_v
        module procedure ccw_vvv
    end interface ccw_v

contains

    function ccw_vvv(s, va, vb, vc) result(det)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: va, vb, vc
        real(dp) :: det
        det = tc_counterclockwise(s%vx(va), s%vy(va), s%vx(vb), s%vy(vb), &
                                  s%vx(vc), s%vy(vc))
    end function ccw_vvv

    function ccw_p(s, va, vb, px, py) result(det)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: va, vb
        real(dp), intent(in) :: px, py
        real(dp) :: det
        det = tc_counterclockwise(s%vx(va), s%vy(va), s%vx(vb), s%vy(vb), &
                                  px, py)
    end function ccw_p

    function incircle_v(s, va, vb, vc, vd) result(det)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: va, vb, vc, vd
        real(dp) :: det
        det = tc_incircle(s%vx(va), s%vy(va), s%vx(vb), s%vy(vb), &
                          s%vx(vc), s%vy(vc), s%vx(vd), s%vy(vd))
    end function incircle_v

end module tc_geom
