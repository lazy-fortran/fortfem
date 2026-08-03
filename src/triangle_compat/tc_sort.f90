module tc_sort
    !! Vertex sorting for divide-and-conquer triangulation: quicksort with
    !! random pivots, median selection, and alternating cut axes.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tc_state, only: tc_state_t, tc_random
    implicit none
    private

    public :: vertexsort, alternateaxes

contains

    recursive subroutine vertexsort(s, a, n)
        type(tc_state_t), intent(inout) :: s
        integer, intent(inout) :: a(:)
        integer, intent(in) :: n
        integer :: pivot, left, right, temp
        real(dp) :: pivotx, pivoty

        if (n == 2) then
            if ((s%vx(a(1)) > s%vx(a(2))) .or. &
                ((s%vx(a(1)) == s%vx(a(2))) .and. (s%vy(a(1)) > s%vy(a(2))))) then
                temp = a(2)
                a(2) = a(1)
                a(1) = temp
            end if
            return
        end if
        pivot = tc_random(s, n)
        pivotx = s%vx(a(pivot + 1))
        pivoty = s%vy(a(pivot + 1))
        left = -1
        right = n
        do while (left < right)
            do
                left = left + 1
                if (left > right) exit
                if (.not. ((s%vx(a(left + 1)) < pivotx) .or. &
                           ((s%vx(a(left + 1)) == pivotx) .and. &
                            (s%vy(a(left + 1)) < pivoty)))) exit
            end do
            do
                right = right - 1
                if (left > right) exit
                if (.not. ((s%vx(a(right + 1)) > pivotx) .or. &
                           ((s%vx(a(right + 1)) == pivotx) .and. &
                            (s%vy(a(right + 1)) > pivoty)))) exit
            end do
            if (left < right) then
                temp = a(left + 1)
                a(left + 1) = a(right + 1)
                a(right + 1) = temp
            end if
        end do
        if (left > 1) call vertexsort(s, a(1:left), left)
        if (right < n - 2) call vertexsort(s, a(right + 2:n), n - right - 1)
    end subroutine vertexsort

    recursive subroutine vertexmedian(s, a, n, median, axis)
        type(tc_state_t), intent(inout) :: s
        integer, intent(inout) :: a(:)
        integer, intent(in) :: n, median, axis
        integer :: pivot, left, right, temp
        real(dp) :: pivot1, pivot2

        if (n == 2) then
            if ((coord(s, a(1), axis) > coord(s, a(2), axis)) .or. &
                ((coord(s, a(1), axis) == coord(s, a(2), axis)) .and. &
                 (coord(s, a(1), 1 - axis) > coord(s, a(2), 1 - axis)))) then
                temp = a(2)
                a(2) = a(1)
                a(1) = temp
            end if
            return
        end if
        pivot = tc_random(s, n)
        pivot1 = coord(s, a(pivot + 1), axis)
        pivot2 = coord(s, a(pivot + 1), 1 - axis)
        left = -1
        right = n
        do while (left < right)
            do
                left = left + 1
                if (left > right) exit
                if (.not. ((coord(s, a(left + 1), axis) < pivot1) .or. &
                           ((coord(s, a(left + 1), axis) == pivot1) .and. &
                            (coord(s, a(left + 1), 1 - axis) < pivot2)))) exit
            end do
            do
                right = right - 1
                if (left > right) exit
                if (.not. ((coord(s, a(right + 1), axis) > pivot1) .or. &
                           ((coord(s, a(right + 1), axis) == pivot1) .and. &
                            (coord(s, a(right + 1), 1 - axis) > pivot2)))) exit
            end do
            if (left < right) then
                temp = a(left + 1)
                a(left + 1) = a(right + 1)
                a(right + 1) = temp
            end if
        end do
        if (left > median) call vertexmedian(s, a(1:left), left, median, axis)
        if (right < median - 1) then
            call vertexmedian(s, a(right + 2:n), n - right - 1, &
                              median - right - 1, axis)
        end if
    end subroutine vertexmedian

    pure function coord(s, v, axis) result(c)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: v, axis
        real(dp) :: c
        if (axis == 0) then
            c = s%vx(v)
        else
            c = s%vy(v)
        end if
    end function coord

    recursive subroutine alternateaxes(s, a, n, axis_in)
        type(tc_state_t), intent(inout) :: s
        integer, intent(inout) :: a(:)
        integer, intent(in) :: n, axis_in
        integer :: divider, axis

        axis = axis_in
        divider = shiftr(n, 1)
        if (n <= 3) axis = 0
        call vertexmedian(s, a, n, divider, axis)
        if (n - divider >= 2) then
            if (divider >= 2) then
                call alternateaxes(s, a(1:divider), divider, 1 - axis)
            end if
            call alternateaxes(s, a(divider + 1:n), n - divider, 1 - axis)
        end if
    end subroutine alternateaxes

end module tc_sort
