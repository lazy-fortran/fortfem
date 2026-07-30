module fortfem_bspline_multipatch
    !! Orientation-aware conforming traces for tensor-product spline patches.
    implicit none
    private

    integer, parameter, public :: BSPLINE_FACE_X_MIN = 1
    integer, parameter, public :: BSPLINE_FACE_X_MAX = 2
    integer, parameter, public :: BSPLINE_FACE_Y_MIN = 3
    integer, parameter, public :: BSPLINE_FACE_Y_MAX = 4

    public :: build_bspline_feec_2d_interface_dofs

contains

    subroutine build_bspline_feec_2d_interface_dofs( &
            nx_left, ny_left, nx_right, ny_right, face_left, face_right, &
            reversed, h1_left, h1_right, hcurl_left, hcurl_right, hcurl_sign, &
            status)
        integer, intent(in) :: nx_left, ny_left, nx_right, ny_right
        integer, intent(in) :: face_left, face_right
        logical, intent(in) :: reversed
        integer, allocatable, intent(out) :: h1_left(:), h1_right(:)
        integer, allocatable, intent(out) :: hcurl_left(:), hcurl_right(:)
        integer, allocatable, intent(out) :: hcurl_sign(:)
        integer, intent(out) :: status

        integer :: trace_count

        status = 1
        if (nx_left < 2 .or. ny_left < 2) return
        if (nx_right < 2 .or. ny_right < 2) return
        trace_count = face_trace_count(nx_left, ny_left, face_left)
        if (trace_count < 2) return
        if (face_trace_count(nx_right, ny_right, face_right) /= &
            trace_count) return
        allocate( &
            h1_left(trace_count), h1_right(trace_count), &
            hcurl_left(trace_count - 1), hcurl_right(trace_count - 1), &
            hcurl_sign(trace_count - 1))
        call build_h1_face_dofs(nx_left, ny_left, face_left, h1_left)
        call build_h1_face_dofs(nx_right, ny_right, face_right, h1_right)
        call build_hcurl_face_dofs( &
            nx_left, ny_left, face_left, hcurl_left)
        call build_hcurl_face_dofs( &
            nx_right, ny_right, face_right, hcurl_right)
        hcurl_sign = 1
        if (reversed) then
            h1_right = h1_right(size(h1_right):1:-1)
            hcurl_right = hcurl_right(size(hcurl_right):1:-1)
            hcurl_sign = -1
        end if
        status = 0
    end subroutine build_bspline_feec_2d_interface_dofs

    pure integer function face_trace_count(nx, ny, face) result(count)
        integer, intent(in) :: nx, ny, face

        select case (face)
        case (BSPLINE_FACE_X_MIN, BSPLINE_FACE_X_MAX)
            count = ny
        case (BSPLINE_FACE_Y_MIN, BSPLINE_FACE_Y_MAX)
            count = nx
        case default
            count = 0
        end select
    end function face_trace_count

    pure subroutine build_h1_face_dofs(nx, ny, face, dofs)
        integer, intent(in) :: nx, ny, face
        integer, intent(out) :: dofs(:)
        integer :: along, fixed

        select case (face)
        case (BSPLINE_FACE_X_MIN, BSPLINE_FACE_X_MAX)
            fixed = 1
            if (face == BSPLINE_FACE_X_MAX) fixed = nx
            do along = 1, ny
                dofs(along) = fixed + (along - 1)*nx
            end do
        case (BSPLINE_FACE_Y_MIN, BSPLINE_FACE_Y_MAX)
            fixed = 1
            if (face == BSPLINE_FACE_Y_MAX) fixed = ny
            do along = 1, nx
                dofs(along) = along + (fixed - 1)*nx
            end do
        end select
    end subroutine build_h1_face_dofs

    pure subroutine build_hcurl_face_dofs(nx, ny, face, dofs)
        integer, intent(in) :: nx, ny, face
        integer, intent(out) :: dofs(:)
        integer :: along, fixed, x_component_count

        x_component_count = (nx - 1)*ny
        select case (face)
        case (BSPLINE_FACE_X_MIN, BSPLINE_FACE_X_MAX)
            fixed = 1
            if (face == BSPLINE_FACE_X_MAX) fixed = nx
            do along = 1, ny - 1
                dofs(along) = &
                    x_component_count + fixed + (along - 1)*nx
            end do
        case (BSPLINE_FACE_Y_MIN, BSPLINE_FACE_Y_MAX)
            fixed = 1
            if (face == BSPLINE_FACE_Y_MAX) fixed = ny
            do along = 1, nx - 1
                dofs(along) = along + (fixed - 1)*(nx - 1)
            end do
        end select
    end subroutine build_hcurl_face_dofs

end module fortfem_bspline_multipatch
