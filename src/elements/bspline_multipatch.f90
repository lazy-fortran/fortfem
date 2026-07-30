module fortfem_bspline_multipatch
    !! Orientation-aware conforming traces for tensor-product spline patches.
    implicit none
    private

    integer, parameter, public :: BSPLINE_FACE_X_MIN = 1
    integer, parameter, public :: BSPLINE_FACE_X_MAX = 2
    integer, parameter, public :: BSPLINE_FACE_Y_MIN = 3
    integer, parameter, public :: BSPLINE_FACE_Y_MAX = 4
    integer, parameter, public :: BSPLINE_FACE_Z_MIN = 5
    integer, parameter, public :: BSPLINE_FACE_Z_MAX = 6

    public :: build_bspline_feec_2d_interface_dofs
    public :: build_bspline_feec_2d_two_patch_maps
    public :: build_bspline_feec_3d_interface_dofs
    public :: build_bspline_feec_3d_two_patch_maps

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

    subroutine build_bspline_feec_2d_two_patch_maps( &
            nx_left, ny_left, nx_right, ny_right, face_left, face_right, &
            reversed, h1_left_map, h1_right_map, hcurl_left_map, &
            hcurl_right_map, hcurl_left_sign, hcurl_right_sign, l2_left_map, &
            l2_right_map, status)
        integer, intent(in) :: nx_left, ny_left, nx_right, ny_right
        integer, intent(in) :: face_left, face_right
        logical, intent(in) :: reversed
        integer, allocatable, intent(out) :: h1_left_map(:), h1_right_map(:)
        integer, allocatable, intent(out) :: hcurl_left_map(:)
        integer, allocatable, intent(out) :: hcurl_right_map(:)
        integer, allocatable, intent(out) :: hcurl_left_sign(:)
        integer, allocatable, intent(out) :: hcurl_right_sign(:)
        integer, allocatable, intent(out) :: l2_left_map(:), l2_right_map(:)
        integer, intent(out) :: status

        integer, allocatable :: h1_left_trace(:), h1_right_trace(:)
        integer, allocatable :: hcurl_left_trace(:), hcurl_right_trace(:)
        integer, allocatable :: trace_sign(:)
        integer :: dof, global_dof, trace

        call build_bspline_feec_2d_interface_dofs( &
            nx_left, ny_left, nx_right, ny_right, face_left, face_right, &
            reversed, h1_left_trace, h1_right_trace, hcurl_left_trace, &
            hcurl_right_trace, trace_sign, status)
        if (status /= 0) return
        allocate( &
            h1_left_map(nx_left*ny_left), h1_right_map(nx_right*ny_right), &
            hcurl_left_map((nx_left - 1)*ny_left + nx_left*(ny_left - 1)), &
            hcurl_right_map( &
            (nx_right - 1)*ny_right + nx_right*(ny_right - 1)), &
            l2_left_map((nx_left - 1)*(ny_left - 1)), &
            l2_right_map((nx_right - 1)*(ny_right - 1)))
        allocate( &
            hcurl_left_sign(size(hcurl_left_map)), &
            hcurl_right_sign(size(hcurl_right_map)))
        h1_left_map = [(dof, dof = 1, size(h1_left_map))]
        h1_right_map = 0
        do trace = 1, size(h1_left_trace)
            h1_right_map(h1_right_trace(trace)) = &
                h1_left_map(h1_left_trace(trace))
        end do
        global_dof = size(h1_left_map)
        do dof = 1, size(h1_right_map)
            if (h1_right_map(dof) /= 0) cycle
            global_dof = global_dof + 1
            h1_right_map(dof) = global_dof
        end do

        hcurl_left_map = [(dof, dof = 1, size(hcurl_left_map))]
        hcurl_left_sign = 1
        hcurl_right_map = 0
        hcurl_right_sign = 1
        do trace = 1, size(hcurl_left_trace)
            hcurl_right_map(hcurl_right_trace(trace)) = &
                hcurl_left_map(hcurl_left_trace(trace))
            hcurl_right_sign(hcurl_right_trace(trace)) = trace_sign(trace)
        end do
        global_dof = size(hcurl_left_map)
        do dof = 1, size(hcurl_right_map)
            if (hcurl_right_map(dof) /= 0) cycle
            global_dof = global_dof + 1
            hcurl_right_map(dof) = global_dof
        end do

        l2_left_map = [(dof, dof = 1, size(l2_left_map))]
        l2_right_map = &
            [(size(l2_left_map) + dof, dof = 1, size(l2_right_map))]
    end subroutine build_bspline_feec_2d_two_patch_maps

    subroutine build_bspline_feec_3d_interface_dofs( &
            nx_left, ny_left, nz_left, nx_right, ny_right, nz_right, &
            face_left, face_right, swap_axes, reverse_u, reverse_v, h1_left, &
            h1_right, hcurl_left, hcurl_right, hcurl_sign, hdiv_left, &
            hdiv_right, hdiv_sign, status)
        integer, intent(in) :: nx_left, ny_left, nz_left
        integer, intent(in) :: nx_right, ny_right, nz_right
        integer, intent(in) :: face_left, face_right
        logical, intent(in) :: swap_axes, reverse_u, reverse_v
        integer, allocatable, intent(out) :: h1_left(:), h1_right(:)
        integer, allocatable, intent(out) :: hcurl_left(:), hcurl_right(:)
        integer, allocatable, intent(out) :: hcurl_sign(:)
        integer, allocatable, intent(out) :: hdiv_left(:), hdiv_right(:)
        integer, allocatable, intent(out) :: hdiv_sign(:)
        integer, intent(out) :: status

        integer :: edge, iu, iv, ju1, ju2, jv1, jv2
        integer :: nu_left, nu_right, nv_left, nv_right
        integer :: node1(3), node2(3), right_node1(3), right_node2(3)
        integer :: orientation

        status = 1
        call face_shape_3d( &
            nx_left, ny_left, nz_left, face_left, nu_left, nv_left)
        call face_shape_3d( &
            nx_right, ny_right, nz_right, face_right, nu_right, nv_right)
        if (min(nu_left, nv_left, nu_right, nv_right) < 2) return
        if (swap_axes) then
            if (nu_left /= nv_right .or. nv_left /= nu_right) return
        else
            if (nu_left /= nu_right .or. nv_left /= nv_right) return
        end if
        allocate( &
            h1_left(nu_left*nv_left), h1_right(nu_left*nv_left), &
            hcurl_left((nu_left - 1)*nv_left + nu_left*(nv_left - 1)), &
            hcurl_right((nu_left - 1)*nv_left + nu_left*(nv_left - 1)), &
            hcurl_sign((nu_left - 1)*nv_left + nu_left*(nv_left - 1)), &
            hdiv_left((nu_left - 1)*(nv_left - 1)), &
            hdiv_right((nu_left - 1)*(nv_left - 1)), &
            hdiv_sign((nu_left - 1)*(nv_left - 1)))
        do iv = 1, nv_left
            do iu = 1, nu_left
                edge = iu + (iv - 1)*nu_left
                call face_node_3d( &
                    nx_left, ny_left, nz_left, face_left, iu, iv, node1)
                call map_face_indices( &
                    iu, iv, nu_right, nv_right, swap_axes, reverse_u, &
                    reverse_v, ju1, jv1)
                call face_node_3d( &
                    nx_right, ny_right, nz_right, face_right, ju1, jv1, &
                    right_node1)
                h1_left(edge) = node_dof_3d( &
                    node1, nx_left, ny_left)
                h1_right(edge) = node_dof_3d( &
                    right_node1, nx_right, ny_right)
            end do
        end do

        edge = 0
        do iv = 1, nv_left
            do iu = 1, nu_left - 1
                edge = edge + 1
                call append_face_edge( &
                    iu, iv, iu + 1, iv, edge)
            end do
        end do
        do iv = 1, nv_left - 1
            do iu = 1, nu_left
                edge = edge + 1
                call append_face_edge( &
                    iu, iv, iu, iv + 1, edge)
            end do
        end do

        orientation = 1
        if (swap_axes) orientation = -orientation
        if (reverse_u) orientation = -orientation
        if (reverse_v) orientation = -orientation
        orientation = orientation*face_form_sign(face_left)* &
            face_form_sign(face_right)
        edge = 0
        do iv = 1, nv_left - 1
            do iu = 1, nu_left - 1
                edge = edge + 1
                call face_cell_dof_3d( &
                    nx_left, ny_left, nz_left, face_left, iu, iv, &
                    hdiv_left(edge))
                call map_face_indices( &
                    iu, iv, nu_right, nv_right, swap_axes, reverse_u, &
                    reverse_v, ju1, jv1)
                call map_face_indices( &
                    iu + 1, iv + 1, nu_right, nv_right, swap_axes, reverse_u, &
                    reverse_v, ju2, jv2)
                call face_cell_dof_3d( &
                    nx_right, ny_right, nz_right, face_right, min(ju1, ju2), &
                    min(jv1, jv2), hdiv_right(edge))
                hdiv_sign(edge) = orientation
            end do
        end do
        status = 0

    contains

        subroutine append_face_edge(iu1, iv1, iu2, iv2, local_edge)
            integer, intent(in) :: iu1, iv1, iu2, iv2, local_edge

            call face_node_3d( &
                nx_left, ny_left, nz_left, face_left, iu1, iv1, node1)
            call face_node_3d( &
                nx_left, ny_left, nz_left, face_left, iu2, iv2, node2)
            call map_face_indices( &
                iu1, iv1, nu_right, nv_right, swap_axes, reverse_u, reverse_v, &
                ju1, jv1)
            call map_face_indices( &
                iu2, iv2, nu_right, nv_right, swap_axes, reverse_u, reverse_v, &
                ju2, jv2)
            call face_node_3d( &
                nx_right, ny_right, nz_right, face_right, ju1, jv1, &
                right_node1)
            call face_node_3d( &
                nx_right, ny_right, nz_right, face_right, ju2, jv2, &
                right_node2)
            call edge_dof_3d( &
                node1, node2, nx_left, ny_left, nz_left, &
                hcurl_left(local_edge), orientation)
            call edge_dof_3d( &
                right_node1, right_node2, nx_right, ny_right, nz_right, &
                hcurl_right(local_edge), hcurl_sign(local_edge))
            hcurl_sign(local_edge) = hcurl_sign(local_edge)*orientation
        end subroutine append_face_edge

    end subroutine build_bspline_feec_3d_interface_dofs

    subroutine build_bspline_feec_3d_two_patch_maps( &
            nx_left, ny_left, nz_left, nx_right, ny_right, nz_right, &
            face_left, face_right, swap_axes, reverse_u, reverse_v, &
            h1_left_map, h1_right_map, hcurl_left_map, hcurl_right_map, &
            hcurl_left_sign, hcurl_right_sign, hdiv_left_map, hdiv_right_map, &
            hdiv_left_sign, hdiv_right_sign, l2_left_map, l2_right_map, status)
        integer, intent(in) :: nx_left, ny_left, nz_left
        integer, intent(in) :: nx_right, ny_right, nz_right
        integer, intent(in) :: face_left, face_right
        logical, intent(in) :: swap_axes, reverse_u, reverse_v
        integer, allocatable, intent(out) :: h1_left_map(:), h1_right_map(:)
        integer, allocatable, intent(out) :: hcurl_left_map(:)
        integer, allocatable, intent(out) :: hcurl_right_map(:)
        integer, allocatable, intent(out) :: hcurl_left_sign(:)
        integer, allocatable, intent(out) :: hcurl_right_sign(:)
        integer, allocatable, intent(out) :: hdiv_left_map(:)
        integer, allocatable, intent(out) :: hdiv_right_map(:)
        integer, allocatable, intent(out) :: hdiv_left_sign(:)
        integer, allocatable, intent(out) :: hdiv_right_sign(:)
        integer, allocatable, intent(out) :: l2_left_map(:), l2_right_map(:)
        integer, intent(out) :: status

        integer, allocatable :: h1_left_trace(:), h1_right_trace(:)
        integer, allocatable :: hcurl_left_trace(:), hcurl_right_trace(:)
        integer, allocatable :: hcurl_trace_sign(:)
        integer, allocatable :: hdiv_left_trace(:), hdiv_right_trace(:)
        integer, allocatable :: hdiv_trace_sign(:)
        integer :: dof, global_dof, trace

        call build_bspline_feec_3d_interface_dofs( &
            nx_left, ny_left, nz_left, nx_right, ny_right, nz_right, &
            face_left, face_right, swap_axes, reverse_u, reverse_v, &
            h1_left_trace, h1_right_trace, hcurl_left_trace, &
            hcurl_right_trace, hcurl_trace_sign, hdiv_left_trace, &
            hdiv_right_trace, hdiv_trace_sign, status)
        if (status /= 0) return
        allocate( &
            h1_left_map(nx_left*ny_left*nz_left), &
            h1_right_map(nx_right*ny_right*nz_right), &
            hcurl_left_map( &
            (nx_left - 1)*ny_left*nz_left + &
            nx_left*(ny_left - 1)*nz_left + &
            nx_left*ny_left*(nz_left - 1)), &
            hcurl_right_map( &
            (nx_right - 1)*ny_right*nz_right + &
            nx_right*(ny_right - 1)*nz_right + &
            nx_right*ny_right*(nz_right - 1)), &
            hdiv_left_map( &
            nx_left*(ny_left - 1)*(nz_left - 1) + &
            (nx_left - 1)*ny_left*(nz_left - 1) + &
            (nx_left - 1)*(ny_left - 1)*nz_left), &
            hdiv_right_map( &
            nx_right*(ny_right - 1)*(nz_right - 1) + &
            (nx_right - 1)*ny_right*(nz_right - 1) + &
            (nx_right - 1)*(ny_right - 1)*nz_right), &
            l2_left_map((nx_left - 1)*(ny_left - 1)*(nz_left - 1)), &
            l2_right_map((nx_right - 1)*(ny_right - 1)*(nz_right - 1)))
        allocate( &
            hcurl_left_sign(size(hcurl_left_map)), &
            hcurl_right_sign(size(hcurl_right_map)), &
            hdiv_left_sign(size(hdiv_left_map)), &
            hdiv_right_sign(size(hdiv_right_map)))
        call merge_unsigned_maps( &
            h1_left_trace, h1_right_trace, h1_left_map, h1_right_map)
        call merge_signed_maps( &
            hcurl_left_trace, hcurl_right_trace, hcurl_trace_sign, &
            hcurl_left_map, hcurl_right_map, hcurl_left_sign, &
            hcurl_right_sign)
        call merge_signed_maps( &
            hdiv_left_trace, hdiv_right_trace, hdiv_trace_sign, &
            hdiv_left_map, hdiv_right_map, hdiv_left_sign, hdiv_right_sign)
        l2_left_map = [(dof, dof = 1, size(l2_left_map))]
        l2_right_map = &
            [(size(l2_left_map) + dof, dof = 1, size(l2_right_map))]

    contains

        subroutine merge_unsigned_maps( &
                left_trace, right_trace, left_map, right_map)
            integer, intent(in) :: left_trace(:), right_trace(:)
            integer, intent(out) :: left_map(:), right_map(:)

            left_map = [(dof, dof = 1, size(left_map))]
            right_map = 0
            do trace = 1, size(left_trace)
                right_map(right_trace(trace)) = left_map(left_trace(trace))
            end do
            global_dof = size(left_map)
            do dof = 1, size(right_map)
                if (right_map(dof) /= 0) cycle
                global_dof = global_dof + 1
                right_map(dof) = global_dof
            end do
        end subroutine merge_unsigned_maps

        subroutine merge_signed_maps( &
                left_trace, right_trace, trace_sign, left_map, right_map, &
                left_sign, right_sign)
            integer, intent(in) :: left_trace(:), right_trace(:), trace_sign(:)
            integer, intent(out) :: left_map(:), right_map(:)
            integer, intent(out) :: left_sign(:), right_sign(:)

            left_map = [(dof, dof = 1, size(left_map))]
            left_sign = 1
            right_map = 0
            right_sign = 1
            do trace = 1, size(left_trace)
                right_map(right_trace(trace)) = left_map(left_trace(trace))
                right_sign(right_trace(trace)) = trace_sign(trace)
            end do
            global_dof = size(left_map)
            do dof = 1, size(right_map)
                if (right_map(dof) /= 0) cycle
                global_dof = global_dof + 1
                right_map(dof) = global_dof
            end do
        end subroutine merge_signed_maps

    end subroutine build_bspline_feec_3d_two_patch_maps

    pure subroutine face_shape_3d(nx, ny, nz, face, nu, nv)
        integer, intent(in) :: nx, ny, nz, face
        integer, intent(out) :: nu, nv

        select case (face)
        case (BSPLINE_FACE_X_MIN, BSPLINE_FACE_X_MAX)
            nu = ny
            nv = nz
        case (BSPLINE_FACE_Y_MIN, BSPLINE_FACE_Y_MAX)
            nu = nx
            nv = nz
        case (BSPLINE_FACE_Z_MIN, BSPLINE_FACE_Z_MAX)
            nu = nx
            nv = ny
        case default
            nu = 0
            nv = 0
        end select
    end subroutine face_shape_3d

    pure subroutine map_face_indices( &
            iu, iv, nu_right, nv_right, swap_axes, reverse_u, reverse_v, ju, jv)
        integer, intent(in) :: iu, iv, nu_right, nv_right
        logical, intent(in) :: swap_axes, reverse_u, reverse_v
        integer, intent(out) :: ju, jv

        if (swap_axes) then
            ju = iv
            jv = iu
        else
            ju = iu
            jv = iv
        end if
        if (reverse_u) ju = nu_right + 1 - ju
        if (reverse_v) jv = nv_right + 1 - jv
    end subroutine map_face_indices

    pure subroutine face_node_3d(nx, ny, nz, face, iu, iv, node)
        integer, intent(in) :: nx, ny, nz, face, iu, iv
        integer, intent(out) :: node(3)

        select case (face)
        case (BSPLINE_FACE_X_MIN, BSPLINE_FACE_X_MAX)
            node = [1, iu, iv]
            if (face == BSPLINE_FACE_X_MAX) node(1) = nx
        case (BSPLINE_FACE_Y_MIN, BSPLINE_FACE_Y_MAX)
            node = [iu, 1, iv]
            if (face == BSPLINE_FACE_Y_MAX) node(2) = ny
        case default
            node = [iu, iv, 1]
            if (face == BSPLINE_FACE_Z_MAX) node(3) = nz
        end select
    end subroutine face_node_3d

    pure integer function node_dof_3d(node, nx, ny) result(dof)
        integer, intent(in) :: node(3), nx, ny

        dof = node(1) + (node(2) - 1)*nx + &
            (node(3) - 1)*nx*ny
    end function node_dof_3d

    pure subroutine edge_dof_3d(node1, node2, nx, ny, nz, dof, sign)
        integer, intent(in) :: node1(3), node2(3), nx, ny, nz
        integer, intent(out) :: dof, sign
        integer :: axis, lower(3)

        axis = maxloc(abs(node2 - node1), dim=1)
        lower = min(node1, node2)
        sign = merge(1, -1, node2(axis) > node1(axis))
        select case (axis)
        case (1)
            dof = lower(1) + (lower(2) - 1)*(nx - 1) + &
                (lower(3) - 1)*(nx - 1)*ny
        case (2)
            dof = (nx - 1)*ny*nz + lower(1) + &
                (lower(2) - 1)*nx + (lower(3) - 1)*nx*(ny - 1)
        case default
            dof = (nx - 1)*ny*nz + nx*(ny - 1)*nz + lower(1) + &
                (lower(2) - 1)*nx + (lower(3) - 1)*nx*ny
        end select
    end subroutine edge_dof_3d

    pure subroutine face_cell_dof_3d( &
            nx, ny, nz, face, iu, iv, dof)
        integer, intent(in) :: nx, ny, nz, face, iu, iv
        integer, intent(out) :: dof
        integer :: bx_count, by_count, fixed

        bx_count = nx*(ny - 1)*(nz - 1)
        by_count = (nx - 1)*ny*(nz - 1)
        select case (face)
        case (BSPLINE_FACE_X_MIN, BSPLINE_FACE_X_MAX)
            fixed = 1
            if (face == BSPLINE_FACE_X_MAX) fixed = nx
            dof = fixed + (iu - 1)*nx + (iv - 1)*nx*(ny - 1)
        case (BSPLINE_FACE_Y_MIN, BSPLINE_FACE_Y_MAX)
            fixed = 1
            if (face == BSPLINE_FACE_Y_MAX) fixed = ny
            dof = bx_count + iu + (fixed - 1)*(nx - 1) + &
                (iv - 1)*(nx - 1)*ny
        case default
            fixed = 1
            if (face == BSPLINE_FACE_Z_MAX) fixed = nz
            dof = bx_count + by_count + iu + (iv - 1)*(nx - 1) + &
                (fixed - 1)*(nx - 1)*(ny - 1)
        end select
    end subroutine face_cell_dof_3d

    pure integer function face_form_sign(face) result(sign)
        integer, intent(in) :: face

        sign = 1
        if (face == BSPLINE_FACE_Y_MIN .or. &
            face == BSPLINE_FACE_Y_MAX) sign = -1
    end function face_form_sign

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
