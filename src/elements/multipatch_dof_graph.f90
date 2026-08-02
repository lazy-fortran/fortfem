module fortfem_multipatch_dof_graph
    !! Fixed-topology signed local-to-global maps for arbitrary patch graphs.
    !!
    !! Each interface relation declares x_left = sign*x_right.  The output
    !! map uses signed global IDs, matching the existing glued FEEC
    !! assemblers.  Geometry, spline bases, trace quadrature, and physical
    !! interface laws remain caller-owned.
    implicit none
    private

    public :: build_multipatch_signed_dof_map

contains

    subroutine build_multipatch_signed_dof_map( &
            patch_offsets, left_patch, left_local, right_patch, right_local, &
            interface_sign, local_to_global, global_count, status)
        integer, intent(in) :: patch_offsets(:)
        integer, intent(in) :: left_patch(:), left_local(:)
        integer, intent(in) :: right_patch(:), right_local(:)
        integer, intent(in) :: interface_sign(:)
        integer, intent(out) :: local_to_global(:)
        integer, intent(out) :: global_count, status

        integer, allocatable :: parent(:), parity(:), root_id(:)
        integer :: local_count, relation_count
        integer :: relation, left, right, left_root, right_root
        integer :: left_sign, right_sign, root, sign

        local_to_global = 0
        global_count = 0
        status = 1
        if (.not. valid_inputs( &
            patch_offsets, left_patch, left_local, right_patch, right_local, &
            interface_sign, local_to_global, local_count)) return
        relation_count = size(left_patch)
        allocate(parent(local_count), parity(local_count), root_id(local_count))
        parent = [(root, root = 1, local_count)]
        parity = 1
        root_id = 0

        do relation = 1, relation_count
            left = patch_offsets(left_patch(relation)) + &
                left_local(relation) - 1
            right = patch_offsets(right_patch(relation)) + &
                right_local(relation) - 1
            call find_signed(left, parent, parity, left_root, left_sign)
            call find_signed(right, parent, parity, right_root, right_sign)
            if (left_root == right_root) then
                if (left_sign /= interface_sign(relation)*right_sign) then
                    status = 2
                    return
                end if
            else
                parent(left_root) = right_root
                parity(left_root) = interface_sign(relation)*right_sign*left_sign
            end if
        end do

        do local_count = 1, size(parent)
            call find_signed(local_count, parent, parity, root, sign)
            if (root_id(root) == 0) then
                global_count = global_count + 1
                root_id(root) = global_count
            end if
            local_to_global(local_count) = sign*root_id(root)
        end do
        status = 0
    end subroutine build_multipatch_signed_dof_map

    recursive subroutine find_signed(node, parent, parity, root, sign)
        integer, intent(in) :: node
        integer, intent(inout) :: parent(:), parity(:)
        integer, intent(out) :: root, sign
        integer :: parent_root, parent_sign

        if (parent(node) == node) then
            root = node
            sign = 1
            return
        end if
        call find_signed(parent(node), parent, parity, parent_root, parent_sign)
        root = parent_root
        sign = parity(node)*parent_sign
        parent(node) = root
        parity(node) = sign
    end subroutine find_signed

    logical function valid_inputs( &
            patch_offsets, left_patch, left_local, right_patch, right_local, &
            interface_sign, local_to_global, local_count) result(valid)
        integer, intent(in) :: patch_offsets(:)
        integer, intent(in) :: left_patch(:), left_local(:)
        integer, intent(in) :: right_patch(:), right_local(:)
        integer, intent(in) :: interface_sign(:), local_to_global(:)
        integer, intent(out) :: local_count
        integer :: patch, relation_count, patch_size

        local_count = 0
        valid = .false.
        if (size(patch_offsets) < 2) return
        if (patch_offsets(1) /= 1) return
        do patch = 1, size(patch_offsets) - 1
            if (patch_offsets(patch + 1) <= patch_offsets(patch)) return
        end do
        local_count = patch_offsets(size(patch_offsets)) - 1
        if (local_count < 2 .or. size(local_to_global) /= local_count) return
        relation_count = size(left_patch)
        if (size(left_local) /= relation_count .or. &
            size(right_patch) /= relation_count .or. &
            size(right_local) /= relation_count .or. &
            size(interface_sign) /= relation_count) return
        do relation_count = 1, size(left_patch)
            if (.not. valid_endpoint( &
                left_patch(relation_count), left_local(relation_count), &
                patch_offsets)) return
            if (.not. valid_endpoint( &
                right_patch(relation_count), right_local(relation_count), &
                patch_offsets)) return
            if (abs(interface_sign(relation_count)) /= 1) return
        end do
        do patch = 1, size(patch_offsets) - 1
            patch_size = patch_offsets(patch + 1) - patch_offsets(patch)
            if (patch_size < 1) return
        end do
        valid = .true.
    end function valid_inputs

    logical function valid_endpoint(patch, local, patch_offsets) result(valid)
        integer, intent(in) :: patch, local, patch_offsets(:)

        valid = patch >= 1 .and. patch < size(patch_offsets)
        if (.not. valid) return
        valid = local >= 1 .and. local <= &
            patch_offsets(patch + 1) - patch_offsets(patch)
    end function valid_endpoint

end module fortfem_multipatch_dof_graph
