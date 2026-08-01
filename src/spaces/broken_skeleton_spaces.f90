module fortfem_broken_skeleton_spaces
    !! Frozen ownership maps for broken volume and skeleton spaces.
    !!
    !! A broken space gives every cell an independent local copy of its
    !! degrees of freedom.  A skeleton space gives every facet one shared
    !! trace block and records the adjacent cells and their orientations.
    !! The module owns topology and numbering only; basis evaluation, metric
    !! maps, and sparse assembly remain caller-owned layers.
    implicit none
    private

    integer, parameter, public :: BROKEN_SPACE_H1 = 1
    integer, parameter, public :: BROKEN_SPACE_HCURL = 2
    integer, parameter, public :: BROKEN_SPACE_HDIV = 3
    integer, parameter, public :: BROKEN_SPACE_L2 = 4

    integer, parameter, public :: SKELETON_SPACE_SCALAR = 1
    integer, parameter, public :: SKELETON_SPACE_TANGENTIAL = 2
    integer, parameter, public :: SKELETON_SPACE_NORMAL = 3

    type, public :: broken_space_layout_t
        private
        integer :: family_kind = 0
        integer :: cell_count = 0
        integer :: local_dof_count = 0
        integer :: global_dof_count = 0
        integer, allocatable :: local_to_global(:, :)
        integer, allocatable :: local_sign(:, :)
    end type broken_space_layout_t

    type, public :: skeleton_space_layout_t
        private
        integer :: trace_kind = 0
        integer :: cell_count = 0
        integer :: facet_count = 0
        integer :: local_dof_count = 0
        integer :: global_dof_count = 0
        integer, allocatable :: facet_to_global(:, :)
        integer, allocatable :: facet_to_cell(:, :)
        integer, allocatable :: facet_sign(:, :)
    end type skeleton_space_layout_t

    public :: initialize_broken_space_layout
    public :: validate_broken_space_layout
    public :: broken_space_layout_maps
    public :: broken_space_layout_global_count
    public :: initialize_skeleton_space_layout
    public :: validate_skeleton_space_layout
    public :: skeleton_space_layout_maps
    public :: skeleton_space_layout_global_count

contains

    subroutine initialize_broken_space_layout( &
            layout, family_kind, cell_count, local_dof_count, local_sign, status)
        type(broken_space_layout_t), intent(inout) :: layout
        integer, intent(in) :: family_kind, cell_count, local_dof_count
        integer, intent(in), optional :: local_sign(:, :)
        integer, intent(out) :: status
        integer :: cell, local

        call clear_broken_layout(layout)
        status = 1
        if (.not. valid_broken_family(family_kind)) return
        if (cell_count < 1 .or. local_dof_count < 1) return
        if (present(local_sign)) then
            if (size(local_sign, 1) /= local_dof_count .or. &
                size(local_sign, 2) /= cell_count) return
            if (any(abs(local_sign) /= 1)) return
            if (family_kind == BROKEN_SPACE_H1 .or. &
                family_kind == BROKEN_SPACE_L2) then
                if (any(local_sign /= 1)) return
            end if
        end if

        layout%family_kind = family_kind
        layout%cell_count = cell_count
        layout%local_dof_count = local_dof_count
        layout%global_dof_count = cell_count*local_dof_count
        allocate(layout%local_to_global(local_dof_count, cell_count))
        allocate(layout%local_sign(local_dof_count, cell_count))
        do cell = 1, cell_count
            do local = 1, local_dof_count
                layout%local_to_global(local, cell) = &
                    (cell - 1)*local_dof_count + local
            end do
        end do
        layout%local_sign = 1
        if (present(local_sign)) layout%local_sign = local_sign
        status = 0
    end subroutine initialize_broken_space_layout

    subroutine validate_broken_space_layout(layout, status)
        type(broken_space_layout_t), intent(in) :: layout
        integer, intent(out) :: status
        integer :: cell, local

        status = 1
        if (.not. valid_broken_family(layout%family_kind)) return
        if (layout%cell_count < 1 .or. layout%local_dof_count < 1) return
        if (layout%global_dof_count /= &
            layout%cell_count*layout%local_dof_count) return
        if (.not. allocated(layout%local_to_global)) return
        if (.not. allocated(layout%local_sign)) return
        if (size(layout%local_to_global, 1) /= layout%local_dof_count .or. &
            size(layout%local_to_global, 2) /= layout%cell_count) return
        if (size(layout%local_sign, 1) /= layout%local_dof_count .or. &
            size(layout%local_sign, 2) /= layout%cell_count) return
        if (any(abs(layout%local_sign) /= 1)) return
        if (layout%family_kind == BROKEN_SPACE_H1 .or. &
            layout%family_kind == BROKEN_SPACE_L2) then
            if (any(layout%local_sign /= 1)) return
        end if
        do cell = 1, layout%cell_count
            do local = 1, layout%local_dof_count
                if (layout%local_to_global(local, cell) /= &
                    (cell - 1)*layout%local_dof_count + local) return
            end do
        end do
        status = 0
    end subroutine validate_broken_space_layout

    subroutine broken_space_layout_maps(layout, local_to_global, local_sign, status)
        type(broken_space_layout_t), intent(in) :: layout
        integer, allocatable, intent(out) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(out) :: status

        if (allocated(local_to_global)) deallocate(local_to_global)
        if (allocated(local_sign)) deallocate(local_sign)
        call validate_broken_space_layout(layout, status)
        if (status /= 0) return
        allocate(local_to_global(size(layout%local_to_global, 1), &
            size(layout%local_to_global, 2)))
        allocate(local_sign(size(layout%local_sign, 1), size(layout%local_sign, 2)))
        local_to_global = layout%local_to_global
        local_sign = layout%local_sign
    end subroutine broken_space_layout_maps

    integer function broken_space_layout_global_count(layout)
        type(broken_space_layout_t), intent(in) :: layout

        broken_space_layout_global_count = layout%global_dof_count
    end function broken_space_layout_global_count

    subroutine initialize_skeleton_space_layout( &
            layout, trace_kind, cell_count, facet_to_cell, facet_sign, &
            local_dof_count, status)
        type(skeleton_space_layout_t), intent(inout) :: layout
        integer, intent(in) :: trace_kind, cell_count
        integer, intent(in) :: facet_to_cell(:, :), facet_sign(:, :)
        integer, intent(in) :: local_dof_count
        integer, intent(out) :: status
        integer :: facet, local, side, facet_count

        call clear_skeleton_layout(layout)
        status = 1
        if (.not. valid_skeleton_kind(trace_kind)) return
        if (cell_count < 1 .or. local_dof_count < 1) return
        if (size(facet_to_cell, 1) /= 2 .or. &
            size(facet_sign, 1) /= 2) return
        if (size(facet_to_cell, 2) < 1 .or. &
            size(facet_sign, 2) /= size(facet_to_cell, 2)) return
        facet_count = size(facet_to_cell, 2)
        if (any(facet_to_cell < 0) .or. &
            any(facet_to_cell > cell_count)) return
        if (any(abs(facet_sign) > 1)) return
        do facet = 1, facet_count
            if (facet_to_cell(1, facet) == 0 .and. &
                facet_to_cell(2, facet) == 0) return
            do side = 1, 2
                if (facet_to_cell(side, facet) == 0) then
                    if (facet_sign(side, facet) /= 0) return
                else
                    if (abs(facet_sign(side, facet)) /= 1) return
                end if
            end do
        end do

        layout%trace_kind = trace_kind
        layout%cell_count = cell_count
        layout%facet_count = facet_count
        layout%local_dof_count = local_dof_count
        layout%global_dof_count = facet_count*local_dof_count
        allocate(layout%facet_to_global(local_dof_count, facet_count))
        allocate(layout%facet_to_cell(2, facet_count))
        allocate(layout%facet_sign(2, facet_count))
        do facet = 1, facet_count
            do local = 1, local_dof_count
                layout%facet_to_global(local, facet) = &
                    (facet - 1)*local_dof_count + local
            end do
        end do
        layout%facet_to_cell = facet_to_cell
        layout%facet_sign = facet_sign
        status = 0
    end subroutine initialize_skeleton_space_layout

    subroutine validate_skeleton_space_layout(layout, status)
        type(skeleton_space_layout_t), intent(in) :: layout
        integer, intent(out) :: status
        integer :: facet, local, side

        status = 1
        if (.not. valid_skeleton_kind(layout%trace_kind)) return
        if (layout%cell_count < 1 .or. layout%facet_count < 1 .or. &
            layout%local_dof_count < 1) return
        if (layout%global_dof_count /= &
            layout%facet_count*layout%local_dof_count) return
        if (.not. allocated(layout%facet_to_global)) return
        if (.not. allocated(layout%facet_to_cell)) return
        if (.not. allocated(layout%facet_sign)) return
        if (size(layout%facet_to_global, 1) /= layout%local_dof_count .or. &
            size(layout%facet_to_global, 2) /= layout%facet_count) return
        if (size(layout%facet_to_cell, 1) /= 2 .or. &
            size(layout%facet_to_cell, 2) /= layout%facet_count) return
        if (size(layout%facet_sign, 1) /= 2 .or. &
            size(layout%facet_sign, 2) /= layout%facet_count) return
        if (any(layout%facet_to_cell < 0) .or. &
            any(layout%facet_to_cell > layout%cell_count)) return
        if (any(abs(layout%facet_sign) > 1)) return
        do facet = 1, layout%facet_count
            if (layout%facet_to_cell(1, facet) == 0 .and. &
                layout%facet_to_cell(2, facet) == 0) return
            do side = 1, 2
                if (layout%facet_to_cell(side, facet) == 0) then
                    if (layout%facet_sign(side, facet) /= 0) return
                else
                    if (abs(layout%facet_sign(side, facet)) /= 1) return
                end if
            end do
            do local = 1, layout%local_dof_count
                if (layout%facet_to_global(local, facet) /= &
                    (facet - 1)*layout%local_dof_count + local) return
            end do
        end do
        status = 0
    end subroutine validate_skeleton_space_layout

    subroutine skeleton_space_layout_maps( &
            layout, facet_to_global, facet_to_cell, facet_sign, status)
        type(skeleton_space_layout_t), intent(in) :: layout
        integer, allocatable, intent(out) :: facet_to_global(:, :)
        integer, allocatable, intent(out) :: facet_to_cell(:, :), facet_sign(:, :)
        integer, intent(out) :: status

        if (allocated(facet_to_global)) deallocate(facet_to_global)
        if (allocated(facet_to_cell)) deallocate(facet_to_cell)
        if (allocated(facet_sign)) deallocate(facet_sign)
        call validate_skeleton_space_layout(layout, status)
        if (status /= 0) return
        allocate(facet_to_global(size(layout%facet_to_global, 1), &
            size(layout%facet_to_global, 2)))
        allocate(facet_to_cell(size(layout%facet_to_cell, 1), &
            size(layout%facet_to_cell, 2)))
        allocate(facet_sign(size(layout%facet_sign, 1), &
            size(layout%facet_sign, 2)))
        facet_to_global = layout%facet_to_global
        facet_to_cell = layout%facet_to_cell
        facet_sign = layout%facet_sign
    end subroutine skeleton_space_layout_maps

    integer function skeleton_space_layout_global_count(layout)
        type(skeleton_space_layout_t), intent(in) :: layout

        skeleton_space_layout_global_count = layout%global_dof_count
    end function skeleton_space_layout_global_count

    pure logical function valid_broken_family(family_kind)
        integer, intent(in) :: family_kind

        valid_broken_family = family_kind >= BROKEN_SPACE_H1 .and. &
            family_kind <= BROKEN_SPACE_L2
    end function valid_broken_family

    pure logical function valid_skeleton_kind(trace_kind)
        integer, intent(in) :: trace_kind

        valid_skeleton_kind = trace_kind >= SKELETON_SPACE_SCALAR .and. &
            trace_kind <= SKELETON_SPACE_NORMAL
    end function valid_skeleton_kind

    subroutine clear_broken_layout(layout)
        type(broken_space_layout_t), intent(inout) :: layout

        layout%family_kind = 0
        layout%cell_count = 0
        layout%local_dof_count = 0
        layout%global_dof_count = 0
        if (allocated(layout%local_to_global)) deallocate(layout%local_to_global)
        if (allocated(layout%local_sign)) deallocate(layout%local_sign)
    end subroutine clear_broken_layout

    subroutine clear_skeleton_layout(layout)
        type(skeleton_space_layout_t), intent(inout) :: layout

        layout%trace_kind = 0
        layout%cell_count = 0
        layout%facet_count = 0
        layout%local_dof_count = 0
        layout%global_dof_count = 0
        if (allocated(layout%facet_to_global)) deallocate(layout%facet_to_global)
        if (allocated(layout%facet_to_cell)) deallocate(layout%facet_to_cell)
        if (allocated(layout%facet_sign)) deallocate(layout%facet_sign)
    end subroutine clear_skeleton_layout

end module fortfem_broken_skeleton_spaces
