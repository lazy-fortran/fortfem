module fortfem_equation_objective_registry
    !! Neutral metadata and deterministic packing for external physics clients.
    !!
    !! The registry deliberately contains no constitutive law or callback.  A
    !! client declares named equation, objective, and constraint blocks; the
    !! registry assigns their offsets in declaration order and validates the
    !! composition.  The packing actions are only the linear bookkeeping
    !! needed to move block-local residual samples to one deterministic vector.
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    integer, parameter, public :: REGISTRY_KIND_EQUATION = 1
    integer, parameter, public :: REGISTRY_KIND_OBJECTIVE = 2
    integer, parameter, public :: REGISTRY_KIND_CONSTRAINT = 3

    type, public :: equation_objective_block_t
        character(len=:), allocatable :: name
        integer :: kind = 0
        integer :: row_offset = 0
        integer :: row_count = 0
        character(len=:), allocatable :: units
        character(len=:), allocatable :: normalization
        character(len=:), allocatable :: provenance
        logical :: active = .false.
        logical :: fixed = .false.
    end type equation_objective_block_t

    type, public :: equation_objective_registry_t
        type(equation_objective_block_t), allocatable :: blocks(:)
        integer :: total_rows = 0
    contains
        procedure, private :: assign_equation_objective_registry
        generic :: assignment(=) => assign_equation_objective_registry
    end type equation_objective_registry_t

    public :: initialize_equation_objective_registry
    public :: validate_equation_objective_registry
    public :: equation_objective_registry_block
    public :: equation_objective_registry_block_count
    public :: equation_objective_registry_total_rows
    public :: pack_equation_objective_values
    public :: pack_equation_objective_values_jvp
    public :: pack_equation_objective_values_vjp
    public :: unpack_equation_objective_values

contains

    subroutine assign_equation_objective_registry(lhs, rhs)
        class(equation_objective_registry_t), intent(out) :: lhs
        type(equation_objective_registry_t), intent(in) :: rhs
        integer :: block

        lhs%total_rows = rhs%total_rows
        if (.not. allocated(rhs%blocks)) return
        allocate(lhs%blocks(size(rhs%blocks)))
        do block = 1, size(rhs%blocks)
            lhs%blocks(block) = rhs%blocks(block)
        end do
    end subroutine assign_equation_objective_registry

    subroutine initialize_equation_objective_registry( &
            registry, names, kinds, row_counts, units, normalization, &
            provenance, active, fixed, status)
        type(equation_objective_registry_t), intent(out) :: registry
        character(len=*), intent(in) :: names(:), units(:), normalization(:)
        character(len=*), intent(in) :: provenance(:)
        integer, intent(in) :: kinds(:), row_counts(:)
        logical, intent(in) :: active(:), fixed(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: block, offset

        call clear_registry(registry)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "equation/objective registry received incompatible arrays")
        if (size(names) < 1 .or. size(kinds) /= size(names) .or. &
            size(row_counts) /= size(names) .or. size(units) /= size(names) .or. &
            size(normalization) /= size(names) .or. &
            size(provenance) /= size(names) .or. size(active) /= size(names) .or. &
            size(fixed) /= size(names)) return

        allocate(registry%blocks(size(names)))
        offset = 1
        do block = 1, size(names)
            registry%blocks(block)%name = trim(names(block))
            registry%blocks(block)%kind = kinds(block)
            registry%blocks(block)%row_offset = offset
            registry%blocks(block)%row_count = row_counts(block)
            registry%blocks(block)%units = trim(units(block))
            registry%blocks(block)%normalization = trim(normalization(block))
            registry%blocks(block)%provenance = trim(provenance(block))
            registry%blocks(block)%active = active(block)
            registry%blocks(block)%fixed = fixed(block)
            if (row_counts(block) >= 0) offset = offset + row_counts(block)
        end do
        registry%total_rows = offset - 1
        if (.not. validate_equation_objective_registry(registry, status)) return
    end subroutine initialize_equation_objective_registry

    logical function validate_equation_objective_registry(registry, status) &
            result(valid)
        type(equation_objective_registry_t), intent(in) :: registry
        type(fortsparse_status_t), intent(out) :: status
        integer :: block, other, expected_offset, expected_rows

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "equation/objective registry has invalid metadata")
        if (.not. allocated(registry%blocks)) return
        if (size(registry%blocks) < 1) return
        expected_offset = 1
        expected_rows = 0
        do block = 1, size(registry%blocks)
            if (.not. valid_kind(registry%blocks(block)%kind)) return
            if (.not. allocated(registry%blocks(block)%name)) return
            if (.not. allocated(registry%blocks(block)%units)) return
            if (.not. allocated(registry%blocks(block)%normalization)) return
            if (.not. allocated(registry%blocks(block)%provenance)) return
            if (len_trim(registry%blocks(block)%name) < 1) return
            if (len_trim(registry%blocks(block)%units) < 1) return
            if (len_trim(registry%blocks(block)%normalization) < 1) return
            if (len_trim(registry%blocks(block)%provenance) < 1) return
            if (registry%blocks(block)%row_count < 0 .or. &
                registry%blocks(block)%row_offset /= expected_offset) return
            do other = block + 1, size(registry%blocks)
                if (registry%blocks(block)%name == &
                    registry%blocks(other)%name) return
            end do
            expected_offset = expected_offset + registry%blocks(block)%row_count
            expected_rows = expected_rows + registry%blocks(block)%row_count
        end do
        if (registry%total_rows /= expected_rows) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_equation_objective_registry

    subroutine equation_objective_registry_block( &
            registry, block_index, block, status)
        type(equation_objective_registry_t), intent(in) :: registry
        integer, intent(in) :: block_index
        type(equation_objective_block_t), intent(out) :: block
        type(fortsparse_status_t), intent(out) :: status

        call clear_block(block)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "equation/objective registry block index is invalid")
        if (.not. validate_equation_objective_registry(registry, status)) return
        if (block_index < 1 .or. block_index > size(registry%blocks)) return
        block = registry%blocks(block_index)
            call status_set(status, FORTSPARSE_OK, "")
        end subroutine equation_objective_registry_block

        integer function equation_objective_registry_block_count(registry)
            type(equation_objective_registry_t), intent(in) :: registry

            if (allocated(registry%blocks)) then
                equation_objective_registry_block_count = size(registry%blocks)
            else
                equation_objective_registry_block_count = 0
            end if
        end function equation_objective_registry_block_count

        integer function equation_objective_registry_total_rows(registry)
            type(equation_objective_registry_t), intent(in) :: registry

            equation_objective_registry_total_rows = registry%total_rows
        end function equation_objective_registry_total_rows

        subroutine pack_equation_objective_values( &
                registry, block_values, packed, status)
            type(equation_objective_registry_t), intent(in) :: registry
            real(dp), intent(in) :: block_values(:, :)
            real(dp), allocatable, intent(out) :: packed(:)
            type(fortsparse_status_t), intent(out) :: status
            integer :: block, row, row_count, max_rows

            if (allocated(packed)) deallocate(packed)
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "equation/objective packing received incompatible arrays")
            if (.not. validate_equation_objective_registry(registry, status)) then
                allocate(packed(0))
                return
            end if
            max_rows = maxval(registry%blocks%row_count)
            if (size(block_values, 2) /= size(registry%blocks) .or. &
                size(block_values, 1) < max_rows) then
                allocate(packed(0))
                return
            end if
            allocate(packed(registry%total_rows))
            packed = 0.0_dp
            do block = 1, size(registry%blocks)
                row_count = registry%blocks(block)%row_count
                do row = 1, row_count
                    packed(registry%blocks(block)%row_offset + row - 1) = &
                        block_values(row, block)
                end do
            end do
            call status_set(status, FORTSPARSE_OK, "")
        end subroutine pack_equation_objective_values

        subroutine pack_equation_objective_values_jvp( &
                registry, block_values_dot, packed_dot, status)
            type(equation_objective_registry_t), intent(in) :: registry
            real(dp), intent(in) :: block_values_dot(:, :)
            real(dp), allocatable, intent(out) :: packed_dot(:)
            type(fortsparse_status_t), intent(out) :: status

            call pack_equation_objective_values( &
                registry, block_values_dot, packed_dot, status)
        end subroutine pack_equation_objective_values_jvp

        subroutine pack_equation_objective_values_vjp( &
                registry, packed_bar, block_values_bar, status)
            type(equation_objective_registry_t), intent(in) :: registry
            real(dp), intent(in) :: packed_bar(:)
            real(dp), allocatable, intent(out) :: block_values_bar(:, :)
            type(fortsparse_status_t), intent(out) :: status
            integer :: block, row, max_rows

            if (allocated(block_values_bar)) deallocate(block_values_bar)
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "equation/objective unpacking received incompatible arrays")
            if (.not. validate_equation_objective_registry(registry, status)) then
                allocate(block_values_bar(0, 0))
                return
            end if
            if (size(packed_bar) /= registry%total_rows) then
                allocate(block_values_bar(0, 0))
                return
            end if
            max_rows = maxval(registry%blocks%row_count)
            allocate(block_values_bar(max_rows, size(registry%blocks)))
            block_values_bar = 0.0_dp
            do block = 1, size(registry%blocks)
                do row = 1, registry%blocks(block)%row_count
                    block_values_bar(row, block) = packed_bar( &
                        registry%blocks(block)%row_offset + row - 1)
                end do
            end do
            call status_set(status, FORTSPARSE_OK, "")
        end subroutine pack_equation_objective_values_vjp

        subroutine unpack_equation_objective_values( &
                registry, packed, block_values, status)
            type(equation_objective_registry_t), intent(in) :: registry
            real(dp), intent(in) :: packed(:)
            real(dp), allocatable, intent(out) :: block_values(:, :)
            type(fortsparse_status_t), intent(out) :: status

            call pack_equation_objective_values_vjp( &
                registry, packed, block_values, status)
        end subroutine unpack_equation_objective_values

        logical function valid_kind(kind_value)
            integer, intent(in) :: kind_value

            valid_kind = kind_value == REGISTRY_KIND_EQUATION .or. &
                kind_value == REGISTRY_KIND_OBJECTIVE .or. &
                kind_value == REGISTRY_KIND_CONSTRAINT
        end function valid_kind

        subroutine clear_registry(registry)
            type(equation_objective_registry_t), intent(inout) :: registry

            if (allocated(registry%blocks)) deallocate(registry%blocks)
            registry%total_rows = 0
        end subroutine clear_registry

        subroutine clear_block(block)
            type(equation_objective_block_t), intent(out) :: block

            if (allocated(block%name)) deallocate(block%name)
            if (allocated(block%units)) deallocate(block%units)
            if (allocated(block%normalization)) deallocate(block%normalization)
            if (allocated(block%provenance)) deallocate(block%provenance)
            block%kind = 0
            block%row_offset = 0
            block%row_count = 0
            block%active = .false.
            block%fixed = .false.
        end subroutine clear_block

    end module fortfem_equation_objective_registry
