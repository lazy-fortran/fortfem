module fortfem_equation_objective_metadata
    !! Neutral numeric metadata for equation, objective, and constraint rows.
    !!
    !! This record is intentionally independent of any physics, optimizer, or
    !! profile convention.  A caller supplies one target/bound/weight/scale
    !! tuple per row and may attach block-level KKT, nullspace, parameter
    !! tangent, and Hessian-vector capability identifiers.  The identifiers
    !! are opaque to FortFEM; their only contract here is deterministic,
    !! non-negative validation.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    integer, parameter, public :: OBJECTIVE_METADATA_KIND_EQUATION = 1
    integer, parameter, public :: OBJECTIVE_METADATA_KIND_OBJECTIVE = 2
    integer, parameter, public :: OBJECTIVE_METADATA_KIND_CONSTRAINT = 3
    integer, parameter, public :: OBJECTIVE_METADATA_UNSET_ID = 0

    type, public :: equation_objective_metadata_t
        character(len=:), allocatable :: name
        character(len=:), allocatable :: units
        character(len=:), allocatable :: provenance
        integer :: kind = 0
        integer :: value_count = 0
        real(dp), allocatable :: target(:)
        real(dp), allocatable :: lower_bound(:)
        real(dp), allocatable :: upper_bound(:)
        real(dp), allocatable :: weight(:)
        real(dp), allocatable :: scale(:)
        logical :: active = .false.
        logical :: fixed = .false.
        integer :: kkt_id = OBJECTIVE_METADATA_UNSET_ID
        integer :: nullspace_id = OBJECTIVE_METADATA_UNSET_ID
        logical :: has_parameter_tangent = .false.
        integer :: parameter_tangent_id = OBJECTIVE_METADATA_UNSET_ID
        logical :: has_hvp = .false.
        integer :: hvp_id = OBJECTIVE_METADATA_UNSET_ID
    contains
        procedure, private :: assign_equation_objective_metadata
        generic :: assignment(=) => assign_equation_objective_metadata
    end type equation_objective_metadata_t

    public :: initialize_equation_objective_metadata
    public :: validate_equation_objective_metadata
    public :: clear_equation_objective_metadata

contains

    subroutine assign_equation_objective_metadata(lhs, rhs)
        class(equation_objective_metadata_t), intent(out) :: lhs
        type(equation_objective_metadata_t), intent(in) :: rhs

        lhs%kind = rhs%kind
        lhs%value_count = rhs%value_count
        lhs%active = rhs%active
        lhs%fixed = rhs%fixed
        lhs%kkt_id = rhs%kkt_id
        lhs%nullspace_id = rhs%nullspace_id
        lhs%has_parameter_tangent = rhs%has_parameter_tangent
        lhs%parameter_tangent_id = rhs%parameter_tangent_id
        lhs%has_hvp = rhs%has_hvp
        lhs%hvp_id = rhs%hvp_id
        if (allocated(rhs%name)) lhs%name = rhs%name
        if (allocated(rhs%units)) lhs%units = rhs%units
        if (allocated(rhs%provenance)) lhs%provenance = rhs%provenance
        if (allocated(rhs%target)) lhs%target = rhs%target
        if (allocated(rhs%lower_bound)) lhs%lower_bound = rhs%lower_bound
        if (allocated(rhs%upper_bound)) lhs%upper_bound = rhs%upper_bound
        if (allocated(rhs%weight)) lhs%weight = rhs%weight
        if (allocated(rhs%scale)) lhs%scale = rhs%scale
    end subroutine assign_equation_objective_metadata

    subroutine clear_equation_objective_metadata(metadata)
        type(equation_objective_metadata_t), intent(out) :: metadata

        metadata%kind = 0
        metadata%value_count = 0
        metadata%active = .false.
        metadata%fixed = .false.
        metadata%kkt_id = OBJECTIVE_METADATA_UNSET_ID
        metadata%nullspace_id = OBJECTIVE_METADATA_UNSET_ID
        metadata%has_parameter_tangent = .false.
        metadata%parameter_tangent_id = OBJECTIVE_METADATA_UNSET_ID
        metadata%has_hvp = .false.
        metadata%hvp_id = OBJECTIVE_METADATA_UNSET_ID
    end subroutine clear_equation_objective_metadata

    subroutine initialize_equation_objective_metadata( &
            metadata, name, kind, target, lower_bound, upper_bound, weight, &
            scale, active, fixed, kkt_id, nullspace_id, status, units, &
            provenance, has_parameter_tangent, parameter_tangent_id, has_hvp, &
            hvp_id)
        type(equation_objective_metadata_t), intent(out) :: metadata
        character(len=*), intent(in) :: name
        integer, intent(in) :: kind
        real(dp), intent(in) :: target(:), lower_bound(:), upper_bound(:)
        real(dp), intent(in) :: weight(:), scale(:)
        logical, intent(in) :: active, fixed
        integer, intent(in) :: kkt_id, nullspace_id
        type(fortsparse_status_t), intent(out) :: status
        character(len=*), intent(in), optional :: units, provenance
        logical, intent(in), optional :: has_parameter_tangent, has_hvp
        integer, intent(in), optional :: parameter_tangent_id, hvp_id

        call clear_equation_objective_metadata(metadata)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "equation/objective metadata received invalid values")
        if (len_trim(name) < 1) return
        if (size(target) < 1) return
        if (size(lower_bound) /= size(target) .or. &
            size(upper_bound) /= size(target) .or. &
            size(weight) /= size(target) .or. size(scale) /= size(target)) return

        metadata%name = trim(name)
        metadata%units = "dimensionless"
        metadata%provenance = "caller"
        if (present(units)) metadata%units = trim(units)
        if (present(provenance)) metadata%provenance = trim(provenance)
        metadata%kind = kind
        metadata%value_count = size(target)
        metadata%target = target
        metadata%lower_bound = lower_bound
        metadata%upper_bound = upper_bound
        metadata%weight = weight
        metadata%scale = scale
        metadata%active = active
        metadata%fixed = fixed
        metadata%kkt_id = kkt_id
        metadata%nullspace_id = nullspace_id
        if (present(has_parameter_tangent)) then
            metadata%has_parameter_tangent = has_parameter_tangent
        end if
        if (present(parameter_tangent_id)) then
            metadata%parameter_tangent_id = parameter_tangent_id
        end if
        if (present(has_hvp)) metadata%has_hvp = has_hvp
        if (present(hvp_id)) metadata%hvp_id = hvp_id
        if (.not. validate_equation_objective_metadata(metadata, status)) return
    end subroutine initialize_equation_objective_metadata

    logical function validate_equation_objective_metadata(metadata, status) &
            result(valid)
        type(equation_objective_metadata_t), intent(in) :: metadata
        type(fortsparse_status_t), intent(out) :: status
        integer :: value

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "equation/objective metadata is invalid")
        if (.not. valid_kind(metadata%kind)) return
        if (.not. allocated(metadata%name) .or. &
            .not. allocated(metadata%units) .or. &
            .not. allocated(metadata%provenance)) return
        if (len_trim(metadata%name) < 1 .or. len_trim(metadata%units) < 1 .or. &
            len_trim(metadata%provenance) < 1) return
        if (metadata%value_count < 1) return
        if (.not. allocated(metadata%target) .or. &
            .not. allocated(metadata%lower_bound) .or. &
            .not. allocated(metadata%upper_bound) .or. &
            .not. allocated(metadata%weight) .or. &
            .not. allocated(metadata%scale)) return
        if (size(metadata%target) /= metadata%value_count .or. &
            size(metadata%lower_bound) /= metadata%value_count .or. &
            size(metadata%upper_bound) /= metadata%value_count .or. &
            size(metadata%weight) /= metadata%value_count .or. &
            size(metadata%scale) /= metadata%value_count) return
        if (metadata%kkt_id < OBJECTIVE_METADATA_UNSET_ID .or. &
            metadata%nullspace_id < OBJECTIVE_METADATA_UNSET_ID) return
        if (metadata%parameter_tangent_id < OBJECTIVE_METADATA_UNSET_ID .or. &
            metadata%hvp_id < OBJECTIVE_METADATA_UNSET_ID) return
        if ((metadata%has_parameter_tangent .neqv. &
            metadata%parameter_tangent_id > OBJECTIVE_METADATA_UNSET_ID) .or. &
            (metadata%has_hvp .neqv. metadata%hvp_id > &
                OBJECTIVE_METADATA_UNSET_ID)) return
        do value = 1, metadata%value_count
            if (.not. ieee_is_finite(metadata%target(value)) .or. &
                .not. ieee_is_finite(metadata%lower_bound(value)) .or. &
                .not. ieee_is_finite(metadata%upper_bound(value)) .or. &
                .not. ieee_is_finite(metadata%weight(value)) .or. &
                .not. ieee_is_finite(metadata%scale(value))) return
            if (metadata%lower_bound(value) > metadata%upper_bound(value)) return
            if (metadata%target(value) < metadata%lower_bound(value) .or. &
                metadata%target(value) > metadata%upper_bound(value)) return
            if (metadata%weight(value) <= 0.0_dp .or. &
                metadata%scale(value) <= 0.0_dp) return
        end do
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_equation_objective_metadata

    logical function valid_kind(kind) result(valid)
        integer, intent(in) :: kind

        valid = kind == OBJECTIVE_METADATA_KIND_EQUATION .or. &
            kind == OBJECTIVE_METADATA_KIND_OBJECTIVE .or. &
            kind == OBJECTIVE_METADATA_KIND_CONSTRAINT
    end function valid_kind

end module fortfem_equation_objective_metadata
