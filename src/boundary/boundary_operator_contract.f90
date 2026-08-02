module fortfem_boundary_operator_contract
    !! Typed metadata boundary for interchangeable exterior operators.
    !!
    !! Numerical actions remain in the FEM, BEM, DtN, PML, Fourier, NESTOR,
    !! BIEST, or virtual-casing client.  This record makes dimensions, equation
    !! space, available derivative/residual actions, units, topology identity,
    !! and provenance inspectable without importing an application format.
    implicit none
    private

    integer, parameter, public :: BOUNDARY_OPERATOR_BACKEND_FEM = 1
    integer, parameter, public :: BOUNDARY_OPERATOR_BACKEND_BEM = 2
    integer, parameter, public :: BOUNDARY_OPERATOR_BACKEND_DTN = 3
    integer, parameter, public :: BOUNDARY_OPERATOR_BACKEND_PML = 4
    integer, parameter, public :: BOUNDARY_OPERATOR_BACKEND_NESTOR = 5
    integer, parameter, public :: BOUNDARY_OPERATOR_BACKEND_BIEST = 6
    integer, parameter, public :: BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING = 7
    integer, parameter, public :: BOUNDARY_OPERATOR_BACKEND_USER = 8

    type, public :: boundary_operator_contract_t
        character(len=32) :: schema_version = "fortfem-boundary-operator-1"
        integer :: backend_kind = 0
        character(len=32) :: backend_name = ""
        character(len=64) :: equation = ""
        character(len=64) :: space = ""
        integer :: row_count = 0
        integer :: column_count = 0
        logical :: complex_valued = .false.
        logical :: matrix_free = .false.
        logical :: assembled = .false.
        logical :: has_jvp = .false.
        logical :: has_vjp = .false.
        logical :: has_residual = .false.
        character(len=32) :: units = "1"
        character(len=64) :: normalization = "1"
        character(len=128) :: provenance = ""
        character(len=64) :: topology_id = ""
    end type boundary_operator_contract_t

    public :: initialize_boundary_operator_contract
    public :: validate_boundary_operator_contract

contains

    subroutine initialize_boundary_operator_contract( &
            contract, backend_kind, equation, space, row_count, column_count, &
            complex_valued, matrix_free, assembled, has_jvp, has_vjp, &
            has_residual, units, normalization, provenance, topology_id, status)
        type(boundary_operator_contract_t), intent(out) :: contract
        integer, intent(in) :: backend_kind, row_count, column_count
        character(len=*), intent(in) :: equation, space, units, normalization
        character(len=*), intent(in) :: provenance, topology_id
        logical, intent(in) :: complex_valued, matrix_free, assembled
        logical, intent(in) :: has_jvp, has_vjp, has_residual
        integer, intent(out) :: status

        contract%backend_kind = backend_kind
        contract%backend_name = backend_name_for_kind(backend_kind)
        contract%equation = ""
        contract%space = ""
        contract%row_count = row_count
        contract%column_count = column_count
        contract%complex_valued = complex_valued
        contract%matrix_free = matrix_free
        contract%assembled = assembled
        contract%has_jvp = has_jvp
        contract%has_vjp = has_vjp
        contract%has_residual = has_residual
        contract%units = ""
        contract%normalization = ""
        contract%provenance = ""
        contract%topology_id = ""
        if (len_trim(equation) <= len(contract%equation)) contract%equation = equation
        if (len_trim(space) <= len(contract%space)) contract%space = space
        if (len_trim(units) <= len(contract%units)) contract%units = units
        if (len_trim(normalization) <= len(contract%normalization)) &
            contract%normalization = normalization
        if (len_trim(provenance) <= len(contract%provenance)) contract%provenance = provenance
        if (len_trim(topology_id) <= len(contract%topology_id)) contract%topology_id = topology_id
        if (len_trim(equation) > len(contract%equation) .or. &
            len_trim(space) > len(contract%space) .or. len_trim(units) > len(contract%units) .or. &
            len_trim(normalization) > len(contract%normalization) .or. &
            len_trim(provenance) > len(contract%provenance) .or. &
            len_trim(topology_id) > len(contract%topology_id)) then
            status = 1
            return
        end if
        if (.not. validate_boundary_operator_contract(contract, status)) then
            contract = boundary_operator_contract_t()
            return
        end if
    end subroutine initialize_boundary_operator_contract

    logical function validate_boundary_operator_contract(contract, status) result(valid)
        type(boundary_operator_contract_t), intent(in) :: contract
        integer, intent(out) :: status

        valid = .false.
        status = 1
        if (contract%schema_version /= "fortfem-boundary-operator-1") return
        if (contract%backend_kind < BOUNDARY_OPERATOR_BACKEND_FEM .or. &
            contract%backend_kind > BOUNDARY_OPERATOR_BACKEND_USER) return
        if (len_trim(contract%backend_name) == 0 .or. &
            len_trim(contract%equation) == 0 .or. len_trim(contract%space) == 0 .or. &
            len_trim(contract%units) == 0 .or. len_trim(contract%normalization) == 0 .or. &
            len_trim(contract%provenance) == 0 .or. len_trim(contract%topology_id) == 0) return
        if (contract%row_count < 1 .or. contract%column_count < 1) return
        if (.not. contract%matrix_free .and. .not. contract%assembled) return
        if (.not. contract%has_jvp .or. .not. contract%has_vjp .or. &
            .not. contract%has_residual) return
        status = 0
        valid = .true.
    end function validate_boundary_operator_contract

    character(len=32) function backend_name_for_kind(backend_kind) result(name)
        integer, intent(in) :: backend_kind

        name = ""
        select case (backend_kind)
        case (BOUNDARY_OPERATOR_BACKEND_FEM); name = "FEM"
        case (BOUNDARY_OPERATOR_BACKEND_BEM); name = "BEM"
        case (BOUNDARY_OPERATOR_BACKEND_DTN); name = "DtN"
        case (BOUNDARY_OPERATOR_BACKEND_PML); name = "PML"
        case (BOUNDARY_OPERATOR_BACKEND_NESTOR); name = "NESTOR"
        case (BOUNDARY_OPERATOR_BACKEND_BIEST); name = "BIEST"
        case (BOUNDARY_OPERATOR_BACKEND_VIRTUAL_CASING); name = "virtual-casing"
        case (BOUNDARY_OPERATOR_BACKEND_USER); name = "user"
        end select
    end function backend_name_for_kind

end module fortfem_boundary_operator_contract
