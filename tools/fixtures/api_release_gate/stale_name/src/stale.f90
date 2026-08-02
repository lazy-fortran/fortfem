module stale_api_fixture
    implicit none
contains
    subroutine stale_call()
        call evaluate_boundary_operator_parity()
        call evaluate_sheet_current_parity()
    end subroutine stale_call
end module stale_api_fixture
