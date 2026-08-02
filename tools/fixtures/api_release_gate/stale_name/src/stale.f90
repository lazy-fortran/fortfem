module stale_api_fixture
    implicit none
contains
    subroutine stale_call()
        call evaluate_boundary_operator_parity()
        call evaluate_sheet_current_parity()
        call evaluate_tree_cotree_iga_parity()
        call evaluate_beltrami_two_region_parity()
        call evaluate_beltrami_shell_parity()
    end subroutine stale_call
end module stale_api_fixture
