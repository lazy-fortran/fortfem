program test_sheet_current_parity
    !! Independent slab oracle for the explicit/resolved sheet-current bridge.
    use check, only: check_condition, check_summary
    use fortfem_api, only: compare_sheet_current_representations
    use fortfem_kinds, only: dp
    use fortfem_interop, only: &
        compare_sheet_current_representations_interop => &
        compare_sheet_current_representations
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: point_count = 401, component_count = 3
    real(dp), parameter :: pi = acos(-1.0_dp), thickness = 0.02_dp
    real(dp), parameter :: half_width = 5.0_dp*thickness
    real(dp), parameter :: surface_measure = 2.5_dp
    real(dp), parameter :: explicit_current(component_count) = &
        [0.0_dp, 1.25_dp, -0.5_dp]
    real(dp) :: distance(point_count), weights(point_count)
    real(dp) :: current(point_count, component_count)
    real(dp) :: regularized(component_count), explicit(component_count)
    real(dp) :: regularized_interop(component_count), explicit_interop(component_count)
    real(dp) :: relative_error, spacing, gaussian, reference(component_count)
    integer :: point
    type(fortsparse_status_t) :: status

    spacing = 2.0_dp*half_width/real(point_count - 1, dp)
    reference = 0.0_dp
    do point = 1, point_count
        distance(point) = -half_width + spacing*real(point - 1, dp)
        weights(point) = spacing
        if (point == 1 .or. point == point_count) weights(point) = 0.5_dp*spacing
        current(point, :) = explicit_current
        gaussian = exp(-(distance(point)/thickness)**2)/(sqrt(pi)*thickness)
        reference = reference + weights(point)*gaussian*surface_measure* &
            explicit_current
    end do
    call compare_sheet_current_representations( &
        distance, weights, current, thickness, surface_measure, explicit_current, &
        regularized, explicit, relative_error, status)
    call check_condition(status%code == 0, &
        "sheet-current parity accepts a resolved slab profile")
    call check_condition(maxval(abs(explicit - surface_measure*explicit_current)) < &
        1.0e-14_dp, "explicit surface ledger equals K times surface measure")
    call check_condition(maxval(abs(regularized - reference)) < 1.0e-14_dp, &
        "resolved ledger matches an independent Gaussian quadrature oracle")
    call check_condition(relative_error < 2.0e-12_dp, &
        "resolved layer approaches the explicit sheet ledger")
    call compare_sheet_current_representations_interop( &
        distance, weights, current, thickness, surface_measure, explicit_current, &
        regularized_interop, explicit_interop, relative_error, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(regularized_interop - regularized)) < 1.0e-14_dp .and. &
        maxval(abs(explicit_interop - explicit)) < 1.0e-14_dp, &
        "interoperability facade preserves both sheet ledgers")

    current(1, 1) = current(1, 1) + 1.0e-6_dp
    call compare_sheet_current_representations( &
        distance, weights, current, thickness, surface_measure, explicit_current, &
        regularized, explicit, relative_error, status)
    call check_condition(status%code /= 0, &
        "parity rejects a normal-varying current that is not one sheet")

    call check_summary("sheet-current parity")
end program test_sheet_current_parity
