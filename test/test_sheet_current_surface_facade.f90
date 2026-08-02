program test_sheet_current_surface_facade
    !! The surface-current representation is reachable without fortfem_api.
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        compare_sheet_current_surface_representations
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    real(dp) :: plus_field(1, 3), minus_field(1, 3), normals(1, 3)
    real(dp) :: surface_weights(1), signed_distance(3), normal_weights(3)
    real(dp) :: fitted(3), regularized(3), relative_error
    type(fortsparse_status_t) :: status

    plus_field = reshape([1.0_dp, 0.0_dp, 0.0_dp], shape(plus_field))
    minus_field = 0.0_dp
    normals = reshape([0.0_dp, 0.0_dp, 1.0_dp], shape(normals))
    surface_weights = 1.0_dp
    signed_distance = [-1.0_dp, 0.0_dp, 1.0_dp]
    normal_weights = [0.25_dp, 0.5_dp, 0.25_dp]

    call compare_sheet_current_surface_representations( &
        plus_field, minus_field, normals, surface_weights, signed_distance, &
        normal_weights, 0.5_dp, fitted, regularized, relative_error, status)
    call check_condition(status%code == 0, &
        "surface-current representation is exported by interop facade")
    call check_condition(all(fitted == fitted) .and. all(regularized == regularized), &
        "facade returns finite fitted and regularized ledgers")
    call check_summary("test_sheet_current_surface_facade")
end program test_sheet_current_surface_facade
