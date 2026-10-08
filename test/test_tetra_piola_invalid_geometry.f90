program test_tetra_piola_invalid_geometry
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant, &
        map_tetra_nedelec_covariant_jvp, map_tetra_nedelec_covariant_vjp
    implicit none
    real(dp) :: jacobian(3, 3), seed(3, 3), values(3, 1), curls(3, 1)
    real(dp) :: mapped_values(3, 1), mapped_curls(3, 1), jacobian_bar(3, 3)
    integer :: geometry, status

    values = 1.0_dp
    curls = 2.0_dp
    seed = 0.1_dp
    do geometry = 1, 3
        jacobian = 0.0_dp
        jacobian(1, 1) = 1.0_dp
        jacobian(2, 2) = 1.0_dp
        select case (geometry)
        case (1)
            jacobian(3, 3) = 0.0_dp ! Collapsed tetrahedron.
        case (2)
            jacobian(3, 3) = -1.0_dp ! Reversed orientation.
        case (3)
            jacobian(3, 3) = epsilon(1.0_dp) ! Below geometry tolerance.
        end select
        call map_tetra_nedelec_covariant(jacobian, values, curls, &
            mapped_values, mapped_curls, status)
        call check_condition(status /= 0, "Primal rejects invalid geometry")
        call check_condition(maxval(abs(mapped_values)) == 0.0_dp .and. &
            maxval(abs(mapped_curls)) == 0.0_dp, "Primal clears invalid output")
        call map_tetra_nedelec_covariant_jvp(jacobian, values, curls, &
            seed, values, curls, mapped_values, mapped_curls, status)
        call check_condition(status /= 0, "JVP rejects invalid geometry")
        call check_condition(maxval(abs(mapped_values)) == 0.0_dp .and. &
            maxval(abs(mapped_curls)) == 0.0_dp, "JVP clears invalid output")
        call map_tetra_nedelec_covariant_vjp(jacobian, values, curls, &
            values, curls, jacobian_bar, mapped_values, mapped_curls, status)
        call check_condition(status /= 0, "VJP rejects invalid geometry")
        call check_condition(maxval(abs(jacobian_bar)) == 0.0_dp .and. &
            maxval(abs(mapped_values)) == 0.0_dp .and. &
            maxval(abs(mapped_curls)) == 0.0_dp, "VJP clears invalid output")
    end do
    call check_summary("Tetrahedral Piola invalid geometry")
end program test_tetra_piola_invalid_geometry
