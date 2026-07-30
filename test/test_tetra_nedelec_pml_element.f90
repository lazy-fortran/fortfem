program test_tetra_nedelec_pml_element
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_tetra_nedelec_curl_mass_element, &
        assemble_tetra_nedelec_pml_element
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wave_number = 1.7_dp
    real(dp), parameter :: vertices(3, 4) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.2_dp, 0.0_dp, 0.0_dp, &
        0.1_dp, 0.9_dp, 0.0_dp, &
        0.0_dp, 0.2_dp, 1.1_dp], [3, 4])
    complex(dp), allocatable :: pml_matrix(:, :)
    real(dp), allocatable :: reference_matrix(:, :)
    complex(dp) :: stretch(3)
    integer :: order, status
    logical :: all_passed

    all_passed = .true.
    stretch = cmplx(1.0_dp, 0.0_dp, dp)
    do order = 1, 4
        call assemble_tetra_nedelec_pml_element( &
            vertices, order, 10, stretch, wave_number, pml_matrix, status)
        if (status /= 0) error stop "Nedelec PML element assembly failed"
        call assemble_tetra_nedelec_curl_mass_element( &
            vertices, order, 10, reference_matrix, status, &
            curl_coefficient=1.0_dp, mass_coefficient=-wave_number**2)
        if (status /= 0) error stop "reference Nedelec assembly failed"
        call record_condition( &
            maxval(abs(pml_matrix - cmplx(reference_matrix, 0.0_dp, dp))) < &
            2.0e-12_dp, &
            "unit-stretch PML reduces to curl-curl minus mass")
        deallocate(pml_matrix, reference_matrix)
    end do

    stretch = [ &
        cmplx(1.0_dp, 0.8_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp)]
    call assemble_tetra_nedelec_pml_element( &
        vertices, 2, 10, stretch, wave_number, pml_matrix, status)
    call record_condition(status == 0, &
        "anisotropic complex stretch assembles for order two")
    call record_condition( &
        maxval(abs(pml_matrix - transpose(pml_matrix))) < 2.0e-12_dp, &
        "reciprocal Maxwell PML element is complex symmetric")
    call record_condition(maxval(abs(aimag(pml_matrix))) > 1.0e-3_dp, &
        "absorbing stretch produces a nonzero imaginary operator")

    call check_summary("Arbitrary-order tetrahedral Nedelec PML element")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_pml_element
