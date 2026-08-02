program test_field_aligned_anisotropic_3d_gallery
    !! Independent constitutive and manufactured-PDE oracle for the gallery.
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_field_aligned_constitutive_tensor
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: k_parallel = 5000.0_dp
    real(dp), parameter :: k_perpendicular = 1.0_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: direction(3) = [1.0_dp, 2.0_dp, 2.0_dp]/3.0_dp
    real(dp), parameter :: point(3) = [0.31_dp, 0.42_dp, 0.63_dp]
    real(dp) :: tensor(3, 3), expected(3, 3), hessian(3, 3)
    real(dp) :: source, independent_source, solution_value
    real(dp) :: direction_outer(3, 3)
    logical :: all_passed
    type(fortsparse_status_t) :: status

    all_passed = .true.
    direction_outer = spread(direction, 2, 3)*spread(direction, 1, 3)
    expected = k_perpendicular*identity_matrix() + &
        (k_parallel - k_perpendicular)*direction_outer
    call evaluate_field_aligned_constitutive_tensor( &
        k_parallel, k_perpendicular, direction, tensor, status)
    call record_condition(status%code == 0, &
        "field-aligned 3-D gallery tensor accepts the unit direction")
    call record_condition(maxval(abs(tensor - expected)) < 1.0e-13_dp, &
        "field-aligned tensor matches the independent projector oracle")
    call record_condition(abs(dot_product(direction, direction) - 1.0_dp) < &
        1.0e-14_dp, "gallery direction is normalized")

    solution_value = sin(pi*point(1))*sin(pi*point(2))*sin(pi*point(3))
    hessian = 0.0_dp
    hessian(1, 1) = -pi*pi*solution_value
    hessian(2, 2) = hessian(1, 1)
    hessian(3, 3) = hessian(1, 1)
    hessian(1, 2) = pi*pi*cos(pi*point(1))*cos(pi*point(2))* &
        sin(pi*point(3))
    hessian(2, 1) = hessian(1, 2)
    hessian(1, 3) = pi*pi*cos(pi*point(1))*sin(pi*point(2))* &
        cos(pi*point(3))
    hessian(3, 1) = hessian(1, 3)
    hessian(2, 3) = pi*pi*sin(pi*point(1))*cos(pi*point(2))* &
        cos(pi*point(3))
    hessian(3, 2) = hessian(2, 3)
    source = -sum(tensor*hessian)
    independent_source = -sum(expected*hessian)
    call record_condition(abs(source - independent_source) < 1.0e-10_dp, &
        "manufactured anisotropic source matches independent Hessian oracle")
    call record_condition(source > 0.0_dp, &
        "manufactured interior source is positive at the oracle point")
    call check_summary("field-aligned anisotropic 3-D gallery")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    pure function identity_matrix() result(matrix)
        real(dp) :: matrix(3, 3)

        matrix = 0.0_dp
        matrix(1, 1) = 1.0_dp
        matrix(2, 2) = 1.0_dp
        matrix(3, 3) = 1.0_dp
    end function identity_matrix

end program test_field_aligned_anisotropic_3d_gallery
