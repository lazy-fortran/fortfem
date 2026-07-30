program test_tetra_rt_face_moment_transforms
    use check, only: check_condition, check_summary
    use fortfem_tetra_face_moment_transforms, only: &
        map_tetra_rt_face_basis_to_local, transform_tetra_rt_face_moments
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    integer, parameter :: permutations(3, 6) = reshape([ &
        1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
    integer :: degree, permutation, status
    real(dp) :: invalid_output(3)
    logical :: all_passed

    all_passed = .true.
    do degree = 0, 5
        do permutation = 1, 6
            call check_transform(degree, permutations(:, permutation))
        end do
    end do
    call transform_tetra_rt_face_moments( &
        2, [1, 1, 3], [1.0_dp, 2.0_dp, 3.0_dp], &
        invalid_output, status)
    call record_condition(status /= 0, &
        "RT face moment transform rejects a repeated vertex")

    call check_summary("Tetrahedral RT face moment transforms")
    if (.not. all_passed) error stop 1

contains

    subroutine check_transform(degree, permutation)
        integer, intent(in) :: degree, permutation(3)

        integer :: dof_count, moment, status
        real(dp), allocatable :: canonical(:), local(:), recovered(:)
        real(dp), allocatable :: coefficients(:)

        dof_count = (degree + 1)*(degree + 2)/2
        allocate( &
            canonical(dof_count), local(dof_count), recovered(dof_count), &
            coefficients(dof_count))
        do moment = 1, dof_count
            coefficients(moment) = &
                real(3*moment - 2*degree - 5, dp)/real(moment + 3, dp)
        end do
        call numerical_moments( &
            degree, [1, 2, 3], coefficients, canonical)
        call numerical_moments(degree, permutation, coefficients, local)
        call transform_tetra_rt_face_moments( &
            degree, permutation, local, recovered, status)
        call record_condition(status == 0 .and. &
            maxval(abs(recovered - canonical)) < 2.0e-12_dp, &
            "FortSym/runtime RT transform matches numerical face moments")

        call map_tetra_rt_face_basis_to_local( &
            degree, permutation, canonical, recovered, status)
        call record_condition(status == 0 .and. &
            maxval(abs(recovered - local)) < 2.0e-12_dp, &
            "FortSym/runtime RT inverse transform recovers local moments")
    end subroutine check_transform

    subroutine numerical_moments( &
            degree, permutation, coefficients, moments)
        integer, intent(in) :: degree, permutation(3)
        real(dp), intent(in) :: coefficients(:)
        real(dp), intent(out) :: moments(:)

        integer, parameter :: quadrature_order = 8
        real(dp) :: affine(2, 2), nodes(quadrature_order)
        real(dp) :: offset(2), point(2), s, t
        real(dp) :: weights(quadrature_order)
        real(dp) :: polynomial_value
        integer :: first, moment, second, total, x_power, y_power

        call permutation_map(permutation, offset, affine)
        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)
        moments = 0.0_dp
        do first = 1, quadrature_order
            s = nodes(first)
            do second = 1, quadrature_order
                t = (1.0_dp - s)*nodes(second)
                point = offset + matmul(affine, [s, t])
                polynomial_value = evaluate_polynomial( &
                    degree, coefficients, point)
                moment = 0
                do total = 0, degree
                    do x_power = 0, total
                        y_power = total - x_power
                        moment = moment + 1
                        moments(moment) = moments(moment) + &
                            weights(first)*weights(second)*(1.0_dp - s)* &
                            polynomial_value*s**x_power*t**y_power
                    end do
                end do
            end do
        end do
    end subroutine numerical_moments

    pure function evaluate_polynomial(degree, coefficients, point) &
            result(value)
        integer, intent(in) :: degree
        real(dp), intent(in) :: coefficients(:), point(2)
        real(dp) :: value

        integer :: basis, total, x_power, y_power

        value = 0.0_dp
        basis = 0
        do total = 0, degree
            do x_power = 0, total
                y_power = total - x_power
                basis = basis + 1
                value = value + coefficients(basis)* &
                    point(1)**x_power*point(2)**y_power
            end do
        end do
    end function evaluate_polynomial

    pure subroutine permutation_map(permutation, offset, affine)
        integer, intent(in) :: permutation(3)
        real(dp), intent(out) :: offset(2), affine(2, 2)

        real(dp), parameter :: vertices(2, 3) = reshape([ &
            0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])

        offset = vertices(:, permutation(1))
        affine(:, 1) = vertices(:, permutation(2)) - offset
        affine(:, 2) = vertices(:, permutation(3)) - offset
    end subroutine permutation_map

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_rt_face_moment_transforms
