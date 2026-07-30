program test_tetra_rt_arbitrary_order
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_tetra_rt, initialize_tetra_rt, &
        tetra_duffy_quadrature, tetra_rt_dof_count, tetra_rt_t
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    integer, parameter :: expected_counts(0:5) = [4, 15, 36, 70, 120, 189]
    type(tetra_rt_t) :: basis
    integer :: degree, status
    logical :: all_passed

    all_passed = .true.
    do degree = 0, 5
        call initialize_tetra_rt(degree, basis, status)
        call record_condition(status == 0 .and. &
            tetra_rt_dof_count(basis) == expected_counts(degree), &
            "Tetrahedral RT dimension matches the exact polynomial space")
        call check_moments(basis, degree)
        call check_divergence(basis)
    end do
    call initialize_tetra_rt(-1, basis, status)
    call record_condition(status /= 0, &
        "Tetrahedral RT rejects a negative degree")
    call initialize_tetra_rt(6, basis, status)
    call record_condition(status /= 0, &
        "Tetrahedral RT rejects an unsupported degree")

    call check_summary("Tetrahedral Raviart-Thomas arbitrary order")
    if (.not. all_passed) error stop 1

contains

    subroutine check_moments(selected_basis, degree)
        type(tetra_rt_t), intent(in) :: selected_basis
        integer, intent(in) :: degree

        integer :: basis_count, component, face, moment
        integer :: total, x_degree, y_degree, z_degree
        real(dp), allocatable :: computed_vector(:), divergences(:)
        real(dp), allocatable :: expected_vector(:), values(:, :)
        real(dp), allocatable :: x(:), y(:), z(:), weights(:)
        real(dp) :: face_error, interior_error, point(3)
        integer :: point_id, local_status

        basis_count = tetra_rt_dof_count(selected_basis)
        allocate(values(3, basis_count), divergences(basis_count))
        allocate(computed_vector(basis_count), expected_vector(basis_count))
        face_error = 0.0_dp
        interior_error = 0.0_dp
        moment = 0
        do face = 1, 4
            do total = 0, degree
                do x_degree = 0, total
                    y_degree = total - x_degree
                    moment = moment + 1
                    call integrate_face_moments( &
                        selected_basis, face, x_degree, y_degree, &
                        computed_vector)
                    expected_vector = 0.0_dp
                    expected_vector(moment) = 1.0_dp
                    face_error = max(face_error, &
                        maxval(abs(computed_vector - expected_vector)))
                end do
            end do
        end do

        call tetra_duffy_quadrature(2*degree + 3, &
            x, y, z, weights, local_status)
        if (local_status /= 0) error stop "RT volume quadrature failed"
        do component = 1, 3
            do total = 0, degree - 1
                do x_degree = 0, total
                    do y_degree = 0, total - x_degree
                        z_degree = total - x_degree - y_degree
                        moment = moment + 1
                        computed_vector = 0.0_dp
                        do point_id = 1, size(weights)
                            point = [x(point_id), y(point_id), z(point_id)]
                            call evaluate_tetra_rt( &
                                selected_basis, point, values, &
                                divergences, local_status)
                            if (local_status /= 0) then
                                error stop "RT evaluation failed"
                            end if
                            computed_vector = computed_vector + &
                                weights(point_id)*values(component, :)* &
                                x(point_id)**x_degree* &
                                y(point_id)**y_degree* &
                                z(point_id)**z_degree
                        end do
                        expected_vector = 0.0_dp
                        expected_vector(moment) = 1.0_dp
                        interior_error = max(interior_error, &
                            maxval(abs(computed_vector - expected_vector)))
                    end do
                end do
            end do
        end do
        if (face_error >= moment_tolerance(degree)) then
            write (*, '(A,I0,A,ES12.4)') &
                "RT degree ", degree, " face error ", face_error
        end if
        if (interior_error >= moment_tolerance(degree)) then
            write (*, '(A,I0,A,ES12.4)') &
                "RT degree ", degree, " interior error ", interior_error
        end if
        call record_condition(face_error < moment_tolerance(degree), &
            "Tetrahedral RT face moment matrix is identity")
        call record_condition(interior_error < moment_tolerance(degree), &
            "Tetrahedral RT interior moment matrix is identity")
        call record_condition(moment == basis_count, &
            "Tetrahedral RT moments have the expected total dimension")
    end subroutine check_moments

    subroutine integrate_face_moments( &
            selected_basis, face, x_power, y_power, integrals)
        type(tetra_rt_t), intent(in) :: selected_basis
        integer, intent(in) :: face, x_power, y_power
        real(dp), intent(out) :: integrals(:)

        integer, parameter :: quadrature_order = 7
        real(dp) :: nodes(quadrature_order), weights(quadrature_order)
        real(dp) :: area_normal(3), affine(3, 2), offset(3), point(3)
        real(dp) :: s, t
        real(dp), allocatable :: divergences(:), values(:, :)
        integer :: first, local_status, second

        allocate(values(3, tetra_rt_dof_count(selected_basis)))
        allocate(divergences(tetra_rt_dof_count(selected_basis)))
        call face_geometry(face, offset, affine, area_normal)
        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)
        integrals = 0.0_dp
        do first = 1, quadrature_order
            s = nodes(first)
            do second = 1, quadrature_order
                t = (1.0_dp - s)*nodes(second)
                point = offset + affine(:, 1)*s + affine(:, 2)*t
                call evaluate_tetra_rt( &
                    selected_basis, point, values, divergences, local_status)
                if (local_status /= 0) error stop "RT face evaluation failed"
                integrals = integrals + weights(first)*weights(second)* &
                    (1.0_dp - s)*matmul(area_normal, values)* &
                    s**x_power*t**y_power
            end do
        end do
    end subroutine integrate_face_moments

    subroutine check_divergence(selected_basis)
        type(tetra_rt_t), intent(in) :: selected_basis

        real(dp), parameter :: step = 1.0e-3_dp
        real(dp) :: point(3), shifted(3)
        real(dp), allocatable :: divergences(:), finite_difference(:)
        real(dp), allocatable :: minus_one(:, :), minus_two(:, :)
        real(dp), allocatable :: plus_one(:, :), plus_two(:, :)
        real(dp) :: divergence_error, divergence_scale
        integer :: component, local_status

        allocate(divergences(tetra_rt_dof_count(selected_basis)))
        allocate(finite_difference(tetra_rt_dof_count(selected_basis)))
        allocate(minus_one(3, tetra_rt_dof_count(selected_basis)))
        allocate(minus_two(3, tetra_rt_dof_count(selected_basis)))
        allocate(plus_one(3, tetra_rt_dof_count(selected_basis)))
        allocate(plus_two(3, tetra_rt_dof_count(selected_basis)))
        point = [0.19_dp, 0.23_dp, 0.17_dp]
        finite_difference = 0.0_dp
        do component = 1, 3
            shifted = point
            shifted(component) = shifted(component) + step
            call evaluate_tetra_rt( &
                selected_basis, shifted, plus_one, divergences, &
                local_status)
            shifted = point
            shifted(component) = shifted(component) + 2.0_dp*step
            call evaluate_tetra_rt( &
                selected_basis, shifted, plus_two, divergences, &
                local_status)
            shifted = point
            shifted(component) = shifted(component) - step
            call evaluate_tetra_rt( &
                selected_basis, shifted, minus_one, divergences, &
                local_status)
            shifted = point
            shifted(component) = shifted(component) - 2.0_dp*step
            call evaluate_tetra_rt( &
                selected_basis, shifted, minus_two, divergences, &
                local_status)
            finite_difference = finite_difference + &
                (-plus_two(component, :) + 8.0_dp*plus_one(component, :) - &
                8.0_dp*minus_one(component, :) + minus_two(component, :))/ &
                (12.0_dp*step)
        end do
        call evaluate_tetra_rt( &
            selected_basis, point, plus_one, divergences, local_status)
        divergence_scale = max(1.0_dp, maxval(abs(divergences)))
        divergence_error = &
            maxval(abs(divergences - finite_difference))/divergence_scale
        if (divergence_error >= 2.0e-9_dp) then
            write (*, '(A,I0,A,ES12.4)') "RT degree ", &
                selected_basis%degree, " divergence error ", divergence_error
        end if
        call record_condition(divergence_error < 2.0e-9_dp, &
            "Tetrahedral RT divergence matches finite differences")
    end subroutine check_divergence

    pure real(dp) function moment_tolerance(degree) result(tolerance)
        integer, intent(in) :: degree

        if (degree <= 4) then
            tolerance = 3.0e-10_dp
        else
            tolerance = 2.0e-7_dp
        end if
    end function moment_tolerance

    pure subroutine face_geometry(face, offset, affine, area_normal)
        integer, intent(in) :: face
        real(dp), intent(out) :: offset(3), affine(3, 2), area_normal(3)

        select case (face)
        case (1)
            offset = [0.0_dp, 0.0_dp, 0.0_dp]
            affine(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
            affine(:, 2) = [0.0_dp, 0.0_dp, 1.0_dp]
            area_normal = [-1.0_dp, 0.0_dp, 0.0_dp]
        case (2)
            offset = [0.0_dp, 0.0_dp, 0.0_dp]
            affine(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
            affine(:, 2) = [0.0_dp, 0.0_dp, 1.0_dp]
            area_normal = [0.0_dp, -1.0_dp, 0.0_dp]
        case (3)
            offset = [0.0_dp, 0.0_dp, 0.0_dp]
            affine(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
            affine(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
            area_normal = [0.0_dp, 0.0_dp, -1.0_dp]
        case default
            offset = [1.0_dp, 0.0_dp, 0.0_dp]
            affine(:, 1) = [-1.0_dp, 1.0_dp, 0.0_dp]
            affine(:, 2) = [-1.0_dp, 0.0_dp, 1.0_dp]
            area_normal = [1.0_dp, 1.0_dp, 1.0_dp]
        end select
    end subroutine face_geometry

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition
end program test_tetra_rt_arbitrary_order
