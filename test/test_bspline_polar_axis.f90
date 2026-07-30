program test_bspline_polar_axis
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_bspline_polar_h1_operator_csc, &
        build_bspline_polar_feec_2d_operators, &
        build_bspline_polar_feec_2d_extractions, &
        build_bspline_polar_h1_extraction, &
        evaluate_periodic_bspline_basis, &
        restrict_bspline_polar_operator_csc
    use fortfem_kinds, only: dp
    use fortsparse, only: &
        csc_from_triplet, csc_matvec, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: azimuth_count = 8
    integer, parameter :: radial_count = 5
    real(dp), parameter :: radial_knots(7) = [ &
        0.0_dp, 0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp, 1.0_dp]
    real(dp), allocatable :: extraction(:, :), polar_points(:, :)
    real(dp), allocatable :: curl(:, :), gradient(:, :), random_state(:)
    real(dp), allocatable :: hcurl_extraction(:, :), l2_extraction(:, :)
    real(dp), allocatable :: periodic_derivatives(:), periodic_values(:)
    real(dp), allocatable :: periodic_left(:), periodic_right(:)
    real(dp), allocatable :: physical_control_points(:, :, :)
    real(dp), allocatable :: physical_weights(:, :)
    real(dp), allocatable :: expected(:), polar_state(:)
    real(dp), allocatable :: tensor_curl(:, :), tensor_gradient(:, :)
    type(csc_t) :: physical_mass, physical_stiffness
    type(csc_t) :: polar_mass, tensor_identity
    type(fortsparse_status_t) :: sparse_status
    real(dp), allocatable :: tensor_points(:, :)
    real(dp) :: angle, polygon_area, radius(azimuth_count)
    integer :: azimuth, radial, status

    call build_bspline_polar_h1_extraction( &
        azimuth_count, radial_count, extraction, status)
    call check_condition(status == 0, &
        "Type-1 polar H1 extraction accepts a periodic radial patch")
    call check_condition(all(shape(extraction) == [ &
        azimuth_count*(radial_count - 2) + 3, &
        azimuth_count*radial_count]), &
        "Polar H1 dimension matches Toshniwal-Hughes equation 74")
    call check_condition(minval(extraction) >= -epsilon(1.0_dp) .and. &
        maxval(extraction) <= 1.0_dp + epsilon(1.0_dp), &
        "Polar extraction is nonnegative")
    call check_condition(maxval(abs(sum(extraction, dim=1) - 1.0_dp)) < &
        2.0e-15_dp, "Polar extraction is column stochastic")
    call check_condition(gram_is_positive_definite(extraction), &
        "Polar extraction has full row rank")

    allocate(polar_points(2, size(extraction, 1)))
    polar_points = 0.0_dp
    polar_points(:, 1) = [0.0_dp, 2.0_dp/sqrt(3.0_dp)]
    polar_points(:, 2) = [-1.0_dp, -1.0_dp/sqrt(3.0_dp)]
    polar_points(:, 3) = [1.0_dp, -1.0_dp/sqrt(3.0_dp)]
    tensor_points = matmul(polar_points, extraction)
    call check_condition(maxval(abs( &
        tensor_points(:, :azimuth_count) - &
        spread(tensor_points(:, 1), 2, azimuth_count))) < 2.0e-15_dp, &
        "Polar extraction collapses the complete periodic edge to one point")
    do azimuth = 1, azimuth_count
        radius(azimuth) = norm2( &
            tensor_points(:, azimuth_count + azimuth) - &
            tensor_points(:, azimuth))
    end do
    call check_condition(minval(radius) > 0.0_dp .and. &
        maxval(abs(radius - radius(1))) < 2.0e-15_dp, &
        "First polar control ring is a nondegenerate uniform fan")
    angle = signed_area( &
        tensor_points(:, azimuth_count + 1), &
        tensor_points(:, azimuth_count + 2), tensor_points(:, 1))
    call check_condition(abs(angle) > 1.0e-2_dp, &
        "Polar control fan has a well-defined orientation")

    call evaluate_periodic_bspline_basis( &
        azimuth_count, 3, 0.0_dp, periodic_values, periodic_derivatives, status)
    call evaluate_periodic_bspline_basis( &
        azimuth_count, 3, 1.0_dp, periodic_right, periodic_left, status)
    call check_condition(maxval(abs(periodic_values - periodic_right)) < &
        2.0e-15_dp, "Periodic B-spline basis closes exactly at the seam")
    call check_condition(abs(sum(periodic_values) - 1.0_dp) < 2.0e-15_dp .and. &
        abs(sum(periodic_derivatives)) < 2.0e-15_dp, &
        "Periodic B-spline basis preserves partition of unity")
    call evaluate_periodic_bspline_basis( &
        azimuth_count, 3, 0.371_dp - 1.0e-7_dp, periodic_left, &
        periodic_derivatives, status)
    call evaluate_periodic_bspline_basis( &
        azimuth_count, 3, 0.371_dp + 1.0e-7_dp, periodic_right, &
        periodic_derivatives, status)
    call evaluate_periodic_bspline_basis( &
        azimuth_count, 3, 0.371_dp, periodic_values, periodic_derivatives, &
        status)
    call check_condition(maxval(abs( &
        (periodic_right - periodic_left)/2.0e-7_dp - &
        periodic_derivatives)) < 2.0e-8_dp, &
        "Periodic B-spline derivatives match a central difference")

    call build_bspline_polar_feec_2d_operators( &
        azimuth_count, radial_count, gradient, curl, status)
    call check_condition(status == 0 .and. &
        all(shape(gradient) == [ &
        2 + 2*azimuth_count*(radial_count - 2), &
        3 + azimuth_count*(radial_count - 2)]) .and. &
        all(shape(curl) == [ &
        azimuth_count*(radial_count - 2), &
        2 + 2*azimuth_count*(radial_count - 2)]), &
        "Polar form dimensions match Toshniwal-Hughes equations 74--76")
    call check_condition(maxval(abs(matmul(curl, gradient))) < 2.0e-15_dp, &
        "Polar differential forms preserve curl(grad)=0")
    call check_condition(maxval(abs(matmul( &
        gradient, [(1.0_dp, azimuth = 1, size(gradient, 2))]))) < &
        2.0e-15_dp, "Polar gradient annihilates constants")
    call check_condition(gram_is_positive_definite(curl), &
        "Polar curl has full row rank")
    call check_condition(gram_is_positive_definite(transpose( &
        gradient(:, 2:))), "Polar gradient has constants as its only kernel")

    call build_bspline_polar_feec_2d_extractions( &
        azimuth_count, radial_count, extraction, hcurl_extraction, &
        l2_extraction, status)
    call build_periodic_tensor_complex( &
        azimuth_count, radial_count, tensor_gradient, tensor_curl)
    call check_condition(maxval(abs( &
        matmul(tensor_gradient, transpose(extraction)) - &
        matmul(transpose(hcurl_extraction), gradient))) < 2.0e-15_dp, &
        "Polar H1 and Hcurl extraction commutes with the gradient")
    call check_condition(maxval(abs( &
        matmul(tensor_curl, transpose(hcurl_extraction)) - &
        matmul(transpose(l2_extraction), curl))) < 2.0e-15_dp, &
        "Polar Hcurl and L2 extraction commutes with the curl")

    call build_sparse_identity(size(extraction, 2), tensor_identity)
    call restrict_bspline_polar_operator_csc( &
        extraction, tensor_identity, polar_mass, sparse_status)
    polar_state = [(sin(real(azimuth, dp)), &
        azimuth = 1, size(extraction, 1))]
    expected = matmul( &
        extraction, matmul(transpose(extraction), polar_state))
    call check_condition(sparse_status%code == 0 .and. maxval(abs( &
        csc_matvec(polar_mass, polar_state) - expected)) < 2.0e-14_dp, &
        "Sparse polar Hodge restriction equals E M E-transpose")

    allocate( &
        physical_control_points(2, azimuth_count, radial_count), &
        physical_weights(azimuth_count, radial_count))
    physical_weights = 1.0_dp
    do radial = 1, radial_count
        do azimuth = 1, azimuth_count
            angle = 2.0_dp*acos(-1.0_dp)* &
                real(azimuth, dp)/real(azimuth_count, dp)
            physical_control_points(:, azimuth, radial) = &
                real(radial - 1, dp)/real(radial_count - 1, dp)* &
                [cos(angle), sin(angle)]
        end do
    end do
    call assemble_bspline_polar_h1_operator_csc( &
        radial_knots, 1, azimuth_count, 1, physical_control_points, &
        physical_weights, 3, physical_mass, sparse_status, &
        stiffness_coefficient=0.0_dp, mass_coefficient=1.0_dp)
    polar_state = 1.0_dp
    expected = csc_matvec(physical_mass, polar_state)
    polygon_area = 0.5_dp*real(azimuth_count, dp)* &
        sin(2.0_dp*acos(-1.0_dp)/real(azimuth_count, dp))
    call check_condition(sparse_status%code == 0 .and. abs( &
        dot_product(polar_state, expected) - polygon_area) < 2.0e-13_dp, &
        "Physical polar mass integrates the exact regular-polygon area")
    call assemble_bspline_polar_h1_operator_csc( &
        radial_knots, 1, azimuth_count, 1, physical_control_points, &
        physical_weights, 3, physical_stiffness, sparse_status, &
        stiffness_coefficient=1.0_dp, mass_coefficient=0.0_dp)
    call check_condition(maxval(abs( &
        csc_matvec(physical_stiffness, polar_state))) < 2.0e-13_dp, &
        "Physical polar stiffness annihilates constants")

    allocate(random_state(size(gradient, 2)))
    call seed_random_numbers()
    do azimuth = 1, 20
        call random_number(random_state)
        call check_condition(maxval(abs( &
            matmul(curl, matmul(gradient, random_state)))) < 5.0e-15_dp, &
            "Random polar exact-sequence trial preserves curl(grad)=0")
    end do

    call check_summary("Isogeometric magnetic-axis polar extraction")

contains

    logical function gram_is_positive_definite(matrix) result(positive)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: factor(size(matrix, 1), size(matrix, 1))
        real(dp) :: gram(size(matrix, 1), size(matrix, 1)), pivot
        integer :: column, row

        gram = matmul(matrix, transpose(matrix))
        factor = 0.0_dp
        positive = .true.
        do column = 1, size(gram, 1)
            pivot = gram(column, column) - &
                dot_product(factor(column, :column - 1), &
                factor(column, :column - 1))
            if (pivot <= 1.0e3_dp*epsilon(1.0_dp)) then
                positive = .false.
                return
            end if
            factor(column, column) = sqrt(pivot)
            do row = column + 1, size(gram, 1)
                factor(row, column) = ( &
                    gram(row, column) - &
                    dot_product(factor(row, :column - 1), &
                    factor(column, :column - 1)))/factor(column, column)
            end do
        end do
    end function gram_is_positive_definite

    pure real(dp) function signed_area(first, second, origin) result(area)
        real(dp), intent(in) :: first(2), second(2), origin(2)

        area = (first(1) - origin(1))*(second(2) - origin(2)) - &
            (first(2) - origin(2))*(second(1) - origin(1))
    end function signed_area

    subroutine seed_random_numbers()
        integer, allocatable :: seed(:)
        integer :: entry

        call random_seed(size=entry)
        allocate(seed(entry))
        seed = [(104729 + 7919*entry, entry = 1, size(seed))]
        call random_seed(put=seed)
    end subroutine seed_random_numbers

    subroutine build_periodic_tensor_complex( &
            na, nr, tensor_gradient, tensor_curl)
        integer, intent(in) :: na, nr
        real(dp), allocatable, intent(out) :: tensor_gradient(:, :)
        real(dp), allocatable, intent(out) :: tensor_curl(:, :)
        integer :: angular, angular_offset, face, next_angular
        integer :: radial, radial_edge, top_angular, bottom_angular

        angular_offset = (nr - 1)*na
        allocate( &
            tensor_gradient(angular_offset + nr*na, nr*na), &
            tensor_curl((nr - 1)*na, angular_offset + nr*na))
        tensor_gradient = 0.0_dp
        tensor_curl = 0.0_dp
        do radial = 1, nr
            do angular = 1, na
                next_angular = modulo(angular, na) + 1
                bottom_angular = angular_offset + &
                    (radial - 1)*na + angular
                tensor_gradient(bottom_angular, &
                    (radial - 1)*na + angular) = -1.0_dp
                tensor_gradient(bottom_angular, &
                    (radial - 1)*na + next_angular) = 1.0_dp
                if (radial == nr) cycle
                radial_edge = (radial - 1)*na + angular
                tensor_gradient(radial_edge, &
                    (radial - 1)*na + angular) = -1.0_dp
                tensor_gradient(radial_edge, radial*na + angular) = 1.0_dp

                face = radial_edge
                top_angular = bottom_angular + na
                tensor_curl(face, radial_edge) = 1.0_dp
                tensor_curl(face, &
                    (radial - 1)*na + next_angular) = -1.0_dp
                tensor_curl(face, top_angular) = 1.0_dp
                tensor_curl(face, bottom_angular) = -1.0_dp
            end do
        end do
    end subroutine build_periodic_tensor_complex

    subroutine build_sparse_identity(dimension, matrix)
        integer, intent(in) :: dimension
        type(csc_t), intent(out) :: matrix
        integer :: indices(dimension), entry
        real(dp) :: values(dimension)
        type(fortsparse_status_t) :: identity_status

        indices = [(entry, entry = 1, dimension)]
        values = 1.0_dp
        call csc_from_triplet( &
            dimension, dimension, indices, indices, values, matrix, &
            identity_status)
        if (identity_status%code /= 0) &
            error stop "Could not construct polar oracle identity"
    end subroutine build_sparse_identity

end program test_bspline_polar_axis
