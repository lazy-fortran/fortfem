program test_generated_tetra_face_moment_transforms
    use check, only: check_condition, check_summary
    use fortfem_generated_tetra_face_moment_transforms, only: &
        map_tetra_face_basis_to_local, transform_tetra_face_moments
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    integer, parameter :: permutations(3, 6) = reshape([ &
        1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
    real(dp), allocatable :: canonical(:), local(:), transformed(:)
    integer :: order, permutation, status
    logical :: all_passed

    all_passed = .true.
    do order = 2, 4
        allocate( &
            canonical(order * (order - 1)), &
            local(order * (order - 1)), &
            transformed(order * (order - 1)))
        call integrate_face_moments( &
            order, permutations(:, 1), canonical, status)
        call record_condition(status == 0, &
            "Canonical face moments integrate successfully")
        do permutation = 1, size(permutations, 2)
            call integrate_face_moments( &
                order, permutations(:, permutation), local, status)
            call transform_tetra_face_moments( &
                order, permutations(:, permutation), local, transformed, &
                status)
            call record_condition(status == 0 .and. &
                maxval(abs(transformed - canonical)) < 3.0e-13_dp, &
                "Generated face transform matches direct moment integration")
            call map_tetra_face_basis_to_local( &
                order, permutations(:, permutation), canonical, transformed, &
                status)
            call record_condition(status == 0 .and. &
                maxval(abs(transformed - local)) < 3.0e-13_dp, &
                "Generated inverse maps canonical basis data to local data")
        end do
        deallocate(canonical, local, transformed)
    end do

    allocate(canonical(2), local(2), transformed(2))
    local = [1.0_dp, -2.0_dp]
    call transform_tetra_face_moments( &
        1, permutations(:, 1), local, transformed, status)
    call record_condition(status /= 0, &
        "Face transform rejects an order without face moments")
    call transform_tetra_face_moments( &
        2, [1, 1, 3], local, transformed, status)
    call record_condition(status /= 0, &
        "Face transform rejects an invalid vertex permutation")

    call check_summary("Generated tetrahedral face moment transforms")
    if (.not. all_passed) error stop 1

contains

    subroutine integrate_face_moments( &
            order, vertex_permutation, moments, status)
        integer, intent(in) :: order, vertex_permutation(3)
        real(dp), intent(out) :: moments(:)
        integer, intent(out) :: status

        real(dp), allocatable :: u(:), v(:), weights(:)
        real(dp) :: canonical_point(2), field(2), local_field(2)
        real(dp) :: affine(2, 2), offset(2)
        integer :: component, moment, node, total_degree, x_degree, y_degree

        moments = 0.0_dp
        status = 1
        if (size(moments) /= order * (order - 1)) return
        call permutation_map(vertex_permutation, offset, affine, status)
        if (status /= 0) return
        call triangle_duffy_quadrature(8, u, v, weights, status)
        if (status /= 0) return

        moment = 0
        do component = 1, 2
            do total_degree = 0, order - 2
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    moment = moment + 1
                    do node = 1, size(u)
                        canonical_point = offset + &
                            matmul(affine, [u(node), v(node)])
                        call polynomial_trace(canonical_point, field)
                        local_field = matmul(transpose(affine), field)
                        moments(moment) = moments(moment) + weights(node) * &
                            u(node)**x_degree * v(node)**y_degree * &
                            local_field(component)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine integrate_face_moments

    pure subroutine permutation_map( &
            vertex_permutation, offset, affine, status)
        integer, intent(in) :: vertex_permutation(3)
        real(dp), intent(out) :: offset(2), affine(2, 2)
        integer, intent(out) :: status

        real(dp), parameter :: vertices(2, 3) = reshape([ &
            0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])

        offset = 0.0_dp
        affine = 0.0_dp
        status = 1
        if (any(vertex_permutation < 1) .or. &
            any(vertex_permutation > 3)) return
        if (vertex_permutation(1) == vertex_permutation(2) .or. &
            vertex_permutation(1) == vertex_permutation(3) .or. &
            vertex_permutation(2) == vertex_permutation(3)) return
        offset = vertices(:, vertex_permutation(1))
        affine(:, 1) = vertices(:, vertex_permutation(2)) - offset
        affine(:, 2) = vertices(:, vertex_permutation(3)) - offset
        status = 0
    end subroutine permutation_map

    pure subroutine polynomial_trace(point, value)
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: value(2)

        value(1) = 2.0_dp - 3.0_dp * point(1) + point(2) + &
            4.0_dp * point(1) * point(2) - 2.0_dp * point(2)**2 + &
            point(1)**3
        value(2) = -1.0_dp + 2.0_dp * point(1) + 5.0_dp * point(2) - &
            3.0_dp * point(1)**2 + point(1) * point(2)**2
    end subroutine polynomial_trace

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_generated_tetra_face_moment_transforms
