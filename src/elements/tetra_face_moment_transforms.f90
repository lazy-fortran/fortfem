module fortfem_tetra_face_moment_transforms
    use fortfem_generated_tetra_face_moment_transforms, only: &
        generated_basis_to_local => map_tetra_face_basis_to_local, &
        generated_transform_moments => transform_tetra_face_moments, &
        generated_rt_basis_to_local => map_tetra_rt_face_basis_to_local, &
        generated_rt_transform_moments => transform_tetra_rt_face_moments
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: dense_solve
    use fortnum_special_jacobi, only: triangle_dubiner
    implicit none
    private

    public :: build_tetra_face_basis_to_local_matrix
    public :: build_tetra_rt_face_basis_to_local_matrix
    public :: map_tetra_face_basis_to_local
    public :: map_tetra_rt_face_basis_to_local
    public :: transform_tetra_face_moments
    public :: transform_tetra_rt_face_moments

contains

    subroutine build_tetra_face_basis_to_local_matrix( &
            order, permutation, transform, status)
        integer, intent(in) :: order, permutation(3)
        real(dp), intent(out) :: transform(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: canonical(:), local(:)
        integer :: column, size_

        transform = 0.0_dp
        size_ = order*(order - 1)
        status = 1
        if (size(transform, 1) /= size_ .or. size(transform, 2) /= size_) &
            return
        if (order > 4) then
            call build_runtime_transform_matrix( &
                order, permutation, .true., transform, status)
            return
        end if
        allocate(canonical(size_), local(size_))
        do column = 1, size_
            canonical = 0.0_dp
            canonical(column) = 1.0_dp
            call generated_basis_to_local( &
                order, permutation, canonical, local, status)
            if (status /= 0) return
            transform(:, column) = local
        end do
        status = 0
    end subroutine build_tetra_face_basis_to_local_matrix

    subroutine map_tetra_face_basis_to_local( &
            order, permutation, canonical, local, status)
        integer, intent(in) :: order, permutation(3)
        real(dp), intent(in) :: canonical(:)
        real(dp), intent(out) :: local(:)
        integer, intent(out) :: status

        if (order <= 4) then
            call generated_basis_to_local( &
                order, permutation, canonical, local, status)
        else
            call apply_runtime_transform( &
                order, permutation, canonical, local, .true., status)
        end if
    end subroutine map_tetra_face_basis_to_local

    subroutine transform_tetra_face_moments( &
            order, permutation, local, canonical, status)
        integer, intent(in) :: order, permutation(3)
        real(dp), intent(in) :: local(:)
        real(dp), intent(out) :: canonical(:)
        integer, intent(out) :: status

        if (order <= 4) then
            call generated_transform_moments( &
                order, permutation, local, canonical, status)
        else
            call apply_runtime_transform( &
                order, permutation, local, canonical, .false., status)
        end if
    end subroutine transform_tetra_face_moments

    subroutine build_tetra_rt_face_basis_to_local_matrix( &
            degree, permutation, transform, status)
        integer, intent(in) :: degree, permutation(3)
        real(dp), intent(out) :: transform(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: canonical(:), local(:)
        integer :: column, size_

        transform = 0.0_dp
        size_ = (degree + 1)*(degree + 2)/2
        status = 1
        if (degree < 0) return
        if (size(transform, 1) /= size_ .or. size(transform, 2) /= size_) &
            return
        if (degree > 4) then
            call build_runtime_rt_transform_matrix( &
                degree, permutation, .true., transform, status)
            return
        end if
        allocate(canonical(size_), local(size_))
        do column = 1, size_
            canonical = 0.0_dp
            canonical(column) = 1.0_dp
            call generated_rt_basis_to_local( &
                degree, permutation, canonical, local, status)
            if (status /= 0) return
            transform(:, column) = local
        end do
        status = 0
    end subroutine build_tetra_rt_face_basis_to_local_matrix

    subroutine map_tetra_rt_face_basis_to_local( &
            degree, permutation, canonical, local, status)
        integer, intent(in) :: degree, permutation(3)
        real(dp), intent(in) :: canonical(:)
        real(dp), intent(out) :: local(:)
        integer, intent(out) :: status

        if (degree <= 4) then
            call generated_rt_basis_to_local( &
                degree, permutation, canonical, local, status)
        else
            call apply_runtime_rt_transform( &
                degree, permutation, canonical, local, .true., status)
        end if
    end subroutine map_tetra_rt_face_basis_to_local

    subroutine transform_tetra_rt_face_moments( &
            degree, permutation, local, canonical, status)
        integer, intent(in) :: degree, permutation(3)
        real(dp), intent(in) :: local(:)
        real(dp), intent(out) :: canonical(:)
        integer, intent(out) :: status

        if (degree <= 4) then
            call generated_rt_transform_moments( &
                degree, permutation, local, canonical, status)
        else
            call apply_runtime_rt_transform( &
                degree, permutation, local, canonical, .false., status)
        end if
    end subroutine transform_tetra_rt_face_moments

    subroutine apply_runtime_rt_transform( &
            degree, permutation, input, output, basis_to_local, status)
        integer, intent(in) :: degree, permutation(3)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        logical, intent(in) :: basis_to_local
        integer, intent(out) :: status

        real(dp), allocatable :: transform(:, :)
        integer :: size_

        output = 0.0_dp
        status = 1
        size_ = (degree + 1)*(degree + 2)/2
        if (degree < 0) return
        if (size(input) /= size_ .or. size(output) /= size_) return
        if (.not. valid_permutation(permutation)) return
        allocate(transform(size_, size_))
        call build_runtime_rt_transform_matrix( &
            degree, permutation, basis_to_local, transform, status)
        if (status /= 0) return
        output = matmul(transform, input)
        status = 0
    end subroutine apply_runtime_rt_transform

    subroutine build_runtime_rt_transform_matrix( &
            degree, permutation, basis_to_local, transform, status)
        integer, intent(in) :: degree, permutation(3)
        logical, intent(in) :: basis_to_local
        real(dp), intent(out) :: transform(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: canonical(:, :), local(:, :)
        real(dp), allocatable :: solution(:, :)
        integer :: info, size_

        transform = 0.0_dp
        status = 1
        size_ = (degree + 1)*(degree + 2)/2
        if (degree < 0) return
        if (size(transform, 1) /= size_ .or. size(transform, 2) /= size_) &
            return
        if (.not. valid_permutation(permutation)) return
        allocate(canonical(size_, size_), local(size_, size_))
        allocate(solution(size_, size_))
        call build_rt_moment_matrix(degree, [1, 2, 3], canonical, status)
        if (status /= 0) return
        call build_rt_moment_matrix(degree, permutation, local, status)
        if (status /= 0) return
        if (basis_to_local) then
            call dense_solve( &
                transpose(canonical), transpose(local), solution, info)
        else
            call dense_solve( &
                transpose(local), transpose(canonical), solution, info)
        end if
        if (info /= 0) return
        transform = transpose(solution)
        status = 0
    end subroutine build_runtime_rt_transform_matrix

    subroutine build_rt_moment_matrix( &
            degree, permutation, matrix, status)
        integer, intent(in) :: degree, permutation(3)
        real(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: u(:), v(:), weights(:)
        integer :: affine(2, 2), offset(2)
        integer :: basis_degree, basis_x, basis_y, column
        integer :: moment_degree, moment_x, moment_y, node, row
        real(dp) :: mapped(2)

        status = 1
        if (.not. valid_permutation(permutation)) return
        call permutation_map(permutation, offset, affine)
        call triangle_duffy_quadrature( &
            2*degree + 4, u, v, weights, status)
        if (status /= 0) return
        matrix = 0.0_dp
        column = 0
        do basis_degree = 0, degree
            do basis_x = 0, basis_degree
                basis_y = basis_degree - basis_x
                column = column + 1
                row = 0
                do moment_degree = 0, degree
                    do moment_x = 0, moment_degree
                        moment_y = moment_degree - moment_x
                        row = row + 1
                        do node = 1, size(u)
                            mapped = real(offset, dp) + matmul( &
                                real(affine, dp), [u(node), v(node)])
                            matrix(row, column) = matrix(row, column) + &
                                weights(node)*rt_face_polynomial( &
                                degree, basis_x, basis_y, &
                                mapped(1), mapped(2))* &
                                rt_face_polynomial( &
                                degree, moment_x, moment_y, &
                                u(node), v(node))
                        end do
                    end do
                end do
            end do
        end do
        status = merge(0, 1, &
            row == size(matrix, 1) .and. column == size(matrix, 2))
    end subroutine build_rt_moment_matrix

    pure real(dp) function rt_face_polynomial( &
            degree, first_degree, second_degree, x, y) result(value)
        integer, intent(in) :: degree, first_degree, second_degree
        real(dp), intent(in) :: x, y

        if (degree <= 5) then
            value = x**first_degree*y**second_degree
        else
            value = triangle_dubiner(first_degree, second_degree, x, y)
        end if
    end function rt_face_polynomial

    subroutine apply_runtime_transform( &
            order, permutation, input, output, basis_to_local, status)
        integer, intent(in) :: order, permutation(3)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        logical, intent(in) :: basis_to_local
        integer, intent(out) :: status

        real(dp), allocatable :: transform(:, :)
        integer :: size_

        output = 0.0_dp
        status = 1
        size_ = order*(order - 1)
        if (order < 2) return
        if (size(input) /= size_ .or. size(output) /= size_) return
        if (.not. valid_permutation(permutation)) return
        allocate(transform(size_, size_))
        call build_runtime_transform_matrix( &
            order, permutation, basis_to_local, transform, status)
        if (status /= 0) return
        output = matmul(transform, input)
        status = 0
    end subroutine apply_runtime_transform

    subroutine build_runtime_transform_matrix( &
            order, permutation, basis_to_local, transform, status)
        integer, intent(in) :: order, permutation(3)
        logical, intent(in) :: basis_to_local
        real(dp), intent(out) :: transform(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: canonical(:, :), local(:, :)
        real(dp), allocatable :: solution(:, :)
        integer :: info, size_

        transform = 0.0_dp
        status = 1
        size_ = order*(order - 1)
        if (order < 2) return
        if (size(transform, 1) /= size_ .or. size(transform, 2) /= size_) &
            return
        if (.not. valid_permutation(permutation)) return
        allocate(canonical(size_, size_), local(size_, size_))
        allocate(solution(size_, size_))
        call build_moment_matrix(order, [1, 2, 3], canonical, status)
        if (status /= 0) return
        call build_moment_matrix(order, permutation, local, status)
        if (status /= 0) return
        if (basis_to_local) then
            call dense_solve( &
                transpose(canonical), transpose(local), solution, info)
        else
            call dense_solve( &
                transpose(local), transpose(canonical), solution, info)
        end if
        if (info /= 0) return
        transform = transpose(solution)
        status = 0
    end subroutine build_runtime_transform_matrix

    subroutine build_moment_matrix(order, permutation, matrix, status)
        integer, intent(in) :: order, permutation(3)
        real(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: u(:), v(:), weights(:)
        integer :: affine(2, 2), offset(2)
        integer :: basis_component, basis_degree, basis_x, basis_y, column
        integer :: moment_component, moment_degree, moment_x, moment_y, node
        integer :: row
        real(dp) :: mapped(2)

        call permutation_map(permutation, offset, affine)
        call triangle_duffy_quadrature(2*order + 2, u, v, weights, status)
        if (status /= 0) return
        matrix = 0.0_dp
        column = 0
        do basis_component = 1, 2
            do basis_degree = 0, order - 2
                do basis_x = 0, basis_degree
                    basis_y = basis_degree - basis_x
                    column = column + 1
                    row = 0
                    do moment_component = 1, 2
                        do moment_degree = 0, order - 2
                            do moment_x = 0, moment_degree
                                moment_y = moment_degree - moment_x
                                row = row + 1
                                do node = 1, size(u)
                                    mapped = real(offset, dp) + matmul( &
                                        real(affine, dp), [u(node), v(node)])
                                    matrix(row, column) = &
                                        matrix(row, column) + weights(node)* &
                                        real(affine( &
                                        basis_component, &
                                        moment_component), dp)* &
                                        face_polynomial( &
                                        order, basis_x, basis_y, &
                                        mapped(1), mapped(2))* &
                                        face_polynomial( &
                                        order, moment_x, moment_y, &
                                        u(node), v(node))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
        status = merge(0, 1, &
            row == size(matrix, 1) .and. column == size(matrix, 2))
    end subroutine build_moment_matrix

    pure real(dp) function face_polynomial( &
            order, first_degree, second_degree, x, y) result(value)
        integer, intent(in) :: order, first_degree, second_degree
        real(dp), intent(in) :: x, y

        if (order <= 5) then
            value = x**first_degree*y**second_degree
        else
            value = triangle_dubiner(first_degree, second_degree, x, y)
        end if
    end function face_polynomial

    pure subroutine permutation_map(permutation, offset, affine)
        integer, intent(in) :: permutation(3)
        integer, intent(out) :: offset(2), affine(2, 2)
        integer, parameter :: vertices(2, 3) = reshape( &
            [0, 0, 1, 0, 0, 1], [2, 3])

        offset = vertices(:, permutation(1))
        affine(:, 1) = vertices(:, permutation(2)) - offset
        affine(:, 2) = vertices(:, permutation(3)) - offset
    end subroutine permutation_map

    pure logical function valid_permutation(permutation) result(valid)
        integer, intent(in) :: permutation(3)

        valid = all(permutation >= 1) .and. all(permutation <= 3) .and. &
            count(permutation == 1) == 1 .and. &
            count(permutation == 2) == 1 .and. &
            count(permutation == 3) == 1
    end function valid_permutation

end module fortfem_tetra_face_moment_transforms
