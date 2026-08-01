module fortfem_assembly_tetra_lagrange_curvilinear_pml_3d
    !! P1 scalar Helmholtz assembly with full curvilinear PML tensors.
    use fortfem_curvilinear_helmholtz_pml, only: &
        curvilinear_scalar_helmholtz_pml_coefficients, &
        curvilinear_scalar_helmholtz_pml_coefficients_jvp, &
        curvilinear_scalar_helmholtz_pml_coefficients_vjp
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3, det3_jvp, det3_vjp, inv3, inv3_jvp, &
        inv3_vjp
    use fortsparse, only: csc_from_triplet, csc_z_t, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_tetra_lagrange_curvilinear_pml_element
    public :: assemble_tetra_lagrange_curvilinear_pml_element_jvp
    public :: assemble_tetra_lagrange_curvilinear_pml_element_vjp
    public :: assemble_tetra_lagrange_curvilinear_pml_csc
    public :: assemble_tetra_lagrange_curvilinear_pml_csc_jvp
    public :: assemble_tetra_lagrange_curvilinear_pml_csc_vjp

contains

    subroutine assemble_tetra_lagrange_curvilinear_pml_element( &
            vertices, degree, quadrature_degree, stretch, wave_number, &
            matrix, status)
        real(dp), intent(in) :: vertices(3, 4)
        integer, intent(in) :: degree, quadrature_degree
        complex(dp), intent(in) :: stretch(3, 3)
        real(dp), intent(in) :: wave_number
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        complex(dp) :: gradient_coefficient(3, 3), mass_coefficient
        real(dp) :: determinant, inverse_jacobian(3, 3), jacobian(3, 3)
        real(dp) :: physical_gradients(3, 4), reference_gradients(3, 4)
        real(dp) :: volume, reference_mass(4, 4)
        integer :: inverse_status, row, column

        status = 1
        if (allocated(matrix)) deallocate(matrix)
        if (.not. valid_element_inputs( &
            vertices, degree, quadrature_degree, stretch, wave_number)) return
        call curvilinear_scalar_helmholtz_pml_coefficients( &
            stretch, gradient_coefficient, mass_coefficient, status)
        if (status /= 0) return
        call tetra_geometry(vertices, jacobian, determinant, inverse_jacobian, &
            inverse_status)
        if (inverse_status /= 0) return
        call p1_reference_data(reference_gradients, reference_mass)
        physical_gradients = matmul(transpose(inverse_jacobian), &
            reference_gradients)
        volume = determinant/6.0_dp
        allocate(matrix(4, 4), source=cmplx(0.0_dp, 0.0_dp, dp))
        do column = 1, 4
            do row = 1, 4
                matrix(row, column) = volume*( &
                    sum(physical_gradients(:, row)*matmul( &
                    gradient_coefficient, physical_gradients(:, column))) - &
                    wave_number**2*mass_coefficient*reference_mass(row, column))
            end do
        end do
        status = 0
    end subroutine assemble_tetra_lagrange_curvilinear_pml_element

    subroutine assemble_tetra_lagrange_curvilinear_pml_element_jvp( &
            vertices, degree, quadrature_degree, stretch, wave_number, &
            vertices_dot, stretch_dot, wave_number_dot, matrix_dot, status)
        real(dp), intent(in) :: vertices(3, 4), vertices_dot(3, 4)
        integer, intent(in) :: degree, quadrature_degree
        complex(dp), intent(in) :: stretch(3, 3), stretch_dot(3, 3)
        real(dp), intent(in) :: wave_number, wave_number_dot
        complex(dp), allocatable, intent(out) :: matrix_dot(:, :)
        integer, intent(out) :: status

        complex(dp) :: gradient_coefficient(3, 3), gradient_dot(3, 3)
        complex(dp) :: mass_coefficient, mass_dot
        complex(dp) :: energy, energy_dot, mass_energy
        real(dp) :: determinant, determinant_dot
        real(dp) :: inverse_jacobian(3, 3), inverse_jacobian_dot(3, 3)
        real(dp) :: jacobian(3, 3), jacobian_dot(3, 3)
        real(dp) :: physical_gradients(3, 4), physical_gradients_dot(3, 4)
        real(dp) :: reference_gradients(3, 4), reference_mass(4, 4)
        real(dp) :: volume, volume_dot
        integer :: row, column, local_status

        status = 1
        if (allocated(matrix_dot)) deallocate(matrix_dot)
        if (.not. valid_element_inputs( &
            vertices, degree, quadrature_degree, stretch, wave_number)) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        call curvilinear_scalar_helmholtz_pml_coefficients( &
            stretch, gradient_coefficient, mass_coefficient, status)
        if (status /= 0) return
        call curvilinear_scalar_helmholtz_pml_coefficients_jvp( &
            stretch, stretch_dot, gradient_dot, mass_dot, status)
        if (status /= 0) return
        call tetra_geometry(vertices, jacobian, determinant, inverse_jacobian, &
            local_status)
        if (local_status /= 0) return
        call tetra_jacobian(vertices_dot, jacobian_dot)
        call det3_jvp(jacobian, jacobian_dot, determinant_dot)
        call inv3_jvp( &
            jacobian, jacobian_dot, inverse_jacobian, inverse_jacobian_dot, &
            local_status)
        if (local_status /= 0) return
        call p1_reference_data(reference_gradients, reference_mass)
        physical_gradients = matmul(transpose(inverse_jacobian), &
            reference_gradients)
        physical_gradients_dot = matmul(transpose(inverse_jacobian_dot), &
            reference_gradients)
        volume = determinant/6.0_dp
        volume_dot = determinant_dot/6.0_dp
        allocate(matrix_dot(4, 4), source=cmplx(0.0_dp, 0.0_dp, dp))
        do column = 1, 4
            do row = 1, 4
                energy = sum(physical_gradients(:, row)*matmul( &
                    gradient_coefficient, physical_gradients(:, column)))
                energy_dot = sum(physical_gradients_dot(:, row)*matmul( &
                    gradient_coefficient, physical_gradients(:, column))) + &
                    sum(physical_gradients(:, row)*matmul( &
                    gradient_coefficient, physical_gradients_dot(:, column))) + &
                    sum(physical_gradients(:, row)*matmul( &
                    gradient_dot, physical_gradients(:, column)))
                mass_energy = mass_coefficient*reference_mass(row, column)
                matrix_dot(row, column) = volume_dot*(energy - &
                    wave_number**2*mass_energy) + volume*(energy_dot - &
                    2.0_dp*wave_number*wave_number_dot*mass_energy - &
                    wave_number**2*mass_dot*reference_mass(row, column))
            end do
        end do
        status = 0
    end subroutine assemble_tetra_lagrange_curvilinear_pml_element_jvp

    subroutine assemble_tetra_lagrange_curvilinear_pml_element_vjp( &
            vertices, degree, quadrature_degree, stretch, wave_number, &
            matrix_bar, vertices_bar, stretch_bar, wave_number_bar, status)
        real(dp), intent(in) :: vertices(3, 4)
        integer, intent(in) :: degree, quadrature_degree
        complex(dp), intent(in) :: stretch(3, 3), matrix_bar(:, :)
        real(dp), intent(in) :: wave_number
        real(dp), intent(out) :: vertices_bar(3, 4), wave_number_bar
        complex(dp), intent(out) :: stretch_bar(3, 3)
        integer, intent(out) :: status

        complex(dp) :: gradient_coefficient(3, 3), gradient_bar(3, 3)
        complex(dp) :: mass_coefficient, mass_bar, seed, energy
        real(dp) :: determinant, determinant_bar
        real(dp) :: inverse_jacobian(3, 3), inverse_bar(3, 3)
        real(dp) :: inverse_jacobian_bar(3, 3), jacobian(3, 3)
        real(dp) :: jacobian_bar(3, 3), determinant_jacobian_bar(3, 3)
        real(dp) :: physical_gradients(3, 4), physical_bar(3, 4)
        real(dp) :: reference_gradients(3, 4), reference_mass(4, 4)
        real(dp) :: volume
        integer :: row, column, local_status, i, j

        vertices_bar = 0.0_dp
        stretch_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wave_number_bar = 0.0_dp
        status = 1
        if (.not. valid_element_inputs( &
            vertices, degree, quadrature_degree, stretch, wave_number)) return
        if (any(shape(matrix_bar) /= [4, 4])) return
        call curvilinear_scalar_helmholtz_pml_coefficients( &
            stretch, gradient_coefficient, mass_coefficient, status)
        if (status /= 0) return
        call tetra_geometry(vertices, jacobian, determinant, inverse_jacobian, &
            local_status)
        if (local_status /= 0) return
        call p1_reference_data(reference_gradients, reference_mass)
        physical_gradients = matmul(transpose(inverse_jacobian), &
            reference_gradients)
        volume = determinant/6.0_dp
        gradient_bar = cmplx(0.0_dp, 0.0_dp, dp)
        mass_bar = cmplx(0.0_dp, 0.0_dp, dp)
        physical_bar = 0.0_dp
        determinant_bar = 0.0_dp
        do column = 1, 4
            do row = 1, 4
                seed = matrix_bar(row, column)
                energy = sum(physical_gradients(:, row)*matmul( &
                    gradient_coefficient, physical_gradients(:, column))) - &
                    wave_number**2*mass_coefficient* &
                    reference_mass(row, column)
                determinant_bar = determinant_bar + &
                    real(conjg(seed)*energy, dp)/6.0_dp
                gradient_bar = gradient_bar + volume*seed*outer_product( &
                    physical_gradients(:, row), physical_gradients(:, column))
                mass_bar = mass_bar - volume*wave_number**2* &
                    reference_mass(row, column)*seed
                wave_number_bar = wave_number_bar + real(conjg(seed)*volume* &
                    (-2.0_dp*wave_number*mass_coefficient* &
                    reference_mass(row, column)), dp)
                physical_bar(:, row) = physical_bar(:, row) + real( &
                    conjg(seed)*volume*matmul(gradient_coefficient, &
                    physical_gradients(:, column)), dp)
                physical_bar(:, column) = physical_bar(:, column) + real( &
                    conjg(seed)*volume*matmul(transpose(gradient_coefficient), &
                    physical_gradients(:, row)), dp)
            end do
        end do
        inverse_bar = 0.0_dp
        do column = 1, 4
            do i = 1, 3
                do j = 1, 3
                    inverse_bar(j, i) = inverse_bar(j, i) + &
                        reference_gradients(j, column)*physical_bar(i, column)
                end do
            end do
        end do
        call inv3_vjp( &
            jacobian, inverse_bar, inverse_jacobian, inverse_jacobian_bar, &
            local_status)
        if (local_status /= 0) return
        call det3_vjp(jacobian, determinant_bar, determinant_jacobian_bar)
        jacobian_bar = inverse_jacobian_bar + determinant_jacobian_bar
        call tetra_jacobian_vjp(jacobian_bar, vertices_bar)
        call curvilinear_scalar_helmholtz_pml_coefficients_vjp( &
            stretch, gradient_bar, mass_bar, stretch_bar, status)
    end subroutine assemble_tetra_lagrange_curvilinear_pml_element_vjp

    subroutine assemble_tetra_lagrange_curvilinear_pml_csc( &
            mesh_vertices, tetrahedra, degree, stretch, wave_number, matrix, &
            status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        complex(dp), intent(in) :: stretch(:, :, :)
        real(dp), intent(in) :: wave_number
        type(csc_z_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: element_matrix(:, :), values(:)
        integer, allocatable :: columns(:), rows(:)
        real(dp) :: vertices(3, 4)
        integer :: column, element, entry, node, row, local_status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Curvilinear tetrahedral scalar PML assembly failed")
        if (.not. valid_mesh_inputs( &
            mesh_vertices, tetrahedra, degree, stretch, wave_number)) return
        allocate(rows(16*size(tetrahedra, 2)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do element = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, element))
            end do
            call assemble_tetra_lagrange_curvilinear_pml_element( &
                vertices, degree, 0, stretch(:, :, element), wave_number, &
                element_matrix, local_status)
            if (local_status /= 0) return
            do column = 1, 4
                do row = 1, 4
                    entry = entry + 1
                    rows(entry) = tetrahedra(row, element)
                    columns(entry) = tetrahedra(column, element)
                    values(entry) = element_matrix(row, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            size(mesh_vertices, 2), size(mesh_vertices, 2), rows, columns, &
            values, matrix, status)
    end subroutine assemble_tetra_lagrange_curvilinear_pml_csc

    subroutine assemble_tetra_lagrange_curvilinear_pml_csc_jvp( &
            mesh_vertices, tetrahedra, degree, stretch, wave_number, &
            mesh_vertices_dot, stretch_dot, wave_number_dot, matrix_dot, &
            status)
        real(dp), intent(in) :: mesh_vertices(:, :), mesh_vertices_dot(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        complex(dp), intent(in) :: stretch(:, :, :), stretch_dot(:, :, :)
        real(dp), intent(in) :: wave_number, wave_number_dot
        type(csc_z_t), intent(out) :: matrix_dot
        type(fortsparse_status_t), intent(out) :: status

        complex(dp), allocatable :: element_dot(:, :), values(:)
        integer, allocatable :: columns(:), rows(:)
        real(dp) :: vertices(3, 4), vertices_dot(3, 4)
        integer :: column, element, entry, node, row, local_status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Curvilinear tetrahedral scalar PML JVP assembly failed")
        if (.not. valid_mesh_inputs( &
            mesh_vertices, tetrahedra, degree, stretch, wave_number)) return
        if (any(shape(mesh_vertices_dot) /= shape(mesh_vertices))) return
        if (any(shape(stretch_dot) /= shape(stretch))) return
        allocate(rows(16*size(tetrahedra, 2)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do element = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, element))
                vertices_dot(:, node) = &
                    mesh_vertices_dot(:, tetrahedra(node, element))
            end do
            call assemble_tetra_lagrange_curvilinear_pml_element_jvp( &
                vertices, degree, 0, stretch(:, :, element), wave_number, &
                vertices_dot, stretch_dot(:, :, element), wave_number_dot, &
                element_dot, local_status)
            if (local_status /= 0) return
            do column = 1, 4
                do row = 1, 4
                    entry = entry + 1
                    rows(entry) = tetrahedra(row, element)
                    columns(entry) = tetrahedra(column, element)
                    values(entry) = element_dot(row, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            size(mesh_vertices, 2), size(mesh_vertices, 2), rows, columns, &
            values, matrix_dot, status)
    end subroutine assemble_tetra_lagrange_curvilinear_pml_csc_jvp

    subroutine assemble_tetra_lagrange_curvilinear_pml_csc_vjp( &
            mesh_vertices, tetrahedra, degree, stretch, wave_number, &
            matrix_values_bar, mesh_vertices_bar, stretch_bar, wave_number_bar, &
            status)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        complex(dp), intent(in) :: stretch(:, :, :)
        real(dp), intent(in) :: wave_number
        complex(dp), intent(in) :: matrix_values_bar(:)
        real(dp), intent(out) :: mesh_vertices_bar(:, :)
        complex(dp), intent(out) :: stretch_bar(:, :, :)
        real(dp), intent(out) :: wave_number_bar
        type(fortsparse_status_t), intent(out) :: status

        type(csc_z_t) :: matrix
        complex(dp), allocatable :: element_bar(:, :)
        real(dp) :: vertices(3, 4), local_vertices_bar(3, 4)
        complex(dp) :: local_stretch_bar(3, 3)
        real(dp) :: local_wave_number_bar
        integer :: column, element, local_status, node, row

        mesh_vertices_bar = 0.0_dp
        stretch_bar = cmplx(0.0_dp, 0.0_dp, dp)
        wave_number_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Curvilinear tetrahedral scalar PML VJP assembly failed")
        if (.not. valid_mesh_inputs( &
            mesh_vertices, tetrahedra, degree, stretch, wave_number)) return
        if (any(shape(mesh_vertices_bar) /= shape(mesh_vertices))) return
        if (any(shape(stretch_bar) /= shape(stretch))) return
        call assemble_tetra_lagrange_curvilinear_pml_csc( &
            mesh_vertices, tetrahedra, degree, stretch, wave_number, matrix, &
            status)
        if (status%code /= 0) return
        if (size(matrix_values_bar) /= matrix%nnz) return
        allocate(element_bar(4, 4))
        do element = 1, size(tetrahedra, 2)
            do node = 1, 4
                vertices(:, node) = &
                    mesh_vertices(:, tetrahedra(node, element))
            end do
            do column = 1, 4
                do row = 1, 4
                    element_bar(row, column) = csc_value_bar_at( &
                        matrix, matrix_values_bar, tetrahedra(row, element), &
                        tetrahedra(column, element))
                end do
            end do
            call assemble_tetra_lagrange_curvilinear_pml_element_vjp( &
                vertices, degree, 0, stretch(:, :, element), wave_number, &
                element_bar, local_vertices_bar, local_stretch_bar, &
                local_wave_number_bar, local_status)
            if (local_status /= 0) return
            do node = 1, 4
                mesh_vertices_bar(:, tetrahedra(node, element)) = &
                    mesh_vertices_bar(:, tetrahedra(node, element)) + &
                    local_vertices_bar(:, node)
            end do
            stretch_bar(:, :, element) = local_stretch_bar
            wave_number_bar = wave_number_bar + local_wave_number_bar
        end do
        call status_set(status, 0, "")
    end subroutine assemble_tetra_lagrange_curvilinear_pml_csc_vjp

    subroutine p1_reference_data(gradients, mass)
        real(dp), intent(out) :: gradients(3, 4), mass(4, 4)
        integer :: index

        gradients(:, 1) = [-1.0_dp, -1.0_dp, -1.0_dp]
        gradients(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        gradients(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        gradients(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
        mass = 1.0_dp/20.0_dp
        do index = 1, 4
            mass(index, index) = 1.0_dp/10.0_dp
        end do
    end subroutine p1_reference_data

    subroutine tetra_geometry(vertices, jacobian, determinant, inverse, status)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3), inverse(3, 3), determinant
        integer, intent(out) :: status

        integer :: inverse_status

        call tetra_jacobian(vertices, jacobian)
        determinant = det3(jacobian)
        status = 1
        if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        call inv3(jacobian, inverse, inverse_status)
        if (inverse_status /= 0) return
        status = 0
    end subroutine tetra_geometry

    subroutine tetra_jacobian(vertices, jacobian)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3)

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
    end subroutine tetra_jacobian

    subroutine tetra_jacobian_vjp(jacobian_bar, vertices_bar)
        real(dp), intent(in) :: jacobian_bar(3, 3)
        real(dp), intent(out) :: vertices_bar(3, 4)

        vertices_bar(:, 1) = -sum(jacobian_bar, dim=2)
        vertices_bar(:, 2) = jacobian_bar(:, 1)
        vertices_bar(:, 3) = jacobian_bar(:, 2)
        vertices_bar(:, 4) = jacobian_bar(:, 3)
    end subroutine tetra_jacobian_vjp

    pure function outer_product(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3, 3)
        integer :: i, j

        do j = 1, 3
            do i = 1, 3
                product(i, j) = left(i)*right(j)
            end do
        end do
    end function outer_product

    pure logical function valid_element_inputs( &
            vertices, degree, quadrature_degree, stretch, wave_number) &
            result(valid)
        real(dp), intent(in) :: vertices(3, 4)
        integer, intent(in) :: degree, quadrature_degree
        complex(dp), intent(in) :: stretch(3, 3)
        real(dp), intent(in) :: wave_number

        valid = degree == 1 .and. quadrature_degree >= 0 .and. &
            wave_number > 0.0_dp .and. all(vertices == vertices) .and. &
            all(real(stretch) == real(stretch)) .and. &
            all(aimag(stretch) == aimag(stretch))
    end function valid_element_inputs

    pure logical function valid_mesh_inputs( &
            vertices, tetrahedra, degree, stretch, wave_number) result(valid)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree
        complex(dp), intent(in) :: stretch(:, :, :)
        real(dp), intent(in) :: wave_number

        valid = size(vertices, 1) == 3 .and. size(vertices, 2) >= 4 .and. &
            size(tetrahedra, 1) == 4 .and. size(tetrahedra, 2) >= 1 .and. &
            all(tetrahedra >= 1) .and. all(tetrahedra <= size(vertices, 2)) .and. &
            degree == 1 .and. size(stretch, 1) == 3 .and. &
            size(stretch, 2) == 3 .and. size(stretch, 3) == size(tetrahedra, 2) .and. &
            wave_number > 0.0_dp
    end function valid_mesh_inputs

    pure function csc_value_bar_at(matrix, values_bar, row, column) result(value)
        type(csc_z_t), intent(in) :: matrix
        complex(dp), intent(in) :: values_bar(:)
        integer, intent(in) :: row, column
        complex(dp) :: value
        integer :: entry

        value = cmplx(0.0_dp, 0.0_dp, dp)
        do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
            if (matrix%row_idx(entry) == row) then
                value = values_bar(entry)
                return
            end if
        end do
    end function csc_value_bar_at

end module fortfem_assembly_tetra_lagrange_curvilinear_pml_3d
