module fortfem_maxwell_virtual_casing
    !! Neutral magnetic-field traces on caller-owned target surfaces.
    !!
    !! The source map is the RWG Biot--Savart/Helmholtz representation. This
    !! wrapper only performs the target-side normal/tangential projection; it
    !! does not impose an equilibrium, wall, or plasma model.
    use fortfem_kinds, only: dp
    use fortfem_maxwell_magnetic_rwg_3d, only: &
        evaluate_maxwell_magnetic_field_rwg_3d_targets, &
        evaluate_maxwell_magnetic_field_rwg_3d_targets_jvp, &
        evaluate_maxwell_magnetic_field_rwg_3d_targets_vjp
    implicit none
    private

    public :: evaluate_maxwell_virtual_casing_rwg_3d
    public :: evaluate_maxwell_virtual_casing_rwg_3d_jvp
    public :: evaluate_maxwell_virtual_casing_rwg_3d_vjp

contains

    subroutine evaluate_maxwell_virtual_casing_rwg_3d( &
            vertices, triangles, coefficients, observations, normals, wave_number, &
            quadrature_degree, magnetic_fields, normal_fields, tangential_fields, &
            status)
        real(dp), intent(in) :: vertices(:, :), observations(:, :), normals(:, :)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: coefficients(:)
        complex(dp), intent(out) :: magnetic_fields(:, :), normal_fields(:)
        complex(dp), intent(out) :: tangential_fields(:, :)
        integer, intent(out) :: status

        integer :: target

        status = 1
        if (size(observations, 1) /= 3 .or. size(normals, 1) /= 3) return
        if (any(shape(normals) /= shape(observations))) return
        if (size(magnetic_fields, 1) /= 3 .or. &
            size(tangential_fields, 1) /= 3) return
        if (any(shape(magnetic_fields) /= shape(observations)) .or. &
            any(shape(tangential_fields) /= shape(observations))) return
        if (size(normal_fields) /= size(observations, 2)) return
        do target = 1, size(observations, 2)
            if (norm2(normals(:, target)) <= tiny(1.0_dp)) return
        end do
        call evaluate_maxwell_magnetic_field_rwg_3d_targets( &
            vertices, triangles, coefficients, observations, wave_number, &
            quadrature_degree, magnetic_fields, status)
        if (status /= 0) return
        do target = 1, size(observations, 2)
            normal_fields(target) = sum( &
                normals(:, target)*magnetic_fields(:, target))
            tangential_fields(:, target) = magnetic_fields(:, target) - &
                normal_fields(target)*normals(:, target)
        end do
        status = 0
    end subroutine evaluate_maxwell_virtual_casing_rwg_3d

    subroutine evaluate_maxwell_virtual_casing_rwg_3d_jvp( &
            vertices, triangles, coefficients, observations, normals, wave_number, &
            quadrature_degree, vertices_dot, coefficients_dot, observations_dot, &
            normals_dot, wave_number_dot, magnetic_fields_dot, normal_fields_dot, &
            tangential_fields_dot, status)
        real(dp), intent(in) :: vertices(:, :), observations(:, :), normals(:, :)
        real(dp), intent(in) :: wave_number, vertices_dot(:, :)
        real(dp), intent(in) :: observations_dot(:, :), normals_dot(:, :)
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: coefficients(:), coefficients_dot(:)
        real(dp), intent(in) :: wave_number_dot
        complex(dp), intent(out) :: magnetic_fields_dot(:, :)
        complex(dp), intent(out) :: normal_fields_dot(:), tangential_fields_dot(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: magnetic_fields(:, :)
        integer :: target

        status = 1
        if (size(observations, 1) /= 3 .or. size(normals, 1) /= 3) return
        if (any(shape(normals) /= shape(observations)) .or. &
            any(shape(observations_dot) /= shape(observations)) .or. &
            any(shape(normals_dot) /= shape(observations))) return
        if (size(magnetic_fields_dot, 1) /= 3 .or. &
            size(tangential_fields_dot, 1) /= 3) return
        if (any(shape(magnetic_fields_dot) /= shape(observations)) .or. &
            any(shape(tangential_fields_dot) /= shape(observations))) return
        if (size(normal_fields_dot) /= size(observations, 2)) return
        allocate(magnetic_fields(3, size(observations, 2)))
        do target = 1, size(observations, 2)
            if (norm2(normals(:, target)) <= tiny(1.0_dp)) return
        end do
        call evaluate_maxwell_magnetic_field_rwg_3d_targets( &
            vertices, triangles, coefficients, observations, wave_number, &
            quadrature_degree, magnetic_fields, status)
        if (status /= 0) return
        call evaluate_maxwell_magnetic_field_rwg_3d_targets_jvp( &
            vertices, triangles, coefficients, observations, wave_number, &
            quadrature_degree, vertices_dot, coefficients_dot, observations_dot, &
            wave_number_dot, magnetic_fields_dot, status)
        if (status /= 0) return
        do target = 1, size(observations, 2)
            normal_fields_dot(target) = sum( &
                normals_dot(:, target)*magnetic_fields(:, target)) + &
                sum(normals(:, target)*magnetic_fields_dot(:, target))
            tangential_fields_dot(:, target) = magnetic_fields_dot(:, target) - &
                normal_fields_dot(target)*normals(:, target) - &
                sum(normals(:, target)*magnetic_fields(:, target))* &
                normals_dot(:, target)
        end do
        status = 0
    end subroutine evaluate_maxwell_virtual_casing_rwg_3d_jvp

    subroutine evaluate_maxwell_virtual_casing_rwg_3d_vjp( &
            vertices, triangles, coefficients, observations, normals, wave_number, &
            quadrature_degree, magnetic_fields_bar, normal_fields_bar, &
            tangential_fields_bar, vertices_bar, coefficients_bar, observations_bar, &
            normals_bar, wave_number_bar, status)
        real(dp), intent(in) :: vertices(:, :), observations(:, :), normals(:, :)
        real(dp), intent(in) :: wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: coefficients(:), magnetic_fields_bar(:, :)
        complex(dp), intent(in) :: normal_fields_bar(:), tangential_fields_bar(:, :)
        real(dp), intent(out) :: vertices_bar(:, :), observations_bar(:, :)
        real(dp), intent(out) :: normals_bar(:, :), wave_number_bar
        complex(dp), allocatable, intent(out) :: coefficients_bar(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: source_fields_bar(:, :), magnetic_fields(:, :)
        complex(dp) :: scalar_bar, scalar_field
        integer :: target

        vertices_bar = 0.0_dp
        observations_bar = 0.0_dp
        normals_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        if (allocated(coefficients_bar)) deallocate(coefficients_bar)
        status = 1
        if (size(observations, 1) /= 3 .or. size(normals, 1) /= 3) return
        if (any(shape(normals) /= shape(observations))) return
        if (any(shape(magnetic_fields_bar) /= shape(observations)) .or. &
            any(shape(tangential_fields_bar) /= shape(observations))) return
        if (size(normal_fields_bar) /= size(observations, 2)) return
        do target = 1, size(observations, 2)
            if (norm2(normals(:, target)) <= tiny(1.0_dp)) return
        end do
        allocate( &
            source_fields_bar(3, size(observations, 2)), &
            magnetic_fields(3, size(observations, 2)))
        call evaluate_maxwell_magnetic_field_rwg_3d_targets( &
            vertices, triangles, coefficients, observations, wave_number, &
            quadrature_degree, magnetic_fields, status)
        if (status /= 0) return
        do target = 1, size(observations, 2)
            scalar_field = sum(normals(:, target)*magnetic_fields(:, target))
            scalar_bar = normal_fields_bar(target) - sum( &
                tangential_fields_bar(:, target)*normals(:, target))
            source_fields_bar(:, target) = magnetic_fields_bar(:, target) + &
                tangential_fields_bar(:, target) + scalar_bar*normals(:, target)
            normals_bar(:, target) = real(conjg(scalar_bar)* &
                magnetic_fields(:, target) - &
                conjg(tangential_fields_bar(:, target))*scalar_field, dp)
        end do
        call evaluate_maxwell_magnetic_field_rwg_3d_targets_vjp( &
            vertices, triangles, coefficients, observations, wave_number, &
            quadrature_degree, source_fields_bar, vertices_bar, coefficients_bar, &
            observations_bar, wave_number_bar, status)
        if (status /= 0) return
        status = 0
    end subroutine evaluate_maxwell_virtual_casing_rwg_3d_vjp

end module fortfem_maxwell_virtual_casing
