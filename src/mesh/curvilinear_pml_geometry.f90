module fortfem_curvilinear_pml_geometry
    !! Geometry-to-stretch map for locally curved PML layers.
    !!
    !! Each cell owns a point on the layer and a unit outward normal.  The
    !! signed distance from the cell centroid to that tangent plane defines a
    !! polynomial attenuation.  The resulting full tensor is
    !!
    !!   S = I + i * (sigma_max/wave_number) * (d/width)**p * n n^T.
    !!
    !! The normal frame is caller-owned.  Consequently this neutral contract
    !! works for curved FEM, IGA, and surface-fitted layers without importing
    !! a particular distance or closest-point algorithm.  Derivatives are
    !! valid for a fixed active set (cell distances strictly away from zero).
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_generated_cartesian_pml_attenuation_jvp, only: &
        generated_cartesian_pml_attenuation_jvp
    use fortfem_generated_cartesian_pml_attenuation_vjp, only: &
        generated_cartesian_pml_attenuation_vjp
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: build_curvilinear_normal_pml_element_stretch
    public :: build_curvilinear_normal_pml_element_stretch_jvp
    public :: build_curvilinear_normal_pml_element_stretch_vjp

contains

    subroutine build_curvilinear_normal_pml_element_stretch( &
            vertices, cells, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree, stretch, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :)
        real(dp), intent(in) :: layer_origins(:, :), layer_normals(:, :)
        real(dp), intent(in) :: layer_width(:)
        real(dp), intent(in) :: wave_number, sigma_max
        integer, intent(in) :: polynomial_degree
        complex(dp), allocatable, intent(out) :: stretch(:, :, :)
        integer, intent(out) :: status

        real(dp) :: centroid(3), distance, attenuation
        integer :: cell, i, j

        status = 1
        if (allocated(stretch)) deallocate(stretch)
        if (.not. valid_inputs( &
            vertices, cells, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree)) return

        allocate(stretch(3, 3, size(cells, 2)))
        stretch = cmplx(0.0_dp, 0.0_dp, dp)
        do cell = 1, size(cells, 2)
            call cell_centroid(vertices, cells(:, cell), centroid)
            distance = dot_product( &
                layer_normals(:, cell), centroid - layer_origins(:, cell))
            do i = 1, 3
                stretch(i, i, cell) = cmplx(1.0_dp, 0.0_dp, dp)
            end do
            if (distance <= 0.0_dp) cycle
            attenuation = sigma_max*(distance/layer_width(cell))** &
                polynomial_degree/wave_number
            do i = 1, 3
                do j = 1, 3
                    stretch(i, j, cell) = cmplx( &
                        merge(1.0_dp, 0.0_dp, i == j), &
                        attenuation*layer_normals(i, cell)* &
                        layer_normals(j, cell), dp)
                end do
            end do
        end do
        status = 0
    end subroutine build_curvilinear_normal_pml_element_stretch

    subroutine build_curvilinear_normal_pml_element_stretch_jvp( &
            vertices, cells, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree, vertices_dot, &
            origins_dot, normals_dot, width_dot, wave_number_dot, &
            sigma_max_dot, stretch_dot, status)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        integer, intent(in) :: cells(:, :)
        real(dp), intent(in) :: layer_origins(:, :), layer_normals(:, :)
        real(dp), intent(in) :: layer_width(:)
        real(dp), intent(in) :: wave_number, sigma_max
        integer, intent(in) :: polynomial_degree
        real(dp), intent(in) :: origins_dot(:, :), normals_dot(:, :)
        real(dp), intent(in) :: width_dot(:)
        real(dp), intent(in) :: wave_number_dot, sigma_max_dot
        complex(dp), allocatable, intent(out) :: stretch_dot(:, :, :)
        integer, intent(out) :: status

        real(dp) :: centroid(3), centroid_dot(3), displacement(3)
        real(dp) :: distance, distance_dot, attenuation, attenuation_dot
        real(dp) :: outer_dot
        integer :: cell, i, j

        status = 1
        if (allocated(stretch_dot)) deallocate(stretch_dot)
        if (.not. valid_inputs( &
            vertices, cells, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree)) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (any(shape(origins_dot) /= shape(layer_origins))) return
        if (any(shape(normals_dot) /= shape(layer_normals))) return
        if (size(width_dot) /= size(layer_width)) return

        allocate(stretch_dot(3, 3, size(cells, 2)))
        stretch_dot = cmplx(0.0_dp, 0.0_dp, dp)
        do cell = 1, size(cells, 2)
            call cell_centroid(vertices, cells(:, cell), centroid)
            call cell_centroid(vertices_dot, cells(:, cell), centroid_dot)
            displacement = centroid - layer_origins(:, cell)
            distance = dot_product(layer_normals(:, cell), displacement)
            if (distance <= 0.0_dp) cycle
            distance_dot = dot_product( &
                normals_dot(:, cell), displacement) + &
                dot_product(layer_normals(:, cell), &
                centroid_dot - origins_dot(:, cell))
            call generated_cartesian_pml_attenuation_jvp( &
                distance, layer_width(cell), wave_number, sigma_max, &
                real(polynomial_degree, dp), distance_dot, width_dot(cell), &
                wave_number_dot, sigma_max_dot, attenuation_dot)
            attenuation = sigma_max*(distance/layer_width(cell))** &
                polynomial_degree/wave_number
            do i = 1, 3
                do j = 1, 3
                    outer_dot = normals_dot(i, cell)*layer_normals(j, cell) + &
                        layer_normals(i, cell)*normals_dot(j, cell)
                    stretch_dot(i, j, cell) = cmplx(0.0_dp, &
                        attenuation_dot*layer_normals(i, cell)* &
                        layer_normals(j, cell) + attenuation*outer_dot, dp)
                end do
            end do
        end do
        status = 0
    end subroutine build_curvilinear_normal_pml_element_stretch_jvp

    subroutine build_curvilinear_normal_pml_element_stretch_vjp( &
            vertices, cells, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree, stretch_bar, &
            vertices_bar, origins_bar, normals_bar, width_bar, &
            wave_number_bar, sigma_max_bar, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :)
        real(dp), intent(in) :: layer_origins(:, :), layer_normals(:, :)
        real(dp), intent(in) :: layer_width(:)
        real(dp), intent(in) :: wave_number, sigma_max
        integer, intent(in) :: polynomial_degree
        complex(dp), intent(in) :: stretch_bar(:, :, :)
        real(dp), intent(out) :: vertices_bar(:, :), origins_bar(:, :)
        real(dp), intent(out) :: normals_bar(:, :), width_bar(:)
        real(dp), intent(out) :: wave_number_bar, sigma_max_bar
        integer, intent(out) :: status

        real(dp) :: centroid(3), displacement(3), normal_bar_local(3)
        real(dp) :: distance, attenuation, attenuation_bar, distance_bar
        real(dp) :: layer_width_bar_local, wave_number_bar_local
        real(dp) :: sigma_max_bar_local, bar_matrix(3, 3), normal(3)
        integer :: cell

        status = 1
        vertices_bar = 0.0_dp
        origins_bar = 0.0_dp
        normals_bar = 0.0_dp
        width_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        sigma_max_bar = 0.0_dp
        if (.not. valid_inputs( &
            vertices, cells, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree)) return
        if (any(shape(stretch_bar) /= [3, 3, size(cells, 2)])) return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (any(shape(origins_bar) /= shape(layer_origins))) return
        if (any(shape(normals_bar) /= shape(layer_normals))) return
        if (size(width_bar) /= size(layer_width)) return
        if (.not. all(ieee_is_finite(real(stretch_bar, dp)))) return
        if (.not. all(ieee_is_finite(aimag(stretch_bar)))) return

        do cell = 1, size(cells, 2)
            call cell_centroid(vertices, cells(:, cell), centroid)
            normal = layer_normals(:, cell)
            displacement = centroid - layer_origins(:, cell)
            distance = dot_product(normal, displacement)
            if (distance <= 0.0_dp) cycle
            attenuation = sigma_max*(distance/layer_width(cell))** &
                polynomial_degree/wave_number
            bar_matrix = aimag(stretch_bar(:, :, cell))
            attenuation_bar = sum(bar_matrix*spread(normal, 2, 3)* &
                spread(normal, 1, 3))
            normal_bar_local = attenuation*(matmul(bar_matrix, normal) + &
                matmul(transpose(bar_matrix), normal))
            call generated_cartesian_pml_attenuation_vjp( &
                distance, layer_width(cell), wave_number, sigma_max, &
                real(polynomial_degree, dp), attenuation_bar, distance_bar, &
                layer_width_bar_local, wave_number_bar_local, &
                sigma_max_bar_local)
            width_bar(cell) = width_bar(cell) + layer_width_bar_local
            wave_number_bar = wave_number_bar + wave_number_bar_local
            sigma_max_bar = sigma_max_bar + sigma_max_bar_local
            normals_bar(:, cell) = normals_bar(:, cell) + normal_bar_local + &
                distance_bar*displacement
            vertices_bar(:, cells(:, cell)) = vertices_bar(:, cells(:, cell)) + &
                spread(distance_bar*normal/real(size(cells, 1), dp), 2, &
                size(cells, 1))
            origins_bar(:, cell) = origins_bar(:, cell) - distance_bar*normal
        end do
        status = 0
    end subroutine build_curvilinear_normal_pml_element_stretch_vjp

    subroutine cell_centroid(vertices, cell, centroid)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cell(:)
        real(dp), intent(out) :: centroid(3)

        centroid = sum(vertices(:, cell), dim=2)/real(size(cell), dp)
    end subroutine cell_centroid

    pure logical function valid_inputs( &
            vertices, cells, layer_origins, layer_normals, layer_width, &
            wave_number, sigma_max, polynomial_degree) result(valid)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :)
        real(dp), intent(in) :: layer_origins(:, :), layer_normals(:, :)
        real(dp), intent(in) :: layer_width(:)
        real(dp), intent(in) :: wave_number, sigma_max
        integer, intent(in) :: polynomial_degree

        integer :: cell
        real(dp) :: normal_length

        valid = .false.
        if (size(vertices, 1) /= 3 .or. size(vertices, 2) < 1) return
        if (size(cells, 1) < 1 .or. size(cells, 2) < 1) return
        if (any(cells < 1) .or. any(cells > size(vertices, 2))) return
        if (any(shape(layer_origins) /= [3, size(cells, 2)])) return
        if (any(shape(layer_normals) /= [3, size(cells, 2)])) return
        if (size(layer_width) /= size(cells, 2)) return
        if (wave_number <= 0.0_dp .or. sigma_max <= 0.0_dp) return
        if (polynomial_degree < 1) return
        if (.not. all(ieee_is_finite(vertices))) return
        if (.not. all(ieee_is_finite(layer_origins))) return
        if (.not. all(ieee_is_finite(layer_normals))) return
        if (.not. all(ieee_is_finite(layer_width))) return
        if (any(layer_width <= 0.0_dp)) return
        do cell = 1, size(cells, 2)
            normal_length = sqrt(dot_product( &
                layer_normals(:, cell), layer_normals(:, cell)))
            if (abs(normal_length - 1.0_dp) > 1.0e-10_dp) return
        end do
        valid = .true.
    end function valid_inputs

end module fortfem_curvilinear_pml_geometry
