module fortfem_maxwell_cfie_regularized_3d
    !! Stable product-algebra regularized CFIE following Scroggs et al.
    use fortfem_kinds, only: dp
    use fortfem_maxwell_bc_surface, only: assemble_maxwell_rwg_rbc_pairing
    use fortfem_maxwell_bc_surface, only: build_maxwell_bc_transformation
    use fortfem_maxwell_efie_bc_3d, only: &
        assemble_maxwell_efie_bc_imaginary_3d
    use fortfem_maxwell_efie_rwg_3d, only: &
        assemble_maxwell_efie_rwg_3d, assemble_maxwell_plane_wave_rhs_rwg_3d
    use fortfem_maxwell_localized_rwg_surface, only: &
        evaluate_maxwell_localized_rwg_basis
    use fortfem_maxwell_mfie_rwg_rbc_3d, only: &
        assemble_maxwell_mfie_rwg_rbc_3d
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private

    public :: assemble_maxwell_regularized_cfie_rwg_3d
    public :: assemble_maxwell_plane_wave_rhs_bc_3d
    public :: solve_maxwell_pec_regularized_cfie_rwg_3d
    public :: solve_maxwell_pec_regularized_cfie_rwg_multiple_3d

    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            complex(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            complex(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine zgesv
    end interface

contains

    subroutine solve_maxwell_pec_regularized_cfie_rwg_3d( &
            vertices, triangles, direction, polarization, wave_number, &
            impedance, quadrature_degree, tolerance, max_depth, density, status)
        real(dp), intent(in) :: vertices(:, :), direction(3), wave_number
        real(dp), intent(in) :: impedance, tolerance
        complex(dp), intent(in) :: polarization(3)
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: density(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: densities(:, :)
        real(dp) :: directions(3, 1)
        complex(dp) :: polarizations(3, 1)

        directions(:, 1) = direction
        polarizations(:, 1) = polarization
        call solve_maxwell_pec_regularized_cfie_rwg_multiple_3d( &
            vertices, triangles, directions, polarizations, wave_number, &
            impedance, quadrature_degree, tolerance, max_depth, densities, &
            status)
        if (status /= 0) return
        allocate(density(size(densities, 1)))
        density = densities(:, 1)
    end subroutine solve_maxwell_pec_regularized_cfie_rwg_3d

    subroutine solve_maxwell_pec_regularized_cfie_rwg_multiple_3d( &
            vertices, triangles, directions, polarizations, wave_number, &
            impedance, quadrature_degree, tolerance, max_depth, densities, &
            status)
        real(dp), intent(in) :: vertices(:, :), directions(:, :)
        real(dp), intent(in) :: wave_number, impedance, tolerance
        complex(dp), intent(in) :: polarizations(:, :)
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        complex(dp), allocatable, intent(out) :: densities(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: bc_rhs(:), cfie(:, :), efie(:, :)
        complex(dp), allocatable :: efie_rhs(:), mapped_rhs(:, :)
        complex(dp), allocatable :: mass(:, :), mfie(:, :), product(:, :)
        complex(dp), allocatable :: regularizer(:, :), right_hand_side(:, :)
        real(dp), allocatable :: real_mass(:, :)
        integer, allocatable :: pivots(:)
        integer :: incidence, incidence_count, info, n

        status = 1
        if (allocated(densities)) deallocate(densities)
        if (size(directions, 1) /= 3) return
        if (size(polarizations, 1) /= 3) return
        incidence_count = size(directions, 2)
        if (incidence_count < 1) return
        if (size(polarizations, 2) /= incidence_count) return
        call assemble_maxwell_regularized_cfie_rwg_3d( &
            vertices, triangles, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, cfie, efie, mfie, regularizer, product, status)
        if (status /= 0) return
        call assemble_maxwell_rwg_rbc_pairing( &
            vertices, triangles, quadrature_degree, real_mass, status)
        if (status /= 0) return
        n = size(real_mass, 1)
        allocate( &
            mass(n, n), mapped_rhs(n, incidence_count), &
            right_hand_side(n, incidence_count), pivots(n), &
            densities(n, incidence_count))
        do incidence = 1, incidence_count
            call assemble_maxwell_plane_wave_rhs_rwg_3d( &
                vertices, triangles, directions(:, incidence), &
                polarizations(:, incidence), wave_number, quadrature_degree, &
                efie_rhs, status)
            if (status /= 0) return
            call assemble_maxwell_plane_wave_rhs_bc_3d( &
                vertices, triangles, directions(:, incidence), &
                polarizations(:, incidence), wave_number, quadrature_degree, &
                bc_rhs, status)
            if (status /= 0) return
            mapped_rhs(:, incidence) = efie_rhs
            right_hand_side(:, incidence) = bc_rhs
        end do
        mass = transpose(cmplx(real_mass, 0.0_dp, dp))
        call zgesv( &
            n, incidence_count, mass, n, pivots, mapped_rhs, n, info)
        if (info /= 0) then
            status = 2
            return
        end if
        right_hand_side = right_hand_side - matmul(regularizer, mapped_rhs)
        call zgesv( &
            n, incidence_count, cfie, n, pivots, right_hand_side, n, info)
        if (info /= 0) then
            status = 3
            return
        end if
        densities = right_hand_side
        status = 0
    end subroutine solve_maxwell_pec_regularized_cfie_rwg_multiple_3d

    subroutine assemble_maxwell_plane_wave_rhs_bc_3d( &
            vertices, triangles, direction, polarization, wave_number, &
            quadrature_degree, right_hand_side, status)
        real(dp), intent(in) :: vertices(:, :), direction(3), wave_number
        complex(dp), intent(in) :: polarization(3)
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: right_hand_side(:)
        integer, intent(out) :: status

        integer, allocatable :: refined_triangles(:, :)
        real(dp), allocatable :: eta(:), refined_vertices(:, :)
        real(dp), allocatable :: transformation(:, :), weights(:), xi(:)
        real(dp) :: divergence, jacobian, local_value(3), point(3)
        complex(dp) :: incident_field(3)
        integer :: basis, local_edge, node, panel, row

        status = 1
        if (wave_number < 0.0_dp) return
        if (abs(norm2(direction) - 1.0_dp) > 128.0_dp*epsilon(1.0_dp)) return
        if (abs(sum(polarization*direction)) > &
            128.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            sqrt(sum(abs(polarization)**2)))) return
        call build_maxwell_bc_transformation( &
            vertices, triangles, refined_vertices, refined_triangles, &
            transformation, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        if (status /= 0) return
        allocate(right_hand_side(size(transformation, 2)))
        right_hand_side = cmplx(0.0_dp, 0.0_dp, dp)
        do panel = 1, size(refined_triangles, 2)
            jacobian = 2.0_dp*triangle_area( &
                refined_vertices(:, refined_triangles(:, panel)))
            do node = 1, size(weights)
                point = triangle_point( &
                    refined_vertices(:, refined_triangles(:, panel)), &
                    xi(node), eta(node))
                incident_field = polarization*exp(cmplx( &
                    0.0_dp, wave_number*dot_product(direction, point), dp))
                do local_edge = 1, 3
                    call evaluate_maxwell_localized_rwg_basis( &
                        refined_vertices, refined_triangles, panel, local_edge, &
                        point, local_value, divergence, status)
                    if (status /= 0) return
                    row = 3*(panel - 1) + local_edge
                    do basis = 1, size(transformation, 2)
                        right_hand_side(basis) = right_hand_side(basis) - &
                            jacobian*weights(node)*transformation(row, basis)* &
                            sum(cmplx(local_value, 0.0_dp, dp)*incident_field)
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine assemble_maxwell_plane_wave_rhs_bc_3d

    subroutine assemble_maxwell_regularized_cfie_rwg_3d( &
            vertices, triangles, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, matrix, efie, mfie, regularizer, &
            regularized_efie, status)
        real(dp), intent(in) :: vertices(:, :), wave_number, impedance
        integer, intent(in) :: triangles(:, :), quadrature_degree, max_depth
        real(dp), intent(in) :: tolerance
        complex(dp), allocatable, intent(out) :: matrix(:, :), efie(:, :)
        complex(dp), allocatable, intent(out) :: mfie(:, :), regularizer(:, :)
        complex(dp), allocatable, intent(out) :: regularized_efie(:, :)
        integer, intent(out) :: status

        complex(dp), allocatable :: mass(:, :), mapped_efie(:, :)
        complex(dp), allocatable :: principal_value(:, :)
        real(dp), allocatable :: real_mass(:, :)
        integer, allocatable :: pivots(:)
        integer :: info, size_system

        status = 1
        if (wave_number <= 0.0_dp .or. impedance <= 0.0_dp) return
        call assemble_maxwell_efie_rwg_3d( &
            vertices, triangles, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, efie, status)
        if (status /= 0) return
        call assemble_maxwell_efie_bc_imaginary_3d( &
            vertices, triangles, wave_number, impedance, quadrature_degree, &
            tolerance, max_depth, regularizer, status)
        if (status /= 0) return
        call assemble_maxwell_mfie_rwg_rbc_3d( &
            vertices, triangles, wave_number, quadrature_degree, tolerance, &
            max_depth, principal_value, mfie, status)
        if (status /= 0) return
        call assemble_maxwell_rwg_rbc_pairing( &
            vertices, triangles, quadrature_degree, real_mass, status)
        if (status /= 0) return
        size_system = size(real_mass, 1)
        if (size(real_mass, 2) /= size_system) return
        if (any(shape(efie) /= [size_system, size_system])) return
        if (any(shape(mfie) /= [size_system, size_system])) return
        if (any(shape(regularizer) /= [size_system, size_system])) return
        allocate( &
            mass(size_system, size_system), &
            mapped_efie(size_system, size_system), pivots(size_system))
        mass = transpose(cmplx(real_mass, 0.0_dp, dp))
        mapped_efie = efie
        call zgesv( &
            size_system, size_system, mass, size_system, pivots, mapped_efie, &
            size_system, info)
        if (info /= 0) then
            status = 2
            return
        end if
        regularized_efie = matmul(regularizer, mapped_efie)
        matrix = mfie - regularized_efie
        status = 0
    end subroutine assemble_maxwell_regularized_cfie_rwg_3d

    pure function triangle_point(vertices, xi, eta) result(point)
        real(dp), intent(in) :: vertices(3, 3), xi, eta
        real(dp) :: point(3)

        point = vertices(:, 1) + xi*(vertices(:, 2) - vertices(:, 1)) + &
            eta*(vertices(:, 3) - vertices(:, 1))
    end function triangle_point

    pure function triangle_area(vertices) result(area)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: area

        area = 0.5_dp*norm2(cross_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1)))
    end function triangle_area

    pure function cross_product(first, second) result(product_result)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product_result(3)

        product_result = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_maxwell_cfie_regularized_3d
