module fortfem_maxwell_cfie_regularized_3d
    !! Stable product-algebra regularized CFIE following Scroggs et al.
    use fortfem_kinds, only: dp
    use fortfem_maxwell_bc_surface, only: assemble_maxwell_rwg_rbc_pairing
    use fortfem_maxwell_efie_bc_3d, only: &
        assemble_maxwell_efie_bc_imaginary_3d
    use fortfem_maxwell_efie_rwg_3d, only: assemble_maxwell_efie_rwg_3d
    use fortfem_maxwell_mfie_rwg_rbc_3d, only: &
        assemble_maxwell_mfie_rwg_rbc_3d
    implicit none
    private

    public :: assemble_maxwell_regularized_cfie_rwg_3d

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

end module fortfem_maxwell_cfie_regularized_3d
