module fortfem_tetra_mixed_poisson_state_3d
    !! End-to-end moving-mesh RT--DG mixed Poisson state products.
    use fortfem_assembly_tetra_rt_arbitrary_order_3d, only: &
        assemble_tetra_rt_div_mass_csc, assemble_tetra_rt_div_mass_csc_jvp, &
        assemble_tetra_rt_div_mass_csc_vjp, &
        assemble_tetra_rt_divergence_csc
    use fortfem_kinds, only: dp
    use fortfem_mixed_rt_system, only: solve_mixed_rt_system, &
        solve_mixed_rt_system_jvp, solve_mixed_rt_system_vjp
    use fortfem_tetra_mixed_poisson_3d, only: &
        assemble_tetra_dg_source_load_samples, &
        assemble_tetra_dg_source_load_samples_jvp, &
        assemble_tetra_dg_source_load_samples_vjp
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none

    private

    public :: solve_tetra_mixed_poisson_state
    public :: solve_tetra_mixed_poisson_state_jvp
    public :: solve_tetra_mixed_poisson_state_vjp
    public :: solve_tetra_mixed_poisson_sampled_state
    public :: solve_tetra_mixed_poisson_sampled_state_jvp
    public :: solve_tetra_mixed_poisson_sampled_state_vjp

contains

    subroutine solve_tetra_mixed_poisson_sampled_state( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, source_values, flux, pressure, status)
        real(dp), intent(in) :: vertices(:, :), mass_coefficient
        real(dp), intent(in) :: source_values(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), allocatable, intent(out) :: flux(:), pressure(:)
        integer, intent(out) :: status

        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: load(:)

        status = 1
        call assemble_tetra_dg_source_load_samples( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            load, sparse_status)
        if (sparse_status%code /= 0) return
        call solve_tetra_mixed_poisson_state( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, load, flux, pressure, status)
    end subroutine solve_tetra_mixed_poisson_sampled_state

    subroutine solve_tetra_mixed_poisson_sampled_state_jvp( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, source_values, source_gradients, vertices_dot, &
            mass_coefficient_dot, source_parameter_dot, flux_dot, &
            pressure_dot, status)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        real(dp), intent(in) :: mass_coefficient, mass_coefficient_dot
        real(dp), intent(in) :: source_values(:, :)
        real(dp), intent(in) :: source_gradients(:, :, :)
        real(dp), intent(in) :: source_parameter_dot(:, :)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), allocatable, intent(out) :: flux_dot(:), pressure_dot(:)
        integer, intent(out) :: status

        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: load(:), load_dot(:)

        status = 1
        call assemble_tetra_dg_source_load_samples( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            load, sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_tetra_dg_source_load_samples_jvp( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            source_gradients, vertices_dot, source_parameter_dot, load_dot, &
            sparse_status)
        if (sparse_status%code /= 0) return
        call solve_tetra_mixed_poisson_state_jvp( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, load, vertices_dot, mass_coefficient_dot, &
            load_dot, flux_dot, pressure_dot, status)
    end subroutine solve_tetra_mixed_poisson_sampled_state_jvp

    subroutine solve_tetra_mixed_poisson_sampled_state_vjp( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, source_values, source_gradients, flux, &
            pressure, flux_bar, pressure_bar, vertices_bar, &
            mass_coefficient_bar, source_values_bar, status)
        real(dp), intent(in) :: vertices(:, :), mass_coefficient
        real(dp), intent(in) :: source_values(:, :)
        real(dp), intent(in) :: source_gradients(:, :, :)
        real(dp), intent(in) :: flux(:), pressure(:), flux_bar(:)
        real(dp), intent(in) :: pressure_bar(:)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), intent(out) :: vertices_bar(:, :)
        real(dp), intent(out) :: mass_coefficient_bar
        real(dp), intent(out) :: source_values_bar(:, :)
        integer, intent(out) :: status

        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: load(:), load_bar(:)
        real(dp) :: load_vertices_bar(size(vertices, 1), size(vertices, 2))

        vertices_bar = 0.0_dp
        mass_coefficient_bar = 0.0_dp
        source_values_bar = 0.0_dp
        status = 1
        call assemble_tetra_dg_source_load_samples( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            load, sparse_status)
        if (sparse_status%code /= 0) return
        allocate(load_bar(size(load)))
        call solve_tetra_mixed_poisson_state_vjp( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, load, flux, pressure, flux_bar, pressure_bar, &
            vertices_bar, mass_coefficient_bar, load_bar, status)
        if (status /= 0) return
        call assemble_tetra_dg_source_load_samples_vjp( &
            vertices, tetrahedra, degree, quadrature_degree, source_values, &
            source_gradients, load_bar, load_vertices_bar, source_values_bar, &
            sparse_status)
        if (sparse_status%code /= 0) then
            status = sparse_status%code
            return
        end if
        vertices_bar = vertices_bar + load_vertices_bar
    end subroutine solve_tetra_mixed_poisson_sampled_state_vjp

    subroutine solve_tetra_mixed_poisson_state( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, load, flux, pressure, status)
        real(dp), intent(in) :: vertices(:, :), mass_coefficient, load(:)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), allocatable, intent(out) :: flux(:), pressure(:)
        integer, intent(out) :: status

        type(csc_t) :: divergence, flux_mass
        type(fortsparse_status_t) :: sparse_status

        status = 1
        if (allocated(flux)) deallocate(flux)
        if (allocated(pressure)) deallocate(pressure)
        if (mass_coefficient <= 0.0_dp) return
        call assemble_operators( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, flux_mass, divergence, sparse_status)
        if (sparse_status%code /= 0) return
        if (size(load) /= divergence%nrow) return
        allocate(flux(flux_mass%nrow), pressure(divergence%nrow))
        call solve_mixed_rt_system( &
            flux_mass, divergence, load, flux, pressure, status)
    end subroutine solve_tetra_mixed_poisson_state

    subroutine solve_tetra_mixed_poisson_state_jvp( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, load, vertices_dot, mass_coefficient_dot, &
            load_dot, flux_dot, pressure_dot, status)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        real(dp), intent(in) :: mass_coefficient, mass_coefficient_dot
        real(dp), intent(in) :: load(:), load_dot(:)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), allocatable, intent(out) :: flux_dot(:), pressure_dot(:)
        integer, intent(out) :: status

        type(csc_t) :: divergence, divergence_dot, flux_mass, flux_mass_dot
        type(fortsparse_status_t) :: sparse_status

        status = 1
        if (allocated(flux_dot)) deallocate(flux_dot)
        if (allocated(pressure_dot)) deallocate(pressure_dot)
        if (mass_coefficient <= 0.0_dp) return
        call assemble_operators( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, flux_mass, divergence, sparse_status)
        if (sparse_status%code /= 0) return
        call assemble_tetra_rt_div_mass_csc_jvp( &
            vertices, tetrahedra, degree, quadrature_degree, 0.0_dp, &
            mass_coefficient, vertices_dot, 0.0_dp, mass_coefficient_dot, &
            flux_mass_dot, sparse_status)
        if (sparse_status%code /= 0) return
        divergence_dot = divergence
        divergence_dot%val = 0.0_dp
        if (size(load) /= divergence%nrow .or. size(load_dot) /= size(load)) &
            return
        allocate(flux_dot(flux_mass%nrow), pressure_dot(divergence%nrow))
        call solve_mixed_rt_system_jvp( &
            flux_mass, divergence, load, flux_mass_dot, divergence_dot, &
            load_dot, flux_dot, pressure_dot, status)
    end subroutine solve_tetra_mixed_poisson_state_jvp

    subroutine solve_tetra_mixed_poisson_state_vjp( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, load, flux, pressure, flux_bar, pressure_bar, &
            vertices_bar, mass_coefficient_bar, load_bar, status)
        real(dp), intent(in) :: vertices(:, :), mass_coefficient, load(:)
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        real(dp), intent(in) :: flux(:), pressure(:), flux_bar(:)
        real(dp), intent(in) :: pressure_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :)
        real(dp), intent(out) :: mass_coefficient_bar, load_bar(:)
        integer, intent(out) :: status

        type(csc_t) :: divergence, flux_mass
        type(fortsparse_status_t) :: sparse_status
        real(dp), allocatable :: divergence_values_bar(:)
        real(dp), allocatable :: flux_mass_values_bar(:)
        real(dp) :: unused_divergence_coefficient_bar

        vertices_bar = 0.0_dp
        mass_coefficient_bar = 0.0_dp
        load_bar = 0.0_dp
        status = 1
        if (mass_coefficient <= 0.0_dp) return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        call assemble_operators( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, flux_mass, divergence, sparse_status)
        if (sparse_status%code /= 0) return
        allocate(flux_mass_values_bar(flux_mass%nnz))
        allocate(divergence_values_bar(divergence%nnz))
        call solve_mixed_rt_system_vjp( &
            flux_mass, divergence, load, flux, pressure, flux_bar, &
            pressure_bar, flux_mass_values_bar, divergence_values_bar, &
            load_bar, status)
        if (status /= 0) return
        call assemble_tetra_rt_div_mass_csc_vjp( &
            vertices, tetrahedra, degree, quadrature_degree, 0.0_dp, &
            mass_coefficient, flux_mass_values_bar, vertices_bar, &
            unused_divergence_coefficient_bar, mass_coefficient_bar, &
            sparse_status)
        status = sparse_status%code
    end subroutine solve_tetra_mixed_poisson_state_vjp

    subroutine assemble_operators( &
            vertices, tetrahedra, degree, quadrature_degree, &
            mass_coefficient, flux_mass, divergence, status)
        real(dp), intent(in) :: vertices(:, :), mass_coefficient
        integer, intent(in) :: tetrahedra(:, :), degree, quadrature_degree
        type(csc_t), intent(out) :: flux_mass, divergence
        type(fortsparse_status_t), intent(out) :: status

        call assemble_tetra_rt_div_mass_csc( &
            vertices, tetrahedra, degree, quadrature_degree, flux_mass, &
            status, 0.0_dp, mass_coefficient)
        if (status%code /= 0) return
        call assemble_tetra_rt_divergence_csc( &
            vertices, tetrahedra, degree, quadrature_degree, divergence, &
            status)
    end subroutine assemble_operators

end module fortfem_tetra_mixed_poisson_state_3d
