module fortfem_elasticity_acoustic_common_2d
    !! Pieces shared by the coupled elasticity-acoustic solvers: the homogeneous
    !! isotropic material from (E, nu, rho_s, omega) through the generated leaf,
    !! the cell group of a uniform material with JVP and VJP, and the Dirichlet
    !! mask of sparse_direct_solve_constrained.
    !!
    !! Cotangents follow sparse_direct_solve_constrained_vjp: a complex input z
    !! gets z_bar = sum conj(dA/dz) A_bar and a real input x gets
    !! Re(sum conj(dA/dx) A_bar).
    use fortfem_block_csc, only: block_csc_pattern_t
    use fortfem_elasticity_p1_csc_2d, only: elasticity_p1_assemble, &
        elasticity_p1_assemble_jvp, elasticity_p1_assemble_vjp
    use fortfem_generated_isotropic_elastic_material, only: &
        generated_isotropic_elastic_material
    use fortfem_generated_isotropic_elastic_material_jvp, only: &
        generated_isotropic_elastic_material_jvp
    use fortfem_generated_isotropic_elastic_material_vjp, only: &
        generated_isotropic_elastic_material_vjp
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: elastic_material_t
    public :: elastic_material
    public :: elastic_material_jvp
    public :: elastic_material_vjp
    public :: uniform_cells_assemble
    public :: uniform_cells_assemble_jvp
    public :: uniform_cells_assemble_vjp
    public :: dirichlet_mask

    type :: elastic_material_t
        complex(dp) :: lambda = (0.0_dp, 0.0_dp)
        complex(dp) :: mu = (0.0_dp, 0.0_dp)
        complex(dp) :: inertia = (0.0_dp, 0.0_dp) ! rho_s omega**2
        complex(dp) :: p_impedance = (0.0_dp, 0.0_dp)
        complex(dp) :: s_impedance = (0.0_dp, 0.0_dp)
    end type elastic_material_t

contains

    pure subroutine elastic_material(young_modulus, poisson_ratio, solid_density, &
            angular_frequency, material)
        complex(dp), intent(in) :: young_modulus
        real(dp), intent(in) :: poisson_ratio, solid_density, angular_frequency
        type(elastic_material_t), intent(out) :: material

        call generated_isotropic_elastic_material(young_modulus, poisson_ratio, &
            solid_density, angular_frequency, material%lambda, material%mu, &
            material%inertia, material%p_impedance, material%s_impedance)
    end subroutine elastic_material

    pure subroutine elastic_material_jvp(young_modulus, poisson_ratio, &
            solid_density, angular_frequency, young_modulus_dot, &
            poisson_ratio_dot, solid_density_dot, angular_frequency_dot, &
            material_dot)
        complex(dp), intent(in) :: young_modulus, young_modulus_dot
        real(dp), intent(in) :: poisson_ratio, solid_density, angular_frequency
        real(dp), intent(in) :: poisson_ratio_dot, solid_density_dot
        real(dp), intent(in) :: angular_frequency_dot
        type(elastic_material_t), intent(out) :: material_dot

        call generated_isotropic_elastic_material_jvp(young_modulus, &
            poisson_ratio, solid_density, angular_frequency, young_modulus_dot, &
            poisson_ratio_dot, solid_density_dot, angular_frequency_dot, &
            material_dot%lambda, material_dot%mu, material_dot%inertia, &
            material_dot%p_impedance, material_dot%s_impedance)
    end subroutine elastic_material_jvp

    pure subroutine elastic_material_vjp(young_modulus, poisson_ratio, &
            solid_density, angular_frequency, material_bar, young_modulus_bar, &
            poisson_ratio_bar, solid_density_bar, angular_frequency_bar)
        !! material_bar holds cotangents in the convention above; the generated
        !! VJP is the plain transpose, so it gets conj(material_bar).
        complex(dp), intent(in) :: young_modulus
        real(dp), intent(in) :: poisson_ratio, solid_density, angular_frequency
        type(elastic_material_t), intent(in) :: material_bar
        complex(dp), intent(out) :: young_modulus_bar
        real(dp), intent(out) :: poisson_ratio_bar, solid_density_bar
        real(dp), intent(out) :: angular_frequency_bar

        complex(dp) :: e_bar, nu_bar, rho_bar, omega_bar

        call generated_isotropic_elastic_material_vjp(young_modulus, &
            poisson_ratio, solid_density, angular_frequency, &
            conjg(material_bar%lambda), conjg(material_bar%mu), &
            conjg(material_bar%inertia), conjg(material_bar%p_impedance), &
            conjg(material_bar%s_impedance), e_bar, nu_bar, rho_bar, omega_bar)
        young_modulus_bar = conjg(e_bar)
        poisson_ratio_bar = real(nu_bar, dp)
        solid_density_bar = real(rho_bar, dp)
        angular_frequency_bar = real(omega_bar, dp)
    end subroutine elastic_material_vjp

    subroutine uniform_cells_assemble(pattern, group, vertices, triangles, &
            material, values, status)
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: group, triangles(:, :)
        real(dp), intent(in) :: vertices(:, :)
        type(elastic_material_t), intent(in) :: material
        complex(dp), intent(inout) :: values(:)
        integer, intent(out) :: status

        integer :: nt

        nt = size(triangles, 2)
        call elasticity_p1_assemble(pattern, group, vertices, triangles, &
            spread(material%lambda, 1, nt), spread(material%mu, 1, nt), &
            spread(material%inertia, 1, nt), values, status)
    end subroutine uniform_cells_assemble

    subroutine uniform_cells_assemble_jvp(pattern, group, vertices, triangles, &
            material, vertices_dot, material_dot, values_dot, status)
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: group, triangles(:, :)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        type(elastic_material_t), intent(in) :: material, material_dot
        complex(dp), intent(inout) :: values_dot(:)
        integer, intent(out) :: status

        integer :: nt

        nt = size(triangles, 2)
        call elasticity_p1_assemble_jvp(pattern, group, vertices, triangles, &
            spread(material%lambda, 1, nt), spread(material%mu, 1, nt), &
            spread(material%inertia, 1, nt), vertices_dot, &
            spread(material_dot%lambda, 1, nt), spread(material_dot%mu, 1, nt), &
            spread(material_dot%inertia, 1, nt), values_dot, status)
    end subroutine uniform_cells_assemble_jvp

    subroutine uniform_cells_assemble_vjp(pattern, group, vertices, triangles, &
            material, values_bar, vertices_bar, material_bar, status)
        !! vertices_bar is overwritten; the Lame and inertia cotangents of all
        !! cells are summed into material_bar (impedances untouched).
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: group, triangles(:, :)
        real(dp), intent(in) :: vertices(:, :)
        type(elastic_material_t), intent(in) :: material
        complex(dp), intent(in) :: values_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :)
        type(elastic_material_t), intent(inout) :: material_bar
        integer, intent(out) :: status

        complex(dp), allocatable :: lambda_bar(:), mu_bar(:), inertia_bar(:)
        integer :: nt

        nt = size(triangles, 2)
        allocate(lambda_bar(nt), mu_bar(nt), inertia_bar(nt))
        call elasticity_p1_assemble_vjp(pattern, group, vertices, triangles, &
            spread(material%lambda, 1, nt), spread(material%mu, 1, nt), &
            spread(material%inertia, 1, nt), values_bar, vertices_bar, &
            lambda_bar, mu_bar, inertia_bar, status)
        if (status /= 0) return
        material_bar%lambda = material_bar%lambda + sum(lambda_bar)
        material_bar%mu = material_bar%mu + sum(mu_bar)
        material_bar%inertia = material_bar%inertia + sum(inertia_bar)
    end subroutine uniform_cells_assemble_vjp

    pure subroutine dirichlet_mask(dof_count, dofs, values, constrained, fixed)
        !! A repeated dof keeps its last value, as in the dense elimination before.
        integer, intent(in) :: dof_count, dofs(:)
        complex(dp), intent(in) :: values(:)
        logical, allocatable, intent(out) :: constrained(:)
        complex(dp), allocatable, intent(out) :: fixed(:)

        integer :: k

        allocate(constrained(dof_count), fixed(dof_count))
        constrained = .false.
        fixed = (0.0_dp, 0.0_dp)
        do k = 1, size(dofs)
            constrained(dofs(k)) = .true.
            fixed(dofs(k)) = values(k)
        end do
    end subroutine dirichlet_mask

end module fortfem_elasticity_acoustic_common_2d
