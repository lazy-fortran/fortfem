module fortfem_dirichlet_constraints
    !! Dirichlet lists (dofs(k), values(k)) as the mask and fixed values of
    !! sparse_direct_solve_constrained, and the cotangent of that map.
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: dirichlet_mask
    public :: dirichlet_mask_vjp

contains

    pure subroutine dirichlet_mask(dof_count, dofs, values, constrained, fixed)
        !! A repeated dof keeps its last value, as a dense row elimination does.
        !! The map is linear in values, so it is also its own JVP.
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

    pure subroutine dirichlet_mask_vjp(dofs, fixed_bar, values_bar)
        !! Adjoint of dirichlet_mask: the last occurrence of a dof receives its
        !! cotangent, earlier (overridden) occurrences zero.
        integer, intent(in) :: dofs(:)
        complex(dp), intent(in) :: fixed_bar(:)
        complex(dp), intent(out) :: values_bar(:)

        logical, allocatable :: seen(:)
        integer :: k

        allocate(seen(size(fixed_bar)))
        seen = .false.
        values_bar = (0.0_dp, 0.0_dp)
        do k = size(dofs), 1, -1
            if (seen(dofs(k))) cycle
            seen(dofs(k)) = .true.
            values_bar(k) = fixed_bar(dofs(k))
        end do
    end subroutine dirichlet_mask_vjp

end module fortfem_dirichlet_constraints
