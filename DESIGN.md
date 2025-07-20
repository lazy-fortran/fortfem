# FortFEM Design Document

## Vision

Create a modern Fortran finite element library that matches the ease-of-use of FreeFEM and FEniCS while leveraging Fortran's performance advantages. The library will provide an intuitive API for defining variational problems and support a progression from simple scalar elements to advanced finite element exterior calculus (FEEC).

## Core Design Principles

1. **User-Friendly API**: Natural mathematical notation for weak forms
2. **Performance**: Pure Fortran implementation with efficient assembly
3. **Modularity**: Clear separation of concerns with well-defined interfaces
4. **Extensibility**: Easy to add new element types and problem formulations
5. **Modern Fortran**: Use Fortran 2018 features for clean, maintainable code

## Architecture Overview

```
fortfem/
├── mesh/           - Mesh data structures and I/O
├── elements/       - Finite element definitions
├── quadrature/     - Integration rules
├── forms/          - Weak form abstractions
├── assembly/       - Matrix and vector assembly
├── solvers/        - Linear solver interfaces
└── problems/       - High-level problem definitions
```

## Key Components

### 1. Mesh Abstraction

```fortran
type :: mesh_t
    integer :: dim                        ! Spatial dimension (2 or 3)
    integer :: nvertices, nedges, ncells
    real(dp), allocatable :: vertices(:,:)
    integer, allocatable :: cells(:,:)    ! Cell connectivity
    integer, allocatable :: cell_type(:)  ! Triangle/quad identifier
    
    ! Topology relations
    integer, allocatable :: edges(:,:)
    integer, allocatable :: cell_edges(:,:)
    integer, allocatable :: edge_cells(:,:)
contains
    procedure :: read_from_file
    procedure :: refine
    procedure :: get_boundary_edges
end type mesh_t
```

### 2. Element Hierarchy

```fortran
! Base finite element type
type, abstract :: finite_element_t
    integer :: degree
    integer :: ndof_per_vertex, ndof_per_edge, ndof_per_cell
    integer :: ndof_total
contains
    procedure(eval_basis_interface), deferred :: eval_basis
    procedure(eval_grad_basis_interface), deferred :: eval_grad_basis
    procedure :: get_quadrature
end type finite_element_t

! Scalar H1-conforming elements
type, extends(finite_element_t) :: lagrange_element_t
    ! P1, P2, etc. Lagrange elements
end type lagrange_element_t

! Vector H(div)-conforming elements  
type, extends(finite_element_t) :: raviart_thomas_element_t
    ! RT0, RT1, etc. for flux approximation
end type raviart_thomas_element_t

! Vector H(curl)-conforming elements
type, extends(finite_element_t) :: nedelec_element_t
    ! For electromagnetic problems
end type nedelec_element_t

! L2-conforming elements
type, extends(finite_element_t) :: discontinuous_element_t
    ! DG elements
end type discontinuous_element_t
```

### 3. Weak Form API

Inspired by FEniCS/UFL but adapted for Fortran:

```fortran
! Function space
type :: function_space_t
    type(mesh_t), pointer :: mesh
    class(finite_element_t), allocatable :: element
contains
    procedure :: create_test_function
    procedure :: create_trial_function
end type function_space_t

! Form components
type :: form_integrand_t
    procedure(integrand_eval), pointer :: eval => null()
end type form_integrand_t

! Bilinear form a(u,v)
type :: bilinear_form_t
    type(function_space_t), pointer :: trial_space
    type(function_space_t), pointer :: test_space
    type(form_integrand_t), allocatable :: integrands(:)
contains
    procedure :: add_domain_integral
    procedure :: add_boundary_integral
    procedure :: assemble => assemble_bilinear_form
end type bilinear_form_t

! Linear form L(v)
type :: linear_form_t
    type(function_space_t), pointer :: test_space
    type(form_integrand_t), allocatable :: integrands(:)
contains
    procedure :: add_domain_integral
    procedure :: add_boundary_integral
    procedure :: assemble => assemble_linear_form
end type linear_form_t
```

### 4. Problem Definition DSL

High-level API for common problems:

```fortran
! Example: Poisson equation -∇²u = f
subroutine setup_poisson_problem(problem, mesh, f)
    type(fem_problem_t), intent(out) :: problem
    type(mesh_t), intent(in) :: mesh
    procedure(scalar_field) :: f
    
    type(function_space_t) :: V
    type(test_function_t) :: v
    type(trial_function_t) :: u
    
    ! Create function space with P1 elements
    call V%init(mesh, lagrange_element(degree=1))
    
    ! Define variational problem
    v = V%test_function()
    u = V%trial_function()
    
    ! Bilinear form: a(u,v) = ∫∇u·∇v dx
    problem%a = integral(dot(grad(u), grad(v)))
    
    ! Linear form: L(v) = ∫fv dx
    problem%L = integral(f * v)
    
    ! Boundary conditions
    call problem%add_dirichlet_bc(boundary_id=1, value=0.0_dp)
end subroutine
```

### 5. Assembly Engine

Efficient element-wise assembly with optimization for different element types:

```fortran
type :: assembler_t
    type(mesh_t), pointer :: mesh
    class(finite_element_t), pointer :: element
    type(quadrature_rule_t) :: quad_rule
contains
    procedure :: assemble_matrix
    procedure :: assemble_vector
    procedure :: apply_boundary_conditions
end type assembler_t
```

### 6. Solver Interface

```fortran
type, abstract :: linear_solver_t
contains
    procedure(solve_interface), deferred :: solve
end type linear_solver_t

type, extends(linear_solver_t) :: lapack_solver_t
    ! Direct solver using LAPACK
end type

type, extends(linear_solver_t) :: iterative_solver_t
    ! Iterative solver (GMRES, CG, etc.)
end type
```

## Implementation Roadmap

### Phase 1: Foundation (Completed)
- [x] Project setup with FPM
- [x] Basic mesh data structures
- [x] P1 Lagrange elements on triangles
- [x] Simple assembly for Poisson equation
- [x] LAPACK solver interface

### Phase 2: Core Features (Partially Complete)
- [ ] P2 and higher order Lagrange elements
- [ ] Quadrilateral elements
- [x] Boundary condition handling
- [x] Mesh generation with boundary curves

### Phase 3: Advanced Elements (In Progress)
- [ ] Raviart-Thomas elements
- [x] Nedelec edge elements
- [ ] Mixed formulations
- [ ] DG elements

### Phase 4: Performance & Features
- [ ] Advanced sparse solvers
- [ ] Parallel assembly
- [ ] Isoparametric elements
- [ ] Error estimation

### Phase 5: FEEC & Beyond
- [ ] Full FEEC complex
- [ ] High-order geometry
- [ ] Advanced solvers
- [ ] Optimization

## Example Usage

```fortran
program heat_equation
    use fortfem
    implicit none
    
    type(mesh_t) :: mesh
    type(fem_problem_t) :: problem
    real(dp), allocatable :: solution(:)
    
    ! Load mesh
    call mesh%read_from_file("domain.mesh")
    
    ! Setup problem: -∇²u + u = 1
    call setup_helmholtz_problem(problem, mesh, &
        source=constant_function(1.0_dp), &
        coefficient=constant_function(1.0_dp))
    
    ! Solve
    call problem%solve(solution)
    
    ! Output/Visualization
    call plot(solution, "solution.png", "Poisson Solution")
end program
```

## Design Decisions

1. **Pure Fortran**: Unlike MEPHIT's mixed approach, we use pure Fortran for maintainability
2. **Object-Oriented**: Modern Fortran OOP for clean interfaces
3. **Functional Style**: Inspired by UFL for natural problem definition
4. **Type Safety**: Strong typing to catch errors at compile time
5. **Modularity**: Each component can be developed and tested independently

## Testing Strategy

- Unit tests for each module using FPM test framework
- Convergence tests for each element type
- Benchmark suite comparing with analytical solutions
- Integration tests for complete problems