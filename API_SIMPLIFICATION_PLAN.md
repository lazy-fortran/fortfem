# FortFEM API Simplification Plan

## Goal
Transform FortFEM into a FEniCS-style API where users can solve PDEs in ~10 lines of code without seeing implementation details.

## Current Issues
1. Manual assembly routines exposed
2. Direct matrix manipulation for BCs
3. Verbose mesh/space creation
4. No high-level BC abstractions
5. Complex DOF management
6. Missing vector function spaces
7. No expression language
8. Solver details exposed to users

## Proposed API Design

### 1. Mesh Creation (Immediate)
```fortran
! Current (verbose):
call mesh%create_rectangular(nx=20, ny=20, x_min=0.0_dp, x_max=1.0_dp, ...)
call mesh%build_connectivity()
call mesh%find_boundary()

! Proposed:
mesh = unit_square_mesh(20, 20)
! or
mesh = rectangle_mesh([0,0], [2,1], nx=40, ny=20)
```

### 2. Function Spaces (Immediate)
```fortran
! Current:
type(function_space_t) :: V
call create_P1_space(mesh, V)

! Proposed:
V = function_space(mesh, "P", 1)         ! P1 Lagrange
W = function_space(mesh, "P", 2)         ! P2 Lagrange  
E = function_space(mesh, "N", 1)         ! Nedelec
RT = function_space(mesh, "RT", 0)       ! Raviart-Thomas
```

### 3. Boundary Conditions (High Priority)
```fortran
! Current (manual):
do i = 1, n
    if (on_boundary(i)) then
        matrix(i,:) = 0; matrix(i,i) = 1; rhs(i) = 0
    end if
end do

! Proposed:
bc = dirichlet_bc(V, 0.0_dp, boundary_marker)
! or with function:
bc = dirichlet_bc(V, g, boundary_marker)  ! where g is a function_t
```

### 4. Weak Forms (High Priority)
```fortran
! Current:
call assemble_stiffness_matrix(...)
call assemble_mass_matrix(...)

! Proposed - FEniCS-style forms in Fortran:
u = trial_function(V)
v = test_function(V)
f = source_function(V)

a = inner(grad(u), grad(v))*dx
L = f*v*dx
```

### 5. Assembly & Solve (Immediate)
```fortran
! Current:
allocate(matrix(n,n), rhs(n))
call assemble_system(...)
call apply_bc(...)
call dgesv(...)

! Proposed:
solve(a == L, u, bc)
! or
A, b = assemble_system(a, L, bc)
solve(A, u, b)
```

### 6. Complete Example - Target API
```fortran
program poisson
    use fortfem
    implicit none
    
    type(mesh_t) :: mesh
    type(function_space_t) :: V
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: uh, f
    type(dirichlet_bc_t) :: bc
    type(bilinear_form_t) :: a
    type(linear_form_t) :: L
    
    ! Define mesh and function space
    mesh = unit_square_mesh(32, 32)
    V = function_space(mesh, "Lagrange", 1)
    
    ! Define variational problem
    u = trial_function(V)
    v = test_function(V)
    f = constant(1.0_dp)
    
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx
    
    ! Define boundary condition
    bc = dirichlet_bc(V, 0.0_dp, boundary)
    
    ! Compute solution
    uh = function(V)
    solve(a == L, uh, bc)
    
    ! Plot solution
    call plot(uh)
    
end program
```

## Implementation Phases

### Phase 1: Core Simplifications (Immediate)
1. Add mesh factory functions: `UnitSquareMesh`, `RectangleMesh`
2. Add `FunctionSpace` factory with string-based element selection
3. Add `DirichletBC` type with automatic application
4. Add simple `solve` interface hiding matrix assembly
5. Add `write_vtk` and `plot` functions for `Function` types

### Phase 2: Form Language (1-2 weeks)
1. Implement form algebra: `inner()`, `dot()`, `cross()`, `outer()`
2. Differential operators: `grad()`, `div()`, `curl()`, `jump()`, `avg()`
3. Measure types: `dx` (cell), `ds` (exterior facet), `dS` (interior facet)
4. Form operations: `+`, `-`, `*` with proper operator overloading
5. Tensor operations for elasticity

### Phase 3: Advanced Features (2-4 weeks)
1. Add `NeumannBC` and `RobinBC` types with form integration
2. Add vector function spaces: `VectorFunctionSpace`
3. Add mixed function spaces for coupled problems
4. Support time derivatives: `dt(u)`
5. Nonlinear forms with automatic linearization

### Phase 4: Performance & Extensions (Future)
1. Parallel assembly and solvers
2. Adaptive mesh refinement
3. High-order elements
4. 3D support

## Breaking Changes
- Remove all manual assembly routines from public API
- Hide DOF management completely
- Replace procedural API with object-oriented design
- Remove direct matrix access

## Migration Guide
Provide examples showing how to convert existing code to new API.

## Benefits
1. **Ease of use**: Solve PDEs in <10 lines
2. **Readability**: Code matches mathematical notation
3. **Safety**: No manual index manipulation
4. **Flexibility**: String expressions allow rapid prototyping
5. **Performance**: Generated code can be optimized

## Testing Strategy
1. Ensure all simplified examples work
2. Benchmark performance vs current implementation
3. Test error messages and user guidance
4. Validate against analytical solutions