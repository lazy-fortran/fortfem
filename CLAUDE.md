# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FortFEM is a modern Fortran finite element library designed for ease of use, inspired by FreeFEM and FEniCS. It provides natural mathematical notation for defining weak forms and supports various element types on triangular and quadrilateral meshes.

## Build and Development Commands

### Build Commands
- **Build the project**: `fpm build`
- **Build with specific compiler**: `fpm build --compiler gfortran`
- **Build with coverage**: `fpm build --flag "--coverage"`
- **Clean build**: `fpm clean --skip && fpm build`
- **Note**: Always use `fpm clean --skip` to avoid interactive prompts

### Run Commands
- **Run the main application**: `fpm run`
- **Run specific example**: `fpm run --example plot_basis`
- **List available examples**: `fpm run --example --list`
- **Use `fpm run --example <name>` to run examples**

### Test Commands
- **Run all tests**: `fpm test`
- **Run specific test**: `fpm test test-name`
- **Run tests with coverage**: `fpm test --flag "--coverage"`
- **List available tests**: `fpm test --list`

### Important Testing Notes
- ALWAYS place test files in the `test/` directory
- ALWAYS use `fpm test` to run tests, not manual compilation
- Test files should be named `test_*.f90`
- Use the `check` module for test assertions

### Development Workflow
1. Follow strict TDD: Write test first, then implementation
2. Source files go in `src/`
3. Tests go in `test/`
4. Examples go in `example/`
5. Documentation goes in `doc/`
6. Check TODO.md for current implementation tasks

## Architecture and Code Structure

### Directory Layout
```
fortfem/
├── src/          # Core library modules
├── test/         # Unit and integration tests
├── example/      # Example programs with visualization
├── app/          # Main applications
└── doc/          # Documentation
```

### Module Organization (Planned)
- **src/fortfem.f90**: Main module exporting public API
- **src/mesh/**: Mesh data structures and I/O
- **src/elements/**: Finite element definitions
- **src/quadrature/**: Integration rules
- **src/forms/**: Weak form abstractions
- **src/assembly/**: Matrix and vector assembly
- **src/solvers/**: Linear solver interfaces

### Dependencies
- **fortplotlib**: For visualization in examples
- **LAPACK**: For linear algebra (initially)
- **SuiteSparse**: For sparse solvers (future)

### Fortran Conventions
- Derived types follow `typename_t` naming convention
- Line limit: 88 characters
- Indentation: 4 spaces
- Arrays extend with: `arr = [arr, new_element]` syntax
- Use associate blocks for unused dummy arguments
- Real precision: `dp = kind(0.0d0)`

### Critical Fortran Design Principles
- **NO POINTERS**: Never use pointers in Fortran. Always use allocatables instead.
- **DEEP COPY ASSIGNMENT**: When a type contains allocatable members, always overload the assignment operator (=) with deep-copy semantics.
- **ALLOCATABLE FUNCTION RETURNS**: Prefer subroutines with intent(out) allocatable arguments over functions returning allocatables. If you do write a function returning an allocatable, ensure that the return value is guaranteed to be allocated.
- **SAFE ALLOCATION**: Always check if allocatable is already allocated before deallocating or reallocating.
- **CLEAN INTERFACES**: Prefer explicit interfaces over implicit ones.
- **INTENT EVERYWHERE**: Always specify intent for all dummy arguments.
- **PURE/ELEMENTAL**: Use pure and elemental procedures where possible for better optimization and safety.
- **STRICT TDD**: You MUST always write tests first. No exceptions.

### Testing Strategy
- Strict Test-Driven Development (TDD)
- Target >90% code coverage
- Tests must pass before implementation
- Use codecov for coverage tracking
- CI/CD runs on every push

### Example Development
- Examples use fortplotlib for visualization
- Each example demonstrates specific features
- Examples should be self-contained
- Include mathematical description in comments

### Continuous Integration
- GitHub Actions workflow in `.github/workflows/ci.yml`
- Tests on Ubuntu and macOS with gfortran
- Automatic codecov upload
- Example compilation verification

## Code Maintenance Principles
- Always delete obsolete parts of the code.