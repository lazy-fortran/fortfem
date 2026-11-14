# Contributing

FortFEM accepts small, independently validated changes. [AGENTS.md](AGENTS.md)
lists the build and test commands and the conventions of the code.

1. Pick one narrowly scoped issue.
2. Make the smallest complete change.
3. Add a test with an independent oracle: a manufactured or exact solution, a
   published value or a separate implementation.
4. Regenerate generated sources with `tools/codegen/generate.sh` when their
   generator changes; the check `tools/codegen/check_generated.sh` must pass.
5. Run `make test` before proposing the change.

Pointwise kernels are generated from symbolic definitions with
[fortsym](https://github.com/lazy-fortran/fortsym); do not edit files under
`src/generated/` by hand.
