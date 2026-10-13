# FortFEM Makefile
# Main targets for building, testing, and benchmarking

.PHONY: all build test test-slow test-all benchmark clean help doc

# Default target
all: build

# Build the project
build:
	@echo "Building FortFEM..."
	fpm build

# Test targets follow the CI split: `*_slow` programs and the `run_*` drivers
# are not part of the standard suite. `fpm test` starts every test program at
# once, so several hundred OpenMP programs oversubscribe the machine and time
# out although each runs in about a second alone; scripts/run_tests.sh builds
# with fpm and runs the programs with bounded parallelism
# (TEST_JOBS, OMP_NUM_THREADS, TEST_TIMEOUT).

# Run the standard (fast) tests
test:
	@echo "Running FortFEM tests..."
	scripts/run_tests.sh fast

# Run only the slow tests
test-slow:
	@echo "Running FortFEM slow tests..."
	scripts/run_tests.sh slow

# Run the fast and the slow tests
test-all:
	scripts/run_tests.sh all

# Run benchmark comparison against FreeFEM
benchmark: build benchmark-fortfem benchmark-freefem compare-results

# Run FortFEM curl-curl benchmark
benchmark-fortfem:
	@echo "Running FortFEM curl-curl benchmark..."
	@fpm test test_curl_curl_system_solver > fortfem_benchmark.log 2>&1 || true
	@echo "FortFEM benchmark completed. Results saved to fortfem_benchmark.log"

# Run FreeFEM benchmark
benchmark-freefem:
	@echo "Running FreeFEM benchmark..."
	@cd benchmark/freefem/curl-curl && make run > ../../../freefem_benchmark.log 2>&1
	@echo "FreeFEM benchmark completed. Results saved to freefem_benchmark.log"

# Compare benchmark results
compare-results:
	@echo ""
	@echo "=== BENCHMARK COMPARISON ==="
	@echo ""
	@echo "FreeFEM Results:"
	@echo "=================="
	@grep -A 10 "=== FreeFEM Convergence Study ===" freefem_benchmark.log || echo "FreeFEM convergence data not found"
	@echo ""
	@echo "FortFEM Results:"
	@echo "=================="
	@grep -A 10 "h        L2_error        H(curl)_error   DOFs" fortfem_benchmark.log || echo "FortFEM convergence data not found"
	@echo ""
	@echo "=== END COMPARISON ==="

# Build documentation with FORD (using fpm.toml configuration)
doc:
	ford README.md
	# Copy artifacts (plots) if they exist - must be done after FORD runs
	if [ -d artifacts/plots ]; then \
		mkdir -p build/doc/artifacts/plots; \
		cp -r artifacts/plots/* build/doc/artifacts/plots/ 2>/dev/null || true; \
	fi

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	fpm clean
	@rm -f *.log
	@rm -f benchmark/freefem/curl-curl/*.dat
	@rm -rf build/doc

# Show help
help:
	@echo "FortFEM Makefile"
	@echo "=================="
	@echo ""
	@echo "Available targets:"
	@echo "  build      - Build the FortFEM project"
	@echo "  test       - Run the standard tests (everything except *_slow)"
	@echo "  test-slow  - Run the *_slow tests"
	@echo "  test-all   - Run the fast and the slow tests"
	@echo "  benchmark  - Run curl-curl benchmark comparison (FortFEM vs FreeFEM)"
	@echo "  doc        - Build documentation with FORD"
	@echo "  clean      - Clean build artifacts and logs"
	@echo "  help       - Show this help message"
	@echo ""
	@echo "Benchmark comparison:"
	@echo "  This will run both FortFEM and FreeFEM curl-curl solvers"
	@echo "  with identical parameters and compare the convergence results."
