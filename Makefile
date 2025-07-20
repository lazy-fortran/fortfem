# FortFEM Makefile
# Main targets for building, testing, and benchmarking

.PHONY: all build test benchmark clean help doc

# Default target
all: build

# Build the project
build:
	@echo "Building FortFEM..."
	fpm build

# Run all tests
test: build
	@echo "Running FortFEM tests..."
	fpm test

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
	# Copy generated example documentation to FORD output
	if [ -d doc/examples ]; then \
		mkdir -p build/doc/page/examples; \
		cp -r doc/examples/* build/doc/page/examples/ 2>/dev/null || true; \
	fi
	# Copy example media files to doc build directory for proper linking
	mkdir -p build/doc/example
	if [ -d build/example ]; then cp -r build/example/* build/doc/example/ 2>/dev/null || true; fi
	# Copy artifacts (plots) if they exist
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
	@echo "  test       - Run all tests"
	@echo "  benchmark  - Run curl-curl benchmark comparison (FortFEM vs FreeFEM)"
	@echo "  doc        - Build documentation with FORD"
	@echo "  clean      - Clean build artifacts and logs"
	@echo "  help       - Show this help message"
	@echo ""
	@echo "Benchmark comparison:"
	@echo "  This will run both FortFEM and FreeFEM curl-curl solvers"
	@echo "  with identical parameters and compare the convergence results."
