# FortFEM Benchmark Comparison System
# Compare FortFEM implementations against FreeFEM reference

.PHONY: all poisson-1d poisson-2d curl-curl clean

# Default target
all: poisson-1d poisson-2d curl-curl

# Poisson 1D benchmark comparison
poisson-1d:
	@echo "=== Poisson 1D Benchmark Comparison ==="
	@echo "Running FortFEM Poisson 1D benchmark..."
	@cd .. && fpm test test_poisson_1d_benchmark > benchmark/fortfem_poisson_1d_results.txt 2>&1
	@echo "Running FreeFEM Poisson 1D benchmark..."
	@cd freefem/poisson-1d && FreeFem++ poisson_1d.edp > freefem_output.log 2>&1
	@echo ""
	@echo "FortFEM Results:"
	@grep -A 20 "h.*L2_error.*H1_error.*DOFs" benchmark/fortfem_poisson_1d_results.txt || echo "No results found"
	@echo ""
	@echo "FreeFEM Results:"
	@cd freefem/poisson-1d && grep -v "^#" freefem_convergence_1d.dat | awk '{printf "%.6f    %.6E    %.6E    %d\n", $$1, $$2, $$3, $$4}' || echo "No results found"

# Poisson 2D benchmark comparison  
poisson-2d:
	@echo "=== Poisson 2D Benchmark Comparison ==="
	@echo "Running FortFEM Poisson 2D benchmark..."
	@cd .. && fpm test test_poisson_2d_benchmark > benchmark/fortfem_poisson_2d_results.txt 2>&1
	@echo "Running FreeFEM Poisson 2D benchmark..."
	@cd freefem/poisson-2d && FreeFem++ poisson_2d.edp > freefem_output.log 2>&1 || true
	@echo ""
	@echo "FortFEM Results:"
	@grep -A 20 "h.*L2_error.*H1_error.*DOFs" fortfem_poisson_2d_results.txt || echo "No results found"
	@echo ""
	@echo "FreeFEM Results:"
	@cd freefem/poisson-2d && grep -v "^#" freefem_convergence_2d.dat | awk '{printf "%.6f    %.6E    %.6E    %d\n", $$1, $$2, $$3, $$4}' || echo "No results found"

# Curl-curl benchmark comparison
curl-curl:
	@echo "=== Curl-Curl Benchmark Comparison ==="
	@echo "Running FortFEM curl-curl benchmark..."
	@cd .. && fpm test test_curl_curl_system_solver > benchmark/fortfem_curl_curl_results.txt 2>&1
	@echo "Running FreeFEM curl-curl benchmark..."
	@cd freefem/curl-curl && FreeFem++ nedelec_benchmark.edp > freefem_output.log 2>&1
	@echo ""
	@echo "FortFEM Results:"
	@grep -A 20 "h.*L2_error.*Hcurl_error.*DOFs" benchmark/fortfem_curl_curl_results.txt || echo "No results found"
	@echo ""
	@echo "FreeFEM Results:"
	@cd freefem/curl-curl && grep -v "^#" freefem_convergence.dat | awk '{printf "%.6f    %.6E    %.6E    %d\n", $$1, $$2, $$3, $$4}' || echo "No results found"

# Clean generated files
clean:
	@echo "Cleaning benchmark output files..."
	@rm -f fortfem_*_results.txt
	@rm -f freefem/*/freefem_output.log
	@rm -f freefem/*/freefem_convergence*.dat
	@echo "Cleanup complete."

# Help target
help:
	@echo "FortFEM Benchmark System"
	@echo "========================"
	@echo "Available targets:"
	@echo "  all        - Run all benchmarks"
	@echo "  poisson-1d - Compare Poisson 1D implementations"
	@echo "  poisson-2d - Compare Poisson 2D implementations"
	@echo "  curl-curl  - Compare curl-curl implementations"
	@echo "  clean      - Remove generated files"
	@echo "  help       - Show this help message"