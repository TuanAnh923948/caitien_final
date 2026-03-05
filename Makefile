# ==============================================================================
# Makefile for SFD Estimator Project - OPTIMIZED VERSION 2.0
# ==============================================================================
#
# This project implements and compares SFD estimation methods:
# - exact_sfd: Exact enumeration (ground truth)
# - mcmc_sfd: Original MCMC from paper (with bug)
# - unified_sfd: Robust adaptive MCMC (original)
# - unified_sfd_v2: OPTIMIZED with OpenMP, Caching, Early Stopping
#
# Author: TuanAnh (Master's Thesis)
# ==============================================================================

# ==============================================================================
# COMPILER SETTINGS
# ==============================================================================

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native

# OpenMP flags
OPENMP_FLAGS = -fopenmp

# ==============================================================================
# INPUT FILES - MODIFY THESE PATHS FOR YOUR DATASET
# ==============================================================================

INPUT1 = /mnt/d/LUAN_VAN_THAC_SI/fix/coauth-DBLP-nverts.txt
INPUT2 = /mnt/d/LUAN_VAN_THAC_SI/fix/coauth-DBLP-simplices.txt

# ==============================================================================
# PARAMETERS
# ==============================================================================

MCMC_SAMPLES = 10000
MCMC_BURNIN = 1000
UNIFIED_SAMPLES = 0

# ==============================================================================
# OUTPUT
# ==============================================================================

RESULTS_DIR = results

# ==============================================================================
# MAIN TARGETS
# ==============================================================================

.PHONY: all clean test help v2

# Build all executables (including v2)
all: exact_sfd mcmc_sfd unified_sfd unified_sfd_v2 compare
	@mkdir -p $(RESULTS_DIR)
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════════════╗"
	@echo "║                    BUILD COMPLETE                                ║"
	@echo "╚══════════════════════════════════════════════════════════════════╝"
	@echo ""
	@echo "Executables:"
	@echo "  ./exact_sfd      - Exact enumeration (ground truth)"
	@echo "  ./mcmc_sfd       - Original MCMC (paper)"
	@echo "  ./unified_sfd    - Robust adaptive MCMC"
	@echo "  ./unified_sfd_v2 - ⭐ OPTIMIZED: OpenMP + Caching + Early Stop"
	@echo "  ./compare        - Comparison tool"
	@echo ""

exact_sfd: exact_sfd.cpp common.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

mcmc_sfd: mcmc_sfd.cpp common.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

unified_sfd: unified_sfd.cpp common.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

# V2 with OpenMP
unified_sfd_v2: unified_sfd_v2.cpp
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $< -o $@

compare: compare.cpp common.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

# ==============================================================================
# OPTIMIZED TEST (RECOMMENDED)
# ==============================================================================

# ⭐ RECOMMENDED: Run optimized v2 with visualization
v2: all run-mcmc run-v2-log compare-mcmc-v2
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════════════╗"
	@echo "║       ✅ OPTIMIZED TEST COMPLETE (v2)                            ║"
	@echo "╚══════════════════════════════════════════════════════════════════╝"
	@echo ""
	@echo "Results: $(RESULTS_DIR)/"
	@echo "Convergence log: $(RESULTS_DIR)/unified_v2.res.convergence.csv"

# Run v2 without logging
run-v2: unified_sfd_v2
	@mkdir -p $(RESULTS_DIR)
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════════════╗"
	@echo "║          RUNNING UNIFIED V2 (OPTIMIZED)                          ║"
	@echo "╚══════════════════════════════════════════════════════════════════╝"
	./unified_sfd_v2 $(INPUT1) $(INPUT2) $(RESULTS_DIR)/unified_v2.res

# Run v2 with convergence logging
run-v2-log: unified_sfd_v2
	@mkdir -p $(RESULTS_DIR)
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════════════╗"
	@echo "║          RUNNING UNIFIED V2 WITH LOGGING                         ║"
	@echo "╚══════════════════════════════════════════════════════════════════╝"
	./unified_sfd_v2 $(INPUT1) $(INPUT2) $(RESULTS_DIR)/unified_v2.res --log-convergence

# Compare MCMC vs V2
compare-mcmc-v2: compare
	@echo ""
	@echo ">>> Comparing MCMC_ORIG vs UNIFIED_V2..."
	@if [ -f $(RESULTS_DIR)/mcmc.res ] && [ -f $(RESULTS_DIR)/unified_v2.res ]; then \
		./compare $(RESULTS_DIR)/mcmc.res $(RESULTS_DIR)/mcmc.res $(RESULTS_DIR)/unified_v2.res --detailed; \
	else \
		echo "Error: Result files not found."; \
	fi

# ==============================================================================
# LEGACY TESTS (original unified_sfd)
# ==============================================================================

test-with-viz: all run-mcmc run-unified-log compare-mcmc
	@echo "✅ Legacy test complete"

run-mcmc: mcmc_sfd
	@mkdir -p $(RESULTS_DIR)
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════════════╗"
	@echo "║              RUNNING MCMC (ORIGINAL PAPER)                       ║"
	@echo "╚══════════════════════════════════════════════════════════════════╝"
	./mcmc_sfd $(INPUT1) $(INPUT2) $(RESULTS_DIR)/mcmc.res $(MCMC_SAMPLES) $(MCMC_BURNIN)

run-unified-log: unified_sfd
	@mkdir -p $(RESULTS_DIR)
	./unified_sfd $(INPUT1) $(INPUT2) $(RESULTS_DIR)/unified.res --log-convergence

compare-mcmc: compare
	@if [ -f $(RESULTS_DIR)/mcmc.res ] && [ -f $(RESULTS_DIR)/unified.res ]; then \
		./compare $(RESULTS_DIR)/mcmc.res $(RESULTS_DIR)/mcmc.res $(RESULTS_DIR)/unified.res --detailed; \
	fi

# ==============================================================================
# QUICK TEST (small data)
# ==============================================================================

testdata:
	@mkdir -p testdata
	@printf "2\n2\n2\n2\n2\n3\n3\n3\n" > testdata/nverts.txt
	@printf "0\n1\n1\n2\n2\n3\n3\n4\n0\n4\n0\n1\n2\n2\n3\n4\n0\n3\n4\n" > testdata/simplices.txt
	@echo "Test data created."

quick-v2: unified_sfd_v2 mcmc_sfd compare testdata
	@mkdir -p $(RESULTS_DIR)
	@echo "Quick test with v2..."
	./mcmc_sfd testdata/nverts.txt testdata/simplices.txt $(RESULTS_DIR)/mcmc.res 5000 500
	./unified_sfd_v2 testdata/nverts.txt testdata/simplices.txt $(RESULTS_DIR)/unified_v2.res --log-convergence
	./compare $(RESULTS_DIR)/mcmc.res $(RESULTS_DIR)/mcmc.res $(RESULTS_DIR)/unified_v2.res

# ==============================================================================
# BENCHMARK: Compare v1 vs v2
# ==============================================================================

benchmark: all
	@mkdir -p $(RESULTS_DIR)/benchmark
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════════════╗"
	@echo "║              BENCHMARK: V1 vs V2                                 ║"
	@echo "╚══════════════════════════════════════════════════════════════════╝"
	@echo ""
	@echo "Running V1 (original)..."
	@time ./unified_sfd $(INPUT1) $(INPUT2) $(RESULTS_DIR)/benchmark/v1.res 2>&1 | tail -20
	@echo ""
	@echo "Running V2 (optimized)..."
	@time ./unified_sfd_v2 $(INPUT1) $(INPUT2) $(RESULTS_DIR)/benchmark/v2.res 2>&1 | tail -20
	@echo ""
	@echo "Benchmark complete. Compare results in $(RESULTS_DIR)/benchmark/"

# ==============================================================================
# UTILITY
# ==============================================================================

clean:
	rm -f exact_sfd mcmc_sfd unified_sfd unified_sfd_v2 compare
	rm -rf $(RESULTS_DIR) testdata

# ==============================================================================
# HELP
# ==============================================================================

help:
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════════════╗"
	@echo "║           SFD Estimator - Makefile Help                          ║"
	@echo "╚══════════════════════════════════════════════════════════════════╝"
	@echo ""
	@echo "⭐ RECOMMENDED (Optimized V2):"
	@echo "  make v2              Run optimized version with logging"
	@echo "  make run-v2          Run optimized version (no logging)"
	@echo "  make quick-v2        Quick test with small data"
	@echo "  make benchmark       Compare V1 vs V2 speed"
	@echo ""
	@echo "LEGACY (Original unified_sfd):"
	@echo "  make test-with-viz   Run original with logging"
	@echo ""
	@echo "OTHER:"
	@echo "  make all             Build all executables"
	@echo "  make clean           Remove build artifacts"
	@echo "  make help            Show this help"
	@echo ""
	@echo "V2 OPTIMIZATIONS:"
	@echo "  1. OpenMP parallel chain execution"
	@echo "  2. O(1) hash-based lookups (unordered_set)"
	@echo "  3. LRU cache for simplex neighbors"
	@echo "  4. Smart early stopping (CI + stability)"
	@echo ""
	@echo "EXPECTED SPEEDUP: 4-10x faster than V1"
	@echo ""