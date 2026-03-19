# ==============================================================================
# Makefile for SFD Estimator - Paper Original vs Improved V2
# ==============================================================================
#
# Author: TuanAnh (Master's Thesis)
# ==============================================================================

# ==============================================================================
# COMPILER
# ==============================================================================

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native
OPENMP_FLAGS = -fopenmp

# ==============================================================================
# DATASET - MODIFY THESE PATHS
# ==============================================================================

INPUT1 = /mnt/d/LUAN_VAN_THAC_SI/fix/coauth-DBLP-full-nverts.txt
INPUT2 = /mnt/d/LUAN_VAN_THAC_SI/fix/coauth-DBLP-full-simplices.txt

# ==============================================================================
# PARAMETERS
# ==============================================================================

N = 2000
MCMC_SAMPLES = 10000

# ==============================================================================
# OUTPUT
# ==============================================================================

R = results

# ==============================================================================
# BUILD
# ==============================================================================

.PHONY: all clean help

all: mcmc_sfd unified_sfd_v2 exact_sfd compare
	@mkdir -p $(R)
	@echo ""
	@echo "BUILD COMPLETE"
	@echo "  ./mcmc_sfd         Original paper (with bugs)"
	@echo "  ./unified_sfd_v2   Improved V2"
	@echo "  ./exact_sfd        Exact enumeration"
	@echo "  ./compare          Comparison tool"
	@echo ""
	@echo "  Run 'make help' to see all commands."
	@echo ""

mcmc_sfd: mcmc_sfd.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

unified_sfd_v2: unified_sfd_v2.cpp
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $< -o $@

exact_sfd: exact_sfd.cpp common.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

compare: compare.cpp common.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

# ==============================================================================
# TEST DATA (small synthetic)
# ==============================================================================

testdata:
	@mkdir -p testdata
	@printf "2\n2\n2\n2\n2\n3\n3\n3\n" > testdata/nverts.txt
	@printf "0\n1\n1\n2\n2\n3\n3\n4\n0\n4\n0\n1\n2\n2\n3\n4\n0\n3\n4\n" > testdata/simplices.txt
	@echo "Test data created."

# ==============================================================================
# QUICK TEST (small data, exact ground truth)
# ==============================================================================

quick: all testdata
	@mkdir -p $(R)
	@echo ""
	@echo "=== QUICK TEST (small data) ==="
	./exact_sfd testdata/nverts.txt testdata/simplices.txt $(R)/exact.res
	@echo ""
	./mcmc_sfd testdata/nverts.txt testdata/simplices.txt $(R)/mcmc.res 0 5000
	@echo ""
	./unified_sfd_v2 testdata/nverts.txt testdata/simplices.txt $(R)/v2.res 0 --log-convergence
	@echo ""
	./compare $(R)/exact.res $(R)/mcmc.res $(R)/v2.res --detailed

# ==============================================================================
# RUN - MCMC (original paper)
# ==============================================================================

run-mcmc: mcmc_sfd
	@mkdir -p $(R)
	@echo ">>> MCMC original, N=$(N), $(MCMC_SAMPLES) samples..."
	./mcmc_sfd $(INPUT1) $(INPUT2) $(R)/mcmc_n$(N).res $(N) $(MCMC_SAMPLES)

run-mcmc-all: mcmc_sfd
	@mkdir -p $(R)
	@echo ">>> MCMC original, full dataset, $(MCMC_SAMPLES) samples..."
	./mcmc_sfd $(INPUT1) $(INPUT2) $(R)/mcmc_all.res 0 $(MCMC_SAMPLES)

# ==============================================================================
# RUN - V2 Paper Mode (size 2-4 only, fair comparison)
# ==============================================================================

run-v2: unified_sfd_v2
	@mkdir -p $(R)
	@echo ">>> V2 paper mode, N=$(N)..."
	./unified_sfd_v2 $(INPUT1) $(INPUT2) $(R)/v2_n$(N).res $(N)

run-v2-log: unified_sfd_v2
	@mkdir -p $(R)
	@echo ">>> V2 paper mode, N=$(N), +log..."
	./unified_sfd_v2 $(INPUT1) $(INPUT2) $(R)/v2_n$(N).res $(N) --log-convergence

run-v2-all: unified_sfd_v2
	@mkdir -p $(R)
	@echo ">>> V2 paper mode, full dataset..."
	./unified_sfd_v2 $(INPUT1) $(INPUT2) $(R)/v2_all.res 0

run-v2-all-log: unified_sfd_v2
	@mkdir -p $(R)
	@echo ">>> V2 paper mode, full dataset, +log..."
	./unified_sfd_v2 $(INPUT1) $(INPUT2) $(R)/v2_all.res 0 --log-convergence

# ==============================================================================
# RUN - V2 Full Decompose (all simplex sizes)
# ==============================================================================

run-v2-full: unified_sfd_v2
	@mkdir -p $(R)
	@echo ">>> V2 full decompose, N=$(N)..."
	./unified_sfd_v2 $(INPUT1) $(INPUT2) $(R)/v2_full_n$(N).res $(N) --full-decompose

run-v2-full-log: unified_sfd_v2
	@mkdir -p $(R)
	@echo ">>> V2 full decompose, N=$(N), +log..."
	./unified_sfd_v2 $(INPUT1) $(INPUT2) $(R)/v2_full_n$(N).res $(N) --full-decompose --log-convergence

run-v2-full-all: unified_sfd_v2
	@mkdir -p $(R)
	@echo ">>> V2 full decompose, full dataset..."
	./unified_sfd_v2 $(INPUT1) $(INPUT2) $(R)/v2_full_all.res 0 --full-decompose

run-v2-full-all-log: unified_sfd_v2
	@mkdir -p $(R)
	@echo ">>> V2 full decompose, full dataset, +log..."
	./unified_sfd_v2 $(INPUT1) $(INPUT2) $(R)/v2_full_all.res 0 --full-decompose --log-convergence

# ==============================================================================
# RUN - Exact Enumeration (ground truth, small datasets only)
# ==============================================================================

run-exact: exact_sfd
	@mkdir -p $(R)
	@echo ">>> Exact enumeration (may be slow!)..."
	./exact_sfd $(INPUT1) $(INPUT2) $(R)/exact.res

# ==============================================================================
# COMPARE
# ==============================================================================

# MCMC vs V2 paper mode, with N
cmp-n: compare
	@echo ">>> Compare: MCMC(N=$(N)) vs V2-paper(N=$(N))..."
	@[ -f $(R)/mcmc_n$(N).res ] && [ -f $(R)/v2_n$(N).res ] && \
		./compare $(R)/mcmc_n$(N).res $(R)/mcmc_n$(N).res $(R)/v2_n$(N).res --detailed || \
		echo "Missing files. Run: make run-mcmc run-v2"

# MCMC vs V2 paper mode, full dataset
cmp-all: compare
	@echo ">>> Compare: MCMC(ALL) vs V2-paper(ALL)..."
	@[ -f $(R)/mcmc_all.res ] && [ -f $(R)/v2_all.res ] && \
		./compare $(R)/mcmc_all.res $(R)/mcmc_all.res $(R)/v2_all.res --detailed || \
		echo "Missing files. Run: make run-mcmc-all run-v2-all"

# MCMC vs V2 full decompose, with N
cmp-full-n: compare
	@echo ">>> Compare: MCMC(N=$(N)) vs V2-full(N=$(N))..."
	@[ -f $(R)/mcmc_n$(N).res ] && [ -f $(R)/v2_full_n$(N).res ] && \
		./compare $(R)/mcmc_n$(N).res $(R)/mcmc_n$(N).res $(R)/v2_full_n$(N).res --detailed || \
		echo "Missing files. Run: make run-mcmc run-v2-full"

# MCMC vs V2 full decompose, full dataset
cmp-full-all: compare
	@echo ">>> Compare: MCMC(ALL) vs V2-full(ALL)..."
	@[ -f $(R)/mcmc_all.res ] && [ -f $(R)/v2_full_all.res ] && \
		./compare $(R)/mcmc_all.res $(R)/mcmc_all.res $(R)/v2_full_all.res --detailed || \
		echo "Missing files. Run: make run-mcmc-all run-v2-full-all"

# V2 paper vs V2 full, same N
cmp-paper-vs-full: compare
	@echo ">>> Compare: V2-paper(N=$(N)) vs V2-full(N=$(N))..."
	@[ -f $(R)/v2_n$(N).res ] && [ -f $(R)/v2_full_n$(N).res ] && \
		./compare $(R)/v2_n$(N).res $(R)/v2_n$(N).res $(R)/v2_full_n$(N).res --detailed || \
		echo "Missing files. Run: make run-v2 run-v2-full"

# V2 paper vs V2 full, full dataset
cmp-paper-vs-full-all: compare
	@echo ">>> Compare: V2-paper(ALL) vs V2-full(ALL)..."
	@[ -f $(R)/v2_all.res ] && [ -f $(R)/v2_full_all.res ] && \
		./compare $(R)/v2_all.res $(R)/v2_all.res $(R)/v2_full_all.res --detailed || \
		echo "Missing files. Run: make run-v2-all run-v2-full-all"

# With exact ground truth
cmp-exact: compare
	@echo ">>> Compare against exact ground truth..."
	@[ -f $(R)/exact.res ] && [ -f $(R)/mcmc_n$(N).res ] && [ -f $(R)/v2_n$(N).res ] && \
		./compare $(R)/exact.res $(R)/mcmc_n$(N).res $(R)/v2_n$(N).res --detailed || \
		echo "Missing files. Run: make run-exact run-mcmc run-v2"

# ==============================================================================
# COMBINED WORKFLOWS
# ==============================================================================

# Fair comparison: MCMC vs V2, paper mode, N vertices
test: all run-mcmc run-v2 cmp-n
	@echo "DONE: MCMC vs V2, N=$(N), paper mode"

# Same + convergence log
test-log: all run-mcmc run-v2-log cmp-n
	@echo "DONE: MCMC vs V2 +log, N=$(N)"

# Full dataset comparison
test-all: all run-mcmc-all run-v2-all cmp-all
	@echo "DONE: MCMC vs V2, full dataset"

# Full dataset + log
test-all-log: all run-mcmc-all run-v2-all-log cmp-all
	@echo "DONE: MCMC vs V2 +log, full dataset"

# Full decompose comparison
test-full: all run-mcmc run-v2-full cmp-full-n
	@echo "DONE: MCMC vs V2-full, N=$(N)"

# Full decompose, full dataset
test-full-all: all run-mcmc-all run-v2-full-all cmp-full-all
	@echo "DONE: MCMC vs V2-full, full dataset"

# Paper vs Full decompose (V2 only)
test-modes: all run-v2 run-v2-full cmp-paper-vs-full
	@echo "DONE: V2-paper vs V2-full, N=$(N)"

# Everything on subgraph
test-everything: all run-mcmc run-v2 run-v2-full cmp-n cmp-full-n cmp-paper-vs-full
	@echo ""
	@echo "ALL DONE (N=$(N)). Files:"
	@echo "  $(R)/mcmc_n$(N).res"
	@echo "  $(R)/v2_n$(N).res"
	@echo "  $(R)/v2_full_n$(N).res"

# ==============================================================================
# BENCHMARK
# ==============================================================================

benchmark: mcmc_sfd unified_sfd_v2
	@mkdir -p $(R)/benchmark
	@echo "=== BENCHMARK: N=$(N) ==="
	@echo ">>> MCMC..."
	@time ./mcmc_sfd $(INPUT1) $(INPUT2) $(R)/benchmark/mcmc.res $(N) $(MCMC_SAMPLES) 2>&1 | tail -5
	@echo ""
	@echo ">>> V2 paper..."
	@time ./unified_sfd_v2 $(INPUT1) $(INPUT2) $(R)/benchmark/v2.res $(N) 2>&1 | tail -5
	@echo ""
	@echo ">>> V2 full..."
	@time ./unified_sfd_v2 $(INPUT1) $(INPUT2) $(R)/benchmark/v2_full.res $(N) --full-decompose 2>&1 | tail -5

# ==============================================================================
# CLEAN
# ==============================================================================

clean:
	rm -f mcmc_sfd unified_sfd_v2 exact_sfd compare
	rm -rf $(R) testdata

# ==============================================================================
# HELP
# ==============================================================================

help:
	@echo ""
	@echo "================================================================"
	@echo "  SFD Estimator - All Commands"
	@echo "================================================================"
	@echo ""
	@echo "BUILD:"
	@echo "  make all                       Build everything"
	@echo "  make clean                     Remove all artifacts"
	@echo ""
	@echo "QUICK TEST (small data):"
	@echo "  make quick                     exact + mcmc + v2 + compare"
	@echo ""
	@echo "RUN MCMC (original paper):"
	@echo "  make run-mcmc                  N=$(N)"
	@echo "  make run-mcmc-all              full dataset"
	@echo ""
	@echo "RUN V2 - PAPER MODE (size 2-4, fair comparison):"
	@echo "  make run-v2                    N=$(N)"
	@echo "  make run-v2-log                N=$(N) + convergence log"
	@echo "  make run-v2-all                full dataset"
	@echo "  make run-v2-all-log            full dataset + convergence log"
	@echo ""
	@echo "RUN V2 - FULL DECOMPOSE (all simplex sizes):"
	@echo "  make run-v2-full               N=$(N)"
	@echo "  make run-v2-full-log           N=$(N) + convergence log"
	@echo "  make run-v2-full-all           full dataset"
	@echo "  make run-v2-full-all-log       full dataset + convergence log"
	@echo ""
	@echo "RUN EXACT:"
	@echo "  make run-exact                 ground truth (slow!)"
	@echo ""
	@echo "COMPARE:"
	@echo "  make cmp-n                     MCMC vs V2-paper, N=$(N)"
	@echo "  make cmp-all                   MCMC vs V2-paper, full dataset"
	@echo "  make cmp-full-n                MCMC vs V2-full, N=$(N)"
	@echo "  make cmp-full-all              MCMC vs V2-full, full dataset"
	@echo "  make cmp-paper-vs-full         V2-paper vs V2-full, N=$(N)"
	@echo "  make cmp-paper-vs-full-all     V2-paper vs V2-full, full"
	@echo "  make cmp-exact                 all vs exact ground truth"
	@echo ""
	@echo "COMBINED WORKFLOWS:"
	@echo "  make test                      mcmc + v2 + cmp, N=$(N)"
	@echo "  make test-log                  same + convergence log"
	@echo "  make test-all                  mcmc + v2 + cmp, full dataset"
	@echo "  make test-all-log              same + convergence log"
	@echo "  make test-full                 mcmc + v2-full + cmp, N=$(N)"
	@echo "  make test-full-all             mcmc + v2-full + cmp, full"
	@echo "  make test-modes                v2-paper vs v2-full, N=$(N)"
	@echo "  make test-everything           all methods + all cmp, N=$(N)"
	@echo ""
	@echo "BENCHMARK:"
	@echo "  make benchmark                 speed comparison"
	@echo ""
	@echo "PARAMETERS (edit Makefile or override on command line):"
	@echo "  N=$(N)                 connected vertices (0=all)"
	@echo "  MCMC_SAMPLES=$(MCMC_SAMPLES)      samples for MCMC"
	@echo ""
	@echo "  Example: make test N=5000 MCMC_SAMPLES=20000"
	@echo ""
	@echo "OUTPUT FILES (in $(R)/):"
	@echo "  mcmc_n$(N).res             MCMC, N=$(N)"
	@echo "  mcmc_all.res               MCMC, full dataset"
	@echo "  v2_n$(N).res               V2 paper, N=$(N)"
	@echo "  v2_all.res                 V2 paper, full dataset"
	@echo "  v2_full_n$(N).res          V2 full decompose, N=$(N)"
	@echo "  v2_full_all.res            V2 full decompose, full"
	@echo "  exact.res                  exact ground truth"
	@echo ""