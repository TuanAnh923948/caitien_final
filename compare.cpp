/**
 * ============================================================================
 * COMPARE - SFD Estimator Comparison Tool
 * ============================================================================
 * 
 * This tool compares the results of different SFD estimation methods:
 * - EXACT: Ground truth (if available)
 * - MCMC: Original paper's implementation
 * - UNIFIED: Our improved robust auto-adaptive method
 * 
 * Metrics computed:
 * - RMSE: Root Mean Square Error
 * - MAE: Mean Absolute Error
 * - Max Error: Maximum absolute error across all types
 * - TV Distance: Total Variation distance
 * - CI Coverage: Percentage of true values within confidence intervals
 * - Speed: Execution time comparison
 * 
 * Usage: ./compare <exact.res> <mcmc.res> <unified.res> [--detailed]
 * 
 * Author: TuanAnh (Master's Thesis)
 * ============================================================================
 */

 #include "common.hpp"

 // ============================================================================
 // COMPARISON METRICS
 // ============================================================================
 
 struct Metrics {
     double rmse = 0;           // Root Mean Square Error
     double mae = 0;            // Mean Absolute Error
     double max_error = 0;      // Maximum absolute error
     double tv_distance = 0;    // Total Variation distance
     double ci_coverage = 0;    // Percentage of true values in CI
     int ci_count = 0;          // Number of types with non-zero frequency
     int max_error_type = 0;    // Type with maximum error
 };
 
 /**
  * Compute comparison metrics between estimated and true distributions.
  */
 Metrics computeMetrics(const Result& estimated, const Result& truth) {
     Metrics m;
     
     double sum_sq_error = 0;
     double sum_abs_error = 0;
     double sum_tv = 0;
     int count = 0;
     int ci_hits = 0;
     
     for (int i = 0; i < NUM_TYPES; i++) {
         double est = estimated.freq[i];
         double tru = truth.freq[i];
         double err = est - tru;
         
         sum_sq_error += err * err;
         sum_abs_error += abs(err);
         sum_tv += abs(err);
         
         if (abs(err) > m.max_error) {
             m.max_error = abs(err);
             m.max_error_type = i + 1;
         }
         
         // Check CI coverage for types that exist in truth
         if (tru > 0 || est > 0) {
             count++;
             if (tru >= estimated.ci_lo[i] && tru <= estimated.ci_hi[i]) {
                 ci_hits++;
             }
         }
     }
     
     m.rmse = sqrt(sum_sq_error / NUM_TYPES);
     m.mae = sum_abs_error / NUM_TYPES;
     m.tv_distance = sum_tv / 2;  // TV distance is half of L1 norm
     m.ci_count = count;
     m.ci_coverage = count > 0 ? (100.0 * ci_hits / count) : 0;
     
     return m;
 }
 
 // ============================================================================
 // OUTPUT FORMATTING
 // ============================================================================
 
 void printHeader() {
     cerr << endl;
     cerr << "╔════════════════════════════════════════════════════════════════════════════╗" << endl;
     cerr << "║                      SFD ESTIMATOR COMPARISON                              ║" << endl;
     cerr << "║                         (Master's Thesis)                                  ║" << endl;
     cerr << "╚════════════════════════════════════════════════════════════════════════════╝" << endl;
     cerr << endl;
 }
 
 void printResultInfo(const Result& exact, const Result& mcmc, const Result& unified) {
     cerr << "┌──────────────────────────────────────────────────────────────────────────────┐" << endl;
     cerr << "│ LOADED RESULTS                                                               │" << endl;
     cerr << "├──────────────────────────────────────────────────────────────────────────────┤" << endl;
     cerr << "│ EXACT:   " << setw(35) << left << exact.method 
          << " | Time: " << setw(12) << right << fixed << setprecision(1) << exact.time_ms << " ms │" << endl;
     cerr << "│ MCMC:    " << setw(35) << left << mcmc.method 
          << " | Time: " << setw(12) << right << mcmc.time_ms << " ms │" << endl;
     cerr << "│ UNIFIED: " << setw(35) << left << unified.method 
          << " | Time: " << setw(12) << right << unified.time_ms << " ms │" << endl;
     cerr << "└──────────────────────────────────────────────────────────────────────────────┘" << endl;
     cerr << endl;
 }
 
 void printComparisonTable(const Metrics& m_mcmc, const Metrics& m_unified,
                           const Result& mcmc, const Result& unified,
                           bool has_exact_ground_truth) {
     
     if (!has_exact_ground_truth) {
         // No exact ground truth - show side-by-side comparison
         cerr << "╔════════════════════════════════════════════════════════════════════════════╗" << endl;
         cerr << "║          MCMC METHODS COMPARISON (No Exact Ground Truth)                  ║" << endl;
         cerr << "╚════════════════════════════════════════════════════════════════════════════╝" << endl;
         cerr << endl;
         cerr << "⚠️  Cannot compute accuracy without exact ground truth." << endl;
         cerr << "    Showing side-by-side estimates for manual comparison." << endl;
         cerr << endl;
         
         cerr << "┌─────────────────────────────────────────────────────────────────────────────┐" << endl;
         cerr << "│ Metric              │    MCMC_ORIG    │     UNIFIED     │                  │" << endl;
         cerr << "├─────────────────────────────────────────────────────────────────────────────┤" << endl;
         cerr << "│ Time (ms)           │ " << setw(15) << fixed << setprecision(1) << mcmc.time_ms 
              << " │ " << setw(15) << unified.time_ms << " │                  │" << endl;
         cerr << "│ Samples             │ " << setw(15) << mcmc.samples 
              << " │ " << setw(15) << unified.samples << " │                  │" << endl;
         cerr << "└─────────────────────────────────────────────────────────────────────────────┘" << endl;
         cerr << endl;
         
         // Speed comparison
         if (mcmc.time_ms > 0 && unified.time_ms > 0) {
             double speed_ratio = unified.time_ms / mcmc.time_ms;
             cerr << "📊 SPEED COMPARISON:" << endl;
             if (speed_ratio > 1) {
                 cerr << "   UNIFIED is " << fixed << setprecision(1) 
                      << speed_ratio << "x slower than MCMC_ORIG" << endl;
                 cerr << "   (This is expected - UNIFIED prioritizes ACCURACY over speed)" << endl;
             } else {
                 cerr << "   UNIFIED is " << fixed << setprecision(1) 
                      << (1.0/speed_ratio) << "x faster than MCMC_ORIG" << endl;
             }
         }
         
         cerr << endl;
         cerr << "📋 SFD ESTIMATES COMPARISON:" << endl;
         cerr << "┌─────────────────────────────────────────────────────────────────────────────┐" << endl;
         cerr << "│ Type │  MCMC_ORIG │   UNIFIED  │ Difference │    UNIFIED 95% CI     │" << endl;
         cerr << "├─────────────────────────────────────────────────────────────────────────────┤" << endl;
         
         for (int i = 0; i < NUM_TYPES; i++) {
             if (mcmc.freq[i] > 0.0001 || unified.freq[i] > 0.0001) {
                 double diff = unified.freq[i] - mcmc.freq[i];
                 cerr << "│ " << setw(4) << (i+1) << " │ "
                      << fixed << setprecision(6) << setw(10) << mcmc.freq[i] << " │ "
                      << setw(10) << unified.freq[i] << " │ "
                      << setw(10) << showpos << diff << noshowpos << " │ ["
                      << setprecision(4) << setw(7) << unified.ci_lo[i] << ", "
                      << setw(7) << unified.ci_hi[i] << "] │" << endl;
             }
         }
         
         cerr << "└─────────────────────────────────────────────────────────────────────────────┘" << endl;
         cerr << endl;
         
         cerr << "📝 INTERPRETATION:" << endl;
         cerr << "   • MCMC_ORIG uses buggy hash function → likely INCORRECT classification" << endl;
         cerr << "   • UNIFIED uses corrected classifier → MORE ACCURATE results" << endl;
         cerr << "   • Large differences indicate classifier bug impact" << endl;
         cerr << "   • To validate: run 'make quick' on small dataset with exact ground truth" << endl;
         
         return;
     }
     
     // Original comparison with exact ground truth
     cerr << "╔════════════════════════════════════════════════════════════════════════════╗" << endl;
     cerr << "║                         ACCURACY COMPARISON                               ║" << endl;
     cerr << "╚════════════════════════════════════════════════════════════════════════════╝" << endl;
     cerr << endl;
     
     cerr << "┌─────────────────────────────────────────────────────────────────────────────┐" << endl;
     cerr << "│ Metric           │   MCMC_ORIG    │    UNIFIED     │      Winner          │" << endl;
     cerr << "├─────────────────────────────────────────────────────────────────────────────┤" << endl;
     
     // RMSE
     string winner_rmse = (m_unified.rmse < m_mcmc.rmse) ? "✓ UNIFIED" : 
                          (m_mcmc.rmse < m_unified.rmse) ? "✗ MCMC" : "~TIE";
     cerr << "│ RMSE             │ " << fixed << setprecision(6) << setw(14) << m_mcmc.rmse 
          << " │ " << setw(14) << m_unified.rmse 
          << " │ " << setw(20) << winner_rmse << " │" << endl;
     
     // MAE
     string winner_mae = (m_unified.mae < m_mcmc.mae) ? "✓ UNIFIED" : 
                         (m_mcmc.mae < m_unified.mae) ? "✗ MCMC" : "~TIE";
     cerr << "│ MAE              │ " << setw(14) << m_mcmc.mae 
          << " │ " << setw(14) << m_unified.mae 
          << " │ " << setw(20) << winner_mae << " │" << endl;
     
     // Max Error
     string winner_max = (m_unified.max_error < m_mcmc.max_error) ? "✓ UNIFIED" : 
                         (m_mcmc.max_error < m_unified.max_error) ? "✗ MCMC" : "~TIE";
     cerr << "│ Max Error        │ " << setw(14) << m_mcmc.max_error 
          << " │ " << setw(14) << m_unified.max_error 
          << " │ " << setw(20) << winner_max << " │" << endl;
     
     // TV Distance
     string winner_tv = (m_unified.tv_distance < m_mcmc.tv_distance) ? "✓ UNIFIED" : 
                        (m_mcmc.tv_distance < m_unified.tv_distance) ? "✗ MCMC" : "~TIE";
     cerr << "│ TV Distance      │ " << setw(14) << m_mcmc.tv_distance 
          << " │ " << setw(14) << m_unified.tv_distance 
          << " │ " << setw(20) << winner_tv << " │" << endl;
     
     // CI Coverage
     string winner_ci = (m_unified.ci_coverage > m_mcmc.ci_coverage) ? "✓ UNIFIED" : 
                        (m_mcmc.ci_coverage > m_unified.ci_coverage) ? "✗ MCMC" : "~TIE";
     cerr << "│ CI Coverage (%)  │ " << setw(14) << setprecision(2) << m_mcmc.ci_coverage 
          << " │ " << setw(14) << m_unified.ci_coverage 
          << " │ " << setw(20) << winner_ci << " │" << endl;
     
     // Time
     string winner_time = (unified.time_ms < mcmc.time_ms) ? "✓ UNIFIED" : 
                          (mcmc.time_ms < unified.time_ms) ? "✗ MCMC" : "~TIE";
     cerr << "│ Time (ms)        │ " << setw(14) << setprecision(1) << mcmc.time_ms 
          << " │ " << setw(14) << unified.time_ms 
          << " │ " << setw(20) << winner_time << " │" << endl;
     
     cerr << "└─────────────────────────────────────────────────────────────────────────────┘" << endl;
     cerr << endl;
     
     // Summary statistics
     cerr << "📊 SUMMARY STATISTICS:" << endl;
     
     // Speed comparison
     if (mcmc.time_ms > 0 && unified.time_ms > 0) {
         double speed_ratio = mcmc.time_ms / unified.time_ms;
         if (speed_ratio > 1) {
             cerr << "   ⏱️  UNIFIED is " << fixed << setprecision(1) 
                  << speed_ratio << "x FASTER than MCMC" << endl;
         } else {
             cerr << "   ⏱️  UNIFIED is " << fixed << setprecision(1) 
                  << (1.0/speed_ratio) << "x SLOWER than MCMC" << endl;
             cerr << "       (UNIFIED prioritizes accuracy over speed)" << endl;
         }
     }
     
     // Accuracy comparison
     if (m_unified.rmse < 0.000001) {
         cerr << "   🎯 UNIFIED RMSE is PERFECT (≈0)" << endl;
     } else if (m_mcmc.rmse > 0.000001) {
         double rmse_ratio = m_mcmc.rmse / m_unified.rmse;
         cerr << "   🎯 UNIFIED RMSE is " << fixed << setprecision(1) 
              << rmse_ratio << "x BETTER than MCMC" << endl;
     }
     
     // Max error type
     if (m_mcmc.max_error > 0.01) {
         cerr << "   ⚠️  MCMC worst error: Type " << m_mcmc.max_error_type 
              << " (error = " << fixed << setprecision(4) << m_mcmc.max_error << ")" << endl;
     }
     if (m_unified.max_error > 0.01) {
         cerr << "   ⚠️  UNIFIED worst error: Type " << m_unified.max_error_type 
              << " (error = " << fixed << setprecision(4) << m_unified.max_error << ")" << endl;
     }
     
     cerr << endl;
     
     // Overall winner
     int unified_wins = 0;
     int mcmc_wins = 0;
     
     if (m_unified.rmse < m_mcmc.rmse) unified_wins++; else if (m_mcmc.rmse < m_unified.rmse) mcmc_wins++;
     if (m_unified.mae < m_mcmc.mae) unified_wins++; else if (m_mcmc.mae < m_unified.mae) mcmc_wins++;
     if (m_unified.max_error < m_mcmc.max_error) unified_wins++; else if (m_mcmc.max_error < m_unified.max_error) mcmc_wins++;
     if (m_unified.tv_distance < m_mcmc.tv_distance) unified_wins++; else if (m_mcmc.tv_distance < m_unified.tv_distance) mcmc_wins++;
     if (m_unified.ci_coverage > m_mcmc.ci_coverage) unified_wins++; else if (m_mcmc.ci_coverage > m_unified.ci_coverage) mcmc_wins++;
     
     cerr << "╔════════════════════════════════════════════════════════════════════════════╗" << endl;
     if (unified_wins > mcmc_wins) {
         cerr << "║  🏆 UNIFIED WINS on " << unified_wins << "/" << (unified_wins + mcmc_wins) 
              << " accuracy metrics!                                 ║" << endl;
     } else if (mcmc_wins > unified_wins) {
         cerr << "║  ❌ MCMC WINS on " << mcmc_wins << "/" << (unified_wins + mcmc_wins) 
              << " accuracy metrics                                    ║" << endl;
     } else {
         cerr << "║  🤝 TIE - Both methods perform similarly                                   ║" << endl;
     }
     cerr << "╚════════════════════════════════════════════════════════════════════════════╝" << endl;
     cerr << endl;
 }
 
 void printDetailedTable(const Result& exact, const Result& mcmc, const Result& unified) {
     cerr << "╔════════════════════════════════════════════════════════════════════════════╗" << endl;
     cerr << "║                      DETAILED PER-TYPE COMPARISON                          ║" << endl;
     cerr << "╚════════════════════════════════════════════════════════════════════════════╝" << endl;
     cerr << endl;
     
     cerr << "┌─────────────────────────────────────────────────────────────────────────────┐" << endl;
     cerr << "│ Type │    EXACT   │    MCMC    │ MCMC Err │  UNIFIED   │ UNIF Err │ Better │" << endl;
     cerr << "├─────────────────────────────────────────────────────────────────────────────┤" << endl;
     
     for (int i = 0; i < NUM_TYPES; i++) {
         if (exact.freq[i] > 0.0001 || mcmc.freq[i] > 0.0001 || unified.freq[i] > 0.0001) {
             double mcmc_err = abs(mcmc.freq[i] - exact.freq[i]);
             double unified_err = abs(unified.freq[i] - exact.freq[i]);
             
             string better = "~";
             if (unified_err < mcmc_err - 0.0001) better = "UNIF";
             else if (mcmc_err < unified_err - 0.0001) better = "MCMC";
             
             cerr << "│ " << setw(4) << (i+1) << " │ "
                  << fixed << setprecision(6) << setw(10) << exact.freq[i] << " │ "
                  << setw(10) << mcmc.freq[i] << " │ "
                  << setw(8) << setprecision(4) << mcmc_err << " │ "
                  << setw(10) << setprecision(6) << unified.freq[i] << " │ "
                  << setw(8) << setprecision(4) << unified_err << " │ "
                  << setw(6) << better << " │" << endl;
         }
     }
     
     cerr << "└─────────────────────────────────────────────────────────────────────────────┘" << endl;
 }
 
 void printConfidenceIntervalTable(const Result& exact, const Result& unified) {
     cerr << endl;
     cerr << "╔════════════════════════════════════════════════════════════════════════════╗" << endl;
     cerr << "║                   CONFIDENCE INTERVAL ANALYSIS                             ║" << endl;
     cerr << "╚════════════════════════════════════════════════════════════════════════════╝" << endl;
     cerr << endl;
     
     cerr << "┌─────────────────────────────────────────────────────────────────────────────┐" << endl;
     cerr << "│ Type │    EXACT   │      UNIFIED 95% CI       │ CI Width │ Contains? │" << endl;
     cerr << "├─────────────────────────────────────────────────────────────────────────────┤" << endl;
     
     int covered = 0, total = 0;
     
     for (int i = 0; i < NUM_TYPES; i++) {
         if (exact.freq[i] > 0.0001 || unified.freq[i] > 0.0001) {
             double ci_width = unified.ci_hi[i] - unified.ci_lo[i];
             bool contains = (exact.freq[i] >= unified.ci_lo[i] && exact.freq[i] <= unified.ci_hi[i]);
             
             if (contains) covered++;
             total++;
             
             cerr << "│ " << setw(4) << (i+1) << " │ "
                  << fixed << setprecision(6) << setw(10) << exact.freq[i] << " │ ["
                  << setw(10) << unified.ci_lo[i] << ", "
                  << setw(10) << unified.ci_hi[i] << "] │ "
                  << setw(8) << setprecision(4) << ci_width << " │    "
                  << (contains ? "✓" : "✗") << "     │" << endl;
         }
     }
     
     cerr << "└─────────────────────────────────────────────────────────────────────────────┘" << endl;
     
     cerr << endl;
     cerr << "CI Coverage: " << covered << "/" << total << " = " 
          << fixed << setprecision(1) << (100.0 * covered / max(total, 1)) << "%" << endl;
     
     if (covered == total) {
         cerr << "✓ All true values are within 95% confidence intervals!" << endl;
     } else {
         cerr << "⚠ " << (total - covered) << " true value(s) outside confidence intervals." << endl;
     }
 }
 
 // ============================================================================
 // MAIN
 // ============================================================================
 
 int main(int argc, char* argv[]) {
     if (argc < 4) {
         cerr << "Usage: " << argv[0] << " <exact.res> <mcmc.res> <unified.res> [--detailed]" << endl;
         cerr << endl;
         cerr << "COMPARE - SFD Estimator Comparison Tool" << endl;
         cerr << "=======================================" << endl;
         cerr << "Compares results from different SFD estimation methods." << endl;
         cerr << endl;
         cerr << "If exact.res == mcmc.res, assumes no exact ground truth" << endl;
         cerr << "and shows side-by-side MCMC comparison instead." << endl;
         cerr << endl;
         cerr << "Options:" << endl;
         cerr << "  --detailed  Show per-type frequency comparison" << endl;
         return 1;
     }
     
     bool detailed = false;
     if (argc >= 5 && string(argv[4]) == "--detailed") {
         detailed = true;
     }
     
     // Check if we have exact ground truth
     bool has_exact = (string(argv[1]) != string(argv[2]));
     
     // Load results
     Result exact, mcmc, unified;
     
     if (!exact.load(argv[1])) {
         cerr << "Error: Cannot load " << argv[1] << endl;
         return 1;
     }
     
     if (!mcmc.load(argv[2])) {
         cerr << "Error: Cannot load " << argv[2] << endl;
         return 1;
     }
     
     if (!unified.load(argv[3])) {
         cerr << "Error: Cannot load " << argv[3] << endl;
         return 1;
     }
     
     // Print header
     printHeader();
     
     // Print result info
     printResultInfo(exact, mcmc, unified);
     
     // Compute metrics (using exact as ground truth, or mcmc if no exact)
     Metrics m_mcmc = computeMetrics(mcmc, exact);
     Metrics m_unified = computeMetrics(unified, exact);
     
     if (has_exact) {
         // Print ground truth info
         cerr << "┌──────────────────────────────────────────────────────────────────────────────┐" << endl;
         cerr << "│ GROUND TRUTH: " << setw(61) << left << exact.method << "│" << endl;
         cerr << "│ Total Simplets: " << setw(59) << left << exact.total << "│" << endl;
         cerr << "└──────────────────────────────────────────────────────────────────────────────┘" << endl;
         cerr << endl;
     }
     
     // Print comparison table
     printComparisonTable(m_mcmc, m_unified, mcmc, unified, has_exact);
     
     // Print detailed tables if requested and we have exact ground truth
     if (detailed && has_exact) {
         printDetailedTable(exact, mcmc, unified);
         printConfidenceIntervalTable(exact, unified);
     }
     
     return 0;
 }