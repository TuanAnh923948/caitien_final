/**
 * COMPARE - SFD Estimator Comparison Tool
 * Usage: ./compare <exact.res> <mcmc.res> <v2.res> [--detailed]
 * If exact.res == mcmc.res, assumes no exact ground truth.
 */

 #include "common.hpp"

 struct Metrics {
     double rmse = 0;
     double mae = 0;
     double max_error = 0;
     double tv_distance = 0;
     double ci_coverage = 0;
     int ci_count = 0;
     int max_error_type = 0;
 };
 
 Metrics computeMetrics(const Result& estimated, const Result& truth) {
     Metrics m;
     double sum_sq = 0, sum_abs = 0, sum_tv = 0;
     int count = 0, ci_hits = 0;
     
     for (int i = 0; i < NUM_TYPES; i++) {
         double est = estimated.freq[i];
         double tru = truth.freq[i];
         double err = est - tru;
         
         sum_sq += err * err;
         sum_abs += abs(err);
         sum_tv += abs(err);
         
         if (abs(err) > m.max_error) {
             m.max_error = abs(err);
             m.max_error_type = i + 1;
         }
         
         if (tru > 0 || est > 0) {
             count++;
             if (tru >= estimated.ci_lo[i] && tru <= estimated.ci_hi[i]) {
                 ci_hits++;
             }
         }
     }
     
     m.rmse = sqrt(sum_sq / NUM_TYPES);
     m.mae = sum_abs / NUM_TYPES;
     m.tv_distance = sum_tv / 2;
     m.ci_count = count;
     m.ci_coverage = count > 0 ? (100.0 * ci_hits / count) : 0;
     
     return m;
 }
 
 void printComparisonNoExact(const Result& mcmc, const Result& v2) {
     cerr << endl;
     cerr << "╔════════════════════════════════════════════════════════════════╗" << endl;
     cerr << "║      MCMC ORIGINAL vs UNIFIED V2 (No Exact Ground Truth)     ║" << endl;
     cerr << "╚════════════════════════════════════════════════════════════════╝" << endl;
     cerr << endl;
     
     cerr << "┌───────────────────────────────────────────────────────────────┐" << endl;
     cerr << "│ Metric         │   MCMC_ORIG    │      V2        │          │" << endl;
     cerr << "├───────────────────────────────────────────────────────────────┤" << endl;
     cerr << "│ Time (ms)      │ " << setw(14) << fixed << setprecision(1) << mcmc.time_ms
          << " │ " << setw(14) << v2.time_ms << " │          │" << endl;
     cerr << "│ Samples        │ " << setw(14) << mcmc.samples
          << " │ " << setw(14) << v2.samples << " │          │" << endl;
     cerr << "└───────────────────────────────────────────────────────────────┘" << endl;
     
     if (mcmc.time_ms > 0 && v2.time_ms > 0) {
         double ratio = v2.time_ms / mcmc.time_ms;
         cerr << endl;
         if (ratio > 1)
             cerr << "V2 is " << fixed << setprecision(1) << ratio << "x slower (prioritizes accuracy)" << endl;
         else
             cerr << "V2 is " << fixed << setprecision(1) << (1.0/ratio) << "x faster" << endl;
     }
     
     cerr << endl;
     cerr << "SFD Estimates:" << endl;
     cerr << "┌───────────────────────────────────────────────────────────────────────┐" << endl;
     cerr << "│ Type │ MCMC_ORIG  │     V2     │ Difference │    V2 95% CI          │" << endl;
     cerr << "├───────────────────────────────────────────────────────────────────────┤" << endl;
     
     for (int i = 0; i < NUM_TYPES; i++) {
         if (mcmc.freq[i] > 0.0001 || v2.freq[i] > 0.0001) {
             double diff = v2.freq[i] - mcmc.freq[i];
             cerr << "│ " << setw(4) << (i+1) << " │ "
                  << fixed << setprecision(6) << setw(10) << mcmc.freq[i] << " │ "
                  << setw(10) << v2.freq[i] << " │ "
                  << setw(10) << showpos << diff << noshowpos << " │ ["
                  << setprecision(4) << setw(7) << v2.ci_lo[i] << ", "
                  << setw(7) << v2.ci_hi[i] << "] │" << endl;
         }
     }
     
     cerr << "└───────────────────────────────────────────────────────────────────────┘" << endl;
     cerr << endl;
     cerr << "Key differences:" << endl;
     cerr << "  MCMC_ORIG: buggy classifier + acceptance always 1" << endl;
     cerr << "  V2:        corrected classifier + proper M-H + auto-adaptive" << endl;
     cerr << "  Large differences indicate classifier bug impact" << endl;
 }
 
 void printComparisonWithExact(const Metrics& m_mcmc, const Metrics& m_v2,
                               const Result& exact, const Result& mcmc, const Result& v2,
                               bool detailed) {
     cerr << endl;
     cerr << "╔════════════════════════════════════════════════════════════════╗" << endl;
     cerr << "║                   ACCURACY COMPARISON                         ║" << endl;
     cerr << "╚════════════════════════════════════════════════════════════════╝" << endl;
     cerr << endl;
     
     cerr << "Ground truth: " << exact.total << " total simplets" << endl;
     cerr << endl;
     
     cerr << "┌───────────────────────────────────────────────────────────────┐" << endl;
     cerr << "│ Metric          │  MCMC_ORIG     │      V2        │ Winner  │" << endl;
     cerr << "├───────────────────────────────────────────────────────────────┤" << endl;
     
     auto winner = [](double a, double b, bool lower_better) -> string {
         if (lower_better) return (b < a) ? "V2" : (a < b) ? "MCMC" : "TIE";
         else return (b > a) ? "V2" : (a > b) ? "MCMC" : "TIE";
     };
     
     cerr << "│ RMSE            │ " << fixed << setprecision(6) << setw(14) << m_mcmc.rmse
          << " │ " << setw(14) << m_v2.rmse
          << " │ " << setw(7) << winner(m_mcmc.rmse, m_v2.rmse, true) << " │" << endl;
     
     cerr << "│ MAE             │ " << setw(14) << m_mcmc.mae
          << " │ " << setw(14) << m_v2.mae
          << " │ " << setw(7) << winner(m_mcmc.mae, m_v2.mae, true) << " │" << endl;
     
     cerr << "│ Max Error       │ " << setw(14) << m_mcmc.max_error
          << " │ " << setw(14) << m_v2.max_error
          << " │ " << setw(7) << winner(m_mcmc.max_error, m_v2.max_error, true) << " │" << endl;
     
     cerr << "│ TV Distance     │ " << setw(14) << m_mcmc.tv_distance
          << " │ " << setw(14) << m_v2.tv_distance
          << " │ " << setw(7) << winner(m_mcmc.tv_distance, m_v2.tv_distance, true) << " │" << endl;
     
     cerr << "│ CI Coverage (%) │ " << setw(14) << setprecision(1) << m_mcmc.ci_coverage
          << " │ " << setw(14) << m_v2.ci_coverage
          << " │ " << setw(7) << winner(m_mcmc.ci_coverage, m_v2.ci_coverage, false) << " │" << endl;
     
     cerr << "│ Time (ms)       │ " << setw(14) << setprecision(1) << mcmc.time_ms
          << " │ " << setw(14) << v2.time_ms
          << " │ " << setw(7) << winner(mcmc.time_ms, v2.time_ms, true) << " │" << endl;
     
     cerr << "└───────────────────────────────────────────────────────────────┘" << endl;
     
     // Score
     int v2_wins = 0, mcmc_wins = 0;
     if (m_v2.rmse < m_mcmc.rmse) v2_wins++; else if (m_mcmc.rmse < m_v2.rmse) mcmc_wins++;
     if (m_v2.mae < m_mcmc.mae) v2_wins++; else if (m_mcmc.mae < m_v2.mae) mcmc_wins++;
     if (m_v2.max_error < m_mcmc.max_error) v2_wins++; else if (m_mcmc.max_error < m_v2.max_error) mcmc_wins++;
     if (m_v2.tv_distance < m_mcmc.tv_distance) v2_wins++; else if (m_mcmc.tv_distance < m_v2.tv_distance) mcmc_wins++;
     if (m_v2.ci_coverage > m_mcmc.ci_coverage) v2_wins++; else if (m_mcmc.ci_coverage > m_v2.ci_coverage) mcmc_wins++;
     
     cerr << endl;
     if (v2_wins > mcmc_wins)
         cerr << "V2 WINS on " << v2_wins << "/" << (v2_wins + mcmc_wins) << " accuracy metrics" << endl;
     else if (mcmc_wins > v2_wins)
         cerr << "MCMC WINS on " << mcmc_wins << "/" << (v2_wins + mcmc_wins) << " accuracy metrics" << endl;
     else
         cerr << "TIE" << endl;
     
     // Detailed per-type table
     if (detailed) {
         cerr << endl;
         cerr << "Per-type comparison:" << endl;
         cerr << "┌───────────────────────────────────────────────────────────────────────┐" << endl;
         cerr << "│ Type │   EXACT    │ MCMC_ORIG  │ MCMC Err │     V2     │ V2 Err   │" << endl;
         cerr << "├───────────────────────────────────────────────────────────────────────┤" << endl;
         
         for (int i = 0; i < NUM_TYPES; i++) {
             if (exact.freq[i] > 0.0001 || mcmc.freq[i] > 0.0001 || v2.freq[i] > 0.0001) {
                 double mcmc_err = abs(mcmc.freq[i] - exact.freq[i]);
                 double v2_err = abs(v2.freq[i] - exact.freq[i]);
                 
                 string better = "~";
                 if (v2_err < mcmc_err - 0.0001) better = "V2";
                 else if (mcmc_err < v2_err - 0.0001) better = "MCMC";
                 
                 cerr << "│ " << setw(4) << (i+1) << " │ "
                      << fixed << setprecision(6) << setw(10) << exact.freq[i] << " │ "
                      << setw(10) << mcmc.freq[i] << " │ "
                      << setw(8) << setprecision(4) << mcmc_err << " │ "
                      << setw(10) << setprecision(6) << v2.freq[i] << " │ "
                      << setw(8) << setprecision(4) << v2_err << " │ "
                      << setw(4) << better << endl;
             }
         }
         
         cerr << "└───────────────────────────────────────────────────────────────────────┘" << endl;
         
         // CI coverage detail
         cerr << endl;
         cerr << "V2 Confidence Interval Coverage:" << endl;
         int covered = 0, total = 0;
         for (int i = 0; i < NUM_TYPES; i++) {
             if (exact.freq[i] > 0.0001 || v2.freq[i] > 0.0001) {
                 bool in_ci = (exact.freq[i] >= v2.ci_lo[i] && exact.freq[i] <= v2.ci_hi[i]);
                 if (in_ci) covered++;
                 total++;
                 
                 cerr << "  Type " << setw(2) << (i+1) << ": exact=" 
                      << fixed << setprecision(4) << exact.freq[i]
                      << " CI=[" << v2.ci_lo[i] << ", " << v2.ci_hi[i] << "] "
                      << (in_ci ? "✓" : "✗") << endl;
             }
         }
         cerr << "Coverage: " << covered << "/" << total << " = "
              << fixed << setprecision(1) << (100.0 * covered / max(total, 1)) << "%" << endl;
     }
 }
 
 int main(int argc, char* argv[]) {
     if (argc < 4) {
         cerr << "Usage: " << argv[0] << " <exact.res> <mcmc.res> <v2.res> [--detailed]" << endl;
         cerr << "If exact.res == mcmc.res, shows side-by-side (no ground truth)." << endl;
         return 1;
     }
     
     bool detailed = (argc >= 5 && string(argv[4]) == "--detailed");
     bool has_exact = (string(argv[1]) != string(argv[2]));
     
     Result exact, mcmc, v2;
     
     if (!exact.load(argv[1])) { cerr << "Cannot load " << argv[1] << endl; return 1; }
     if (!mcmc.load(argv[2]))  { cerr << "Cannot load " << argv[2] << endl; return 1; }
     if (!v2.load(argv[3]))    { cerr << "Cannot load " << argv[3] << endl; return 1; }
     
     cerr << endl;
     cerr << "Loaded:" << endl;
     cerr << "  [1] " << exact.method << " (time: " << fixed << setprecision(1) << exact.time_ms << " ms)" << endl;
     cerr << "  [2] " << mcmc.method << " (time: " << mcmc.time_ms << " ms)" << endl;
     cerr << "  [3] " << v2.method << " (time: " << v2.time_ms << " ms)" << endl;
     
     if (has_exact) {
         Metrics m_mcmc = computeMetrics(mcmc, exact);
         Metrics m_v2 = computeMetrics(v2, exact);
         printComparisonWithExact(m_mcmc, m_v2, exact, mcmc, v2, detailed);
     } else {
         printComparisonNoExact(mcmc, v2);
     }
     
     return 0;
 }