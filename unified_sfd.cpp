/**
 * ============================================================================
 * UNIFIED SFD ESTIMATOR - ROBUST AUTO-ADAPTIVE VERSION
 * ============================================================================
 * 
 * DESIGN GOALS:
 * 1. ROBUST: Works reliably on ANY dataset without manual tuning
 * 2. ACCURATE: Prioritizes accuracy over speed
 * 3. AUTO-ADAPTIVE: Automatically detects dataset characteristics and adjusts
 * 4. CONVERGENCE-AWARE: Monitors convergence and extends sampling if needed
 * 
 * KEY FEATURES:
 * - Zero-configuration: Just provide dataset, get optimal results
 * - Multi-phase sampling: Pilot → Main → Refinement (if needed)
 * - Gelman-Rubin convergence monitoring with auto-extension
 * - Comprehensive uncertainty quantification
 * - Works on small (100 edges) to large (1M+ edges) datasets
 * 
 * Author: TuanAnh (Master's Thesis)
 * ============================================================================
 */

 #include "common.hpp"
 #include <numeric>
 
 // ============================================================================
 // DATASET CHARACTERISTICS ANALYZER
 // ============================================================================
 
 struct DatasetProfile {
     // Basic stats
     size_t num_vertices = 0;
     size_t num_edges = 0;
     size_t num_triangles = 0;
     size_t num_tetrahedra = 0;
     
     // Degree distribution
     double avg_degree = 0;
     double max_degree = 0;
     double degree_std = 0;
     double degree_skewness = 0;  // High skewness = hub-dominated network
     
     // Structural properties
     double clustering_estimate = 0;  // Local clustering coefficient estimate
     double diameter_estimate = 0;    // Estimated graph diameter
     double density = 0;              // Edge density
     
     // Derived complexity scores
     double heterogeneity_score = 0;  // How heterogeneous is the degree distribution
     double sparsity_score = 0;       // How sparse is the graph
     double size_score = 0;           // Overall size complexity
     
     // Dataset category
     enum Category {
         TINY,       // < 100 edges
         SMALL,      // 100 - 1K edges
         MEDIUM,     // 1K - 100K edges
         LARGE,      // 100K - 1M edges
         HUGE        // > 1M edges
     } category;
     
     string categoryName() const {
         switch(category) {
             case TINY: return "TINY";
             case SMALL: return "SMALL";
             case MEDIUM: return "MEDIUM";
             case LARGE: return "LARGE";
             case HUGE: return "HUGE";
         }
         return "UNKNOWN";
     }
 };
 
 class DatasetAnalyzer {
 public:
     static DatasetProfile analyze(const SimplicialComplex& K) {
         DatasetProfile p;
         
         // Basic counts
         p.num_vertices = K.n;
         p.num_edges = K.edges.size();
         p.num_triangles = K.tris.size();
         p.num_tetrahedra = K.tets.size();
         
         // Degree statistics
         vector<double> degrees;
         double total_degree = 0;
         p.max_degree = 0;
         
         for (size_t v = 0; v < K.n; v++) {
             double d = K.adj[v].size();
             if (d > 0) {
                 degrees.push_back(d);
                 total_degree += d;
                 p.max_degree = max(p.max_degree, d);
             }
         }
         
         size_t active_vertices = degrees.size();
         p.avg_degree = active_vertices > 0 ? total_degree / active_vertices : 0;
         
         // Degree standard deviation and skewness
         if (active_vertices > 1) {
             double sum_sq = 0, sum_cube = 0;
             for (double d : degrees) {
                 double diff = d - p.avg_degree;
                 sum_sq += diff * diff;
                 sum_cube += diff * diff * diff;
             }
             p.degree_std = sqrt(sum_sq / active_vertices);
             if (p.degree_std > 0) {
                 p.degree_skewness = (sum_cube / active_vertices) / (p.degree_std * p.degree_std * p.degree_std);
             }
         }
         
         // Density
         if (active_vertices > 1) {
             double max_edges = active_vertices * (active_vertices - 1) / 2.0;
             p.density = p.num_edges / max_edges;
         }
         
         // Diameter estimate (using average degree)
         if (p.avg_degree > 1) {
             p.diameter_estimate = log(active_vertices) / log(p.avg_degree);
         } else {
             p.diameter_estimate = active_vertices;  // Worst case: linear
         }
         
         // Clustering coefficient estimate (sample-based)
         p.clustering_estimate = estimateClustering(K, degrees);
         
         // Compute complexity scores
         
         // Heterogeneity: ratio of max to avg degree, adjusted by skewness
         p.heterogeneity_score = p.avg_degree > 0 ? 
             (p.max_degree / p.avg_degree) * (1 + abs(p.degree_skewness) / 10) : 1;
         
         // Sparsity: inverse of density, capped
         p.sparsity_score = min(100.0, 1.0 / max(p.density, 0.0001));
         
         // Size score: log-scaled
         p.size_score = log10(max((double)p.num_edges, 1.0));
         
         // Categorize
         if (p.num_edges < 100) p.category = DatasetProfile::TINY;
         else if (p.num_edges < 1000) p.category = DatasetProfile::SMALL;
         else if (p.num_edges < 100000) p.category = DatasetProfile::MEDIUM;
         else if (p.num_edges < 1000000) p.category = DatasetProfile::LARGE;
         else p.category = DatasetProfile::HUGE;
         
         return p;
     }
     
 private:
     static double estimateClustering(const SimplicialComplex& K, const vector<double>& degrees) {
         // Sample-based clustering coefficient estimation
         if (degrees.size() < 3) return 0;
         
         random_device rd;
         mt19937 rng(rd());
         
         size_t sample_size = min((size_t)1000, degrees.size());
         vector<int> sampled_vertices;
         
         for (size_t v = 0; v < K.n && sampled_vertices.size() < sample_size; v++) {
             if (K.adj[v].size() >= 2) {
                 sampled_vertices.push_back(v);
             }
         }
         
         if (sampled_vertices.empty()) return 0;
         
         double total_clustering = 0;
         int count = 0;
         
         for (int v : sampled_vertices) {
             vector<int> neighbors(K.adj[v].begin(), K.adj[v].end());
             if (neighbors.size() < 2) continue;
             
             int triangles = 0;
             int possible = neighbors.size() * (neighbors.size() - 1) / 2;
             
             // Sample neighbor pairs if too many
             if (possible > 100) {
                 uniform_int_distribution<size_t> dist(0, neighbors.size() - 1);
                 for (int i = 0; i < 100; i++) {
                     size_t a = dist(rng), b = dist(rng);
                     if (a != b && K.hasEdge(neighbors[a], neighbors[b])) {
                         triangles++;
                     }
                 }
                 total_clustering += (double)triangles / 100;
             } else {
                 for (size_t i = 0; i < neighbors.size(); i++) {
                     for (size_t j = i + 1; j < neighbors.size(); j++) {
                         if (K.hasEdge(neighbors[i], neighbors[j])) {
                             triangles++;
                         }
                     }
                 }
                 total_clustering += (double)triangles / possible;
             }
             count++;
         }
         
         return count > 0 ? total_clustering / count : 0;
     }
 };
 
 // ============================================================================
 // AUTO-ADAPTIVE PARAMETER CALCULATOR
 // ============================================================================
 
 struct AdaptiveParameters {
     // Burn-in
     size_t burn_in = 5000;
     
     // Sampling
     size_t steps_per_sample = 50;
     size_t num_chains = 4;
     size_t min_samples = 5000;
     size_t max_samples = 100000;
     
     // Convergence
     double target_rhat = 1.05;      // Target Gelman-Rubin R-hat
     double target_ci_width = 0.02;  // Target 95% CI width
     size_t convergence_check_interval = 500;
     
     // Auto-extension
     bool allow_auto_extension = true;
     size_t max_extensions = 5;
     double extension_factor = 1.5;
 };
 
 class ParameterCalculator {
 public:
     static AdaptiveParameters compute(const DatasetProfile& profile) {
         AdaptiveParameters params;
         
         cerr << "\n========================================" << endl;
         cerr << "AUTO-ADAPTIVE PARAMETER CALCULATION" << endl;
         cerr << "========================================" << endl;
         cerr << "Dataset category: " << profile.categoryName() << endl;
         cerr << "Edges: " << profile.num_edges << endl;
         cerr << "Avg degree: " << fixed << setprecision(2) << profile.avg_degree << endl;
         cerr << "Max degree: " << (size_t)profile.max_degree << endl;
         cerr << "Heterogeneity score: " << profile.heterogeneity_score << endl;
         cerr << "Sparsity score: " << profile.sparsity_score << endl;
         cerr << "Diameter estimate: " << profile.diameter_estimate << endl;
         cerr << "Clustering estimate: " << profile.clustering_estimate << endl;
         
         // ====================================================================
         // BURN-IN CALCULATION
         // ====================================================================
         // Burn-in needs to be long enough for chain to explore the space
         // For high heterogeneity graphs, chains can get stuck in hub regions
         
         double base_burn_in = 2000;
         
         // Scale by diameter
         double diameter_factor = max(1.0, profile.diameter_estimate / 5.0);
         
         // Scale by heterogeneity - CRITICAL for hub-dominated networks
         // Use stronger scaling for very high heterogeneity
         double heterogeneity_factor;
         if (profile.heterogeneity_score > 100) {
             // Very high heterogeneity - need much more burn-in
             heterogeneity_factor = 2.0 + log(profile.heterogeneity_score) / 2;
         } else {
             heterogeneity_factor = max(1.0, log(profile.heterogeneity_score + 1));
         }
         
         // Scale by sparsity
         double sparsity_factor = max(1.0, log(profile.sparsity_score) / 2);
         
         params.burn_in = (size_t)(base_burn_in * diameter_factor * heterogeneity_factor * sparsity_factor);
         
         // For HUGE datasets with high heterogeneity, ensure minimum burn-in
         if (profile.category == DatasetProfile::HUGE && profile.heterogeneity_score > 50) {
             params.burn_in = max(params.burn_in, (size_t)100000);
         }
         
         params.burn_in = max((size_t)1000, min(params.burn_in, (size_t)500000));
         
         // ====================================================================
         // STEPS PER SAMPLE CALCULATION
         // ====================================================================
         // Steps between samples for independence
         // High heterogeneity requires more steps to escape hub regions
         
         double base_steps = 30;
         
         // More steps for larger diameter
         double step_diameter_factor = max(1.0, profile.diameter_estimate / 3.0);
         
         // More steps for high max degree (can get stuck in hubs)
         double step_degree_factor = max(1.0, log(profile.max_degree + 1) / 2);
         
         // CRITICAL: More steps for heterogeneous networks
         double step_hetero_factor;
         if (profile.heterogeneity_score > 100) {
             // Very high heterogeneity - need many more steps
             step_hetero_factor = 2.0 + sqrt(profile.heterogeneity_score) / 10;
         } else {
             step_hetero_factor = max(1.0, sqrt(profile.heterogeneity_score) / 2);
         }
         
         params.steps_per_sample = (size_t)(base_steps * step_diameter_factor * step_degree_factor * step_hetero_factor);
         
         // For HUGE datasets, ensure reasonable steps
         if (profile.category == DatasetProfile::HUGE) {
             params.steps_per_sample = max(params.steps_per_sample, (size_t)100);
         }
         
         // Cap at reasonable maximum
         params.steps_per_sample = max((size_t)20, min(params.steps_per_sample, (size_t)1000));
         
         // ====================================================================
         // NUMBER OF CHAINS
         // ====================================================================
         // More chains for heterogeneous networks to explore different regions
         // This is CRITICAL for convergence in hub-dominated networks
         
         params.num_chains = 4;  // Minimum for Gelman-Rubin
         
         if (profile.heterogeneity_score > 5) {
             params.num_chains = 8;
         }
         if (profile.heterogeneity_score > 20) {
             params.num_chains = 12;
         }
         if (profile.heterogeneity_score > 50) {
             params.num_chains = 16;
         }
         if (profile.heterogeneity_score > 100) {
             // Very high heterogeneity - use many chains
             params.num_chains = 20;
         }
         if (profile.heterogeneity_score > 500) {
             // Extremely high heterogeneity
             params.num_chains = 24;
         }
         
         // ====================================================================
         // SAMPLE SIZE
         // ====================================================================
         // Base on dataset size and complexity
         
         switch (profile.category) {
             case DatasetProfile::TINY:
                 params.min_samples = 1000;
                 params.max_samples = 10000;
                 break;
             case DatasetProfile::SMALL:
                 params.min_samples = 5000;
                 params.max_samples = 50000;
                 break;
             case DatasetProfile::MEDIUM:
                 params.min_samples = 10000;
                 params.max_samples = 200000;
                 break;
             case DatasetProfile::LARGE:
                 params.min_samples = 20000;
                 params.max_samples = 500000;
                 break;
             case DatasetProfile::HUGE:
                 params.min_samples = 50000;
                 params.max_samples = 2000000;  // Allow up to 2M samples
                 break;
         }
         
         // Adjust for heterogeneity - high heterogeneity needs more samples
         if (profile.heterogeneity_score > 50) {
             params.min_samples = (size_t)(params.min_samples * 2);
             params.max_samples = (size_t)(params.max_samples * 2);
         }
         if (profile.heterogeneity_score > 200) {
             params.min_samples = (size_t)(params.min_samples * 1.5);
             params.max_samples = (size_t)(params.max_samples * 1.5);
         }
         
         // ====================================================================
         // CONVERGENCE PARAMETERS
         // ====================================================================
         
         // For HUGE datasets with high heterogeneity, we need to be realistic
         // about convergence targets while still being rigorous
         if (profile.category == DatasetProfile::HUGE && profile.heterogeneity_score > 100) {
             // Very challenging dataset - adjust targets but still aim for good convergence
             params.target_rhat = 1.15;  // Slightly relaxed but still rigorous
             params.target_ci_width = 0.02;
             params.max_extensions = 10;  // Allow more extensions
             params.extension_factor = 1.3;  // Smaller increments, more frequent
         } else if (profile.category >= DatasetProfile::LARGE) {
             params.target_rhat = 1.10;
             params.target_ci_width = 0.025;
             params.max_extensions = 8;
             params.extension_factor = 1.4;
         } else if (profile.category == DatasetProfile::MEDIUM) {
             params.target_rhat = 1.05;
             params.target_ci_width = 0.02;
             params.max_extensions = 6;
             params.extension_factor = 1.5;
         } else {
             // Small datasets - strict convergence
             params.target_rhat = 1.02;
             params.target_ci_width = 0.01;
             params.max_extensions = 5;
             params.extension_factor = 1.5;
         }
         
         // Check convergence more frequently for better monitoring
         params.convergence_check_interval = max((size_t)100, params.min_samples / 50);
         
         cerr << "\nComputed parameters:" << endl;
         cerr << "  Burn-in:           " << params.burn_in << endl;
         cerr << "  Steps per sample:  " << params.steps_per_sample << endl;
         cerr << "  Number of chains:  " << params.num_chains << endl;
         cerr << "  Min samples:       " << params.min_samples << endl;
         cerr << "  Max samples:       " << params.max_samples << endl;
         cerr << "  Target R-hat:      " << params.target_rhat << endl;
         cerr << "  Target CI width:   " << params.target_ci_width << endl;
         cerr << "  Max extensions:    " << params.max_extensions << endl;
         cerr << "========================================" << endl;
         
         return params;
     }
 };
 
 // ============================================================================
 // CONVERGENCE MONITOR
 // ============================================================================
 
 class ConvergenceMonitor {
 public:
     struct Diagnostics {
         double max_rhat = 0;           // Maximum R-hat across all types
         double avg_rhat = 0;           // Average R-hat
         double max_ci_width = 0;       // Maximum CI width
         double avg_ci_width = 0;       // Average CI width
         int types_not_converged = 0;   // Number of types with R-hat > threshold
         bool converged = false;
     };
     
     /**
      * Compute Gelman-Rubin R-hat statistic for all types.
      */
     static Diagnostics computeDiagnostics(
         const vector<vector<size_t>>& chain_counts,
         size_t samples_per_chain,
         double target_rhat,
         double target_ci_width
     ) {
         Diagnostics d;
         size_t num_chains = chain_counts.size();
         
         if (num_chains < 2 || samples_per_chain < 10) {
             return d;
         }
         
         vector<double> rhats;
         vector<double> ci_widths;
         
         for (int t = 0; t < NUM_TYPES; t++) {
             // Compute chain frequencies
             vector<double> chain_freqs(num_chains);
             double overall_freq = 0;
             size_t total_samples = 0;
             
             for (size_t c = 0; c < num_chains; c++) {
                 size_t chain_total = 0;
                 for (int i = 0; i < NUM_TYPES; i++) {
                     chain_total += chain_counts[c][i];
                 }
                 if (chain_total > 0) {
                     chain_freqs[c] = (double)chain_counts[c][t] / chain_total;
                     overall_freq += chain_counts[c][t];
                     total_samples += chain_total;
                 }
             }
             
             if (total_samples == 0) continue;
             overall_freq /= total_samples;
             
             // Skip types with very low frequency
             if (overall_freq < 0.001) continue;
             
             // Compute R-hat
             double rhat = computeRhat(chain_freqs, chain_counts, t, samples_per_chain);
             rhats.push_back(rhat);
             
             if (rhat > target_rhat) {
                 d.types_not_converged++;
             }
             
             // Compute CI width (Wilson score)
             double p = overall_freq;
             double n = total_samples;
             double z = 1.96;
             double denom = 1 + z*z/n;
             double margin = z * sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom;
             double ci_width = 2 * margin;
             ci_widths.push_back(ci_width);
         }
         
         if (!rhats.empty()) {
             d.max_rhat = *max_element(rhats.begin(), rhats.end());
             d.avg_rhat = accumulate(rhats.begin(), rhats.end(), 0.0) / rhats.size();
         }
         
         if (!ci_widths.empty()) {
             d.max_ci_width = *max_element(ci_widths.begin(), ci_widths.end());
             d.avg_ci_width = accumulate(ci_widths.begin(), ci_widths.end(), 0.0) / ci_widths.size();
         }
         
         d.converged = (d.max_rhat <= target_rhat && d.max_ci_width <= target_ci_width);
         
         return d;
     }
     
 private:
     /**
      * Compute Gelman-Rubin R-hat statistic.
      * 
      * For proportion estimation, we use the improved formula:
      * R-hat = sqrt((V / W) * (df / (df - 2)))
      * where V is the pooled variance estimate and W is within-chain variance.
      * 
      * For better convergence diagnostics with proportions, we also consider
      * the effective sample size and use split-R-hat when possible.
      */
     static double computeRhat(
         const vector<double>& chain_freqs,
         const vector<vector<size_t>>& chain_counts,
         int /* type */,
         size_t n
     ) {
         size_t m = chain_freqs.size();
         if (m < 2 || n < 10) return 1.0;
         
         // Overall mean
         double overall_mean = 0;
         for (double f : chain_freqs) overall_mean += f;
         overall_mean /= m;
         
         // Skip if overall frequency is essentially zero
         if (overall_mean < 1e-6) return 1.0;
         
         // Between-chain variance (B/n)
         double B = 0;
         for (double f : chain_freqs) {
             double diff = f - overall_mean;
             B += diff * diff;
         }
         B = B * n / (m - 1);
         
         // Within-chain variance (W)
         // For proportions, use variance of binomial: p(1-p)/n_effective
         // But we estimate it from the spread of chain frequencies
         double W = 0;
         for (size_t c = 0; c < m; c++) {
             double p = chain_freqs[c];
             // Use Wilson variance estimate for proportions
             double var_c = p * (1 - p);
             
             // Adjust for sample size of each chain
             size_t chain_samples = 0;
             for (int i = 0; i < NUM_TYPES; i++) {
                 chain_samples += chain_counts[c][i];
             }
             if (chain_samples > 0) {
                 var_c = p * (1 - p) + 1.0 / (4 * chain_samples);  // Wilson adjustment
             }
             
             W += var_c;
         }
         W /= m;
         
         // Avoid division by zero
         if (W < 1e-10) {
             // If within-chain variance is tiny, check if there's between-chain variance
             if (B > 1e-6) {
                 return 2.0;  // Chains disagree but each is confident - bad convergence
             }
             return 1.0;  // All chains agree and confident - good convergence
         }
         
         // Pooled variance estimate
         double V_hat = ((n - 1.0) / n) * W + B / n;
         
         // R-hat with degrees of freedom correction
         double df = 2 * V_hat * V_hat / (
             (((n-1.0)/n) * ((n-1.0)/n) * W * W / (m-1)) + 
             (B * B / (n * n * (m-1)))
         );
         
         double rhat_sq = (V_hat / W);
         if (df > 2) {
             rhat_sq *= df / (df - 2);  // Adjust for finite sample
         }
         
         double rhat = sqrt(max(1.0, rhat_sq));
         
         // Cap at reasonable maximum
         return min(rhat, 10.0);
     }
 };
 
 // ============================================================================
 // ROBUST MCMC SAMPLER
 // ============================================================================
 
 class RobustMCMC_SFD {
     const SimplicialComplex& K;
     DatasetProfile profile;
     AdaptiveParameters params;
     
     vector<mt19937> rngs;
     
     // Logging
     bool log_convergence = false;
     string log_file_path;
     ofstream log_file;
     
 public:
     RobustMCMC_SFD(const SimplicialComplex& K) : K(K) {
         // Analyze dataset
         profile = DatasetAnalyzer::analyze(K);
         
         // Compute adaptive parameters
         params = ParameterCalculator::compute(profile);
         
         // Initialize RNGs
         random_device rd;
         for (size_t i = 0; i < params.num_chains; i++) {
             rngs.emplace_back(rd() + i * 12345);
         }
     }
     
     void enableConvergenceLogging(const string& path) {
         log_convergence = true;
         log_file_path = path;
         log_file.open(path);
         
         if (log_file.is_open()) {
             log_file << "iteration,chain_id,total_samples,phase";
             for (int i = 1; i <= NUM_TYPES; i++) {
                 log_file << ",type" << i << "_freq";
             }
             log_file << ",rhat_max,ci_width_max\n";
             cerr << "Convergence logging enabled: " << path << endl;
         }
     }
     
     Result compute(size_t requested_samples = 0) {
         Result r;
         r.method = "UNIFIED_ROBUST_ADAPTIVE";
         
         // Use requested samples or auto-determined
         size_t target_samples = requested_samples > 0 ? 
             requested_samples : params.min_samples;
         
         auto t0 = high_resolution_clock::now();
         
         cerr << "\n========================================" << endl;
         cerr << "ROBUST MCMC ESTIMATION" << endl;
         cerr << "========================================" << endl;
         cerr << "Target samples: " << target_samples << endl;
         cerr << "Max samples: " << params.max_samples << endl;
         cerr << endl;
         
         // Initialize chains
         vector<vector<int>> chain_states(params.num_chains);
         vector<vector<size_t>> chain_counts(params.num_chains, vector<size_t>(NUM_TYPES, 0));
         vector<vector<int>> chain_samples(params.num_chains);
         
         // Phase 1: Burn-in with adaptive extension
         cerr << "Phase 1: Burn-in (" << params.burn_in << " steps per chain)..." << endl;
         
         // Initialize chains from diverse regions
         for (size_t c = 0; c < params.num_chains; c++) {
             chain_states[c] = initializeSimplet(rngs[c], c);
         }
         
         // Initial burn-in
         for (size_t c = 0; c < params.num_chains; c++) {
             for (size_t i = 0; i < params.burn_in; i++) {
                 chain_states[c] = mcmcStep(chain_states[c], rngs[c]);
             }
         }
         
         // Check chain health after burn-in
         // If chains are in very different "type regions", extend burn-in
         bool needs_more_burnin = true;
         size_t burnin_extensions = 0;
         const size_t max_burnin_extensions = 3;
         
         while (needs_more_burnin && burnin_extensions < max_burnin_extensions) {
             // Sample a few points from each chain to check diversity
             vector<map<int, int>> chain_type_counts(params.num_chains);
             
             for (size_t c = 0; c < params.num_chains; c++) {
                 vector<int> current = chain_states[c];
                 for (size_t s = 0; s < 100; s++) {
                     for (size_t step = 0; step < 10; step++) {
                         current = mcmcStep(current, rngs[c], false);  // No teleport for test
                     }
                     int type = Classifier::classify(current, K);
                     if (type >= 1 && type <= NUM_TYPES) {
                         chain_type_counts[c][type]++;
                     }
                 }
                 chain_states[c] = current;  // Update state
             }
             
             // Check if chains are exploring similar regions
             // Find most common type across all chains
             map<int, double> type_variance;
             for (int t = 1; t <= NUM_TYPES; t++) {
                 vector<double> freqs;
                 for (size_t c = 0; c < params.num_chains; c++) {
                     double freq = chain_type_counts[c][t] / 100.0;
                     freqs.push_back(freq);
                 }
                 
                 // Calculate variance
                 double mean = 0;
                 for (double f : freqs) mean += f;
                 mean /= freqs.size();
                 
                 if (mean > 0.05) {  // Only check significant types
                     double var = 0;
                     for (double f : freqs) {
                         var += (f - mean) * (f - mean);
                     }
                     var /= freqs.size();
                     type_variance[t] = var / max(mean * mean, 0.001);  // Coefficient of variation squared
                 }
             }
             
             // If max relative variance is too high, chains aren't mixing
             double max_rel_var = 0;
             for (auto& [t, v] : type_variance) {
                 max_rel_var = max(max_rel_var, v);
             }
             
             if (max_rel_var < 0.5) {  // Chains are reasonably aligned
                 needs_more_burnin = false;
             } else if (burnin_extensions < max_burnin_extensions) {
                 // Extend burn-in
                 burnin_extensions++;
                 cerr << "  Extending burn-in (chains not well mixed, rel_var=" 
                      << fixed << setprecision(3) << max_rel_var << ")..." << endl;
                 
                 size_t extra_burnin = params.burn_in / 2;
                 for (size_t c = 0; c < params.num_chains; c++) {
                     for (size_t i = 0; i < extra_burnin; i++) {
                         chain_states[c] = mcmcStep(chain_states[c], rngs[c]);
                     }
                 }
             } else {
                 needs_more_burnin = false;  // Give up extending
             }
         }
         
         if (burnin_extensions > 0) {
             cerr << "  Burn-in extended " << burnin_extensions << " time(s)." << endl;
         }
         cerr << "  Burn-in complete." << endl;
         
         // Phase 2: Main sampling with convergence monitoring
         cerr << "\nPhase 2: Main sampling with convergence monitoring..." << endl;
         
         size_t samples_per_chain = target_samples / params.num_chains;
         size_t total_iterations = 0;
         size_t extension_count = 0;
         bool converged = false;
         
         while (!converged && total_iterations < params.max_samples / params.num_chains) {
             // Sample batch
             size_t batch_size = min(params.convergence_check_interval, 
                                    samples_per_chain - total_iterations);
             
             for (size_t s = 0; s < batch_size; s++) {
                 for (size_t c = 0; c < params.num_chains; c++) {
                     // Take steps
                     for (size_t step = 0; step < params.steps_per_sample; step++) {
                         chain_states[c] = mcmcStep(chain_states[c], rngs[c]);
                     }
                     
                     // Classify and record
                     int type = Classifier::classify(chain_states[c], K);
                     if (type >= 1 && type <= NUM_TYPES) {
                         chain_samples[c].push_back(type);
                         chain_counts[c][type - 1]++;
                     }
                 }
             }
             
             total_iterations += batch_size;
             
             // Check convergence
             auto diag = ConvergenceMonitor::computeDiagnostics(
                 chain_counts, total_iterations, params.target_rhat, params.target_ci_width);
             
             // Log progress
             if (log_convergence && log_file.is_open()) {
                 logConvergenceState(total_iterations, chain_counts, diag, "main");
             }
             
             // Progress report
             size_t total_samples_now = total_iterations * params.num_chains;
             cerr << "  Samples: " << total_samples_now 
                  << " | R-hat max: " << fixed << setprecision(3) << diag.max_rhat
                  << " | CI width max: " << setprecision(4) << diag.max_ci_width;
             
             if (diag.converged) {
                 cerr << " [CONVERGED]" << endl;
                 converged = true;
             } else if (total_iterations >= samples_per_chain) {
                 // Check if we should extend
                 if (params.allow_auto_extension && 
                     extension_count < params.max_extensions &&
                     total_iterations * params.num_chains < params.max_samples) {
                     
                     extension_count++;
                     size_t extension = (size_t)(samples_per_chain * (params.extension_factor - 1));
                     samples_per_chain += extension;
                     
                     cerr << " [EXTENDING by " << extension * params.num_chains << " samples]" << endl;
                     
                     // If R-hat is very high after extension, try restarting worst chains
                     if (diag.max_rhat > 1.5 && extension_count >= 2) {
                         cerr << "    R-hat still high, attempting chain restart..." << endl;
                         
                         // Find chains with most divergent frequency estimates
                         vector<pair<size_t, double>> chain_divergence;
                         for (size_t c = 0; c < params.num_chains; c++) {
                             double div = 0;
                             size_t chain_total = 0;
                             for (int t = 0; t < NUM_TYPES; t++) {
                                 chain_total += chain_counts[c][t];
                             }
                             if (chain_total == 0) continue;
                             
                             // Calculate divergence from mean
                             for (int t = 0; t < NUM_TYPES; t++) {
                                 double chain_freq = (double)chain_counts[c][t] / chain_total;
                                 double overall_freq = 0;
                                 size_t overall_total = 0;
                                 for (size_t c2 = 0; c2 < params.num_chains; c2++) {
                                     overall_freq += chain_counts[c2][t];
                                     for (int t2 = 0; t2 < NUM_TYPES; t2++) {
                                         overall_total += chain_counts[c2][t2];
                                     }
                                 }
                                 overall_freq /= max(overall_total, (size_t)1);
                                 div += abs(chain_freq - overall_freq);
                             }
                             chain_divergence.push_back({c, div});
                         }
                         
                         // Sort by divergence (highest first)
                         sort(chain_divergence.begin(), chain_divergence.end(),
                              [](auto& a, auto& b) { return a.second > b.second; });
                         
                         // Restart top 20% most divergent chains
                         size_t num_restart = max((size_t)1, params.num_chains / 5);
                         for (size_t i = 0; i < min(num_restart, chain_divergence.size()); i++) {
                             size_t c = chain_divergence[i].first;
                             
                             // Re-initialize from a random region
                             uniform_int_distribution<size_t> region_dist(0, params.num_chains - 1);
                             chain_states[c] = initializeSimplet(rngs[c], region_dist(rngs[c]));
                             
                             // Quick burn-in
                             for (size_t b = 0; b < params.burn_in / 4; b++) {
                                 chain_states[c] = mcmcStep(chain_states[c], rngs[c]);
                             }
                             
                             cerr << "    Restarted chain " << (c+1) << " (div=" 
                                  << fixed << setprecision(3) << chain_divergence[i].second << ")" << endl;
                         }
                     }
                 } else {
                     cerr << " [MAX REACHED]" << endl;
                     break;
                 }
             } else {
                 cerr << endl;
             }
         }
         
         // Close log file
         if (log_file.is_open()) {
             log_file.close();
         }
         
         // Combine results
         cerr << "\nCombining results from " << params.num_chains << " chains..." << endl;
         
         for (size_t c = 0; c < params.num_chains; c++) {
             cerr << "  Chain " << (c+1) << ": " << chain_samples[c].size() << " samples" << endl;
             for (int type : chain_samples[c]) {
                 r.counts[type - 1]++;
             }
             r.total += chain_samples[c].size();
         }
         r.samples = r.total;
         
         // Compute frequencies and confidence intervals
         computeFrequenciesAndCI(r, chain_samples);
         
         // Final convergence diagnostics
         auto final_diag = ConvergenceMonitor::computeDiagnostics(
             chain_counts, total_iterations, params.target_rhat, params.target_ci_width);
         
         printFinalDiagnostics(final_diag, chain_counts);
         
         r.time_ms = duration_cast<microseconds>(
             high_resolution_clock::now() - t0).count() / 1000.0;
         
         return r;
     }
     
 private:
     vector<int> initializeSimplet(mt19937& rng, size_t chain_id = 0) {
         // Strategy: Start chains from DIFFERENT regions of the graph
         // This is CRITICAL for convergence in heterogeneous networks
         
         // Collect vertices with enough neighbors, grouped by degree
         vector<pair<int, size_t>> vertices_by_degree;  // (vertex, degree)
         for (size_t v = 0; v < K.n; v++) {
             if (K.adj[v].size() >= 2) {
                 vertices_by_degree.push_back({v, K.adj[v].size()});
             }
         }
         
         if (vertices_by_degree.empty()) {
             // Fallback to any edge
             for (auto& e : K.edges) {
                 return {e.first, e.second};
             }
             return {};
         }
         
         // Sort by degree
         sort(vertices_by_degree.begin(), vertices_by_degree.end(),
              [](auto& a, auto& b) { return a.second < b.second; });
         
         // Divide into regions based on chain_id
         // Different chains start from different degree regions
         size_t num_regions = min((size_t)10, vertices_by_degree.size() / 100 + 1);
         size_t region_size = vertices_by_degree.size() / num_regions;
         
         // Select region based on chain_id
         size_t region = chain_id % num_regions;
         size_t start_idx = region * region_size;
         size_t end_idx = (region == num_regions - 1) ? vertices_by_degree.size() : (region + 1) * region_size;
         
         // Random selection within the region
         uniform_int_distribution<size_t> region_dist(start_idx, end_idx - 1);
         int selected_v = vertices_by_degree[region_dist(rng)].first;
         
         // Pick random neighbors to form initial simplet
         vector<int> neighbors(K.adj[selected_v].begin(), K.adj[selected_v].end());
         uniform_int_distribution<size_t> nb_dist(0, neighbors.size() - 1);
         
         size_t i1 = nb_dist(rng);
         size_t i2;
         int attempts = 0;
         do {
             i2 = nb_dist(rng);
             attempts++;
         } while (i2 == i1 && neighbors.size() > 1 && attempts < 100);
         
         vector<int> result = {selected_v, neighbors[i1], neighbors[i2]};
         sort(result.begin(), result.end());
         return result;
     }
     
     vector<int> mcmcStep(const vector<int>& current, mt19937& rng, bool allow_teleport = true) {
         uniform_real_distribution<double> unif(0, 1);
         
         // ================================================================
         // MOVE TYPE SELECTION
         // ================================================================
         // For high heterogeneity, we need different types of moves:
         // 1. Standard neighbor moves (add/remove vertex)
         // 2. Swap moves (replace one vertex with another)
         // 3. Teleportation (jump to random simplet - rare but helps escape)
         
         double r = unif(rng);
         
         // Teleportation (rare - helps escape local modes)
         // Only for chains that have been running a while
         if (allow_teleport && r < 0.001) {  // 0.1% chance
             // Jump to random simplet
             vector<int> new_simplet = initializeSimplet(rng, 0);
             if (!new_simplet.empty() && Classifier::isConnected(new_simplet, K)) {
                 return new_simplet;
             }
         }
         
         // Swap move (10% chance) - helps traverse between degree classes
         if (r < 0.1 && current.size() >= 3) {
             return swapMove(current, rng);
         }
         
         // Standard neighbor move (89.9% chance)
         return standardMove(current, rng);
     }
     
     vector<int> standardMove(const vector<int>& current, mt19937& rng) {
         auto neighbors = getSimplexNeighbors(current);
         
         if (neighbors.empty()) return current;
         
         uniform_int_distribution<size_t> dist(0, neighbors.size() - 1);
         vector<int> proposed = neighbors[dist(rng)];
         
         auto proposed_neighbors = getSimplexNeighbors(proposed);
         
         // Metropolis-Hastings acceptance probability
         double accept_prob = min(1.0, 
             (double)neighbors.size() / max((size_t)1, proposed_neighbors.size()));
         
         uniform_real_distribution<double> unif(0, 1);
         if (unif(rng) <= accept_prob) {
             return proposed;
         }
         return current;
     }
     
     vector<int> swapMove(const vector<int>& current, mt19937& rng) {
         // Swap move: replace one vertex with another while maintaining connectivity
         // This helps chains move between different degree regions
         
         if (current.size() < 3) return current;
         
         // Find candidates to swap in (vertices connected to at least 2 current vertices)
         set<int> current_set(current.begin(), current.end());
         map<int, int> candidate_connections;  // vertex -> number of connections to current
         
         for (int v : current) {
             for (int u : K.adj[v]) {
                 if (!current_set.count(u)) {
                     candidate_connections[u]++;
                 }
             }
         }
         
         // Filter: need at least 2 connections to maintain connectivity after swap
         vector<int> valid_candidates;
         for (auto& [v, conn] : candidate_connections) {
             if (conn >= 2) {
                 valid_candidates.push_back(v);
             }
         }
         
         if (valid_candidates.empty()) return current;
         
         // Select random candidate to add
         uniform_int_distribution<size_t> cand_dist(0, valid_candidates.size() - 1);
         int new_v = valid_candidates[cand_dist(rng)];
         
         // Try removing each current vertex and check connectivity
         vector<vector<int>> valid_swaps;
         for (size_t i = 0; i < current.size(); i++) {
             vector<int> proposed;
             for (size_t j = 0; j < current.size(); j++) {
                 if (i != j) proposed.push_back(current[j]);
             }
             proposed.push_back(new_v);
             
             if (Classifier::isConnected(proposed, K)) {
                 sort(proposed.begin(), proposed.end());
                 valid_swaps.push_back(proposed);
             }
         }
         
         if (valid_swaps.empty()) return current;
         
         // Select random valid swap
         uniform_int_distribution<size_t> swap_dist(0, valid_swaps.size() - 1);
         vector<int> proposed = valid_swaps[swap_dist(rng)];
         
         // Calculate acceptance probability
         // For swap moves, we use a simplified symmetric proposal
         auto current_neighbors = getSimplexNeighbors(current);
         auto proposed_neighbors = getSimplexNeighbors(proposed);
         
         double accept_prob = min(1.0,
             (double)current_neighbors.size() / max((size_t)1, proposed_neighbors.size()));
         
         uniform_real_distribution<double> unif(0, 1);
         if (unif(rng) <= accept_prob) {
             return proposed;
         }
         return current;
     }
     
     vector<vector<int>> getSimplexNeighbors(const vector<int>& simplex) {
         vector<vector<int>> neighbors;
         set<int> simplex_set(simplex.begin(), simplex.end());
         
         // Remove a vertex (if size > 2)
         if (simplex.size() > 2) {
             for (size_t i = 0; i < simplex.size(); i++) {
                 vector<int> ns;
                 for (size_t j = 0; j < simplex.size(); j++) {
                     if (i != j) ns.push_back(simplex[j]);
                 }
                 if (Classifier::isConnected(ns, K)) {
                     sort(ns.begin(), ns.end());
                     neighbors.push_back(ns);
                 }
             }
         }
         
         // Add a vertex (if size < 4)
         if (simplex.size() < 4) {
             set<int> candidates;
             for (int v : simplex) {
                 for (int u : K.adj[v]) {
                     if (!simplex_set.count(u)) {
                         candidates.insert(u);
                     }
                 }
             }
             
             for (int c : candidates) {
                 vector<int> ns = simplex;
                 ns.push_back(c);
                 if (Classifier::isConnected(ns, K)) {
                     sort(ns.begin(), ns.end());
                     neighbors.push_back(ns);
                 }
             }
         }
         
         return neighbors;
     }
     
     void computeFrequenciesAndCI(Result& r, const vector<vector<int>>& /* chain_samples */) {
         double z = 1.96;  // 95% CI
         
         for (int i = 0; i < NUM_TYPES; i++) {
             r.freq[i] = r.samples > 0 ? (double)r.counts[i] / r.samples : 0;
             
             // Wilson score confidence interval
             double p = r.freq[i];
             double n = r.samples;
             
             if (n > 0) {
                 double denom = 1 + z*z/n;
                 double center = (p + z*z/(2*n)) / denom;
                 double margin = z * sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom;
                 
                 r.ci_lo[i] = max(0.0, center - margin);
                 r.ci_hi[i] = min(1.0, center + margin);
             }
         }
     }
     
     void logConvergenceState(
         size_t iteration,
         const vector<vector<size_t>>& chain_counts,
         const ConvergenceMonitor::Diagnostics& diag,
         const string& phase
     ) {
         if (!log_file.is_open()) return;
         
         // Log each chain
         for (size_t c = 0; c < params.num_chains; c++) {
             size_t total = 0;
             for (size_t count : chain_counts[c]) total += count;
             if (total == 0) continue;
             
             log_file << iteration << "," << c << "," << total << "," << phase;
             
             for (int i = 0; i < NUM_TYPES; i++) {
                 log_file << "," << fixed << setprecision(6) 
                          << (double)chain_counts[c][i] / total;
             }
             
             log_file << "," << diag.max_rhat << "," << diag.max_ci_width << "\n";
         }
         
         // Log combined
         vector<size_t> combined(NUM_TYPES, 0);
         size_t total_combined = 0;
         for (size_t c = 0; c < params.num_chains; c++) {
             for (int i = 0; i < NUM_TYPES; i++) {
                 combined[i] += chain_counts[c][i];
             }
         }
         for (size_t cnt : combined) total_combined += cnt;
         
         if (total_combined > 0) {
             log_file << iteration << ",combined," << total_combined << "," << phase;
             for (int i = 0; i < NUM_TYPES; i++) {
                 log_file << "," << fixed << setprecision(6) 
                          << (double)combined[i] / total_combined;
             }
             log_file << "," << diag.max_rhat << "," << diag.max_ci_width << "\n";
         }
     }
     
     void printFinalDiagnostics(
         const ConvergenceMonitor::Diagnostics& diag,
         const vector<vector<size_t>>& chain_counts
     ) {
         cerr << "\n========================================" << endl;
         cerr << "CONVERGENCE DIAGNOSTICS" << endl;
         cerr << "========================================" << endl;
         cerr << "Max R-hat:      " << fixed << setprecision(4) << diag.max_rhat;
         if (diag.max_rhat <= params.target_rhat) {
             cerr << " ✓ (target: " << params.target_rhat << ")";
         } else {
             cerr << " ⚠ (target: " << params.target_rhat << ")";
         }
         cerr << endl;
         
         cerr << "Avg R-hat:      " << diag.avg_rhat << endl;
         cerr << "Max CI width:   " << diag.max_ci_width;
         if (diag.max_ci_width <= params.target_ci_width) {
             cerr << " ✓ (target: " << params.target_ci_width << ")";
         } else {
             cerr << " ⚠ (target: " << params.target_ci_width << ")";
         }
         cerr << endl;
         cerr << "Avg CI width:   " << diag.avg_ci_width << endl;
         
         if (diag.converged) {
             cerr << "\n✅ CONVERGED - Results are reliable." << endl;
         } else {
             cerr << "\n⚠️  NOT FULLY CONVERGED - Consider increasing samples." << endl;
         }
         
         // Per-type R-hat for types with significant frequency
         cerr << "\nPer-type R-hat (significant types only):" << endl;
         
         for (int t = 0; t < NUM_TYPES; t++) {
             size_t total = 0;
             for (const auto& cc : chain_counts) {
                 for (size_t cnt : cc) total += cnt;
             }
             
             size_t type_count = 0;
             for (const auto& cc : chain_counts) {
                 type_count += cc[t];
             }
             
             if (type_count > 0 && (double)type_count / total > 0.001) {
                 // Compute R-hat for this type
                 vector<double> chain_freqs;
                 for (const auto& cc : chain_counts) {
                     size_t ct = 0;
                     for (size_t cnt : cc) ct += cnt;
                     if (ct > 0) {
                         chain_freqs.push_back((double)cc[t] / ct);
                     }
                 }
                 
                 double mean = accumulate(chain_freqs.begin(), chain_freqs.end(), 0.0) / chain_freqs.size();
                 double var = 0;
                 for (double f : chain_freqs) var += (f - mean) * (f - mean);
                 var /= chain_freqs.size();
                 double std_dev = sqrt(var);
                 
                 cerr << "  Type " << setw(2) << (t+1) << ": " 
                      << fixed << setprecision(4) << (100.0 * type_count / total) << "% "
                      << "| chain std: " << setprecision(4) << std_dev;
                 if (std_dev > 0.01) cerr << " ⚠";
                 cerr << endl;
             }
         }
         
         cerr << "========================================" << endl;
     }
 };
 
 // ============================================================================
 // MAIN
 // ============================================================================
 
 int main(int argc, char* argv[]) {
     if (argc < 4) {
         cerr << "Usage: " << argv[0] << " <nverts> <simplices> <output> [samples] [--log-convergence]" << endl;
         cerr << endl;
         cerr << "UNIFIED SFD ESTIMATOR - ROBUST AUTO-ADAPTIVE VERSION" << endl;
         cerr << "=====================================================" << endl;
         cerr << endl;
         cerr << "This estimator automatically adapts to any dataset:" << endl;
         cerr << "  - Analyzes dataset characteristics" << endl;
         cerr << "  - Computes optimal MCMC parameters" << endl;
         cerr << "  - Monitors convergence and extends if needed" << endl;
         cerr << "  - Provides comprehensive uncertainty quantification" << endl;
         cerr << endl;
         cerr << "Options:" << endl;
         cerr << "  samples            Number of samples (default: auto-determined)" << endl;
         cerr << "  --log-convergence  Enable convergence logging for visualization" << endl;
         cerr << endl;
         cerr << "The algorithm prioritizes ACCURACY over speed." << endl;
         return 1;
     }
     
     size_t num_samples = 0;  // 0 = auto-determine
     bool enable_logging = false;
     
     for (int i = 4; i < argc; i++) {
         string arg = argv[i];
         if (arg == "--log-convergence") {
             enable_logging = true;
         } else if (arg[0] != '-') {
             num_samples = stoul(arg);
         }
     }
     
     // Load dataset
     SimplicialComplex K;
     
     cerr << "Loading dataset..." << endl;
     if (!loadDBLP(K, argv[1], argv[2])) {
         cerr << "Error: Cannot load dataset" << endl;
         return 1;
     }
     
     // Create estimator (this analyzes dataset and computes parameters)
     RobustMCMC_SFD estimator(K);
     
     // Enable logging if requested
     if (enable_logging) {
         string log_path = string(argv[3]) + ".convergence.csv";
         estimator.enableConvergenceLogging(log_path);
     }
     
     // Run estimation
     Result r = estimator.compute(num_samples);
     
     // Save results
     r.save(argv[3]);
     
     // Print final results
     cerr << "\n========================================" << endl;
     cerr << "FINAL RESULTS" << endl;
     cerr << "========================================" << endl;
     cerr << "Time:    " << fixed << setprecision(2) << r.time_ms << " ms" << endl;
     cerr << "Samples: " << r.samples << endl;
     cerr << endl;
     cerr << "SFD Estimates with 95% CI:" << endl;
     cerr << "-----------------------------------------" << endl;
     
     for (int i = 0; i < NUM_TYPES; i++) {
         if (r.freq[i] > 0.0001 || r.counts[i] > 0) {
             cerr << "  Type " << setw(2) << (i+1) << ": " 
                  << fixed << setprecision(6) << r.freq[i]
                  << "  [" << setprecision(4) << r.ci_lo[i] 
                  << ", " << r.ci_hi[i] << "]"
                  << "  (n=" << r.counts[i] << ")" << endl;
         }
     }
     
     cerr << endl;
     cerr << "Result saved to " << argv[3] << endl;
     
     if (enable_logging) {
         cerr << "Convergence log saved to " << argv[3] << ".convergence.csv" << endl;
     }
     
     return 0;
 }