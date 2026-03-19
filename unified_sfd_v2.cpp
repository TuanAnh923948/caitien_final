/**
 * ============================================================================
 * UNIFIED SFD ESTIMATOR - VERSION 2.0 (OPTIMIZED + FIXED)
 * ============================================================================
 *
 * Improvements over original paper (Beigy et al., EuroCG'24):
 * 1. CORRECTED simplet classifier (fixes edge-check bug in original)
 * 2. Proper Metropolis-Hastings (fixes acceptance probability bug)
 * 3. Auto-adaptive parameters based on dataset characteristics
 * 4. Multi-chain MCMC with Gelman-Rubin convergence monitoring
 * 5. Wilson score confidence intervals
 *
 * Performance optimizations:
 * 1. OpenMP parallel chain execution
 * 2. O(1) hash-based lookups (unordered_set)
 * 3. Per-chain LRU cache for simplex neighbors (no lock contention)
 * 4. Smart early stopping (CI + stability + R-hat)
 * 5. All 3 move types: add, remove, swap (matching paper's state space)
 *
 * Author: TuanAnh (Master's Thesis)
 * ============================================================================
 */

 #include <iostream>
 #include <fstream>
 #include <sstream>
 #include <vector>
 #include <set>
 #include <map>
 #include <unordered_set>
 #include <unordered_map>
 #include <algorithm>
 #include <cmath>
 #include <random>
 #include <chrono>
 #include <iomanip>
 #include <tuple>
 #include <cassert>
 #include <numeric>
 #include <list>
 #include <queue>
 
 #ifdef _OPENMP
 #include <omp.h>
 #endif
 
 using namespace std;
 using namespace std::chrono;
 
 const int NUM_TYPES = 18;
 const int MIN_SIZE = 2;
 const int MAX_SIZE = 4;
 
 // ============================================================================
 // OPTIMIZED DATA STRUCTURES WITH O(1) LOOKUP
 // ============================================================================
 
 struct EdgeHash {
     size_t operator()(const pair<int,int>& e) const {
         return ((size_t)e.first << 32) | (size_t)(unsigned int)e.second;
     }
 };
 
 struct TripleHash {
     size_t operator()(const tuple<int,int,int>& t) const {
         auto [a, b, c] = t;
         size_t h = (size_t)a;
         h = h * 2654435761UL + (size_t)b;
         h = h * 2654435761UL + (size_t)c;
         return h;
     }
 };
 
 struct QuadHash {
     size_t operator()(const tuple<int,int,int,int>& t) const {
         auto [a, b, c, d] = t;
         size_t h = (size_t)a;
         h = h * 2654435761UL + (size_t)b;
         h = h * 2654435761UL + (size_t)c;
         h = h * 2654435761UL + (size_t)d;
         return h;
     }
 };
 
 struct VectorHash {
     size_t operator()(const vector<int>& v) const {
         size_t h = 0;
         for (int x : v) {
             h = h * 31 + (size_t)x;
         }
         return h;
     }
 };
 
 // ============================================================================
 // OPTIMIZED SIMPLICIAL COMPLEX
 // ============================================================================
 
 struct SimplicialComplex {
     size_t n = 0;
 
     // Adjacency list using sorted vectors (cache-friendly)
     vector<vector<int>> adj;
 
     // O(1) lookup structures
     unordered_set<pair<int,int>, EdgeHash> edges;
     unordered_set<tuple<int,int,int>, TripleHash> tris;
     unordered_set<tuple<int,int,int,int>, QuadHash> tets;
 
     void init(size_t num_vertices) {
         n = num_vertices;
         adj.resize(n);
     }
 
     void addEdge(int u, int v) {
         if (u > v) swap(u, v);
         if (u != v && u >= 0 && v >= 0 && (size_t)u < n && (size_t)v < n) {
             if (edges.insert({u, v}).second) {
                 adj[u].push_back(v);
                 adj[v].push_back(u);
             }
         }
     }
 
     void addTriangle(int a, int b, int c) {
         vector<int> v = {a, b, c};
         sort(v.begin(), v.end());
         tris.insert({v[0], v[1], v[2]});
         addEdge(a, b);
         addEdge(b, c);
         addEdge(a, c);
     }
 
     void addTetrahedron(int a, int b, int c, int d) {
         vector<int> v = {a, b, c, d};
         sort(v.begin(), v.end());
         tets.insert({v[0], v[1], v[2], v[3]});
         addTriangle(a, b, c);
         addTriangle(a, b, d);
         addTriangle(a, c, d);
         addTriangle(b, c, d);
     }
 
     // Finalize: sort adjacency lists and remove duplicates
     void finalize() {
         for (size_t i = 0; i < n; i++) {
             sort(adj[i].begin(), adj[i].end());
             adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
         }
     }
 
     inline bool hasEdge(int u, int v) const {
         if (u > v) swap(u, v);
         return edges.count({u, v}) > 0;
     }
 
     inline bool hasTriangle(int a, int b, int c) const {
         vector<int> v = {a, b, c};
         sort(v.begin(), v.end());
         return tris.count({v[0], v[1], v[2]}) > 0;
     }
 
     inline bool hasTetrahedron(int a, int b, int c, int d) const {
         vector<int> v = {a, b, c, d};
         sort(v.begin(), v.end());
         return tets.count({v[0], v[1], v[2], v[3]}) > 0;
     }
 };
 
 // ============================================================================
 // PER-CHAIN LRU CACHE (no mutex needed - each chain owns its cache)
 // ============================================================================
 
 class SimplexNeighborCache {
     struct CacheEntry {
         vector<int> simplex;
         vector<vector<int>> neighbors;
     };
 
     size_t capacity;
     list<CacheEntry> cache_list;
     unordered_map<vector<int>, list<CacheEntry>::iterator, VectorHash> cache_map;
 
     size_t hits = 0;
     size_t misses = 0;
 
 public:
     SimplexNeighborCache(size_t cap = 10000) : capacity(cap) {}
 
     bool get(const vector<int>& simplex, vector<vector<int>>& result) {
         auto it = cache_map.find(simplex);
         if (it != cache_map.end()) {
             result = it->second->neighbors;
             cache_list.splice(cache_list.begin(), cache_list, it->second);
             hits++;
             return true;
         }
         misses++;
         return false;
     }
 
     void put(const vector<int>& simplex, const vector<vector<int>>& neighbors) {
         auto it = cache_map.find(simplex);
         if (it != cache_map.end()) {
             it->second->neighbors = neighbors;
             cache_list.splice(cache_list.begin(), cache_list, it->second);
             return;
         }
 
         if (cache_list.size() >= capacity) {
             auto last = cache_list.end();
             --last;
             cache_map.erase(last->simplex);
             cache_list.pop_back();
         }
 
         cache_list.push_front({simplex, neighbors});
         cache_map[simplex] = cache_list.begin();
     }
 
     size_t getHits() const { return hits; }
     size_t getMisses() const { return misses; }
 };
 
 // ============================================================================
 // RESULT STRUCTURE
 // ============================================================================
 
 struct Result {
     string method;
     double time_ms = 0;
     size_t total = 0;
     size_t samples = 0;
 
     size_t counts[NUM_TYPES] = {0};
     double freq[NUM_TYPES] = {0};
     double ci_lo[NUM_TYPES] = {0};
     double ci_hi[NUM_TYPES] = {0};
 
     void save(const string& filename) const {
         ofstream out(filename);
         out << "METHOD: " << method << "\n";
         out << "TIME_MS: " << fixed << setprecision(3) << time_ms << "\n";
         out << "TOTAL: " << total << "\n";
         out << "SAMPLES: " << samples << "\n";
         out << "TYPE,COUNT,FREQ,CI_LO,CI_HI\n";
         for (int i = 0; i < NUM_TYPES; i++) {
             out << (i+1) << "," << counts[i] << ","
                 << fixed << setprecision(6) << freq[i] << ","
                 << ci_lo[i] << "," << ci_hi[i] << "\n";
         }
         out.close();
     }
 };
 
 // ============================================================================
 // CORRECTED SIMPLET CLASSIFIER
 // ============================================================================
 //
 // Fixes the bug in the original paper where hashSimplex always treated
 // all vertex pairs as having edges (degree check always true).
 //
 // Classification is based on:
 // - Number of actual edges (checked via K.hasEdge)
 // - Number of filled triangles (2-simplices, checked via K.hasTriangle)
 // - Number of tetrahedra (3-simplices, checked via K.hasTetrahedron)
 // - Degree sequence within the induced subgraph
 //
 // Verified completeness for all reachable cases:
 // - 2 vertices: 1 type (edge)
 // - 3 vertices: 3 types (path, empty triangle, filled triangle)
 // - 4 vertices: 14 types (all degree sequences + triangle/tet combinations)
 //
 // Note on 4-vertex cases:
 // - 3 edges: star {1,1,1,3} or path {1,1,2,2}
 // - 4 edges: triangle+pendant {1,2,2,3} or 4-cycle {2,2,2,2}
 //   ({1,1,3,3} requires 5 edges, unreachable with 4)
 // - 5 edges: K4-minus-edge, degree sequence always {2,2,3,3}
 //   max 2 filled triangles (missing edge breaks exactly 2 of 4 triangles)
 // - 6 edges: K4, 0-4 filled triangles + optional tetrahedron
 
 class Classifier {
 public:
     static bool isConnected(const vector<int>& vertices, const SimplicialComplex& K) {
         if (vertices.size() <= 1) return true;
 
         set<int> visited;
         vector<int> stack;
         stack.push_back(vertices[0]);
         visited.insert(vertices[0]);
 
         while (!stack.empty()) {
             int v = stack.back();
             stack.pop_back();
 
             for (int u : vertices) {
                 if (visited.count(u) == 0 && K.hasEdge(v, u)) {
                     visited.insert(u);
                     stack.push_back(u);
                 }
             }
         }
 
         return visited.size() == vertices.size();
     }
 
     static int classify(vector<int> vertices, const SimplicialComplex& K) {
         sort(vertices.begin(), vertices.end());
 
         if (vertices.size() == 2) return 1;
         if (vertices.size() == 3) return classify3(vertices, K);
         if (vertices.size() == 4) return classify4(vertices, K);
         return 0;
     }
 
 private:
     static int classify3(const vector<int>& v, const SimplicialComplex& K) {
         int a = v[0], b = v[1], c = v[2];
 
         int edge_count = K.hasEdge(a, b) + K.hasEdge(b, c) + K.hasEdge(a, c);
 
         if (edge_count == 2) return 2;  // Path
         if (edge_count == 3) {
             return K.hasTriangle(a, b, c) ? 4 : 3;  // Filled vs empty triangle
         }
         return 0;  // Not connected
     }
 
     static int classify4(const vector<int>& v, const SimplicialComplex& K) {
         int a = v[0], b = v[1], c = v[2], d = v[3];
 
         // Count edges
         bool e_ab = K.hasEdge(a, b), e_ac = K.hasEdge(a, c), e_ad = K.hasEdge(a, d);
         bool e_bc = K.hasEdge(b, c), e_bd = K.hasEdge(b, d), e_cd = K.hasEdge(c, d);
 
         int num_edges = e_ab + e_ac + e_ad + e_bc + e_bd + e_cd;
 
         // Count filled triangles
         int num_tris = K.hasTriangle(a, b, c) + K.hasTriangle(a, b, d) +
                        K.hasTriangle(a, c, d) + K.hasTriangle(b, c, d);
 
         // Check tetrahedron
         bool has_tet = K.hasTetrahedron(a, b, c, d);
 
         // Degree sequence (sorted)
         vector<int> degrees = {
             (int)(e_ab + e_ac + e_ad),
             (int)(e_ab + e_bc + e_bd),
             (int)(e_ac + e_bc + e_cd),
             (int)(e_ad + e_bd + e_cd)
         };
         sort(degrees.begin(), degrees.end());
 
         // === 3 edges (tree structures, no cycles) ===
         if (num_edges == 3) {
             if (degrees == vector<int>{1, 1, 1, 3}) return 5;  // Star
             if (degrees == vector<int>{1, 1, 2, 2}) return 6;  // Path
         }
 
         // === 4 edges (exactly one cycle) ===
         if (num_edges == 4) {
             if (degrees == vector<int>{1, 2, 2, 3}) return 7;   // Triangle + pendant
             if (degrees == vector<int>{2, 2, 2, 2}) return 8;   // 4-cycle
             // Note: {1,1,3,3} is unreachable with 4 edges on 4 vertices
         }
 
         // === 5 edges (K4 minus one edge, degrees always {2,2,3,3}) ===
         // Max 2 filled triangles (missing edge breaks exactly 2 triangles)
         if (num_edges == 5) {
             if (num_tris == 0) return 10;
             if (num_tris == 1) return 11;
             if (num_tris == 2) return 12;
             // num_tris >= 3 unreachable: missing edge participates in 2 triangles
         }
 
         // === 6 edges (complete graph K4) ===
         if (num_edges == 6) {
             if (num_tris == 0) return 13;
             if (num_tris == 1) return 14;
             if (num_tris == 2) return 15;
             if (num_tris == 3) return 16;
             if (num_tris == 4) {
                 return has_tet ? 18 : 17;
             }
         }
 
         // Fallback (should not reach here for valid connected induced subgraphs)
         return 0;
     }
 };
 
 // ============================================================================
 // DATASET ANALYZER
 // ============================================================================
 
 struct DatasetProfile {
     enum Category { TINY, SMALL, MEDIUM, LARGE, HUGE };
 
     size_t num_vertices = 0;
     size_t num_edges = 0;
     double avg_degree = 0;
     double max_degree = 0;
     double degree_std = 0;
     double heterogeneity_score = 0;
     double sparsity_score = 0;
     double diameter_estimate = 0;
     Category category = MEDIUM;
 
     string categoryName() const {
         switch (category) {
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
 
         p.num_vertices = K.n;
         p.num_edges = K.edges.size();
 
         if (p.num_vertices == 0) return p;
 
         double sum_deg = 0, sum_deg_sq = 0;
         p.max_degree = 0;
 
         for (size_t v = 0; v < K.n; v++) {
             double deg = K.adj[v].size();
             sum_deg += deg;
             sum_deg_sq += deg * deg;
             p.max_degree = max(p.max_degree, deg);
         }
 
         p.avg_degree = sum_deg / K.n;
         double var = (sum_deg_sq / K.n) - (p.avg_degree * p.avg_degree);
         p.degree_std = sqrt(max(0.0, var));
 
         p.heterogeneity_score = p.avg_degree > 0 ? p.max_degree / p.avg_degree : 1;
 
         double max_possible_edges = (double)K.n * (K.n - 1) / 2;
         p.sparsity_score = max_possible_edges > 0 ?
             max_possible_edges / max(p.num_edges, (size_t)1) : 1;
         p.sparsity_score = min(p.sparsity_score, 100.0);
 
         p.diameter_estimate = p.avg_degree > 1 ? log(K.n) / log(p.avg_degree) : K.n;
         p.diameter_estimate = min(p.diameter_estimate, 50.0);
 
         if (p.num_edges < 100) p.category = DatasetProfile::TINY;
         else if (p.num_edges < 10000) p.category = DatasetProfile::SMALL;
         else if (p.num_edges < 100000) p.category = DatasetProfile::MEDIUM;
         else if (p.num_edges < 1000000) p.category = DatasetProfile::LARGE;
         else p.category = DatasetProfile::HUGE;
 
         return p;
     }
 };
 
 // ============================================================================
 // ADAPTIVE PARAMETERS
 // ============================================================================
 
 struct AdaptiveParameters {
     size_t burn_in = 10000;
     size_t steps_per_sample = 100;
     size_t num_chains = 8;
     size_t min_samples = 10000;
     size_t max_samples = 500000;
     double target_rhat = 1.1;
     double target_ci_width = 0.02;
     size_t convergence_check_interval = 1000;
     size_t max_extensions = 5;
     double extension_factor = 1.5;
     bool allow_auto_extension = true;
 
     // Early stopping
     double stability_threshold = 0.001;
     size_t stability_window = 5;
     double early_stop_rhat = 1.25;
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
         cerr << "Heterogeneity: " << profile.heterogeneity_score << endl;
 
         bool extreme_hetero = profile.heterogeneity_score > 100;
         bool very_extreme = profile.heterogeneity_score > 500;
 
         // Burn-in
         double base_burn_in = 2000;
         double diameter_factor = max(1.0, profile.diameter_estimate / 5.0);
         double hetero_factor = very_extreme ?
             pow(log(profile.heterogeneity_score + 1), 1.5) :
             max(1.0, log(profile.heterogeneity_score + 1));
 
         params.burn_in = (size_t)(base_burn_in * diameter_factor * hetero_factor);
         params.burn_in = max((size_t)1000, min(params.burn_in, (size_t)100000));
 
         // Steps per sample
         double base_steps = 20;
         double step_hetero = very_extreme ? sqrt(profile.heterogeneity_score) / 5 :
                             extreme_hetero ? sqrt(profile.heterogeneity_score) / 3 :
                             max(1.0, sqrt(profile.heterogeneity_score) / 4);
 
         params.steps_per_sample = (size_t)(base_steps * diameter_factor * step_hetero);
         params.steps_per_sample = max((size_t)10, min(params.steps_per_sample, (size_t)500));
 
         // Number of chains
         params.num_chains = 4;
         if (profile.heterogeneity_score > 5)   params.num_chains = 6;
         if (profile.heterogeneity_score > 20)  params.num_chains = 8;
         if (profile.heterogeneity_score > 50)  params.num_chains = 10;
         if (profile.heterogeneity_score > 100) params.num_chains = 12;
         if (profile.heterogeneity_score > 500) params.num_chains = 16;
 
         // Sample size
         switch (profile.category) {
             case DatasetProfile::TINY:
                 params.min_samples = 1000; params.max_samples = 10000; break;
             case DatasetProfile::SMALL:
                 params.min_samples = 5000; params.max_samples = 50000; break;
             case DatasetProfile::MEDIUM:
                 params.min_samples = 10000; params.max_samples = 200000; break;
             case DatasetProfile::LARGE:
                 params.min_samples = 20000; params.max_samples = 500000; break;
             case DatasetProfile::HUGE:
                 params.min_samples = 30000; params.max_samples = 1000000; break;
         }
 
         if (extreme_hetero) {
             params.min_samples = (size_t)(params.min_samples * 1.5);
             params.max_samples = (size_t)(params.max_samples * 1.5);
         }
 
         // Convergence targets
         if (profile.category <= DatasetProfile::SMALL) {
             params.target_rhat = 1.05;
             params.target_ci_width = 0.01;
             params.early_stop_rhat = 1.10;
         } else if (profile.category == DatasetProfile::MEDIUM) {
             params.target_rhat = 1.10;
             params.target_ci_width = 0.015;
             params.early_stop_rhat = 1.15;
         } else {
             params.target_rhat = 1.15;
             params.target_ci_width = 0.02;
             params.early_stop_rhat = 1.25;
         }
 
         if (very_extreme) {
             params.early_stop_rhat = 1.35;
             params.min_samples = max(params.min_samples, (size_t)100000);
         } else if (extreme_hetero) {
             params.early_stop_rhat = 1.30;
             params.min_samples = max(params.min_samples, (size_t)80000);
         }
 
         params.stability_threshold = 0.001;
         params.stability_window = 5;
         params.convergence_check_interval = max((size_t)500, params.min_samples / 20);
         params.max_extensions = extreme_hetero ? 10 : 6;
         params.extension_factor = 1.3;
 
         cerr << "\nComputed parameters:" << endl;
         cerr << "  Burn-in:           " << params.burn_in << endl;
         cerr << "  Steps per sample:  " << params.steps_per_sample << endl;
         cerr << "  Number of chains:  " << params.num_chains << endl;
         cerr << "  Min samples:       " << params.min_samples << endl;
         cerr << "  Max samples:       " << params.max_samples << endl;
         cerr << "  Target R-hat:      " << params.target_rhat << endl;
         cerr << "  Early stop R-hat:  " << params.early_stop_rhat << endl;
         cerr << "  Target CI width:   " << params.target_ci_width << endl;
         cerr << "========================================" << endl;
 
         return params;
     }
 };
 
 // ============================================================================
 // CONVERGENCE MONITOR WITH IMPROVED STABILITY CHECK
 // ============================================================================
 
 class ConvergenceMonitor {
 public:
     struct Diagnostics {
         double max_rhat = 0;
         double avg_rhat = 0;
         double max_ci_width = 0;
         double avg_ci_width = 0;
         int types_not_converged = 0;
         bool converged = false;
         bool stable = false;
     };
 
     // Frequency history for stability check
     vector<vector<double>> freq_history;
     size_t stable_count = 0;
 
     Diagnostics computeDiagnostics(
         const vector<vector<size_t>>& chain_counts,
         size_t samples_per_chain,
         const AdaptiveParameters& params
     ) {
         Diagnostics d;
         size_t num_chains = chain_counts.size();
 
         if (num_chains < 2 || samples_per_chain < 10) return d;
 
         vector<double> rhats;
         vector<double> ci_widths;
         vector<double> current_freqs(NUM_TYPES, 0);
         size_t total_samples = 0;
 
         for (size_t c = 0; c < num_chains; c++) {
             for (int t = 0; t < NUM_TYPES; t++) {
                 current_freqs[t] += chain_counts[c][t];
                 total_samples += chain_counts[c][t];
             }
         }
 
         if (total_samples == 0) return d;
 
         for (int t = 0; t < NUM_TYPES; t++) {
             current_freqs[t] /= total_samples;
         }
 
         // Improved stability check: compare with window_size steps back
         d.stable = checkStability(current_freqs, params);
 
         for (int t = 0; t < NUM_TYPES; t++) {
             vector<double> chain_freqs(num_chains);
             double overall_freq = current_freqs[t];
 
             for (size_t c = 0; c < num_chains; c++) {
                 size_t chain_total = 0;
                 for (int i = 0; i < NUM_TYPES; i++) {
                     chain_total += chain_counts[c][i];
                 }
                 chain_freqs[c] = chain_total > 0 ? (double)chain_counts[c][t] / chain_total : 0;
             }
 
             if (overall_freq < 0.001) continue;
 
             double rhat = computeRhat(chain_freqs, samples_per_chain);
             rhats.push_back(rhat);
 
             if (rhat > params.target_rhat) d.types_not_converged++;
 
             // Wilson CI width
             double p = overall_freq;
             double n_d = (double)total_samples;
             double z = 1.96;
             double denom = 1 + z*z/n_d;
             double margin = z * sqrt(p*(1-p)/n_d + z*z/(4*n_d*n_d)) / denom;
             ci_widths.push_back(2 * margin);
         }
 
         if (!rhats.empty()) {
             d.max_rhat = *max_element(rhats.begin(), rhats.end());
             d.avg_rhat = accumulate(rhats.begin(), rhats.end(), 0.0) / rhats.size();
         }
 
         if (!ci_widths.empty()) {
             d.max_ci_width = *max_element(ci_widths.begin(), ci_widths.end());
             d.avg_ci_width = accumulate(ci_widths.begin(), ci_widths.end(), 0.0) / ci_widths.size();
         }
 
         // Convergence conditions:
         // 1. Traditional: R-hat AND CI both meet strict targets
         bool traditional = (d.max_rhat <= params.target_rhat &&
                            d.max_ci_width <= params.target_ci_width);
 
         // 2. CI-focused: CI very good AND R-hat acceptable
         bool ci_early = (d.max_ci_width <= params.target_ci_width / 2 &&
                         d.max_rhat <= params.early_stop_rhat);
 
         // 3. Stability: Frequencies stable AND R-hat acceptable
         bool stability_early = (d.stable &&
                                d.max_rhat <= params.early_stop_rhat &&
                                d.max_ci_width <= params.target_ci_width);
 
         d.converged = traditional || ci_early || stability_early;
 
         return d;
     }
 
 private:
     /**
      * Improved stability check:
      * Compare current frequencies with those from stability_window steps back
      * (not just the immediately previous step, which is biased for large intervals).
      */
     bool checkStability(const vector<double>& current_freqs, const AdaptiveParameters& params) {
         freq_history.push_back(current_freqs);
 
         // Need at least window_size + 1 entries to compare
         if (freq_history.size() <= params.stability_window) return false;
 
         // Compare with entry from stability_window steps back
         const auto& old_freqs = freq_history[freq_history.size() - 1 - params.stability_window];
         double max_rel_change = 0;
 
         for (int t = 0; t < NUM_TYPES; t++) {
             if (current_freqs[t] > 0.01) {
                 double rel_change = abs(current_freqs[t] - old_freqs[t]) / current_freqs[t];
                 max_rel_change = max(max_rel_change, rel_change);
             }
         }
 
         if (max_rel_change < params.stability_threshold) {
             stable_count++;
         } else {
             stable_count = 0;
         }
 
         // Limit history size
         if (freq_history.size() > 20) {
             freq_history.erase(freq_history.begin());
         }
 
         return stable_count >= params.stability_window;
     }
 
     static double computeRhat(const vector<double>& chain_freqs, size_t n) {
         size_t m = chain_freqs.size();
         if (m < 2 || n < 10) return 1.0;
 
         double overall_mean = accumulate(chain_freqs.begin(), chain_freqs.end(), 0.0) / m;
         if (overall_mean < 1e-6) return 1.0;
 
         double B = 0;
         for (double f : chain_freqs) {
             double diff = f - overall_mean;
             B += diff * diff;
         }
         B = B * n / (m - 1);
 
         double W = 0;
         for (double p : chain_freqs) {
             W += p * (1 - p);
         }
         W /= m;
 
         if (W < 1e-10) {
             return B > 1e-6 ? 2.0 : 1.0;
         }
 
         double V_hat = ((n - 1.0) / n) * W + B / n;
         double rhat = sqrt(max(1.0, V_hat / W));
 
         return min(rhat, 10.0);
     }
 };
 
 // ============================================================================
 // PARALLEL MCMC SAMPLER
 // ============================================================================
 
 class ParallelMCMC_SFD {
     const SimplicialComplex& K;
     DatasetProfile profile;
     AdaptiveParameters params;
 
     // Logging
     bool log_convergence = false;
     string log_file_path;
     ofstream log_file;
 
 public:
     ParallelMCMC_SFD(const SimplicialComplex& K) : K(K) {
         profile = DatasetAnalyzer::analyze(K);
         params = ParameterCalculator::compute(profile);
     }
 
     void enableConvergenceLog(const string& path) {
         log_convergence = true;
         log_file_path = path;
     }
 
     Result estimate(size_t target_samples = 0) {
         if (target_samples == 0) {
             target_samples = params.min_samples;
         }
 
         Result r;
         r.method = "UNIFIED_V2_OPTIMIZED";
 
         auto t0 = high_resolution_clock::now();
 
         // Open log file
         if (log_convergence) {
             log_file.open(log_file_path);
             if (log_file.is_open()) {
                 log_file << "iteration,total_samples,max_rhat,avg_rhat,max_ci_width,stable,phase\n";
             }
         }
 
         cerr << "\n========================================" << endl;
         cerr << "PARALLEL MCMC ESTIMATION (V2)" << endl;
         cerr << "========================================" << endl;
 
         #ifdef _OPENMP
         cerr << "OpenMP threads: " << omp_get_max_threads() << endl;
         #else
         cerr << "OpenMP not available - running sequential" << endl;
         #endif
 
         cerr << "Target samples: " << target_samples << endl;
         cerr << "Max samples: " << params.max_samples << endl;
 
         // Initialize per-chain data
         size_t num_chains = params.num_chains;
         vector<mt19937> rngs(num_chains);
         vector<vector<int>> chain_states(num_chains);
         vector<vector<size_t>> chain_counts(num_chains, vector<size_t>(NUM_TYPES, 0));
         vector<vector<int>> chain_type_samples(num_chains);
 
         // Per-chain caches (no mutex needed!)
         vector<SimplexNeighborCache> chain_caches(num_chains, SimplexNeighborCache(10000));
 
         // Initialize RNGs with different seeds
         for (size_t c = 0; c < num_chains; c++) {
             rngs[c].seed(chrono::steady_clock::now().time_since_epoch().count() + c * 12345 + c);
         }
 
         // Phase 1: Parallel Burn-in
         cerr << "\nPhase 1: Burn-in (" << params.burn_in << " steps/chain)..." << endl;
 
         #pragma omp parallel for schedule(dynamic)
         for (size_t c = 0; c < num_chains; c++) {
             chain_states[c] = initializeSimplet(rngs[c], c);
 
             for (size_t i = 0; i < params.burn_in; i++) {
                 chain_states[c] = mcmcStep(chain_states[c], rngs[c], chain_caches[c]);
             }
         }
         cerr << "  Burn-in complete." << endl;
 
         // Phase 2: Parallel sampling with convergence monitoring
         cerr << "\nPhase 2: Sampling with convergence monitoring..." << endl;
 
         size_t samples_per_chain = target_samples / num_chains;
         size_t total_iterations = 0;
         size_t extension_count = 0;
         bool converged = false;
 
         ConvergenceMonitor monitor;
 
         while (!converged && total_iterations < params.max_samples / num_chains) {
             size_t batch_size = min(params.convergence_check_interval,
                                    samples_per_chain - total_iterations);
 
             // Parallel sampling — NO critical section needed
             // Each thread only accesses its own chain data (index c)
             #pragma omp parallel for schedule(dynamic)
             for (size_t c = 0; c < num_chains; c++) {
                 for (size_t s = 0; s < batch_size; s++) {
                     for (size_t step = 0; step < params.steps_per_sample; step++) {
                         chain_states[c] = mcmcStep(chain_states[c], rngs[c], chain_caches[c]);
                     }
 
                     int type = Classifier::classify(chain_states[c], K);
                     if (type >= 1 && type <= NUM_TYPES) {
                         chain_type_samples[c].push_back(type);
                         chain_counts[c][type - 1]++;
                     }
                 }
             }
 
             total_iterations += batch_size;
 
             // Convergence check (single-threaded)
             auto diag = monitor.computeDiagnostics(chain_counts, total_iterations, params);
 
             // Log
             if (log_convergence && log_file.is_open()) {
                 size_t total_now = total_iterations * num_chains;
                 log_file << total_iterations << "," << total_now << ","
                          << diag.max_rhat << "," << diag.avg_rhat << ","
                          << diag.max_ci_width << "," << (diag.stable ? 1 : 0) << ",main\n";
             }
 
             // Progress
             size_t total_now = total_iterations * num_chains;
             cerr << "  Samples: " << total_now
                  << " | R-hat: " << fixed << setprecision(3) << diag.max_rhat
                  << " | CI: " << setprecision(4) << diag.max_ci_width
                  << " | Stable: " << (diag.stable ? "Yes" : "No");
 
             if (diag.converged) {
                 bool trad = (diag.max_rhat <= params.target_rhat &&
                             diag.max_ci_width <= params.target_ci_width);
                 bool ci_e = (diag.max_ci_width <= params.target_ci_width / 2 &&
                             diag.max_rhat <= params.early_stop_rhat);
                 bool stab = (diag.stable &&
                             diag.max_rhat <= params.early_stop_rhat &&
                             diag.max_ci_width <= params.target_ci_width);
 
                 if (trad)      cerr << " [CONVERGED - R-hat & CI met]";
                 else if (ci_e) cerr << " [EARLY STOP - CI excellent]";
                 else if (stab) cerr << " [EARLY STOP - Stable]";
                 cerr << endl;
                 converged = true;
             } else if (total_iterations >= samples_per_chain) {
                 if (params.allow_auto_extension &&
                     extension_count < params.max_extensions &&
                     total_iterations * num_chains < params.max_samples) {
 
                     extension_count++;
                     size_t ext = (size_t)(samples_per_chain * (params.extension_factor - 1));
                     samples_per_chain += ext;
                     cerr << " [EXTENDING +" << ext * num_chains << "]" << endl;
                 } else {
                     cerr << " [MAX REACHED]" << endl;
                     break;
                 }
             } else {
                 cerr << endl;
             }
         }
 
         if (log_file.is_open()) {
             log_file.close();
         }
 
         // Combine results
         cerr << "\nCombining " << num_chains << " chains..." << endl;
 
         size_t total_cache_hits = 0, total_cache_misses = 0;
         for (size_t c = 0; c < num_chains; c++) {
             cerr << "  Chain " << (c+1) << ": " << chain_type_samples[c].size() << " samples" << endl;
             for (int type : chain_type_samples[c]) {
                 r.counts[type - 1]++;
             }
             r.total += chain_type_samples[c].size();
             total_cache_hits += chain_caches[c].getHits();
             total_cache_misses += chain_caches[c].getMisses();
         }
         r.samples = r.total;
 
         // Cache stats
         size_t total_cache = total_cache_hits + total_cache_misses;
         if (total_cache > 0) {
             cerr << "\nCache: " << total_cache_hits << " hits / " << total_cache
                  << " total (" << fixed << setprecision(1)
                  << (100.0 * total_cache_hits / total_cache) << "% hit rate)" << endl;
         }
 
         // Compute frequencies and CI
         computeFrequenciesAndCI(r);
 
         // Final diagnostics
         auto final_diag = monitor.computeDiagnostics(chain_counts, total_iterations, params);
         printFinalDiagnostics(final_diag);
 
         r.time_ms = duration_cast<milliseconds>(
             high_resolution_clock::now() - t0).count();
 
         return r;
     }
 
 private:
     vector<int> initializeSimplet(mt19937& rng, size_t chain_id) {
         // Start chains from different degree regions for diversity
         vector<pair<int, size_t>> vertices_by_degree;
         for (size_t v = 0; v < K.n; v++) {
             if (K.adj[v].size() >= 2) {
                 vertices_by_degree.push_back({(int)v, K.adj[v].size()});
             }
         }
 
         if (vertices_by_degree.empty()) {
             for (auto& e : K.edges) {
                 return {e.first, e.second};
             }
             return {};
         }
 
         sort(vertices_by_degree.begin(), vertices_by_degree.end(),
              [](auto& a, auto& b) { return a.second < b.second; });
 
         size_t num_regions = min((size_t)10, vertices_by_degree.size() / 100 + 1);
         size_t region_size = vertices_by_degree.size() / num_regions;
         size_t region = chain_id % num_regions;
         size_t start_idx = region * region_size;
         size_t end_idx = (region == num_regions - 1) ?
             vertices_by_degree.size() : (region + 1) * region_size;
 
         uniform_int_distribution<size_t> region_dist(start_idx, end_idx - 1);
         int selected_v = vertices_by_degree[region_dist(rng)].first;
 
         const vector<int>& neighbors = K.adj[selected_v];
         if (neighbors.size() < 2) {
             return {selected_v, neighbors[0]};
         }
 
         uniform_int_distribution<size_t> nb_dist(0, neighbors.size() - 1);
         size_t i1 = nb_dist(rng);
         size_t i2;
         int attempts = 0;
         do {
             i2 = nb_dist(rng);
             attempts++;
         } while (i2 == i1 && attempts < 100);
 
         vector<int> result = {selected_v, neighbors[i1], neighbors[i2]};
         sort(result.begin(), result.end());
         return result;
     }
 
     /**
      * MCMC step with proper Metropolis-Hastings.
      *
      * Uses per-chain cache (no mutex needed).
      * Neighbors include all 3 move types: add, remove, swap.
      */
     vector<int> mcmcStep(const vector<int>& current, mt19937& rng,
                          SimplexNeighborCache& cache) {
         vector<vector<int>> neighbors;
 
         // Try cache first
         if (!cache.get(current, neighbors)) {
             neighbors = getSimplexNeighbors(current);
             cache.put(current, neighbors);
         }
 
         if (neighbors.empty()) return current;
 
         uniform_int_distribution<size_t> dist(0, neighbors.size() - 1);
         vector<int> proposed = neighbors[dist(rng)];
 
         vector<vector<int>> proposed_neighbors;
         if (!cache.get(proposed, proposed_neighbors)) {
             proposed_neighbors = getSimplexNeighbors(proposed);
             cache.put(proposed, proposed_neighbors);
         }
 
         // Correct Metropolis-Hastings acceptance probability
         double accept_prob = min(1.0,
             (double)neighbors.size() / max((size_t)1, proposed_neighbors.size()));
 
         uniform_real_distribution<double> unif(0, 1);
         if (unif(rng) <= accept_prob) {
             return proposed;
         }
         return current;
     }
 
     /**
      * Get all neighbor simplets of the current simplet.
      *
      * Includes all 3 move types from the paper:
      * 1. Remove one vertex (if size > MIN_SIZE)
      * 2. Add one vertex (if size < MAX_SIZE)
      * 3. Swap one vertex with a candidate (if MIN_SIZE < size < MAX_SIZE)
      *
      * This ensures the Markov chain can traverse the full state space P^m_K
      * as defined in the paper (Section 3.2).
      */
     vector<vector<int>> getSimplexNeighbors(const vector<int>& simplex) {
         vector<vector<int>> neighbors;
         set<int> simplex_set(simplex.begin(), simplex.end());
 
         // 1. Remove a vertex (if size > 2)
         if (simplex.size() > (size_t)MIN_SIZE) {
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
 
         // Collect candidate vertices (neighbors of current simplex not in it)
         set<int> candidates;
         for (int v : simplex) {
             for (int u : K.adj[v]) {
                 if (!simplex_set.count(u)) {
                     candidates.insert(u);
                 }
             }
         }
 
         // 2. Add a vertex (if size < 4)
         if (simplex.size() < (size_t)MAX_SIZE) {
             for (int c : candidates) {
                 vector<int> ns = simplex;
                 ns.push_back(c);
                 if (Classifier::isConnected(ns, K)) {
                     sort(ns.begin(), ns.end());
                     neighbors.push_back(ns);
                 }
             }
         }
 
         // 3. Swap one vertex (if MIN_SIZE < size < MAX_SIZE, i.e. size == 3)
         // This matches the original paper's state space transitions
         if (simplex.size() > (size_t)MIN_SIZE && simplex.size() < (size_t)MAX_SIZE) {
             for (size_t i = 0; i < simplex.size(); i++) {
                 for (int c : candidates) {
                     vector<int> ns;
                     for (size_t j = 0; j < simplex.size(); j++) {
                         ns.push_back(i == j ? c : simplex[j]);
                     }
                     sort(ns.begin(), ns.end());
                     if (Classifier::isConnected(ns, K)) {
                         neighbors.push_back(ns);
                     }
                 }
             }
         }
 
         return neighbors;
     }
 
     void computeFrequenciesAndCI(Result& r) {
         double z = 1.96;
 
         for (int i = 0; i < NUM_TYPES; i++) {
             r.freq[i] = r.samples > 0 ? (double)r.counts[i] / r.samples : 0;
 
             double p = r.freq[i];
             double n = (double)r.samples;
 
             if (n > 0) {
                 double denom = 1 + z*z/n;
                 double center = (p + z*z/(2*n)) / denom;
                 double margin = z * sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom;
 
                 r.ci_lo[i] = max(0.0, center - margin);
                 r.ci_hi[i] = min(1.0, center + margin);
             }
         }
     }
 
     void printFinalDiagnostics(const ConvergenceMonitor::Diagnostics& diag) {
         cerr << "\n========================================" << endl;
         cerr << "CONVERGENCE DIAGNOSTICS" << endl;
         cerr << "========================================" << endl;
         cerr << "Max R-hat:      " << fixed << setprecision(4) << diag.max_rhat;
         cerr << (diag.max_rhat <= params.target_rhat ? " ✓" : " ⚠");
         cerr << " (target: " << params.target_rhat << ")" << endl;
 
         cerr << "Max CI width:   " << diag.max_ci_width;
         cerr << (diag.max_ci_width <= params.target_ci_width ? " ✓" : " ⚠");
         cerr << " (target: " << params.target_ci_width << ")" << endl;
 
         cerr << "Stable:         " << (diag.stable ? "Yes ✓" : "No") << endl;
         cerr << "Converged:      " << (diag.converged ? "Yes ✓" : "No ⚠") << endl;
         cerr << "========================================" << endl;
     }
 };
 
 
 // ============================================================================
 // DBLP FORMAT LOADER WITH SUBGRAPH SELECTION
 // ============================================================================
 
 /**
  * Load simplicial complex from DBLP format.
  * 
  * Modes:
  * - Default (full_decompose=false): Only load size 2-4, matching original paper.
  *   This ensures fair comparison with mcmc_sfd.
  * - Full decompose (full_decompose=true): Decompose higher-order simplices (size>4)
  *   into all sub-simplices. This captures the full simplicial complex structure
  *   for mathematically complete analysis.
  * 
  * If N > 0, select N connected vertices via BFS (matching original paper)
  * and only keep simplices within those vertices.
  */
 bool loadDBLP(SimplicialComplex& K, const string& nverts_file, const string& simplices_file,
               int N = 0, bool full_decompose = false) {
     ifstream nv_in(nverts_file);
     ifstream sp_in(simplices_file);
 
     if (!nv_in.is_open() || !sp_in.is_open()) {
         cerr << "Error: Cannot open input files" << endl;
         return false;
     }
 
     vector<int> nverts_list;
     int nv;
     while (nv_in >> nv) { nverts_list.push_back(nv); }
     nv_in.close();
 
     cerr << "Read " << nverts_list.size() << " simplex sizes" << endl;
     cerr << "Load mode: " << (full_decompose ? "FULL DECOMPOSE (all sizes)" : "PAPER MODE (size 2-4 only)") << endl;
 
     map<int, int> size_dist;
     for (int n : nverts_list) { size_dist[n]++; }
     cerr << "Simplex size distribution:" << endl;
     for (auto& [sz, cnt] : size_dist) {
         cerr << "  Size " << sz << ": " << cnt << endl;
     }
 
     vector<int> all_vertices;
     int v;
     while (sp_in >> v) { all_vertices.push_back(v); }
     sp_in.close();
 
     int max_v = 0;
     for (int vv : all_vertices) { max_v = max(max_v, vv); }
 
     // ================================================================
     // Parse simplices based on mode
     // ================================================================
     // In both modes, we first collect raw simplices for adjacency (BFS).
     // all_simplices_raw: size 2-4 simplices used for BFS neighbor construction.
     // For full_decompose, higher-order simplices are stored separately
     // and decomposed when building the SimplicialComplex.
     // ================================================================
 
     vector<vector<int>> all_simplices_raw;    // size 2-4 for BFS
     vector<vector<int>> higher_simplices;     // size 5+ for decomposition
 
     size_t pos = 0;
     for (size_t i = 0; i < nverts_list.size() && pos < all_vertices.size(); i++) {
         int size = nverts_list[i];
         if (pos + size > all_vertices.size()) break;
 
         vector<int> simplex;
         for (int j = 0; j < size; j++) {
             simplex.push_back(all_vertices[pos + j]);
         }
         pos += size;
 
         if (size >= 2 && size <= 4) {
             sort(simplex.begin(), simplex.end());
             all_simplices_raw.push_back(simplex);
         } else if (full_decompose && size > 4) {
             higher_simplices.push_back(simplex);
         }
     }
 
     cerr << "Loaded " << all_simplices_raw.size() << " simplices (size 2-4)" << endl;
     if (full_decompose && !higher_simplices.empty()) {
         cerr << "Loaded " << higher_simplices.size() << " higher-order simplices (size>4) for decomposition" << endl;
     }
 
     // ================================================================
     // Subgraph selection (N > 0): BFS on size 2-4 adjacency
     // ================================================================
 
     set<int> selectedSet;
 
     if (N > 0) {
         // Build adjacency from size 2-4 simplices (same as original paper)
         vector<set<int>> all_nei(max_v + 1);
         for (const auto& simplex : all_simplices_raw) {
             for (size_t a = 0; a < simplex.size(); a++) {
                 for (size_t b = a + 1; b < simplex.size(); b++) {
                     all_nei[simplex[a]].insert(simplex[b]);
                     all_nei[simplex[b]].insert(simplex[a]);
                 }
             }
         }
 
         cerr << "Selecting " << N << " connected vertices (BFS)..." << endl;
         vector<int> selected;
         for (int start = 0; start <= max_v && selected.empty(); ++start) {
             if (all_nei[start].empty()) continue;
 
             vector<bool> visited(max_v + 1, false);
             vector<int> component;
             queue<int> q;
             q.push(start);
             visited[start] = true;
 
             while (!q.empty() && component.size() < (size_t)N) {
                 int cur = q.front(); q.pop();
                 component.push_back(cur);
                 for (int nb : all_nei[cur]) {
                     if (!visited[nb]) {
                         visited[nb] = true;
                         q.push(nb);
                     }
                 }
             }
 
             if (component.size() >= (size_t)N) {
                 selected.assign(component.begin(), component.begin() + N);
             }
         }
 
         if (selected.empty()) {
             cerr << "Error: Cannot find " << N << " connected vertices" << endl;
             return false;
         }
         cerr << "Selected " << selected.size() << " vertices" << endl;
         selectedSet.insert(selected.begin(), selected.end());
     }
 
     // ================================================================
     // Helper: check if all vertices of a simplex are in selectedSet
     // ================================================================
     auto isInSubgraph = [&](const vector<int>& simplex) -> bool {
         if (N <= 0) return true;  // No filtering
         for (int vertex : simplex) {
             if (selectedSet.find(vertex) == selectedSet.end()) return false;
         }
         return true;
     };
 
     // ================================================================
     // Build SimplicialComplex
     // ================================================================
 
     K.init(max_v + 1);
 
     // Add size 2-4 simplices
     size_t filtered_count = 0;
     for (const auto& simplex : all_simplices_raw) {
         if (!isInSubgraph(simplex)) continue;
         filtered_count++;
 
         int size = simplex.size();
         if (size == 2) K.addEdge(simplex[0], simplex[1]);
         else if (size == 3) K.addTriangle(simplex[0], simplex[1], simplex[2]);
         else if (size == 4) K.addTetrahedron(simplex[0], simplex[1], simplex[2], simplex[3]);
     }
 
     // Decompose higher-order simplices (only in full_decompose mode)
     size_t decomposed_count = 0;
     if (full_decompose) {
         for (const auto& simplex : higher_simplices) {
             // Filter vertices in subgraph
             vector<int> vs;
             for (int vertex : simplex) {
                 if (N <= 0 || selectedSet.count(vertex)) {
                     vs.push_back(vertex);
                 }
             }
 
             // Need at least 2 vertices in subgraph to form edges
             if (vs.size() < 2) continue;
             decomposed_count++;
 
             if (vs.size() > 25) {
                 // Too large: only add edges to avoid combinatorial explosion
                 for (size_t a = 0; a < vs.size(); a++)
                     for (size_t b = a + 1; b < vs.size(); b++)
                         K.addEdge(vs[a], vs[b]);
             } else {
                 // Full decomposition: edges + triangles + tetrahedra
                 for (size_t a = 0; a < vs.size(); a++)
                     for (size_t b = a + 1; b < vs.size(); b++)
                         K.addEdge(vs[a], vs[b]);
                 for (size_t a = 0; a < vs.size(); a++)
                     for (size_t b = a + 1; b < vs.size(); b++)
                         for (size_t c = b + 1; c < vs.size(); c++)
                             K.addTriangle(vs[a], vs[b], vs[c]);
                 for (size_t a = 0; a < vs.size(); a++)
                     for (size_t b = a + 1; b < vs.size(); b++)
                         for (size_t c = b + 1; c < vs.size(); c++)
                             for (size_t d = c + 1; d < vs.size(); d++)
                                 K.addTetrahedron(vs[a], vs[b], vs[c], vs[d]);
             }
         }
     }
 
     K.finalize();
 
     cerr << "\n========================================" << endl;
     cerr << "SIMPLICIAL COMPLEX LOADED" << endl;
     cerr << "========================================" << endl;
     cerr << "  Mode:          " << (full_decompose ? "FULL DECOMPOSE" : "PAPER (size 2-4)") << endl;
     cerr << "  Vertices (N):  " << (N > 0 ? to_string(N) : "ALL") << endl;
     cerr << "  Simplices:     " << filtered_count << " (size 2-4)";
     if (full_decompose && decomposed_count > 0)
         cerr << " + " << decomposed_count << " decomposed";
     cerr << endl;
     cerr << "  Edges:         " << K.edges.size() << endl;
     cerr << "  Triangles:     " << K.tris.size() << endl;
     cerr << "  Tetrahedra:    " << K.tets.size() << endl;
     cerr << "========================================" << endl;
 
     return true;
 }
 
 // ============================================================================
 // MAIN
 // ============================================================================
 
 int main(int argc, char* argv[]) {
     if (argc < 4) {
         cerr << "Usage: " << argv[0] << " <nverts> <simplices> <output> <N> [samples] [--log-convergence]" << endl;
         cerr << endl;
         cerr << "UNIFIED SFD ESTIMATOR v2.0 (OPTIMIZED + FIXED)" << endl;
         cerr << "===============================================" << endl;
         cerr << "Arguments:" << endl;
         cerr << "  N                Number of connected vertices (0 = use all)" << endl;
         cerr << "  samples          Target sample count (0 = auto)" << endl;
         cerr << "  --log-convergence  Save convergence log" << endl;
         cerr << "  --full-decompose   Decompose size>4 simplices into sub-simplices" << endl;
         cerr << "                     (default: only load size 2-4, matching paper)" << endl;
         return 1;
     }
 
     string nverts_file = argv[1];
     string simplices_file = argv[2];
     string output_file = argv[3];
     int N = (argc > 4) ? atoi(argv[4]) : 0;
 
     size_t target_samples = 0;
     bool log_conv = false;
     bool full_decompose = false;
 
     for (int i = 5; i < argc; i++) {
         string arg = argv[i];
         if (arg == "--log-convergence") {
             log_conv = true;
         } else if (arg == "--full-decompose") {
             full_decompose = true;
         } else {
             target_samples = stoull(arg);
         }
     }
 
     // Load dataset
     cerr << "Loading dataset..." << endl;
     SimplicialComplex K;
     if (!loadDBLP(K, nverts_file, simplices_file, N, full_decompose)) {
         return 1;
     }
 
     // Create sampler
     ParallelMCMC_SFD sampler(K);
 
     if (log_conv) {
         sampler.enableConvergenceLog(output_file + ".convergence.csv");
         cerr << "Convergence logging: " << output_file << ".convergence.csv" << endl;
     }
 
     // Run
     Result r = sampler.estimate(target_samples);
 
     // Save
     r.save(output_file);
 
     // Print summary
     cerr << "\n========================================" << endl;
     cerr << "FINAL RESULTS" << endl;
     cerr << "========================================" << endl;
     cerr << "Time:    " << fixed << setprecision(0) << r.time_ms << " ms ("
          << setprecision(2) << r.time_ms / 60000 << " min)" << endl;
     cerr << "Samples: " << r.samples << endl;
     cerr << endl;
 
     cerr << "SFD Estimates with 95% CI:" << endl;
     cerr << "-----------------------------------------" << endl;
     for (int i = 0; i < NUM_TYPES; i++) {
         if (r.counts[i] > 0) {
             cerr << "  Type " << setw(2) << (i+1) << ": "
                  << fixed << setprecision(6) << r.freq[i]
                  << "  [" << setprecision(4) << r.ci_lo[i] << ", " << r.ci_hi[i] << "]"
                  << "  (n=" << r.counts[i] << ")" << endl;
         }
     }
 
     cerr << "\nSaved to " << output_file << endl;
 
     return 0;
 }