/**
 * ============================================================================
 * EXACT SFD ENUMERATION
 * ============================================================================
 * 
 * This program computes the EXACT Simplet Frequency Distribution by
 * enumerating ALL simplets in the simplicial complex.
 * 
 * Purpose: Generate GROUND TRUTH for comparing approximate methods.
 * 
 * Complexity: O(n * d^3) where d is max degree
 * - Suitable for small/medium datasets
 * - Not scalable for large datasets (use for validation only)
 * 
 * Usage: ./exact_sfd <nverts_file> <simplices_file> <output_file>
 * 
 * Author: TuanAnh (Master's Thesis)
 * ============================================================================
 */

 #include "common.hpp"

 class ExactSFD {
     const SimplicialComplex& K;
     
 public:
     ExactSFD(const SimplicialComplex& K) : K(K) {}
     
     Result compute() {
         Result r;
         r.method = "EXACT";
         
         auto t0 = high_resolution_clock::now();
         
         // ====================================================================
         // TYPE 1: All edges (2-simplets)
         // ====================================================================
         r.counts[0] = K.edges.size();
         
         // ====================================================================
         // TYPES 2-4: All 3-simplets
         // ====================================================================
         // Enumerate all connected 3-vertex subsets
         // Use vertex-centric approach with deduplication
         
         for (size_t v = 0; v < K.n; v++) {
             vector<int> neighbors(K.adj[v].begin(), K.adj[v].end());
             
             for (size_t i = 0; i < neighbors.size(); i++) {
                 for (size_t j = i + 1; j < neighbors.size(); j++) {
                     int u = neighbors[i];
                     int w = neighbors[j];
                     
                     // Check if u-w edge exists
                     bool has_uw = K.hasEdge(u, w);
                     
                     if (!has_uw) {
                         // PATH: v is center vertex (v-u and v-w edges, no u-w)
                         // Count once per path (v is always the center)
                         r.counts[1]++;  // Type 2: Path
                     } else {
                         // TRIANGLE: All three edges exist
                         // Count only when v is the smallest vertex (deduplication)
                         if ((int)v < u && (int)v < w) {
                             if (K.hasTriangle(v, u, w)) {
                                 r.counts[3]++;  // Type 4: Filled triangle
                             } else {
                                 r.counts[2]++;  // Type 3: Empty triangle
                             }
                         }
                     }
                 }
             }
         }
         
         // ====================================================================
         // TYPES 5-18: All 4-simplets
         // ====================================================================
         // Use hash-based deduplication for efficiency
         
         auto hashSimplet = [](vector<int> s) -> uint64_t {
             sort(s.begin(), s.end());
             return ((uint64_t)s[0] << 48) | ((uint64_t)s[1] << 32) | 
                    ((uint64_t)s[2] << 16) | s[3];
         };
         
         set<uint64_t> seen4;
         
         for (size_t v = 0; v < K.n; v++) {
             vector<int> neighbors(K.adj[v].begin(), K.adj[v].end());
             
             // Case 1: v + 3 neighbors
             for (size_t i = 0; i < neighbors.size(); i++) {
                 for (size_t j = i + 1; j < neighbors.size(); j++) {
                     for (size_t k = j + 1; k < neighbors.size(); k++) {
                         vector<int> simplet = {(int)v, neighbors[i], neighbors[j], neighbors[k]};
                         uint64_t h = hashSimplet(simplet);
                         
                         if (!seen4.count(h) && Classifier::isConnected(simplet, K)) {
                             seen4.insert(h);
                             int type = Classifier::classify(simplet, K);
                             if (type >= 5 && type <= 18) {
                                 r.counts[type - 1]++;
                             }
                         }
                     }
                 }
             }
             
             // Case 2: v + 2 neighbors + 1 vertex at distance 2
             for (int u : neighbors) {
                 for (int w : K.adj[u]) {
                     if (w != (int)v && K.adj[v].count(w) == 0) {
                         // w is at distance 2 from v
                         for (int x : neighbors) {
                             if (x != u) {
                                 vector<int> simplet = {(int)v, u, w, x};
                                 uint64_t h = hashSimplet(simplet);
                                 
                                 if (!seen4.count(h) && Classifier::isConnected(simplet, K)) {
                                     seen4.insert(h);
                                     int type = Classifier::classify(simplet, K);
                                     if (type >= 5 && type <= 18) {
                                         r.counts[type - 1]++;
                                     }
                                 }
                             }
                         }
                     }
                 }
             }
         }
         
         // ====================================================================
         // COMPUTE FREQUENCIES
         // ====================================================================
         r.total = 0;
         for (int i = 0; i < NUM_TYPES; i++) {
             r.total += r.counts[i];
         }
         
         for (int i = 0; i < NUM_TYPES; i++) {
             r.freq[i] = r.total > 0 ? (double)r.counts[i] / r.total : 0;
             // Exact method has no confidence intervals (perfect accuracy)
             r.ci_lo[i] = r.freq[i];
             r.ci_hi[i] = r.freq[i];
         }
         
         r.samples = r.total;  // All simplets enumerated
         
         r.time_ms = duration_cast<microseconds>(
             high_resolution_clock::now() - t0).count() / 1000.0;
         
         return r;
     }
 };
 
 // ============================================================================
 // MAIN
 // ============================================================================
 
 int main(int argc, char* argv[]) {
     if (argc < 4) {
         cerr << "Usage: " << argv[0] << " <nverts_file> <simplices_file> <output_file>" << endl;
         cerr << endl;
         cerr << "EXACT SFD ENUMERATION" << endl;
         cerr << "=====================" << endl;
         cerr << "Computes exact simplet frequency distribution by" << endl;
         cerr << "enumerating all simplets. Use as ground truth." << endl;
         cerr << endl;
         cerr << "Warning: May be slow for large datasets!" << endl;
         return 1;
     }
     
     SimplicialComplex K;
     
     cerr << "Loading dataset..." << endl;
     if (!loadDBLP(K, argv[1], argv[2])) {
         cerr << "Error: Cannot load dataset" << endl;
         return 1;
     }
     
     cerr << "Running EXACT enumeration..." << endl;
     
     ExactSFD exact(K);
     Result r = exact.compute();
     
     r.save(argv[3]);
     
     cerr << "Result saved to " << argv[3] << endl;
     cerr << "Time: " << fixed << setprecision(3) << r.time_ms << " ms" << endl;
     cerr << "Total simplets: " << r.total << endl;
     cerr << endl;
     
     // Print summary
     cerr << "SFD (Simplet Frequency Distribution):" << endl;
     cerr << "-------------------------------------" << endl;
     for (int i = 0; i < NUM_TYPES; i++) {
         if (r.counts[i] > 0) {
             cerr << "  Type " << setw(2) << (i+1) << ": " 
                  << setw(8) << r.counts[i] << " ("
                  << fixed << setprecision(4) << (r.freq[i] * 100) << "%)" << endl;
         }
     }
     
     return 0;
 }