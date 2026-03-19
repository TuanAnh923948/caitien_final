/**
 * ============================================================================
 * MCMC SFD Estimator - FAITHFUL REPRODUCTION OF ORIGINAL PAPER
 * ============================================================================
 * 
 * This is a faithful reimplementation of the MCMC algorithm from:
 * "Approximating Simplet Frequency Distribution for Simplicial Complexes"
 * (Beigy et al., EuroCG'24)
 * 
 * ALL BUGS FROM THE ORIGINAL CODE (complex.cpp) ARE PRESERVED INTENTIONALLY
 * for fair comparison purposes:
 * 
 * BUG 1 (hashSimplex): Edge existence check always returns true.
 *   The original checks "simplex.count(nodes[i]) && simplex.count(nodes[j])"
 *   which is always true since nodes[i] and nodes[j] are already in the simplex.
 *   This causes all degree sequences to be identical (complete graph),
 *   leading to incorrect classification of simplet types 5-18.
 * 
 * BUG 2 (randomWalk): nextNeighbors computed from currentSimplex instead of
 *   nextSimplex. This makes acceptance_probability = min(|N(s)|/|N(s)|, 1) = 1,
 *   so every proposed move is accepted, violating Metropolis-Hastings.
 * 
 * STRUCTURE: Matches original complex.cpp as closely as possible:
 * - Same global data structures
 * - Same hash function and simplet index map
 * - Same getNeighbors (add + remove + swap moves)
 * - Same isConnectedFast
 * - Same randomWalk logic
 * - Fixed 20 steps per sample (as in original)
 * - No burn-in by default (original had none)
 * 
 * Only adaptation: DBLP format loading for compatibility with comparison tools.
 * 
 * Usage: ./mcmc_sfd <nverts.txt> <simplices.txt> <output.res> [samples] [burnin]
 * ============================================================================
 */

 #include <iostream>
 #include <fstream>
 #include <vector>
 #include <set>
 #include <map>
 #include <queue>
 #include <algorithm>
 #include <random>
 #include <chrono>
 #include <iomanip>
 #include <cmath>
 #include <sstream>
 
 using namespace std;
 using namespace std::chrono;
 
 // ============================================================================
 // CONSTANTS (exactly from original)
 // ============================================================================
 
 const int NUM_TYPES = 18;
 const int MIN_SIZE = 2;
 const int MAX_SIZE = 4;
 const int PRIME = 7;
 
 // ============================================================================
 // GLOBAL DATA STRUCTURES (exactly from original complex.cpp)
 // ============================================================================
 
 vector<vector<int>> all_simplices;  // All loaded simplices (before filtering)
 vector<set<int>> all_nei;           // Adjacency list for all vertices
 
 vector<vector<int>> simplices;      // Filtered simplices (within N vertices)
 vector<set<int>> nei;               // Adjacency list (filtered)
 set<int> tri;                       // Triangle hashes
 set<int> quad;                      // Tetrahedron hashes
 int maxVertex = 0;
 
 // ============================================================================
 // SIMPLET INDEX MAP (exactly from original)
 // ============================================================================
 
 map<int, int> simpletIndexMap = {
     {2005, 1}, {3015, 2}, {3020, 3}, {3031, 4}, {4029, 5},
     {4039, 6}, {4050, 7}, {4046, 8}, {4057, 9}, {4068, 10},
     {4031, 11}, {4034, 12}, {4051, 13}, {4062, 14}, {4073, 15},
     {4084, 16}, {4095, 17}, {4108, 18}
 };
 
 int findSimpletIndex(int number) {
     auto it = simpletIndexMap.find(number);
     if (it != simpletIndexMap.end()) {
         return it->second;
     } else {
         return -1;
     }
 }
 
 // ============================================================================
 // HASH FUNCTIONS (exactly from original - WITH BUG)
 // ============================================================================
 
 int computeHash(const vector<int>& vertices) {
     int hashValue = 0;
     int power = 1;
     for (int v : vertices) {
         hashValue += v * power;
         power *= PRIME;
     }
     return hashValue;
 }
 
 /**
  * BUG 1: hashSimplex - edge check always returns true
  * 
  * The line "if (simplex.count(nodes[i]) && simplex.count(nodes[j]))"
  * checks if nodes[i] and nodes[j] are IN the simplex, which is always true
  * since they were extracted FROM the simplex.
  * 
  * It SHOULD check if there is an EDGE between nodes[i] and nodes[j] in the
  * adjacency list, i.e. "if (nei[nodes[i]].count(nodes[j]))"
  * 
  * Result: all degree sequences are identical (complete graph degrees),
  * so classification depends ONLY on triangle/tetrahedron counts.
  */
 int hashSimplex(const set<int>& simplex) {
     vector<int> nodes(simplex.begin(), simplex.end());
     int n = nodes.size();
     vector<int> degrees(n, 0);
 
     // BUG: This always returns true - checks membership, not adjacency
     for (int i = 0; i < n; ++i) {
         for (int j = i + 1; j < n; ++j) {
             if (simplex.count(nodes[i]) && simplex.count(nodes[j])) {
                 degrees[i]++;
                 degrees[j]++;
             }
         }
     }
 
     sort(degrees.begin(), degrees.end());
 
     const vector<int> primes = {2, 3, 5, 7};
     int hash = 0;
     for (int i = 0; i < n; ++i) {
         hash += degrees[i] * primes[i];
     }
 
     // Count filled triangles
     int triangle_count = 0;
     if (n >= 3) {
         for (int i = 0; i < n; ++i) {
             for (int j = i + 1; j < n; ++j) {
                 for (int k = j + 1; k < n; ++k) {
                     vector<int> triangle = {nodes[i], nodes[j], nodes[k]};
                     sort(triangle.begin(), triangle.end());
                     int tri_hash = computeHash(triangle);
                     if (tri.count(tri_hash)) {
                         triangle_count++;
                     }
                 }
             }
         }
     }
 
     // Count tetrahedra
     int tetrahedron_count = 0;
     if (n == 4) {
         vector<int> tetrahedron = {nodes[0], nodes[1], nodes[2], nodes[3]};
         sort(tetrahedron.begin(), tetrahedron.end());
         int quad_hash = computeHash(tetrahedron);
         if (quad.count(quad_hash)) {
             tetrahedron_count++;
         }
     }
 
     hash += triangle_count * 11;
     hash += tetrahedron_count * 13;
     hash += 1000 * n;
 
     return hash;
 }
 
 // ============================================================================
 // LOAD SIMPLICIAL COMPLEX (matching original complex.cpp flow)
 // ============================================================================
 
 /**
  * Step 1: Load all simplices (size 2-4) into all_simplices.
  * Matches original: "if (size <= 4) simplex.push_back(vertex)"
  */
 bool loadSimplices(const string& nverts_file, const string& simplices_file) {
     ifstream nv_in(nverts_file);
     ifstream sp_in(simplices_file);
 
     if (!nv_in.is_open() || !sp_in.is_open()) {
         cerr << "Error: Cannot open input files" << endl;
         return false;
     }
 
     vector<int> nverts_list;
     int nv;
     while (nv_in >> nv) {
         nverts_list.push_back(nv);
     }
     nv_in.close();
 
     cerr << "Read " << nverts_list.size() << " simplex sizes" << endl;
 
     map<int, int> size_dist;
     for (int n : nverts_list) { size_dist[n]++; }
     cerr << "Simplex size distribution:" << endl;
     for (auto& [sz, cnt] : size_dist) {
         cerr << "  Size " << sz << ": " << cnt << endl;
     }
 
     vector<int> all_vertices;
     int v;
     while (sp_in >> v) {
         all_vertices.push_back(v);
         maxVertex = max(maxVertex, v);
     }
     sp_in.close();
 
     cerr << "Read " << all_vertices.size() << " vertex indices" << endl;
     cerr << "Max vertex ID: " << maxVertex << endl;
 
     size_t pos = 0;
     for (size_t i = 0; i < nverts_list.size() && pos < all_vertices.size(); i++) {
         int size = nverts_list[i];
         if (pos + size > all_vertices.size()) break;
 
         if (size >= 2 && size <= 4) {
             vector<int> simplex;
             for (int j = 0; j < size; j++) {
                 simplex.push_back(all_vertices[pos + j]);
             }
             sort(simplex.begin(), simplex.end());
             all_simplices.push_back(simplex);
         }
         pos += size;
     }
 
     cerr << "Loaded " << all_simplices.size() << " simplices (size 2-4)" << endl;
     return true;
 }
 
 /**
  * Step 2: Build adjacency list from ALL simplices.
  * Exactly from original: construct_all_neighbors()
  */
 void constructAllNeighbors() {
     all_nei.resize(maxVertex + 1);
 
     for (const auto& simplex : all_simplices) {
         int size = simplex.size();
         for (int i = 0; i < size; i++) {
             for (int j = i + 1; j < size; j++) {
                 all_nei[simplex[i]].insert(simplex[j]);
                 all_nei[simplex[j]].insert(simplex[i]);
             }
         }
     }
 }
 
 /**
  * Step 3: BFS to select N connected vertices.
  * Exactly from original: selectConnectedVertices()
  */
 vector<int> selectConnectedVertices(int k) {
     int n = all_nei.size();
     for (int start = 0; start < n; ++start) {
         vector<bool> visited(n, false);
         vector<int> component;
         queue<int> q;
 
         q.push(start);
         visited[start] = true;
 
         while (!q.empty() && component.size() < (size_t)k) {
             int cur = q.front();
             q.pop();
             component.push_back(cur);
 
             for (int neighbor : all_nei[cur]) {
                 if (neighbor >= 0 && neighbor < n && !visited[neighbor]) {
                     visited[neighbor] = true;
                     q.push(neighbor);
                 }
             }
         }
 
         if (component.size() >= (size_t)k) {
             return vector<int>(component.begin(), component.begin() + k);
         }
     }
     return {};
 }
 
 /**
  * Step 4: Filter simplices to keep only those within selected vertices.
  * Exactly from original: filterSimplices()
  */
 void filterSimplices(const vector<int>& connectedVertices) {
     set<int> selectedVertices(connectedVertices.begin(), connectedVertices.end());
 
     for (const auto& simplex : all_simplices) {
         bool valid = true;
         for (int vertex : simplex) {
             if (selectedVertices.find(vertex) == selectedVertices.end()) {
                 valid = false;
                 break;
             }
         }
         if (valid) {
             simplices.push_back(simplex);
         }
     }
 }
 
 // ============================================================================
 // CONSTRUCT COMPLEX (exactly from original)
 // ============================================================================
 
 void constructComplex() {
     nei.resize(maxVertex + 1);
 
     for (const auto& simplex : simplices) {
         int size = simplex.size();
 
         // Process edges (pairs)
         for (int i = 0; i < size; i++) {
             for (int j = i + 1; j < size; j++) {
                 nei[simplex[i]].insert(simplex[j]);
                 nei[simplex[j]].insert(simplex[i]);
             }
         }
 
         // Process triangles (triples)
         if (size >= 3) {
             for (int i = 0; i < size; i++) {
                 for (int j = i + 1; j < size; j++) {
                     for (int k = j + 1; k < size; k++) {
                         vector<int> triVertices = {simplex[i], simplex[j], simplex[k]};
                         sort(triVertices.begin(), triVertices.end());
                         int triHash = computeHash(triVertices);
                         tri.insert(triHash);
                     }
                 }
             }
         }
 
         // Process tetrahedrons (quadruples)
         if (size >= 4) {
             for (int i = 0; i < size; i++) {
                 for (int j = i + 1; j < size; j++) {
                     for (int k = j + 1; k < size; k++) {
                         for (int l = k + 1; l < size; l++) {
                             vector<int> quadVertices = {simplex[i], simplex[j], simplex[k], simplex[l]};
                             sort(quadVertices.begin(), quadVertices.end());
                             int quadHash = computeHash(quadVertices);
                             quad.insert(quadHash);
                         }
                     }
                 }
             }
         }
     }
 
     // Count edges
     size_t num_edges = 0;
     for (size_t v = 0; v < nei.size(); v++) {
         num_edges += nei[v].size();
     }
     num_edges /= 2;
 
     cerr << "Complex constructed:" << endl;
     cerr << "  Vertices: " << nei.size() << endl;
     cerr << "  Simplices: " << simplices.size() << endl;
     cerr << "  Edges: " << num_edges << endl;
     cerr << "  Triangles: " << tri.size() << endl;
     cerr << "  Tetrahedra: " << quad.size() << endl;
 }
 
 // ============================================================================
 // RANDOM WALK (exactly from original complex.cpp, WITH BUGS)
 // ============================================================================
 
 int getRandomIndex(int maxIndex, mt19937& gen) {
     uniform_int_distribution<int> dist(0, maxIndex - 1);
     return dist(gen);
 }
 
 // Connectivity check (exactly from original)
 bool isConnectedFast(const vector<int>& simplex) {
     int sz = simplex.size();
     if (sz < 2) return false;
     if (sz == 2) {
         return nei[simplex[0]].count(simplex[1]);
     }
 
     int count = 0;
     for (int i = 0; i < sz; ++i) {
         for (int j = i + 1; j < sz; ++j) {
             if (nei[simplex[i]].count(simplex[j]))
                 count++;
         }
     }
 
     return count >= sz - 1;
 }
 
 /**
  * getNeighbors - exactly from original complex.cpp
  * 
  * Includes all 3 move types:
  * 1. Remove one vertex (if size > MIN_SIZE)
  * 2. Add one vertex (if size < MAX_SIZE)
  * 3. Swap one vertex (if MIN_SIZE < size < MAX_SIZE, i.e. size == 3)
  */
 vector<vector<int>> getNeighbors(const set<int>& currentSimplex) {
     vector<vector<int>> neighbors;
 
     vector<int> simplexVec(currentSimplex.begin(), currentSimplex.end());
 
     // 1. Remove one vertex (if size > 2)
     if (currentSimplex.size() > (size_t)MIN_SIZE) {
         for (int v : simplexVec) {
             vector<int> newSimplex = simplexVec;
             newSimplex.erase(remove(newSimplex.begin(), newSimplex.end(), v), newSimplex.end());
             if (isConnectedFast(newSimplex)) {
                 neighbors.push_back(newSimplex);
             }
         }
     }
 
     set<int> candidateVertices;
 
     // 2. Add one vertex (if size < 4)
     if (currentSimplex.size() < (size_t)MAX_SIZE) {
         for (int v : simplexVec) {
             for (int neighbor : nei[v]) {
                 if (currentSimplex.count(neighbor) == 0) {
                     candidateVertices.insert(neighbor);
                 }
             }
         }
 
         for (int newV : candidateVertices) {
             vector<int> newSimplex = simplexVec;
             newSimplex.push_back(newV);
             sort(newSimplex.begin(), newSimplex.end());
             if (isConnectedFast(newSimplex)) {
                 neighbors.push_back(newSimplex);
             }
         }
     }
 
     // 3. Swap one vertex (exactly from original: size > MIN_SIZE && size < MAX_SIZE)
     if (currentSimplex.size() > (size_t)MIN_SIZE && currentSimplex.size() < (size_t)MAX_SIZE) {
         for (int v : simplexVec) {
             for (int newV : candidateVertices) {
                 vector<int> newSimplex = simplexVec;
                 replace(newSimplex.begin(), newSimplex.end(), v, newV);
                 sort(newSimplex.begin(), newSimplex.end());
                 if (isConnectedFast(newSimplex)) {
                     neighbors.push_back(newSimplex);
                 }
             }
         }
     }
 
     return neighbors;
 }
 
 /**
  * randomWalk - exactly from original complex.cpp
  * 
  * BUG 2: Line "getNeighbors(currentSimplex)" should be "getNeighbors(nextSimplex)".
  * This makes acceptance_probability always = 1 (since neighbors.size() == nextNeighbors.size()),
  * so every proposed move is accepted unconditionally.
  */
 set<int> randomWalk(int steps, mt19937& gen) {
     set<int> currentSimplex;
 
     if (simplices.empty()) {
         cerr << "Error: No simplices available!" << endl;
         return currentSimplex;
     }
 
     // Pick a random starting simplex (size >= 2)
     int currentIndex;
     do {
         currentIndex = getRandomIndex(simplices.size(), gen);
     } while (simplices[currentIndex].size() <= 1);
     currentSimplex = set<int>(simplices[currentIndex].begin(), simplices[currentIndex].end());
 
     uniform_real_distribution<> dis(0.0, 1.0);
 
     vector<vector<int>> neighbors = getNeighbors(currentSimplex);
     for (int step = 0; step < steps; step++) {
 
         if (neighbors.empty()) {
             break;
         }
 
         int random = getRandomIndex(neighbors.size(), gen);
         set<int> nextSimplex = set<int>(neighbors[random].begin(), neighbors[random].end());
 
         // BUG: Uses currentSimplex instead of nextSimplex
         // Original: vector<vector<int>> nextNeighbors = getNeighbors(currentSimplex);
         vector<vector<int>> nextNeighbors = getNeighbors(currentSimplex);
 
         double acceptance_probability = min(static_cast<double>(neighbors.size()) /
                                             static_cast<double>(nextNeighbors.size()), 1.0);
         if (dis(gen) <= acceptance_probability) {
             currentSimplex = nextSimplex;
             neighbors = nextNeighbors;
         }
     }
 
     return currentSimplex;
 }
 
 // ============================================================================
 // RESULT STRUCTURE (compatible with comparison tool)
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
 
     bool save(const string& filename) {
         ofstream out(filename);
         if (!out.is_open()) return false;
 
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
 
         return true;
     }
 };
 
 // ============================================================================
 // APPROXIMATE SFD (matching original algorithm)
 // ============================================================================
 
 Result approximateSFD(size_t num_samples, size_t burn_in) {
     Result r;
     r.method = "MCMC_ORIGINAL_PAPER";
 
     auto t0 = high_resolution_clock::now();
 
     mt19937 gen(chrono::steady_clock::now().time_since_epoch().count());
 
     const int steps = 20;  // Fixed 20 steps (exactly from original)
 
     // Optional burn-in (original had none, but allow for testing)
     if (burn_in > 0) {
         cerr << "Burn-in: " << burn_in << " steps..." << endl;
         for (size_t i = 0; i < burn_in; i++) {
             randomWalk(steps, gen);
         }
     }
 
     // Debug: Track hash statistics
     map<int, int> hash_counts;
     size_t unmatched = 0;
 
     // Sampling phase
     for (size_t i = 0; i < num_samples; i++) {
         set<int> simplex = randomWalk(steps, gen);
         int hash = hashSimplex(simplex);
         hash_counts[hash]++;
 
         int idx = findSimpletIndex(hash);
         if (idx >= 1 && idx <= NUM_TYPES) {
             r.counts[idx - 1]++;
         } else {
             unmatched++;
         }
 
         // Progress report
         if ((i + 1) % 5000 == 0) {
             cerr << "  Sampled " << (i + 1) << "/" << num_samples << endl;
         }
     }
 
     // Debug output
     cerr << "\nHash statistics:" << endl;
     cerr << "  Total samples: " << num_samples << endl;
     cerr << "  Matched: " << (num_samples - unmatched) << endl;
     cerr << "  Unmatched: " << unmatched << " ("
          << fixed << setprecision(2) << (100.0 * unmatched / num_samples) << "%)" << endl;
 
     if (unmatched > 0) {
         cerr << "  Top unmatched hash values:" << endl;
         vector<pair<int,int>> sorted_hashes(hash_counts.begin(), hash_counts.end());
         sort(sorted_hashes.begin(), sorted_hashes.end(),
              [](auto& a, auto& b) { return a.second > b.second; });
 
         int shown = 0;
         for (auto& [hash, count] : sorted_hashes) {
             if (findSimpletIndex(hash) == -1) {
                 cerr << "    Hash " << hash << ": " << count << " times" << endl;
                 if (++shown >= 10) break;
             }
         }
     }
 
     // Compute frequencies
     r.samples = num_samples;
     r.total = num_samples;
 
     for (int i = 0; i < NUM_TYPES; i++) {
         r.freq[i] = (double)r.counts[i] / num_samples;
         // Original MCMC has no confidence intervals
         r.ci_lo[i] = r.freq[i];
         r.ci_hi[i] = r.freq[i];
     }
 
     r.time_ms = duration_cast<microseconds>(
         high_resolution_clock::now() - t0).count() / 1000.0;
 
     return r;
 }
 
 // ============================================================================
 // MAIN
 // ============================================================================
 
 int main(int argc, char* argv[]) {
     if (argc < 4) {
         cerr << "Usage: " << argv[0] << " <nverts.txt> <simplices.txt> <output.res> <N> [samples]" << endl;
         cerr << endl;
         cerr << "MCMC_ORIGINAL - Faithful reproduction of Beigy et al. (EuroCG'24)" << endl;
         cerr << "=================================================================" << endl;
         cerr << "All bugs from the original code are preserved for fair comparison." << endl;
         cerr << endl;
         cerr << "Arguments:" << endl;
         cerr << "  N         Number of connected vertices to select (0 = use all)" << endl;
         cerr << "  samples   Number of MCMC samples (default: 10000)" << endl;
         cerr << endl;
         cerr << "Known bugs preserved:" << endl;
         cerr << "  1. hashSimplex: edge check always true (incorrect classification)" << endl;
         cerr << "  2. randomWalk: acceptance always 1 (not proper M-H)" << endl;
         return 1;
     }
 
     string nverts_file = argv[1];
     string simplices_file = argv[2];
     string output_file = argv[3];
     int N = (argc > 4) ? atoi(argv[4]) : 0;  // 0 = use all vertices
     size_t num_samples = (argc > 5) ? atoi(argv[5]) : 10000;
 
     cerr << "========================================" << endl;
     cerr << "MCMC ORIGINAL PAPER IMPLEMENTATION" << endl;
     cerr << "========================================" << endl;
     cerr << "N (subgraph vertices): " << (N > 0 ? to_string(N) : "ALL") << endl;
     cerr << "Samples: " << num_samples << endl;
     cerr << "Steps per sample: 20 (fixed)" << endl;
     cerr << "========================================" << endl;
     cerr << endl;
 
     // Step 1: Load all simplices
     cerr << "Loading dataset..." << endl;
     if (!loadSimplices(nverts_file, simplices_file)) {
         cerr << "Error: Cannot load dataset" << endl;
         return 1;
     }
 
     // Step 2: Build full adjacency list
     cerr << "Building adjacency list..." << endl;
     constructAllNeighbors();
 
     // Step 3: Select N connected vertices (or use all)
     if (N > 0) {
         cerr << "Selecting " << N << " connected vertices (BFS)..." << endl;
         vector<int> connectedVertices = selectConnectedVertices(N);
         if (connectedVertices.empty()) {
             cerr << "Error: Cannot find " << N << " connected vertices" << endl;
             return 1;
         }
         cerr << "Selected " << connectedVertices.size() << " vertices" << endl;
 
         // Step 4: Filter simplices
         cerr << "Filtering simplices..." << endl;
         filterSimplices(connectedVertices);
         cerr << "Filtered: " << simplices.size() << " simplices" << endl;
     } else {
         // Use all simplices
         simplices = all_simplices;
         cerr << "Using all " << simplices.size() << " simplices" << endl;
     }
 
     // Step 5: Construct complex (edges, triangles, tetrahedra)
     cerr << "Constructing complex..." << endl;
     constructComplex();
 
     // Step 6: Run MCMC
     cerr << "\nRunning MCMC (" << num_samples << " samples)..." << endl;
     Result result = approximateSFD(num_samples, 0);  // 0 burn-in (matching original)
 
     // Save result
     result.save(output_file);
 
     // Print summary
     cerr << "\n========================================" << endl;
     cerr << "RESULTS" << endl;
     cerr << "========================================" << endl;
     cerr << "Time: " << fixed << setprecision(3) << result.time_ms << " ms" << endl;
     cerr << endl;
     cerr << "SFD (Simplet Frequency Distribution):" << endl;
     for (int i = 0; i < NUM_TYPES; i++) {
         if (result.counts[i] > 0) {
             cerr << "  Type " << setw(2) << (i+1) << ": "
                  << fixed << setprecision(6) << result.freq[i]
                  << " (n=" << result.counts[i] << ")" << endl;
         }
     }
     cerr << "\nResult saved to " << output_file << endl;
 
     return 0;
 }