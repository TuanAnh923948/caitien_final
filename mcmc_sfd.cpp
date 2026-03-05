/**
 * MCMC SFD Estimator - ORIGINAL IMPLEMENTATION FROM PAPER
 * ========================================================
 * 
 * This is a faithful reimplementation of the MCMC algorithm from:
 * "Simplet Frequency Distribution Estimation" (Beigy et al., EuroCG'24)
 * 
 * THE CODE IS KEPT AS CLOSE AS POSSIBLE TO THE ORIGINAL TO ENSURE
 * FAIR COMPARISON. ALL BUGS ARE PRESERVED INTENTIONALLY.
 * 
 * Known issues in original:
 * 1. Edge existence check always returns true (line: if simplex.count(nodes[i]) && simplex.count(nodes[j]))
 * 2. This causes incorrect degree counting within simplex
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
 
 using namespace std;
 using namespace std::chrono;
 
 const int NUM_TYPES = 18;
 const int MIN_SIZE = 2;
 const int MAX_SIZE = 4;
 const int PRIME = 7;
 
 // ============================================================================
 // GLOBAL DATA STRUCTURES (matching original code structure)
 // ============================================================================
 
 vector<vector<int>> simplices;      // All simplices
 vector<set<int>> nei;               // Adjacency list
 set<int> tri;                       // Triangle hashes
 set<int> quad;                      // Tetrahedron hashes
 int maxVertex = 0;
 
 // ============================================================================
 // SIMPLET INDEX MAP (EXACTLY FROM ORIGINAL CODE)
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
         return -1;  // not found
     }
 }
 
 // ============================================================================
 // HASH FUNCTIONS (EXACTLY FROM ORIGINAL CODE)
 // ============================================================================
 
 // Function to compute a hash for a sorted set of vertices
 int computeHash(const vector<int>& vertices) {
     int hashValue = 0;
     int power = 1;
     for (int v : vertices) {
         hashValue += v * power;
         power *= PRIME;
     }
     return hashValue;
 }
 
 // Updated function to hash the simplex considering both degree counts and 3-simplices (triangles)
 // NOTE: THIS FUNCTION HAS A BUG - the edge check always returns true
 int hashSimplex(const set<int>& simplex) {
     vector<int> nodes(simplex.begin(), simplex.end());
     int n = nodes.size();
     vector<int> degrees(n, 0);
 
     // Count degrees within the simplex
     // BUG: This check always returns true because both nodes are in the simplex by definition
     for (int i = 0; i < n; ++i) {
         for (int j = i + 1; j < n; ++j) {
             // If there's an edge between two vertices, increment their degrees
             if (simplex.count(nodes[i]) && simplex.count(nodes[j])) {
                 degrees[i]++;
                 degrees[j]++;
             }
         }
     }
 
     // Sort degrees to make the order irrelevant
     sort(degrees.begin(), degrees.end());
 
     // Multiply each degree by a prime number
     const vector<int> primes = {2, 3, 5, 7}; // For up to 4 nodes
     int hash = 0;
     for (int i = 0; i < n; ++i) {
         hash += degrees[i] * primes[i];
     }
 
     // Count 3-simplices (triangles) and multiply by a prime number
     int triangle_count = 0;
     if (n >= 3) {
         for (int i = 0; i < n; ++i) {
             for (int j = i + 1; j < n; ++j) {
                 for (int k = j + 1; k < n; ++k) {
                     vector<int> triangle = {nodes[i], nodes[j], nodes[k]};
                     sort(triangle.begin(), triangle.end());
                     int tri_hash = computeHash(triangle);
                     if (tri.count(tri_hash)) {
                         triangle_count++; // Count the triangle
                     }
                 }
             }
         }
     }
 
     // Count 4-simplices (tetrahedra) and multiply by a different prime number
     int tetrahedron_count = 0;
     if (n == 4) {
         vector<int> tetrahedron = {nodes[0], nodes[1], nodes[2], nodes[3]};
         sort(tetrahedron.begin(), tetrahedron.end());
         int quad_hash = computeHash(tetrahedron);
         if (quad.count(quad_hash)) {
             tetrahedron_count++; // Count the tetrahedron
         }
     }
 
     // Multiply the triangle count by a prime number (e.g., 11) and tetrahedron count by another prime number (e.g., 13)
     hash += triangle_count * 11;  // Multiply by prime for triangles
     hash += tetrahedron_count * 13;  // Multiply by prime for tetrahedra
 
     // Add the size of the simplex as an offset to separate different simplex sizes
     hash += 1000 * n;
 
     return hash;
 }
 
 // ============================================================================
 // LOAD AND CONSTRUCT COMPLEX (matching original structure)
 // ============================================================================
 
 bool loadSimplices(const string& nverts_file, const string& simplices_file) {
     ifstream nv_in(nverts_file);
     ifstream sp_in(simplices_file);
     
     if (!nv_in.is_open() || !sp_in.is_open()) {
         cerr << "Error: Cannot open input files" << endl;
         return false;
     }
     
     // Read all nverts values
     vector<int> nverts_list;
     int nv;
     while (nv_in >> nv) {
         nverts_list.push_back(nv);
     }
     nv_in.close();
     
     cerr << "Read " << nverts_list.size() << " simplex sizes" << endl;
     
     // Debug: Show distribution
     map<int, int> size_dist;
     for (int n : nverts_list) {
         size_dist[n]++;
     }
     cerr << "Simplex size distribution:" << endl;
     for (auto& [sz, cnt] : size_dist) {
         cerr << "  Size " << sz << ": " << cnt << endl;
     }
     
     // Read ALL vertex indices sequentially
     vector<int> all_vertices;
     int v;
     while (sp_in >> v) {
         all_vertices.push_back(v);
         maxVertex = max(maxVertex, v);
     }
     sp_in.close();
     
     cerr << "Read " << all_vertices.size() << " vertex indices" << endl;
     
     // Parse simplices using nverts as boundaries
     size_t pos = 0;
     for (size_t i = 0; i < nverts_list.size() && pos < all_vertices.size(); i++) {
         int size = nverts_list[i];
         
         if (pos + size > all_vertices.size()) break;
         
         // Only keep simplices with size 2-4 (matching original)
         if (size >= 2 && size <= 4) {
             vector<int> simplex;
             for (int j = 0; j < size; j++) {
                 simplex.push_back(all_vertices[pos + j]);
             }
             sort(simplex.begin(), simplex.end());
             simplices.push_back(simplex);
         }
         pos += size;
     }
     
     cerr << "Loaded " << simplices.size() << " simplices (size 2-4)" << endl;
     
     return true;
 }
 
 void constructComplex() {
     // Resize adjacency list
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
     
     cerr << "Constructing complex..." << endl;
     cerr << "  Vertices: " << nei.size() << endl;
     cerr << "  Simplices: " << simplices.size() << endl;
     cerr << "  Edges: " << num_edges << endl;
     cerr << "  Triangles: " << tri.size() << endl;
     cerr << "  Tetrahedra: " << quad.size() << endl;
 }
 
 // ============================================================================
 // RANDOM WALK (EXACTLY FROM ORIGINAL CODE)
 // ============================================================================
 
 int getRandomIndex(int maxIndex, mt19937& gen) {
     uniform_int_distribution<int> dist(0, maxIndex - 1);
     return dist(gen);
 }
 
 bool isConnectedFast(const vector<int>& simplex) {
     int sz = simplex.size();
     if (sz < 2) return false;
     if (sz == 2) {
         // Two vertices: check if there's an edge
         return nei[simplex[0]].count(simplex[1]);
     }
 
     // Build subgraph adjacency between vertices in simplex
     int count = 0;
     for (int i = 0; i < sz; ++i) {
         for (int j = i + 1; j < sz; ++j) {
             if (nei[simplex[i]].count(simplex[j]))
                 count++;
         }
     }
 
     // For k nodes, minimum k - 1 edges needed to be connected
     return count >= sz - 1;
 }
 
 // Generate all neighboring simplices while ensuring connectivity
 vector<vector<int>> getNeighbors(const set<int>& currentSimplex) {
     vector<vector<int>> neighbors;
 
     // Convert set to sorted vector
     vector<int> simplexVec(currentSimplex.begin(), currentSimplex.end());
 
     // 1️⃣ **Remove one vertex (if size > 2)**
     if (currentSimplex.size() > MIN_SIZE) {
         for (int v : simplexVec) {
             vector<int> newSimplex = simplexVec;
             newSimplex.erase(remove(newSimplex.begin(), newSimplex.end(), v), newSimplex.end());
             if (isConnectedFast(newSimplex)) {
                 neighbors.push_back(newSimplex);
             }
         }
     }
 
     set<int> candidateVertices;
 
     // 2️⃣ **Add one vertex (if size < 4)**
     if (currentSimplex.size() < MAX_SIZE) {
         // create candidate vertices
         for (int v : simplexVec) {
             for (int neighbor : nei[v]) {
                 if (currentSimplex.count(neighbor) == 0) { // Ensure uniqueness
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
 
     // 3️⃣ **Swap one vertex (remove one, add another)**
     if (currentSimplex.size() > MIN_SIZE && currentSimplex.size() < MAX_SIZE) {
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
 
 set<int> randomWalk(int steps, mt19937& gen) {
     set<int> currentSimplex;
 
     if (simplices.empty()) {
         cerr << "Error: No simplices available!" << endl;
         return currentSimplex;
     }
 
     // Pick a random starting simplex
     int currentIndex;
     do {
         currentIndex = getRandomIndex(simplices.size(), gen);
     } while (simplices[currentIndex].size() <= 1);
     currentSimplex = set<int>(simplices[currentIndex].begin(), simplices[currentIndex].end());
 
     // Create a random distribution for acceptance probability
     uniform_real_distribution<> dis(0.0, 1.0);
 
     vector<vector<int>> neighbors = getNeighbors(currentSimplex);
     for (int step = 0; step < steps; step++) {
 
         if (neighbors.empty()) {
             break;
         }
 
         int random = getRandomIndex(neighbors.size(), gen);
         set<int> nextSimplex = set<int>(neighbors[random].begin(), neighbors[random].end());
         vector<vector<int>> nextNeighbors = getNeighbors(nextSimplex);
 
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
 // RESULT STRUCTURE
 // ============================================================================
 
 struct Result {
     string method;
     double time_ms;
     long long total;
     long long samples;
     long long counts[NUM_TYPES] = {0};
     double freq[NUM_TYPES] = {0};
     double ci_lo[NUM_TYPES] = {0};
     double ci_hi[NUM_TYPES] = {0};
     
     bool save(const string& filename) {
         ofstream out(filename);
         if (!out.is_open()) return false;
         
         out << "METHOD: " << method << "\n";
         out << "TIME_MS: " << time_ms << "\n";
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
     r.method = "MCMC_ORIGINAL";
     
     auto t0 = high_resolution_clock::now();
     
     mt19937 gen(chrono::steady_clock::now().time_since_epoch().count());
     
     const int steps = 20;  // Fixed steps as in original
     
     vector<int> sfd(NUM_TYPES + 1, 0);
     
     // Burn-in phase (not in original, but good practice)
     for (size_t i = 0; i < burn_in; i++) {
         randomWalk(steps, gen);
     }
     
     // Debug: Track unmatched hashes
     map<int, int> hash_counts;
     int unmatched = 0;
     
     // Sampling phase
     for (size_t i = 0; i < num_samples; i++) {
         set<int> simplex = randomWalk(steps, gen);
         int hash = hashSimplex(simplex);
         hash_counts[hash]++;
         
         int idx = findSimpletIndex(hash);
         if (idx >= 1 && idx <= NUM_TYPES) {
             sfd[idx]++;
             r.counts[idx - 1]++;
         } else {
             unmatched++;
         }
     }
     
     // Debug output
     cerr << "DEBUG: Hash statistics:" << endl;
     cerr << "  Total samples: " << num_samples << endl;
     cerr << "  Matched hashes: " << (num_samples - unmatched) << endl;
     cerr << "  Unmatched hashes: " << unmatched << " (" 
          << fixed << setprecision(2) << (100.0 * unmatched / num_samples) << "%)" << endl;
     
     if (unmatched > 0) {
         cerr << "  Top 10 unmatched hash values:" << endl;
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
         cerr << "Usage: " << argv[0] << " <nverts.txt> <simplices.txt> <output.res> [samples] [burnin]" << endl;
         cerr << endl;
         cerr << "MCMC_ORIGINAL - Original MCMC from Beigy et al. (EuroCG'24)" << endl;
         cerr << "=======================================================" << endl;
         cerr << "This implementation matches the original paper's code exactly," << endl;
         cerr << "including all known bugs, for fair comparison purposes." << endl;
         return 1;
     }
     
     string nverts_file = argv[1];
     string simplices_file = argv[2];
     string output_file = argv[3];
     size_t num_samples = (argc > 4) ? atoi(argv[4]) : 10000;
     size_t burn_in = (argc > 5) ? atoi(argv[5]) : 1000;
     
     // Load dataset
     cerr << "Loading dataset..." << endl;
     if (!loadSimplices(nverts_file, simplices_file)) {
         cerr << "Error: Cannot load dataset" << endl;
         return 1;
     }
     
     // Construct complex
     constructComplex();
     
     // Run MCMC
     cerr << "Running MCMC (" << num_samples << " samples, " << burn_in << " burn-in)..." << endl;
     Result result = approximateSFD(num_samples, burn_in);
     
     // Save result
     result.save(output_file);
     cerr << "Result saved to " << output_file << endl;
     cerr << "Time: " << fixed << setprecision(3) << result.time_ms << " ms" << endl;
     
     return 0;
 }