/**
 * ============================================================================
 * COMMON.HPP - Shared Classes and Utilities for SFD Estimation
 * ============================================================================
 * 
 * This file contains all shared data structures and utilities used by:
 * - exact_sfd.cpp: Exact enumeration (ground truth)
 * - mcmc_sfd.cpp: Original MCMC from paper (with hash function bug)
 * - unified_sfd.cpp: Our improved MCMC (Robust Auto-Adaptive)
 * - compare.cpp: Comparison tool
 * 
 * Author: TuanAnh (Master's Thesis)
 * Based on: "Approximating Simplet Frequency Distribution for Simplicial Complexes"
 *           (Beigy et al., EuroCG'24)
 * ============================================================================
 */

 #ifndef COMMON_HPP
 #define COMMON_HPP
 
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
 #include <numeric>  // For std::accumulate
 
 using namespace std;
 using namespace std::chrono;
 
 // ============================================================================
 // CONSTANTS
 // ============================================================================
 
 const int NUM_TYPES = 18;  // 18 simplet types (1-18)
 
 // ============================================================================
 // SIMPLICIAL COMPLEX DATA STRUCTURE
 // ============================================================================
 
 /**
  * SimplicialComplex represents a simplicial complex with:
  * - Vertices (0 to n-1)
  * - Edges (2-simplices)
  * - Triangles (3-simplices, "filled" triangles)
  * - Tetrahedra (4-simplices)
  */
 struct SimplicialComplex {
     size_t n = 0;                           // Number of vertices
     vector<set<int>> adj;                   // Adjacency list (edges)
     set<pair<int,int>> edges;               // All edges as pairs
     set<tuple<int,int,int>> tris;           // Filled triangles (3-simplices)
     set<tuple<int,int,int,int>> tets;       // Tetrahedra (4-simplices)
     
     void init(size_t num_vertices) {
         n = num_vertices;
         adj.resize(n);
     }
     
     void addEdge(int u, int v) {
         if (u > v) swap(u, v);
         if (u != v && u >= 0 && v >= 0 && (size_t)u < n && (size_t)v < n) {
             adj[u].insert(v);
             adj[v].insert(u);
             edges.insert({u, v});
         }
     }
     
     void addTriangle(int a, int b, int c) {
         vector<int> v = {a, b, c};
         sort(v.begin(), v.end());
         tris.insert({v[0], v[1], v[2]});
         // Also add edges
         addEdge(a, b);
         addEdge(b, c);
         addEdge(a, c);
     }
     
     void addTetrahedron(int a, int b, int c, int d) {
         vector<int> v = {a, b, c, d};
         sort(v.begin(), v.end());
         tets.insert({v[0], v[1], v[2], v[3]});
         // Also add triangles and edges
         addTriangle(a, b, c);
         addTriangle(a, b, d);
         addTriangle(a, c, d);
         addTriangle(b, c, d);
     }
     
     bool hasEdge(int u, int v) const {
         if (u > v) swap(u, v);
         return edges.count({u, v}) > 0;
     }
     
     bool hasTriangle(int a, int b, int c) const {
         vector<int> v = {a, b, c};
         sort(v.begin(), v.end());
         return tris.count({v[0], v[1], v[2]}) > 0;
     }
     
     bool hasTetrahedron(int a, int b, int c, int d) const {
         vector<int> v = {a, b, c, d};
         sort(v.begin(), v.end());
         return tets.count({v[0], v[1], v[2], v[3]}) > 0;
     }
 };
 
 // ============================================================================
 // RESULT STRUCTURE
 // ============================================================================
 
 /**
  * Result stores the output of an SFD estimation method.
  */
 struct Result {
     string method;
     double time_ms = 0;
     size_t total = 0;      // Total simplets (population or estimate)
     size_t samples = 0;    // Number of samples taken
     
     size_t counts[NUM_TYPES] = {0};    // Count of each type
     double freq[NUM_TYPES] = {0};      // Frequency of each type
     double ci_lo[NUM_TYPES] = {0};     // 95% CI lower bound
     double ci_hi[NUM_TYPES] = {0};     // 95% CI upper bound
     
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
     
     bool load(const string& filename) {
         ifstream in(filename);
         if (!in.is_open()) return false;
         
         string line;
         while (getline(in, line)) {
             if (line.substr(0, 7) == "METHOD:") {
                 method = line.substr(8);
             } else if (line.substr(0, 8) == "TIME_MS:") {
                 time_ms = stod(line.substr(9));
             } else if (line.substr(0, 6) == "TOTAL:") {
                 total = stoull(line.substr(7));
             } else if (line.substr(0, 8) == "SAMPLES:") {
                 samples = stoull(line.substr(9));
             } else if (line[0] >= '0' && line[0] <= '9') {
                 // Parse data line: TYPE,COUNT,FREQ,CI_LO,CI_HI
                 stringstream ss(line);
                 string token;
                 vector<string> tokens;
                 while (getline(ss, token, ',')) {
                     tokens.push_back(token);
                 }
                 if (tokens.size() >= 3) {
                     int type = stoi(tokens[0]) - 1;
                     if (type >= 0 && type < NUM_TYPES) {
                         counts[type] = stoull(tokens[1]);
                         freq[type] = stod(tokens[2]);
                         if (tokens.size() >= 5) {
                             ci_lo[type] = stod(tokens[3]);
                             ci_hi[type] = stod(tokens[4]);
                         }
                     }
                 }
             }
         }
         in.close();
         return true;
     }
 };
 
 // ============================================================================
 // CORRECTED SIMPLET CLASSIFIER
 // ============================================================================
 
 /**
  * Classifier provides CORRECT simplet classification.
  * 
  * This fixes the bug in the original paper where the hash function
  * always returned true for edge existence checks.
  * 
  * Simplet Types:
  * - Type 1: Edge (2 vertices)
  * - Types 2-4: 3-simplets (path, empty triangle, filled triangle)
  * - Types 5-18: 4-simplets (14 different topological structures)
  */
 class Classifier {
 public:
     /**
      * Check if a set of vertices forms a connected subgraph.
      */
     static bool isConnected(const vector<int>& vertices, const SimplicialComplex& K) {
         if (vertices.size() <= 1) return true;
         
         set<int> visited;
         vector<int> queue;
         queue.push_back(vertices[0]);
         visited.insert(vertices[0]);
         
         set<int> vertex_set(vertices.begin(), vertices.end());
         
         while (!queue.empty()) {
             int v = queue.back();
             queue.pop_back();
             
             for (int u : vertices) {
                 if (visited.count(u) == 0 && K.hasEdge(v, u)) {
                     visited.insert(u);
                     queue.push_back(u);
                 }
             }
         }
         
         return visited.size() == vertices.size();
     }
     
     /**
      * Classify a simplet (2, 3, or 4 vertices) into one of 18 types.
      */
     static int classify(vector<int> vertices, const SimplicialComplex& K) {
         sort(vertices.begin(), vertices.end());
         
         if (vertices.size() == 2) {
             // Type 1: Edge
             return 1;
         }
         
         if (vertices.size() == 3) {
             return classify3(vertices, K);
         }
         
         if (vertices.size() == 4) {
             return classify4(vertices, K);
         }
         
         return 0;  // Invalid
     }
     
 private:
     /**
      * Classify 3-simplet into types 2, 3, or 4.
      */
     static int classify3(const vector<int>& v, const SimplicialComplex& K) {
         int a = v[0], b = v[1], c = v[2];
         
         // Count edges
         int edge_count = 0;
         if (K.hasEdge(a, b)) edge_count++;
         if (K.hasEdge(b, c)) edge_count++;
         if (K.hasEdge(a, c)) edge_count++;
         
         if (edge_count == 2) {
             // Type 2: Path (2 edges, no triangle possible)
             return 2;
         }
         
         if (edge_count == 3) {
             // 3 edges = triangle structure
             if (K.hasTriangle(a, b, c)) {
                 // Type 4: Filled triangle (has 2-simplex)
                 return 4;
             } else {
                 // Type 3: Empty triangle (edges only, no 2-simplex)
                 return 3;
             }
         }
         
         return 0;  // Not connected or invalid
     }
     
     /**
      * Classify 4-simplet into types 5-18.
      * 
      * Classification is based on:
      * - Number of edges (3-6)
      * - Number of triangles (0-4)
      * - Number of tetrahedra (0-1)
      * - Structural properties (cycles, degrees)
      */
     static int classify4(const vector<int>& v, const SimplicialComplex& K) {
         int a = v[0], b = v[1], c = v[2], d = v[3];
         
         // Count edges
         bool e_ab = K.hasEdge(a, b);
         bool e_ac = K.hasEdge(a, c);
         bool e_ad = K.hasEdge(a, d);
         bool e_bc = K.hasEdge(b, c);
         bool e_bd = K.hasEdge(b, d);
         bool e_cd = K.hasEdge(c, d);
         
         int num_edges = e_ab + e_ac + e_ad + e_bc + e_bd + e_cd;
         
         // Count filled triangles
         bool t_abc = K.hasTriangle(a, b, c);
         bool t_abd = K.hasTriangle(a, b, d);
         bool t_acd = K.hasTriangle(a, c, d);
         bool t_bcd = K.hasTriangle(b, c, d);
         
         int num_tris = t_abc + t_abd + t_acd + t_bcd;
         
         // Check tetrahedron
         bool has_tet = K.hasTetrahedron(a, b, c, d);
         
         // Compute degree sequence
         vector<int> degrees(4);
         degrees[0] = e_ab + e_ac + e_ad;
         degrees[1] = e_ab + e_bc + e_bd;
         degrees[2] = e_ac + e_bc + e_cd;
         degrees[3] = e_ad + e_bd + e_cd;
         sort(degrees.begin(), degrees.end());
         
         // Classify based on (num_edges, num_tris, has_tet, degree_sequence)
         
         // 3 edges (tree structures)
         if (num_edges == 3) {
             if (degrees == vector<int>{1, 1, 1, 3}) {
                 return 5;  // Star (one center connected to all)
             }
             if (degrees == vector<int>{1, 1, 2, 2}) {
                 return 6;  // Path
             }
         }
         
         // 4 edges
         if (num_edges == 4) {
             if (degrees == vector<int>{1, 2, 2, 3}) {
                 return 7;  // Triangle + pendant
             }
             if (degrees == vector<int>{2, 2, 2, 2}) {
                 return 8;  // 4-cycle
             }
             if (degrees == vector<int>{1, 1, 3, 3}) {
                 return 9;  // Two triangles sharing edge? Check structure
             }
         }
         
         // 5 edges
         if (num_edges == 5) {
             if (num_tris == 0) {
                 return 10;  // No triangles, 5 edges
             }
             if (num_tris == 1) {
                 return 11;  // One filled triangle
             }
             if (num_tris == 2) {
                 return 12;  // Two filled triangles
             }
         }
         
         // 6 edges (complete graph K4)
         if (num_edges == 6) {
             if (num_tris == 0) {
                 return 13;  // K4 with no filled triangles
             }
             if (num_tris == 1) {
                 return 14;  // K4 with 1 filled triangle
             }
             if (num_tris == 2) {
                 return 15;  // K4 with 2 filled triangles
             }
             if (num_tris == 3) {
                 return 16;  // K4 with 3 filled triangles
             }
             if (num_tris == 4) {
                 if (has_tet) {
                     return 18;  // Complete tetrahedron
                 } else {
                     return 17;  // K4 with 4 triangles but no tetrahedron
                 }
             }
         }
         
         // Fallback: use edge count for partial classification
         if (num_edges == 4) {
             // Distinguish between structures with 4 edges
             // Check for specific patterns
             if (num_tris >= 1) return 7;
             
             // Check if it's a cycle
             bool is_cycle = (degrees == vector<int>{2, 2, 2, 2});
             if (is_cycle) return 8;
             
             return 9;
         }
         
         return 5;  // Default for tree-like structures
     }
 };
 
 // ============================================================================
 // DBLP FORMAT LOADER
 // ============================================================================
 
 /**
  * Load a simplicial complex from DBLP format files.
  * 
  * DBLP Format:
  * - nverts.txt: Each line contains ONE integer = number of vertices in that simplex
  * - simplices.txt: Each line contains ONE vertex ID (all vertices listed sequentially)
  * 
  * Example:
  *   nverts.txt:     simplices.txt:
  *   1               1
  *   2               1
  *   3               2
  *                   3
  *                   4
  *                   5
  *   
  *   Meaning:
  *   - Simplex 1: vertex 1 (isolated vertex, size=1)
  *   - Simplex 2: vertices 1,2 (edge, size=2)
  *   - Simplex 3: vertices 3,4,5 (triangle, size=3)
  */
 bool loadDBLP(SimplicialComplex& K, const string& nverts_file, const string& simplices_file) {
     ifstream nv_in(nverts_file);
     ifstream sp_in(simplices_file);
     
     if (!nv_in.is_open()) {
         cerr << "Error: Cannot open " << nverts_file << endl;
         return false;
     }
     if (!sp_in.is_open()) {
         cerr << "Error: Cannot open " << simplices_file << endl;
         return false;
     }
     
     // Step 1: Read all nverts values
     vector<int> nverts_list;
     int nv;
     while (nv_in >> nv) {
         nverts_list.push_back(nv);
     }
     nv_in.close();
     
     cerr << "Read " << nverts_list.size() << " simplex sizes from nverts file" << endl;
     
     // Debug: Show distribution of simplex sizes
     map<int, int> size_dist;
     for (int n : nverts_list) {
         size_dist[n]++;
     }
     cerr << "Simplex size distribution:" << endl;
     for (auto& [sz, cnt] : size_dist) {
         cerr << "  Size " << sz << ": " << cnt << " simplices" << endl;
     }
     
     // Step 2: Read ALL vertex indices from simplices file
     vector<int> all_vertices;
     int v;
     while (sp_in >> v) {
         all_vertices.push_back(v);
     }
     sp_in.close();
     
     cerr << "Read " << all_vertices.size() << " vertex indices from simplices file" << endl;
     
     // Verify total count
     size_t expected_total = 0;
     for (int n : nverts_list) {
         expected_total += n;
     }
     
     cerr << "Expected total vertices: " << expected_total << endl;
     
     if (expected_total != all_vertices.size()) {
         cerr << "WARNING: Mismatch! Expected " << expected_total 
              << " but got " << all_vertices.size() << endl;
     }
     
     // Step 3: Find max vertex ID
     int max_v = 0;
     for (int v : all_vertices) {
         max_v = max(max_v, v);
     }
     
     cerr << "Max vertex ID: " << max_v << endl;
     
     // Step 4: Initialize complex
     K.init(max_v + 1);
     
     // Step 5: Parse simplices using nverts as boundaries
     size_t pos = 0;
     size_t num_isolated = 0, num_edges = 0, num_tris = 0, num_tets = 0, num_higher = 0;
     
     for (size_t i = 0; i < nverts_list.size() && pos < all_vertices.size(); i++) {
         int size = nverts_list[i];
         
         if (pos + size > all_vertices.size()) {
             cerr << "Warning: Not enough vertices for simplex " << i 
                  << " (need " << size << ", have " << (all_vertices.size() - pos) << ")" << endl;
             break;
         }
         
         // Extract vertices for this simplex
         vector<int> vertices;
         for (int j = 0; j < size; j++) {
             vertices.push_back(all_vertices[pos + j]);
         }
         pos += size;
         
         // Add simplex based on size
         if (size == 1) {
             // Isolated vertex - skip (no edges)
             num_isolated++;
         }
         else if (size == 2) {
             K.addEdge(vertices[0], vertices[1]);
             num_edges++;
         } 
         else if (size == 3) {
             K.addTriangle(vertices[0], vertices[1], vertices[2]);
             num_tris++;
         }
         else if (size == 4) {
             K.addTetrahedron(vertices[0], vertices[1], vertices[2], vertices[3]);
             num_tets++;
         }
         else if (size > 4) {
             // Higher-order simplex: add all sub-simplices
             // LIMIT: Only decompose if size <= 25 to prevent memory explosion
             // For size=25: C(25,4) = 12,650 tetrahedra (manageable)
             // For size=280: C(280,4) = 251 million tetrahedra (OOM!)
             
             if (size > 25) {
                 // Too large - only add direct edges (co-authorship)
                 // This is a reasonable interpretation: vertices in same simplex are connected
                 for (size_t a = 0; a < vertices.size(); a++) {
                     for (size_t b = a + 1; b < vertices.size(); b++) {
                         K.addEdge(vertices[a], vertices[b]);
                     }
                 }
                 num_higher++;
             } else {
                 num_higher++;
                 
                 // Add all edges
                 for (size_t a = 0; a < vertices.size(); a++) {
                     for (size_t b = a + 1; b < vertices.size(); b++) {
                         K.addEdge(vertices[a], vertices[b]);
                     }
                 }
                 // Add all triangles
                 for (size_t a = 0; a < vertices.size(); a++) {
                     for (size_t b = a + 1; b < vertices.size(); b++) {
                         for (size_t c = b + 1; c < vertices.size(); c++) {
                             K.addTriangle(vertices[a], vertices[b], vertices[c]);
                         }
                     }
                 }
                 // Add all tetrahedra
                 for (size_t a = 0; a < vertices.size(); a++) {
                     for (size_t b = a + 1; b < vertices.size(); b++) {
                         for (size_t c = b + 1; c < vertices.size(); c++) {
                             for (size_t d = c + 1; d < vertices.size(); d++) {
                                 K.addTetrahedron(vertices[a], vertices[b], vertices[c], vertices[d]);
                             }
                         }
                     }
                 }
             }
         }
         
         // Progress report every 500K simplices
         if ((i + 1) % 500000 == 0) {
             cerr << "  Processed " << (i + 1) << "/" << nverts_list.size() << " simplices..." << endl;
         }
     }
     
     cerr << endl;
     cerr << "========================================" << endl;
     cerr << "SIMPLICIAL COMPLEX LOADED" << endl;
     cerr << "========================================" << endl;
     cerr << "Graph structure:" << endl;
     cerr << "  Vertices:    " << K.n << endl;
     cerr << "  Edges:       " << K.edges.size() << endl;
     cerr << "  Triangles:   " << K.tris.size() << endl;
     cerr << "  Tetrahedra:  " << K.tets.size() << endl;
     cerr << endl;
     cerr << "Input simplices: " << nverts_list.size() << endl;
     cerr << "  Isolated vertices (size=1): " << num_isolated << endl;
     cerr << "  2-simplices (edges):        " << num_edges << endl;
     cerr << "  3-simplices (triangles):    " << num_tris << endl;
     cerr << "  4-simplices (tetrahedra):   " << num_tets << endl;
     cerr << "  Higher-order (size>4):      " << num_higher << endl;
     cerr << "========================================" << endl;
     
     // Large dataset warning
     if (K.edges.size() > 1000000) {
         cerr << endl;
         cerr << "⚠️  LARGE DATASET DETECTED (" << K.edges.size() << " edges)" << endl;
         cerr << "    EXACT enumeration may be slow or run out of memory." << endl;
         cerr << "    Consider using only MCMC methods: make test-mcmc-only" << endl;
     }
     
     if (K.edges.size() == 0) {
         cerr << endl;
         cerr << "WARNING: No edges found! Check if:" << endl;
         cerr << "  1. File format is correct (each line = one number)" << endl;
         cerr << "  2. nverts.txt contains values >= 2 for edges" << endl;
         cerr << "  3. Files don't have Windows line endings (\\r\\n)" << endl;
     }
     
     return true;
 }
 
 #endif // COMMON_HPP