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
 #include <numeric>
 
 using namespace std;
 using namespace std::chrono;
 
 const int NUM_TYPES = 18;
 
 
 struct SimplicialComplex {
     size_t n = 0;
     vector<set<int>> adj;
     set<pair<int,int>> edges;
     set<tuple<int,int,int>> tris;
     set<tuple<int,int,int,int>> tets;
     
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
 
 
 
 class Classifier {
 public:
     static bool isConnected(const vector<int>& vertices, const SimplicialComplex& K) {
         if (vertices.size() <= 1) return true;
         
         set<int> visited;
         vector<int> queue;
         queue.push_back(vertices[0]);
         visited.insert(vertices[0]);
         
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
         
         if (edge_count == 2) return 2;
         if (edge_count == 3) {
             return K.hasTriangle(a, b, c) ? 4 : 3;
         }
         return 0;
     }
     
     static int classify4(const vector<int>& v, const SimplicialComplex& K) {
         int a = v[0], b = v[1], c = v[2], d = v[3];
         
         bool e_ab = K.hasEdge(a, b), e_ac = K.hasEdge(a, c), e_ad = K.hasEdge(a, d);
         bool e_bc = K.hasEdge(b, c), e_bd = K.hasEdge(b, d), e_cd = K.hasEdge(c, d);
         
         int num_edges = e_ab + e_ac + e_ad + e_bc + e_bd + e_cd;
         
         int num_tris = K.hasTriangle(a, b, c) + K.hasTriangle(a, b, d) + 
                        K.hasTriangle(a, c, d) + K.hasTriangle(b, c, d);
         
         bool has_tet = K.hasTetrahedron(a, b, c, d);
         
         vector<int> degrees = {
             (int)(e_ab + e_ac + e_ad),
             (int)(e_ab + e_bc + e_bd),
             (int)(e_ac + e_bc + e_cd),
             (int)(e_ad + e_bd + e_cd)
         };
         sort(degrees.begin(), degrees.end());
         
         if (num_edges == 3) {
             if (degrees == vector<int>{1, 1, 1, 3}) return 5;
             if (degrees == vector<int>{1, 1, 2, 2}) return 6;
         }
         
         if (num_edges == 4) {
            if (degrees == vector<int>{1, 2, 2, 3}) return num_tris == 1 ? 9 : 7;
            if (degrees == vector<int>{2, 2, 2, 2}) return 8;
        }
         
         if (num_edges == 5) {
             if (num_tris == 0) return 10;
             if (num_tris == 1) return 11;
             if (num_tris == 2) return 12;
         }
         
         if (num_edges == 6) {
             if (num_tris == 0) return 13;
             if (num_tris == 1) return 14;
             if (num_tris == 2) return 15;
             if (num_tris == 3) return 16;
             if (num_tris == 4) return has_tet ? 18 : 17;
         }
         
         return 0;
     }
 };
 
 
 
 bool loadDBLP(SimplicialComplex& K, const string& nverts_file, const string& simplices_file) {
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
     
     vector<int> all_vertices;
     int v;
     while (sp_in >> v) { all_vertices.push_back(v); }
     sp_in.close();
     
     int max_v = 0;
     for (int vv : all_vertices) { max_v = max(max_v, vv); }
     
     K.init(max_v + 1);
     
     size_t pos = 0;
     for (size_t i = 0; i < nverts_list.size() && pos < all_vertices.size(); i++) {
         int size = nverts_list[i];
         if (pos + size > all_vertices.size()) break;
         
         vector<int> vertices;
         for (int j = 0; j < size; j++) {
             vertices.push_back(all_vertices[pos + j]);
         }
         pos += size;
         
         if (size == 1) continue;
         else if (size == 2) K.addEdge(vertices[0], vertices[1]);
         else if (size == 3) K.addTriangle(vertices[0], vertices[1], vertices[2]);
         else if (size == 4) K.addTetrahedron(vertices[0], vertices[1], vertices[2], vertices[3]);
         else if (size <= 25) {
             for (size_t a = 0; a < vertices.size(); a++)
                 for (size_t b = a + 1; b < vertices.size(); b++)
                     K.addEdge(vertices[a], vertices[b]);
             for (size_t a = 0; a < vertices.size(); a++)
                 for (size_t b = a + 1; b < vertices.size(); b++)
                     for (size_t c = b + 1; c < vertices.size(); c++)
                         K.addTriangle(vertices[a], vertices[b], vertices[c]);
             for (size_t a = 0; a < vertices.size(); a++)
                 for (size_t b = a + 1; b < vertices.size(); b++)
                     for (size_t c = b + 1; c < vertices.size(); c++)
                         for (size_t d = c + 1; d < vertices.size(); d++)
                             K.addTetrahedron(vertices[a], vertices[b], vertices[c], vertices[d]);
         } else {
             for (size_t a = 0; a < vertices.size(); a++)
                 for (size_t b = a + 1; b < vertices.size(); b++)
                     K.addEdge(vertices[a], vertices[b]);
         }
     }
     
     cerr << "Loaded: " << K.edges.size() << " edges, " 
          << K.tris.size() << " triangles, " << K.tets.size() << " tetrahedra" << endl;
     
     return true;
 }
 
 #endif 