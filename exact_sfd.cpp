

 #include "common.hpp"

 class ExactSFD {
     const SimplicialComplex& K;
     
 public:
     ExactSFD(const SimplicialComplex& K) : K(K) {}
     
     Result compute() {
         Result r;
         r.method = "EXACT";
         
         auto t0 = high_resolution_clock::now();
         
         
         r.counts[0] = K.edges.size();
         
         
         for (size_t v = 0; v < K.n; v++) {
             vector<int> neighbors(K.adj[v].begin(), K.adj[v].end());
             
             for (size_t i = 0; i < neighbors.size(); i++) {
                 for (size_t j = i + 1; j < neighbors.size(); j++) {
                     int u = neighbors[i];
                     int w = neighbors[j];
                     
                     bool has_uw = K.hasEdge(u, w);
                     
                     if (!has_uw) {
                         r.counts[1]++;  
                     } else {
                         if ((int)v < u && (int)v < w) {
                             if (K.hasTriangle(v, u, w)) {
                                 r.counts[3]++;  
                             } else {
                                 r.counts[2]++;  
                             }
                         }
                     }
                 }
             }
         }
         
         
         auto hashSimplet = [](vector<int> s) -> uint64_t {
             sort(s.begin(), s.end());
             return ((uint64_t)s[0] << 48) | ((uint64_t)s[1] << 32) | 
                    ((uint64_t)s[2] << 16) | s[3];
         };
         
         set<uint64_t> seen4;
         
         for (size_t v = 0; v < K.n; v++) {
             vector<int> neighbors(K.adj[v].begin(), K.adj[v].end());
             
             
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
             
             
             for (int u : neighbors) {
                 for (int w : K.adj[u]) {
                     if (w != (int)v && K.adj[v].count(w) == 0) {
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
         
         
         r.total = 0;
         for (int i = 0; i < NUM_TYPES; i++) {
             r.total += r.counts[i];
         }
         
         for (int i = 0; i < NUM_TYPES; i++) {
             r.freq[i] = r.total > 0 ? (double)r.counts[i] / r.total : 0;
             r.ci_lo[i] = r.freq[i];
             r.ci_hi[i] = r.freq[i];
         }
         
         r.samples = r.total;
         
         r.time_ms = duration_cast<microseconds>(
             high_resolution_clock::now() - t0).count() / 1000.0;
         
         return r;
     }
 };
 
 int main(int argc, char* argv[]) {
     if (argc < 4) {
         cerr << "Usage: " << argv[0] << " <nverts_file> <simplices_file> <output_file>" << endl;
         return 1;
     }
     
     SimplicialComplex K;
     
     cerr << "Loading dataset..." << endl;
     if (!loadDBLP(K, argv[1], argv[2])) {
         return 1;
     }
     
     cerr << "Running EXACT enumeration..." << endl;
     
     ExactSFD exact(K);
     Result r = exact.compute();
     
     r.save(argv[3]);
     
     cerr << "Time: " << fixed << setprecision(3) << r.time_ms << " ms" << endl;
     cerr << "Total simplets: " << r.total << endl;
     
     for (int i = 0; i < NUM_TYPES; i++) {
         if (r.counts[i] > 0) {
             cerr << "  Type " << setw(2) << (i+1) << ": " 
                  << setw(8) << r.counts[i] << " ("
                  << fixed << setprecision(4) << (r.freq[i] * 100) << "%)" << endl;
         }
     }
     
     cerr << "Saved to " << argv[3] << endl;
     return 0;
 }