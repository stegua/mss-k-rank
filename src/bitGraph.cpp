#include "BitGraph.h"

// Read a dimacs instance
#include <fstream>
void read_graph(const std::string& filename, BitGraph& G) {
   std::ifstream infile(filename);

   // p edges 60 89
   int n = 0;
   int m = 0;
   char buffer[1024];
   infile >> buffer >> buffer >> n >> m;

   G.init(n);

   int i = 0;
   int j = 0;
   while (infile >> buffer >> i >> j) {
      G.add_edge(i - 1, j - 1);
   }
}

// Select best node to explore
size_t select_best_candidate(const BitGraph* G, const vertex_set_t& Vs, int n_vs) {
   size_t i_best = n_vs - 1;
   for (size_t i = 0, i_max = n_vs - 1; i < i_max; ++i)
      if (G->node_weight(Vs[i]) > G->node_weight(Vs[i_best]))
         i_best = i;
      else {
         if (fabs(G->node_weight(Vs[i]) - G->node_weight(Vs[i_best])) < INT_TOL
               && G->node_degree(Vs[i]) > G->node_degree(Vs[i_best]))
            i_best = i;
      }
   return i_best;
}