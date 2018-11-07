#pragma once

#include "Utils.h"

#include "VertexTuples.h"

/*
* @brief
*/
void find_best_3_rank(vertex_set_t& Vs,
                      int n_vs,
                      double vs_cost,
                      vertex_set_t& Rs,
                      double rs_cost,
                      vertex_set_t& Best,
                      double& best_cost,
                      const BitGraph* G,
                      VertexTuples& Es,
                      bitset_t& Cs,
                      double CUT_VIOL,
                      uint64_t& nodes,
                      tpoint start_time) {
   // Number of nodes of B&B
   nodes++;
   // Select heavier clique
   if (n_vs == 0) {
      if (rs_cost > best_cost) {
         best_cost = rs_cost;
         Best = Rs;
      } else {
         // Take for maximal set
         if (fabs(rs_cost - best_cost) < INT_TOL && Rs.size() > Best.size()) {
            best_cost = rs_cost;
            Best = Rs;
         }
      }
      return;
   }

   // Naive pruning
   if (rs_cost + vs_cost <= std::max(3.0 + CUT_VIOL, best_cost))
      return;

   if (best_cost > 3.0 + CUT_VIOL) {
      if (limit == std::numeric_limits<uint64_t>::max()) {
         nodes = 1;
         limit = 5000;
      }

      if (nodes > limit)
         return;
   }

   // TIMELIMIT
   if (std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - start_time).count() > TIMEOUT * 1000)
      return;

   size_t i_best = select_best_candidate(G, Vs, n_vs);
   std::swap(Vs[i_best], Vs[n_vs - 1]);

   vertex_t head = Vs[n_vs - 1];
   int n_ns = 1;
   double ns_cost = G->node_weight(head);

   int i = 0;
   while (i < n_vs - n_ns) {
      if (TestTuples(Vs[i], head, Cs, Es, G, 2)) {
         i++;
      } else {
         ns_cost += G->node_weight(Vs[i]);
         std::swap(Vs[n_vs - n_ns - 1], Vs[i]);
         ++n_ns;
      }
   }

   // LEFT BRANCH
   int added = 0;
   for (vertex_t v : Rs)
      if (!G->is_edge(head, v)) {
         Es.add_pair(head, v);
         added++;
      }

   Rs.push_back(head);
   Cs.set(head);
   find_best_3_rank(Vs, n_vs - n_ns, vs_cost - ns_cost,  // Candidates vertices
                    Rs, rs_cost + G->node_weight(head),   // Current partial solution
                    Best, best_cost, G, Es, Cs, CUT_VIOL, nodes, start_time);             // Best solution + Graph

   Rs.pop_back();
   Cs.reset(head);
   Es.remove_last(added);

   // RIGHT BRANCH
   find_best_3_rank(Vs, n_vs - 1, vs_cost - G->node_weight(head),
                    Rs, rs_cost,
                    Best, best_cost, G, Es, Cs, CUT_VIOL, nodes, start_time);
}

/*
* @brief
*/
std::pair<double, vertex_set_t>
best_3_rank(const BitGraph* G, double CUT_VIOL, tpoint start_time) {
   vertex_set_t Vs;
   for (vertex_t v = 0; v < G->num_vertices(); ++v)
      Vs.push_back(v);

   vertex_set_t Rs;
   Rs.reserve(G->num_vertices());
   bitset_t Cs = G->get_candidates();

   // Non-edge which leads to potential quadrpuple
   VertexTuples Es;
   Es.init(G->num_vertices()*G->num_vertices());

   vertex_set_t Best;
   Best.reserve(G->num_vertices());

   double vs_cost = 0.0;
   for (vertex_t v = 0; v < G->num_vertices(); ++v)
      vs_cost += G->node_weight(v);

   // Reset node limit
   limit = std::numeric_limits<uint64_t>::max();

   double best_cost = 3.0;

   // node limit
   uint64_t nodes = 0;

   // When it finshed Best is the best set, and best_cost its cost
   std::sort(Vs.begin(), Vs.end(), [&G](int i, int j) {
      return G->node_weight(i) / G->node_degree(i) < G->node_weight(j) / G->node_degree(j);
   });

   find_best_3_rank(Vs, Vs.size(), vs_cost, Rs, 0.0, Best, best_cost, G, Es, Cs, CUT_VIOL, nodes, start_time);

   return std::make_pair(best_cost, Best);
}