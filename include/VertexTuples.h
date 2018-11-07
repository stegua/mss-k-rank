#pragma once

#include "BitGraph.h"

class VertexTuples {
 public:

   void init(int _n) {
      data.reserve(_n);
   }

   // Setters
   void add_pair(vertex_t i, vertex_t j) {
      data.emplace_back();
      data.back().reset();  // NON MI PIACE QUESTO (!)
      data.back().set(i);
      data.back().set(j);
   }

   void add_triple(vertex_t i, vertex_t j, vertex_t h) {
      data.emplace_back();
      data.back().reset();  // NON MI PIACE QUESTO (!)
      data.back().set(i);
      data.back().set(j);
      data.back().set(h);
   }

   void add_quadruple(vertex_t i, vertex_t j, vertex_t h, vertex_t k) {
      data.emplace_back();
      data.back().reset();  // NON MI PIACE QUESTO (!)
      data.back().set(i);
      data.back().set(j);
      data.back().set(h);
      data.back().set(k);
   }

   void remove_last(size_t i) {
      data.erase(data.end() - i, data.end());
      data.shrink_to_fit();
   }

   void clear() {
      data.clear();
   }

   size_t size() const {
      return data.size();
   }

   set_bitset_t data;
};


// Test for potential complete anti-clique (clique on the complement)
bool TestTuples(vertex_t v, vertex_t h, const bitset_t& Cs, const VertexTuples& Es, const BitGraph* G, int k) {
   if (G->is_edge(h, v))
      return true;

   return has_tuple(G->anti_n(v), G->anti_n(h), Cs, Es.data, k);
}