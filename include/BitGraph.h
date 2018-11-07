/*
*  Main authors:
*     Stefano Gualandi <stefano.gualandi@gmail.com>
*
*  Copyright:
*     Stefano Gualandi, 2017
*
*  Last update: November, 2018
*/
#pragma once

#include "VectorBitSet.h"

#include "Utils.h"

// Basic graph data structure
class BitGraph {
 public:
   // std c'tor
   BitGraph() : n(0) {}

   BitGraph(const BitGraph& o)
      : n(o.n), matrix(o.matrix), weights(o.weights), degrees(o.degrees),  anti_neigh(o.anti_neigh)
   {}

   // Init the internal memory
   void init(int _n) {
      if ( _n > std::min<int>(N_MAX, std::numeric_limits<vertex_t>::max()))
         throw std::invalid_argument("This datastructure support at most 512 vertices");

      n = _n;
      matrix.resize(n * n, false);
      for (int i = 0; i < n; ++i )
         matrix[i * n + i] = true;

      weights.resize(n, 0.0);
      degrees.resize(n, 0);

      anti_neigh.resize(n);
   }

   // Add the edge {i,j}
   void add_edge(int i, int j) {
      matrix[i * n + j] = true;
      matrix[j * n + i] = true;
      degrees[i]++;
      degrees[j]++;

      // Set the bit set for anti-neigh
      anti_neigh[i].reset(j);
      anti_neigh[j].reset(i);
   }

   // Check if {i,j} is an edge
   bool is_edge(int i, int j) const {
      return matrix[i * n + j];
   }

   // Number of vertices
   int num_vertices() const {
      return n;
   }

   // Get node weight
   double node_weight(vertex_t v) const {
      return weights[v];
   }

   // Get node weight
   double node_degree(vertex_t v) const {
      return degrees[v];
   }

   // Get candidates vertices for max stable set (reduce bit set size)
   inline bitset_t get_candidates() const {
      bitset_t Cs;
      Cs.set();
      for (size_t i = n; i < N_MAX; i++)
         Cs.reset(i);
      return Cs;
   }

   // Get bit set for anti neighbour
   inline bitset_t anti_n(size_t i) const {
      return anti_neigh[i];
   }

   // Set a node weight
   void set_node_weight(vertex_t v, double weight)  {
      weights[v] = weight;
   }

   // Print adjacent list
   void print() {
      for (int i = 0; i < n; ++i) {
         fprintf(stdout, "%d#%.3f: ", i, weights[i]);
         for (int j = 0; j < n; ++j)
            if (is_edge(i, j))
               fprintf(stdout, "%d ", j);
         fprintf(stdout, "\n");
      }
   }

 private:
   // Number of vertices
   int n;
   // Matrix in a vector format (compact in memory)
   vector<bool>    matrix;
   // Vector of weights
   vector<double>  weights;
   // Vector of vertex degree
   vertex_set_t    degrees;
   // Anti neihgbors
   set_bitset_t    anti_neigh;
};

size_t select_best_candidate(const BitGraph* G, const vertex_set_t& Vs, int n_vs);

// Read a dimacs instance
void read_graph(const std::string& filename, BitGraph& G);