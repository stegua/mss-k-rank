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

#include <bitset>
#include <array>
#include <vector>

//#include "intrin.h"

#define POPCNT __builtin_popcountll

const int N_MAX = 512;

class VectorBitSet {
 public:
   // Standard C'tor
   VectorBitSet(size_t n = 512) {
      _n = n;
      _words = n / 64 + 1;
      set();
      data.shrink_to_fit();
   }

   // Standard C'tor
   void set() {
      data.assign(_words, (uint64_t)~0);
   }

   // Standard C'tor
   void reset() {
      data.assign(_words, (uint64_t)0);
   }

   // Set a given bit to 1
   void set(size_t pos) noexcept {
      data[pos / 64] |= (uint64_t)1 << pos % 64;
   }

   // Set a given bit to 0
   void reset(size_t pos) noexcept {
      data[pos / 64] &= ~((uint64_t)1 << pos % 64);
   }

   // Test a given bit
   bool test(size_t pos) const noexcept {
      return ((data[pos / 64] & ((uint64_t)1 << pos % 64)) != 0);
   }

   // Size of bitset
   size_t size() const {
      return _n;
   }

   // Intersect three other sets
   void triple_intersect(const VectorBitSet& fst,
                         const VectorBitSet& snd,
                         const VectorBitSet& Cs) {
      for (size_t pos = 0; pos < _words; ++pos)
         data[pos] = (fst.get_word(pos) & (snd.get_word(pos) & Cs.get_word(pos)));
   }

   void intersect(const VectorBitSet& a, const VectorBitSet& b) {
      for (size_t pos = 0; pos < _words; ++pos)
         data[pos] = a.get_word(pos) & b.get_word(pos);
   }

   uint64_t get_word(size_t pos) const {
      return data[pos];
   }

   int count() const {
      uint64_t c = 0;
      for (size_t pos = 0; pos < _words; ++pos)
         c += POPCNT(data[pos]);
      return int(c);
   }

   bool more_than_k(const VectorBitSet& b, uint64_t k) const {
      uint64_t c = 0;
      for (size_t pos = 0; pos < _words; ++pos) {
         c += POPCNT(data[pos] & b.get_word(pos));
         if (c >= k)
            return true;
      }
      return false;
   }

 private:
   size_t _n;
   size_t _words;
   std::vector<uint64_t>  data;
};

typedef VectorBitSet            bitset_t;
typedef std::vector<bitset_t>   set_bitset_t;


// Test intersection
// Ts: set of "k" complete subgraph in the complement graph
//
// Check if cardinality(r=(fst intersect snd) intersect Cs) < k  ====> no problem
//       if any T in Ts has cardinality(r intersect T) == k      ====> problem
//       else no problem
bool has_tuple(const VectorBitSet& fst,
               const VectorBitSet& snd,
               const VectorBitSet& Cs,
               const set_bitset_t& Ts,
               int k);
