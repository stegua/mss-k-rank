/*
*  Main authors:
*     Stefano Gualandi <stefano.gualandi@gmail.com>
*
*  Copyright:
*     Stefano Gualandi, 2017
*
*  Last update: November, 2020
*/
#pragma once

#include <numeric>
#include <algorithm>
#include <cassert>
#include <random>

#include <vector>
using std::vector;

#include <array>
using std::array;

#include <set>
using std::set;

#include "BitGraph.h"

// Typedef for graph types
typedef int      node_t;
typedef int      arc_t;
typedef double   cost_t;


// Largest possible value
const node_t      node_last = std::numeric_limits<node_t>::max();
const arc_t       arc_last = std::numeric_limits<arc_t>::max();
const cost_t      cost_max = std::numeric_limits<cost_t>::max();

// Support for labeling algorithm: Labels dijkstra algorithm
enum Label { UNREACHED, LABELED, SCANNED };


// Basic Heap implementation
#define UP(i) i >> 1

class BinaryHeap {

public:
	// Standard c'tor
	BinaryHeap(size_t n) {
		size = 1;
		max_size = n;
		data.resize(n + 1, std::pair<cost_t, node_t>(std::numeric_limits<cost_t>::max(), node_t(n)));
		pos.resize(n + 1, n);
		data[0].first = 0;
	}

	// Rule of five: Move constructor
	// Deleted constructor (this is a singleton object)
	BinaryHeap(BinaryHeap&& o) = delete;
	BinaryHeap(const BinaryHeap& o) = delete;
	BinaryHeap& operator=(const BinaryHeap&) = delete;
	BinaryHeap& operator=(BinaryHeap&&) = delete;

	// Check if empty
	bool empty() const {
		return size == 1;
	}

	// get distance of given node
	cost_t get_top_distance() const {
		return data[1].first;
	}

	// get distance of given node
	cost_t get_distance(node_t v) const {
		return data[pos[v]].first;
	}

	// insert a new item
	void insert(cost_t t, node_t i) {
		// Percolate up
		int64_t hole = size;
		++size;
		for (; hole > 1 && t < data[UP(hole)].first; hole = UP(hole)) {
			data[hole] = data[UP(hole)];
			pos[data[hole].second] = hole;
		}

		data[hole].first = t;
		data[hole].second = i;
		pos[i] = hole;
	}

	// take the min
	auto extract_min() {
		auto p = data[1];
		auto hole = 1;
		size--;
		data[hole] = data[size];
		// Sentile value at the last spot
		data[size].first = std::numeric_limits<cost_t>::max();
		// Percolate down
		auto tmp = data[1];
		auto child = 0;
		//for (; 2 * hole <= size; hole = child) {
		for (; (hole << 1) <= size; hole = child) {
			child = hole << 1;  // hole *2
			//if (child != size && data[child + 1].first < data[child].first)
			//if (data[child + 1].first < data[child].first)
			//   child++;
			child += (data[child + 1].first < data[child].first);

			if (data[child].first >= tmp.first)
				break;

			data[hole] = data[child];
			pos[data[hole].second] = hole;
		}

		data[hole] = tmp;
		pos[data[hole].second] = hole;
		// Return minimum element
		return p;
	}

	// if Dv is a smaller distance, then update key of node v
	void increase_key(node_t v, cost_t t) {
		size_t hole = pos[v];

		for (; hole > 1 && t < data[UP(hole)].first; hole = UP(hole)) {
			data[hole] = data[UP(hole)];
			pos[data[hole].second] = hole;
		}

		data[hole].first = t;
		data[hole].second = v;
		pos[v] = hole;
	}

private:
	// current and maximum size
	size_t size;
	size_t max_size;
	// Data array
	vector<std::pair<cost_t, node_t>>  data;
	// Vector of nodes of iterators
	vector<size_t> pos;
};


// Simple Arc class: store tuple (i,j,c)
class BiArc {
public:
	node_t      v;  // Source node
	node_t      w;  // Target node
	cost_t      l;  // Time to traverse the arc
	// Standard constructor
	BiArc(node_t _v, node_t _w, cost_t _l)
		: v(_v), w(_w), l(_l)
	{}

	// Get source node
	node_t getSource() const {
		return v;
	}

	// Get destination node
	node_t getTarget() const {
		return w;
	}

	// Get destination node
	cost_t getTime() const {
		return l;
	}
};


// Directed graph implementation
class AsDigraph {
	// Forward star
	typedef vector<arc_t>            ArcList;
	typedef ArcList::iterator        ArcIter;
	typedef ArcList::const_iterator  CArcIter;

public:
	// Standard c'tor
	AsDigraph() {}

	// Rule of five: Move constructor
	AsDigraph(AsDigraph&& o) : n(o.n), m(o.m), Fs(o.Fs), Bs(o.Bs), Ac(o.Ac)
	{}
	// Deleted constructor (this is a singleton object)
	AsDigraph(const AsDigraph& o) = delete;
	AsDigraph& operator=(const AsDigraph&) = delete;
	AsDigraph& operator=(AsDigraph&&) = delete;

	// get number of nodes
	node_t num_nodes() const {
		return n;
	}

	// get number of arcs
	arc_t num_arcs() const {
		return m;
	}

	size_t get_node_scans() const {
		return node_scans;
	}

	size_t get_arcs_scans() const {
		return arcs_scans;
	}

	// Init the graph and reserve memory
	void init(node_t _n, arc_t _m) {
		n = _n;
		m = _m;
		assert(n < node_last&& m < arc_last);
		// Reserve memory
		Fs.reserve(n);
		Bs.reserve(n);
		Ac.reserve(m);
		// Reserve memory for node list
		for (node_t i = 0; i < n; ++i) {
			Fs.emplace_back(ArcList());
			Bs.emplace_back(ArcList());
		}
	}

	// Add an arc to the current digraph
	void addArc(node_t i, node_t j, cost_t c) {
		arc_t idx = static_cast<arc_t>(Ac.size());
		Ac.emplace_back(i, j, c);
		Fs[i].push_back(idx);
		Bs[j].push_back(idx);
	}

	// Shrink the memory used
	void shrink_to_fit() {
		Fs.shrink_to_fit();
		Bs.shrink_to_fit();
		Ac.shrink_to_fit();
		for (auto& fs : Fs)
			fs.shrink_to_fit();
		for (auto& bs : Fs)
			bs.shrink_to_fit();
	}

	// Print this graph
	void print() const {
		fprintf(stdout, "p sp %d %d\n", int(n), int(m));

		for (size_t i = 0, i_max = Fs.size(); i < i_max; ++i)
			for (CArcIter it = Fs[i].begin(), it_end = Fs[i].end(); it != it_end; ++it)
				fprintf(stdout, "a %d %d %f\n", int(i), Ac[*it].getTarget(), int(Ac[*it].getTime()));
	}

	//--------------------------------------------------
	// Select randomly the distance from a single landmark
	void landmarkSetting() {
		// random naive test
		std::mt19937 gen(13);
		std::uniform_int_distribution<> Uniform0N(0, n - 1);

		auto landmark = Uniform0N(gen);
		LandmarkFS = spp_tree(landmark);
		LandmarkBS = spp_reverse_tree(landmark);
	}

	// compute the lower bound usingthe available landmarks
	cost_t getLowerBound(node_t v, node_t T) const {
		return std::max<cost_t>(LandmarkBS[v] - LandmarkBS[T], LandmarkFS[T] - LandmarkFS[v]);
	}

	//--------------------------------------------------
	// Shortest Path for a graph with positive weights
	//std::pair<cost_t, > 
	cost_t spp(node_t S, node_t T, vector<node_t>& sol) {
		// Priority queue
		BinaryHeap  H(n);
		// Predecessor vector
		vector<node_t> P(n, 0);
		P.shrink_to_fit();
		// Vector of labels
		vector<Label>  Q(n, UNREACHED);
		Q.shrink_to_fit();
		// Initialize the source distance
		size_t scans = 0;
		H.insert(getLowerBound(S, T), S);
		// Measure efficiency
		node_scans = 0;
		arcs_scans = 0;
		while (!H.empty()) {
			auto p = H.extract_min();
			node_t u = p.second;
			Q[u] = SCANNED;
			cost_t Du = p.first;
			node_scans++;
			if (u == T) {
				// Return the full path 
				sol.clear();
				// Quando separo non serve
				sol.push_back(T);
				node_t z = P[u];
				while (z != S) {
					sol.push_back(z);
					z = P[z];
				}
				//sol.push_back(S);

				return Du;
			}
			// for all edges (u, v) \in E
			for (ArcIter it = Fs[u].begin(), it_end = Fs[u].end(); it != it_end; ++it) {
				node_t v = Ac[*it].getTarget();
				arcs_scans++;
				if (Q[v] != SCANNED) {
					cost_t Duv = Ac[*it].getTime() - getLowerBound(u, T) + getLowerBound(v, T);
					cost_t Dv = Du + Duv;
					if (Q[v] == UNREACHED) {
						P[v] = u;  // predecessor
						Q[v] = LABELED;
						H.insert(Dv, v);
					}
					else {
						if (H.get_distance(v) > Dv) {
							P[v] = u;
							H.increase_key(v, Dv);
						}
					}
				}
			}
		}
		//fprintf(stdout, "not reacheble\n");

		// return time at destination
		return cost_max;
	}

	// Shortest Path Distance to node Q from every other node i
	vector<cost_t> spp_tree(node_t S) {
		// Priority queue
		BinaryHeap  H(n);
		// Predecessor vector
		vector<node_t> P(n, 0);
		P.shrink_to_fit();
		// Distance vector
		vector<cost_t> D(n, cost_max);
		D.shrink_to_fit();
		D[S] = 0;
		// Vector of labels
		vector<Label>  Q(n, UNREACHED);
		Q.shrink_to_fit();
		// Initialize the source distance
		H.insert(0, S);
		while (!H.empty()) {
			auto p = H.extract_min();
			node_t u = p.second;
			Q[u] = SCANNED;
			cost_t Du = p.first;
			// for all edges (u, v) \in E
			for (ArcIter it = Fs[u].begin(), it_end = Fs[u].end(); it != it_end; ++it) {
				node_t v = Ac[*it].getTarget();
				if (Q[v] != SCANNED) {
					cost_t Duv = Ac[*it].getTime();
					cost_t Dv = Du + Duv;
					if (Q[v] == UNREACHED) {
						P[v] = u;  // predecessor
						Q[v] = LABELED;
						D[v] = Dv;
						H.insert(Dv, v);
					}
					else {
						if (H.get_distance(v) > Dv) {
							P[v] = u;
							D[v] = Dv;
							H.increase_key(v, Dv);
						}
					}
				}
			}
		}

		// return time at destination
		return D;
	}

	// Shortest Path Distance to node Q from every other node i
	vector<cost_t> spp_reverse_tree(node_t S) {
		// Priority queue
		BinaryHeap  H(n);
		// Predecessor vector
		vector<node_t> P(n, 0);
		P.shrink_to_fit();
		// Distance vector
		vector<cost_t> D(n, cost_max);
		D.shrink_to_fit();
		D[S] = 0;
		// Vector of labels
		vector<Label>  Q(n, UNREACHED);
		Q.shrink_to_fit();
		// Initialize the source distance
		H.insert(0, S);
		while (!H.empty()) {
			auto p = H.extract_min();
			node_t u = p.second;
			Q[u] = SCANNED;
			cost_t Du = p.first;
			// for all edges (u, v) \in E
			for (ArcIter it = Bs[u].begin(), it_end = Bs[u].end(); it != it_end; ++it) {
				node_t v = Ac[*it].getSource();
				if (Q[v] != SCANNED) {
					cost_t Duv = Ac[*it].getTime();
					cost_t Dv = Du + Duv;
					if (Q[v] == UNREACHED) {
						P[v] = u;  // predecessor
						Q[v] = LABELED;
						D[v] = Dv;
						H.insert(Dv, v);
					}
					else {
						if (H.get_distance(v) > Dv) {
							P[v] = u;
							D[v] = Dv;
							H.increase_key(v, Dv);
						}
					}
				}
			}
		}

		// return time at destination
		return D;
	}

	// Find max arc cost
	cost_t findMaxTime() const {
		return std::max_element(Ac.begin(), Ac.end(), [](const auto& a, const auto& b) {
			return a.getTime() < b.getTime();
			})->getTime();
	}

	// Iterate over outgoing arc
	auto out_begin(node_t v) {
		return Fs[v].begin();
	}
	auto out_end(node_t v) {
		return Fs[v].end();
	}
	auto getTarget(ArcList::iterator it) const {
		return Ac[*it].getTarget();
	}
	auto getTime(ArcList::iterator it) const {
		return Ac[*it].getTime();
	}

private:
	node_t  n = 0;
	arc_t   m = 0;

	vector<ArcList>  Fs;   // Nodes container
	vector<ArcList>  Bs;   // Nodes container
	vector<BiArc>    Ac;   // Arcs container

	// Measure efficiency
	size_t node_scans = 0;
	size_t arcs_scans = 0;

	// Landmarks
	node_t landmark;
	vector<cost_t> LandmarkFS;
	vector<cost_t> LandmarkBS;
};
///---------------------------------------------------
class SmallClique {
public:
	SmallClique(const BitGraph& H) : n(H.num_vertices()), H0(H) {
		Fs.resize(n);
		for (int u = 0; u < n; u++) {
			Fs.reserve(n);
			for (int v = u + 1; v < n; v++)
				if (H.is_edge(u, v))
					Fs[u].push_back(v);
		}
	}

	array<int, 3> separateK3(double* xbar) const {
		array<int, 3> r = { -1, -1, -1 };
		double viol = 1.0 + INT_TOL;

		for (int u = 0; u < n; u++) {
			for (int v: Fs[u])
				for (int w: Fs[v])
					if (H0.is_edge(u, w) && xbar[u] + xbar[v] + xbar[w] > viol) {
						viol = xbar[u] + xbar[v] + xbar[w];
						r = { v, u, w };
					}
		}

		return r;
	}

	array<int, 4> separateK4(double* xbar) const {
		array<int, 4> r = { -1, -1, -1, -1 };
		double viol = 1.0 + INT_TOL;

		for (int u = 0; u < n; ++u) {
			for (int v : Fs[u])
				for (int w : Fs[v])
					if (H0.is_edge(u, w)) {
						for (int z : Fs[w])
							if (H0.is_edge(u, z) && H0.is_edge(v, z) && xbar[u] + xbar[v] + xbar[w] + xbar[z] > viol) {
								viol = xbar[u] + xbar[v] + xbar[w] + xbar[z];
								r = { u, v, w, z };
							}
					}
		}

		return r;
	}

private:
	const int n;
	vector<vector<int>> Fs;

	BitGraph H0;
};

///---------------------------------------------------
class SeparatorOddCycle {
public:
	SeparatorOddCycle(const BitGraph& H) : cutter_smallclique(H) {
		n = H.num_vertices();
		V1.resize(n, 0);
		V2.resize(n, 0);
		B.resize(n * 2, 0);

		for (int v = 0; v < n; v++) {
			V1[v] = v;
			V2[v] = n + v;
			B[v] = v;
			B[n + v] = v;
		}

	}

	// Separate a xbar returning an error code
	int separate(const BitGraph& H, double* x_bar, GRBmodel* master,
		const std::string& msg) {
		auto start_time = std::chrono::steady_clock::now();

		// First try to separate triangles
		auto tri = cutter_smallclique.separateK3(x_bar);
		if (tri[0] != -1) {
			int* ind = (int*)malloc(n * sizeof(int));
			double* val = (double*)malloc(n * sizeof(double));
			int nz = 0;
			for (node_t i : tri) {
				val[nz] = 1.0;
				ind[nz] = int(i);
				nz++;
			}

			POST("add triangle", GRBaddconstr(master, nz, ind, val, GRB_LESS_EQUAL, double(nz - 1) / 2, (const char*)NULL));

			free(ind);
			free(val);

			auto end = std::chrono::steady_clock::now();
			auto elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start_time).count()) / 1000;
			fprintf(stdout, "%s Violation: %.3f |U| = %d Time %.3f\n",
				msg.c_str(), 0.0, (int)tri.size(), elapsed);

			return 1;
		}

		// Separate with general bipartite graph
		AsDigraph G;
		G.init(2 * n, 4 * n * n);
		for (int u = 0; u < n; u++)
			for (int v = u + 1; v < n; v++)
				if (H.is_edge(u, v)) {
					double l = std::max<double>(1e-09, 1.0 - x_bar[u] - x_bar[v]);

					G.addArc(V1[u], V2[v], l);
					G.addArc(V1[v], V2[u], l);

					G.addArc(V2[v], V1[u], l);
					G.addArc(V2[u], V1[v], l);
				}
		G.landmarkSetting();

		double bv = 1.0 - CUT_VIOL;
		size_t best_l = n + 1;
		vector<node_t> best_cut;
		for (size_t v = 0; v < n; v++) {
			vector<node_t> tc;
			cost_t p = G.spp(V1[v], V2[v], tc);
			if (p <= bv && tc.size() < best_l) {
				bv = p;
				best_l = tc.size();
				best_cut.clear();
				for (node_t u : tc) {
					best_cut.push_back(B[u]);
				}
			}
		}

		int ret_code = 0;

		if (!best_cut.empty()) {
			int* ind = (int*)malloc(n * sizeof(int));
			double* val = (double*)malloc(n * sizeof(double));
			int nz = 0;
			for (node_t i : best_cut) {
				val[nz] = 1.0;
				ind[nz] = int(i);
				nz++;
			}

			POST("add cut", GRBaddconstr(master, nz, ind, val, GRB_LESS_EQUAL, double(nz - 1) / 2, (const char*)NULL));

			free(ind);
			free(val);

			ret_code = 1;

			auto end = std::chrono::steady_clock::now();
			auto elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start_time).count()) / 1000;
			fprintf(stdout, "%s Violation: %.3f |U| = %d Time %.3f\n",
				msg.c_str(), bv, (int)best_cut.size(), elapsed);
		}


		return ret_code;
	}

private:
	size_t n;

	vector<int> V1;
	vector<int> V2;
	vector<int> B;

	SmallClique cutter_smallclique;
};

///---------------------------------------------------
class SeparatorOddWheel {
public:
	SeparatorOddWheel(const BitGraph& H) : cutter_smallclique(H) {
		n = H.num_vertices();
	}

	// Separate a xbar returning an error code
	int separate(const BitGraph& H, double* x_bar, GRBmodel* master,
		const std::string& msg) {
		auto start_time = std::chrono::steady_clock::now();

		// First try to separate K4
		auto tri = cutter_smallclique.separateK4(x_bar);

		if (tri[0] != -1) {
			int* ind = (int*)malloc(n * sizeof(int));
			double* val = (double*)malloc(n * sizeof(double));
			int nz = 0;
			for (node_t i : tri) {
				val[nz] = 1.0;
				ind[nz] = int(i);
				nz++;
			}

			POST("add quadri", GRBaddconstr(master, nz, ind, val, GRB_LESS_EQUAL, 1.0, (const char*)NULL));

			free(ind);
			free(val);

			auto end = std::chrono::steady_clock::now();
			auto elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start_time).count()) / 1000;
			fprintf(stdout, "%s K4 Violation: %.3f |U| = %d Time %.3f\n",
				msg.c_str(), 0.0, (int)tri.size(), elapsed);

			return 1;
		}

		// Separate with general bipartite graph
		int ret_code = 0;

		std::set<int> U;
		vector<int> V;
		vector<int> W;
		V.reserve(n);
		W.reserve(n);

		int bw = -1;
		double bvv = 1.0 - CUT_VIOL;
		size_t best_l = n + 1;
		vector<node_t> best_cut;

		for (int w = 0; w < n; w++) {
			U.clear();
			V.clear();
			W.clear();

			for (int u = 0; u < n; u++)
				if (u != w && H.is_edge(u, w))
					for (int v = 0; v < n; v++)
						if (v != w && v != u && H.is_edge(w, v) && H.is_edge(u, v)) {
							U.insert(u);
							U.insert(v);
							V.push_back(u);
							W.push_back(v);
						}

			int t = U.size();
			vector<int> V1(n, 0);
			vector<int> V2(n, 0);
			vector<int> B(2 * t, 0);

			int c = 0;
			for (int u : U) {
				V1[u] = c;
				V2[u] = c + 1;
				B[c] = u;
				B[c + 1] = u;
				c = c + 2;
			}

			if (U.size() >= 3 && V.size() >= 3) {
				AsDigraph G;

				G.init(2 * t, 4 * W.size());
				
				for (int i = 0, i_max = W.size(); i < i_max; ++i) {
					double l = std::max<double>(
						1e-09,
						1.0 - x_bar[V[i]] - x_bar[W[i]] - x_bar[w]
						);

					G.addArc(V1[V[i]], V2[W[i]], l);
					G.addArc(V1[W[i]], V2[V[i]], l);

					G.addArc(V2[W[i]], V1[V[i]], l);
					G.addArc(V2[V[i]], V1[W[i]], l);
				}

				G.landmarkSetting();

				double bv = 1.0 - x_bar[w] - CUT_VIOL;

				for (int v : U) {
					vector<node_t> tc;
					cost_t p = G.spp(V1[v], V2[v], tc);
					if (p <= bv && p < bvv && tc.size() < best_l && !tc.empty()) {
						bv = p;
						bvv = p;
						bw = w;
						best_l = tc.size();
						best_cut.clear();
						for (node_t u : tc) {
							best_cut.push_back(B[u]);
						}
					}
				}
			}
		}

		if (best_cut.size() > 2) {
			int* ind = (int*)malloc(n * sizeof(int));
			double* val = (double*)malloc(n * sizeof(double));
			int nz = 0;
			for (node_t i : best_cut) {
				val[nz] = 1.0;
				ind[nz] = int(i);
				nz++;
			}
			double b = 0.5 * (double(nz) - 1.0);
			val[nz] = b;
			ind[nz] = int(bw);
			nz++;

			POST("add odd wheel", GRBaddconstr(master, nz, ind, val, GRB_LESS_EQUAL, b, (const char*)NULL));

			free(ind);
			free(val);

			ret_code = 1;

			auto end = std::chrono::steady_clock::now();
			auto elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start_time).count()) / 1000;
			fprintf(stdout, "%s Wheel Violation: %.3f |U| = %d Time %.3f\n",
				msg.c_str(), bvv, nz, elapsed);

			return 1;
		}

		return ret_code;
	}

private:
	int n;

	SmallClique cutter_smallclique;
};