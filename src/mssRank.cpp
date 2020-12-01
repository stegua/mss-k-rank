/*
*  Main authors:
*     Stefano Gualandi <stefano.gualandi@gmail.com>
*
*  Copyright:
*     Stefano Gualandi, 2017
*
*  Last update: November, 2020
*/

#include <string>
#include <cstdlib>

#include "Utils.h"

#include "separators.h"
#include "DigraphSpp.h"

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
	int k) {
	VectorBitSet r(fst.size());
	r.triple_intersect(fst, snd, Cs);

	if (r.count() < k)
		return true;

	for (const bitset_t& d : Ts)
		if (r.more_than_k(d, k))
			return false;

	return true;
}

/**
* @brief Solution is integer
*/
bool isInteger(const double* x_bar, int n) {
	for (int i = 0; i < n; i++)
		if (x_bar[i] > INT_TOL&& x_bar[i] < 1.0 - INT_TOL)
			return false;
	return true;
}

/**
* @brief Solve LP relaxation (master) of k-rank relaxations
*/
void solveLP(BitGraph& G, double alpha, bool BnC) {
	// Random generator
	std::mt19937 gen(13);  // Seed = 13
	std::uniform_real_distribution<> perturba(1e-07, 1e-05); // white noise

	// Copy graph
	graph_t* g;
	graph_t* h;
	int n = G.num_vertices();

	// Cliquer graph
	g = graph_new(n);
	h = graph_new(n);
	for (int i = 0; i < n - 1; ++i)
		for (int j = i + 1; j < n; ++j) {
			if (G.is_edge(i, j))
				GRAPH_ADD_EDGE(g, i, j);
			else
				GRAPH_ADD_EDGE(h, i, j);
		}

	int      status;
	int      n_iter = 0;
	double   alpha_lp;

	// Pointers for Gurobi
	GRBenv* env_master;
	GRBmodel* master;
	// Arrays for passing data to Gurobi
	char* le = (char*)malloc(sizeof(char) * (n + 1));
	int* ind = (int*)malloc(sizeof(int) * (n + 1));
	double* val = (double*)malloc(sizeof(double) * (n + 1));
	double* ub = (double*)malloc(sizeof(double) * n);
	double* xbar = (double*)malloc(sizeof(double) * n);
	double* oj = (double*)malloc(sizeof(double) * n);
	double* pert = (double*)malloc(sizeof(double) * n);

	for (int i = 0; i < n; ++i) {
		ub[i] = 1.0;
		pert[i] = 1;
		oj[i] = 1.0;
	}

	// Build master problem
	GRBloadenv(&env_master, NULL);
	GRBsetintparam(env_master, GRB_INT_PAR_OUTPUTFLAG, 0);
	GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
	GRBsetintparam(env_master, GRB_INT_PAR_THREADS, 1);
	GRBnewmodel(env_master, &master, "mss", n, oj, NULL, ub, NULL, NULL);

	GRBsetintattr(master, "ModelSense", GRB_MAXIMIZE);

	// Solve the problem
	int ccs = 0;
	int U_card = 0;
	double rhs = 0.0;
	double violation = 0.0;
	auto start = std::chrono::steady_clock::now();
	auto end = std::chrono::steady_clock::now();
	double elapsed;

	/* Edge constraints */
	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			if (G.is_edge(i, j)) {
				ind[0] = i;
				val[0] = 1.0;
				ind[1] = j;
				val[1] = 1.0;

				POST("add edge constraint", GRBaddconstr(master, 2, ind, val, GRB_LESS_EQUAL, 1, ""));
			}
		}
	}

	double alpha0 = 1.0;
	vector<double> cutvs = { 0.5, 0.25, 0.1, 0.05, 0.01 };
	size_t idx = 0;
	elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	bool isInt = false;
	double nl = 0;

	while (elapsed < TIMEOUT && alpha0 <= alpha) {
		GRBoptimize(master);

		GRBgetintattr(master, "Status", &status);
		if (status == GRB_UNBOUNDED || status == GRB_INFEASIBLE) {
			fprintf(stdout, "Unbounded or Infeasible\n");
			goto QUIT;
		}
		// Take the current LP decision vector
		POST("get X", GRBgetdblattrarray(master, "X", 0, n, xbar));
		POST("get obj", GRBgetdblattr(master, "ObjVal", &alpha_lp));

		// Check if the solution is integer
		if (isInteger(xbar, n)) {
			isInt = true;
			break;
		}

		end = std::chrono::steady_clock::now();
		std::string msg = "LogLocal: it " + std::to_string(n_iter++) + " Bound " + std::to_string(alpha_lp);
		// Incremental timeout
		double t0 = TIMEOUT - elapsed;
		if (t0 < 1e-03) break;

		for (int i = 0; i < n; ++i)
			xbar[i] += perturba(gen);  // Diversity a random in [1e-05 ... 1e0-7]
				
				int sss = 0;
				while (alpha0 <= alpha) {
					if (BnC)
						sss = separatorBnC(g, h, xbar, master, t0, 0, alpha0, cutvs[idx], start, msg);
					else
						sss = separatorBnB(&G, xbar, master, alpha0, cutvs[idx], start, msg);

					if (sss > 0) {
						if (sss == 2)
							if (idx < 4)
								idx++;
						break;
					}

					if (idx < 4)
						idx++;
					else {
						if (alpha0 <= alpha0) {
							alpha0++;
							idx = 0;
							if (alpha0 == 3)
								idx = 2;
							if (alpha0 >= 4)
								idx = 3;
						}
						else
							break;
					}
				}			
	}

	end = std::chrono::steady_clock::now();
	elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	printf("LogGlobal: it %d Bound %.4f Time %.3f alpha: %.f integer: %d\n",
		n_iter++, alpha_lp, elapsed, std::min(alpha0, alpha), isInt);

	/*for (int i = 0; i < n; ++i)
		POST("set to binary", GRBsetcharattrelement(master, GRB_CHAR_ATTR_VTYPE, i, GRB_BINARY));
	GRBupdatemodel(master);

	GRBoptimize(master);

	end = std::chrono::steady_clock::now();
	elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

	POST("get obj", GRBgetdblattr(master, "ObjVal", &alpha_lp));
	POST("get obj p", GRBgetdblattr(master, "NodeCount", &nl));
	fprintf(stdout, "ott: %f time: %f nodes: %.f\n", alpha_lp, elapsed, nl);*/
QUIT:
	free(ind);
	free(val);
	free(xbar);
	free(oj);
	free(le);
	free(ub);
	free(pert);

	POST("free", GRBfreemodel(master));
	GRBfreeenv(env_master);
}


/**
* @brief Solve LP relaxation (master) of k-rank relaxations
*/
void solveLP_oddcycle(BitGraph& G, bool odd_wheel = false) {
	// Random generator
	std::mt19937 gen(13);  // Seed = 13
	std::uniform_real_distribution<> perturba(1e-07, 1e-05); // white noise

	// Copy graph
	int n = G.num_vertices();

	int      status;
	int      n_iter = 0;
	double   alpha_lp;

	// Pointers for Gurobi
	GRBenv* env_master;
	GRBmodel* master;
	// Arrays for passing data to Gurobi
	char* le = (char*)malloc(sizeof(char) * (n + 1));
	int* ind = (int*)malloc(sizeof(int) * (n + 1));
	double* val = (double*)malloc(sizeof(double) * (n + 1));
	double* ub = (double*)malloc(sizeof(double) * n);
	double* xbar = (double*)malloc(sizeof(double) * n);
	double* oj = (double*)malloc(sizeof(double) * n);
	double* pert = (double*)malloc(sizeof(double) * n);

	for (int i = 0; i < n; ++i) {
		ub[i] = 1.0;
		pert[i] = 1;
		oj[i] = 1.0;
	}

	// Build master problem
	GRBloadenv(&env_master, NULL);
	GRBsetintparam(env_master, GRB_INT_PAR_OUTPUTFLAG, 0);
	GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
	GRBsetintparam(env_master, GRB_INT_PAR_THREADS, 1);

	GRBnewmodel(env_master, &master, "mss", n, oj, NULL, ub, NULL, NULL);

	GRBsetintattr(master, "ModelSense", GRB_MAXIMIZE);

	// Solve the problem
	int ccs = 0;
	int U_card = 0;
	double rhs = 0.0;
	double violation = 0.0;
	auto start = std::chrono::steady_clock::now();
	auto end = std::chrono::steady_clock::now();
	double elapsed;

	/* Edge constraints */
	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			if (G.is_edge(i, j)) {
				ind[0] = i;
				val[0] = 1.0;
				ind[1] = j;
				val[1] = 1.0;

				POST("add edge constraint", GRBaddconstr(master, 2, ind, val, GRB_LESS_EQUAL, 1, ""));
			}
		}
	}

	size_t idx = 0;
	elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	bool isInt = false;
	double nl = 0;

	SeparatorOddCycle cutter(G);
	SeparatorOddWheel cutter_wheel(G);

	while (elapsed < TIMEOUT) {
		GRBoptimize(master);

		GRBgetintattr(master, "Status", &status);
		if (status == GRB_UNBOUNDED || status == GRB_INFEASIBLE) {
			fprintf(stdout, "Unbounded or Infeasible\n");
			goto QUIT;
		}
		// Take the current LP decision vector
		POST("get X", GRBgetdblattrarray(master, "X", 0, n, xbar));
		POST("get obj", GRBgetdblattr(master, "ObjVal", &alpha_lp));

		// Check if the solution is integer
		if (isInteger(xbar, n)) {
			isInt = true;
			break;
		}

		end = std::chrono::steady_clock::now();
		elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
		
		std::string msg = "LogLocal: it " + std::to_string(n_iter++) + " Bound " + std::to_string(alpha_lp);
		
		// Incremental timeout
		double t0 = TIMEOUT - elapsed;
		if (t0 < 1e-03) break;

		for (int i = 0; i < n; ++i)
			xbar[i] += perturba(gen);  // Diversity a random in [1e-05 ... 1e0-7]

		int res_code = cutter.separate(G, xbar, master, msg);

		if (res_code == 0 && odd_wheel == false)
			break;

		if (res_code == 0) {
			res_code = cutter_wheel.separate(G, xbar, master, msg);
			if (res_code == 0)
				break;
		}
	}

	end = std::chrono::steady_clock::now();
	elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	printf("LogGlobal: it %d Bound %.4f Time %.3f integer: %d\n",
		n_iter++, alpha_lp, elapsed, isInt);

QUIT:
	free(ind);
	free(val);
	free(xbar);
	free(oj);
	free(le);
	free(ub);
	free(pert);

	POST("free", GRBfreemodel(master));
	GRBfreeenv(env_master);
}


/**
* @brief Solve LP relaxation (master) of k-rank relaxations
*/
void solveLP_zero_half(BitGraph& G, bool odd_wheel = false) {
	// Random generator
	std::mt19937 gen(13);  // Seed = 13
	std::uniform_real_distribution<> perturba(1e-07, 1e-05); // white noise

	// Copy graph
	int n = G.num_vertices();

	int      status;
	int      n_iter = 0;
	double   alpha_lp;

	// Pointers for Gurobi
	GRBenv* env_master;
	GRBmodel* master;
	// Arrays for passing data to Gurobi
	char* le = (char*)malloc(sizeof(char) * (n + 1));
	int* ind = (int*)malloc(sizeof(int) * (n + 1));
	double* val = (double*)malloc(sizeof(double) * (n + 1));
	double* ub = (double*)malloc(sizeof(double) * n);
	double* xbar = (double*)malloc(sizeof(double) * n);
	double* oj = (double*)malloc(sizeof(double) * n);
	double* pert = (double*)malloc(sizeof(double) * n);
	char* vtype = (char*)malloc(sizeof(char) * n);

	for (int i = 0; i < n; ++i) {
		ub[i] = 1.0;
		pert[i] = 1;
		oj[i] = 1.0;
		vtype[i] = GRB_BINARY;
	}

	// Build master problem
	GRBloadenv(&env_master, NULL);
	GRBsetintparam(env_master, GRB_INT_PAR_OUTPUTFLAG, 0);
	GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
	GRBsetdblparam(env_master, GRB_DBL_PAR_NODELIMIT, 1);
	GRBsetintparam(env_master, GRB_INT_PAR_THREADS, 1);
	GRBsetintparam(env_master, GRB_INT_PAR_PRESOLVE, 0);
	GRBsetintparam(env_master, GRB_INT_PAR_DUALREDUCTIONS, 0);
	GRBsetdblparam(env_master, GRB_DBL_PAR_HEURISTICS, 0);
	GRBsetintparam(env_master, GRB_INT_PAR_CUTS, 0);
	GRBsetintparam(env_master, GRB_INT_PAR_ZEROHALFCUTS, 2);

	GRBnewmodel(env_master, &master, "mss", n, oj, NULL, ub, vtype, NULL);

	GRBsetintattr(master, "ModelSense", GRB_MAXIMIZE);

	// Solve the problem
	auto start = std::chrono::steady_clock::now();
	auto end = std::chrono::steady_clock::now();
	double elapsed;

	/* Edge constraints */
	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			if (G.is_edge(i, j)) {
				ind[0] = i;
				val[0] = 1.0;
				ind[1] = j;
				val[1] = 1.0;

				POST("add edge constraint", GRBaddconstr(master, 2, ind, val, GRB_LESS_EQUAL, 1, ""));
			}
		}
	}

	size_t idx = 0;
	elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	bool isInt = false;
	double nl = 0;

	GRBoptimize(master);

	GRBgetintattr(master, "Status", &status);
	if (status == GRB_UNBOUNDED || status == GRB_INFEASIBLE) {
		fprintf(stdout, "Unbounded or Infeasible\n");
		goto QUIT;
	}
	// Take the current LP decision vector
	POST("get obj", GRBgetdblattr(master, "ObjBoundC", &alpha_lp));

	// Check if the solution is integer
	if (isInteger(xbar, n)) {
		isInt = true;
	}

	end = std::chrono::steady_clock::now();
	elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	printf("LogGlobal: it %d Bound %.4f Time %.3f integer: %d\n",
		n_iter++, alpha_lp, elapsed, isInt);

QUIT:
	free(ind);
	free(val);
	free(xbar);
	free(oj);
	free(le);
	free(ub);
	free(pert);

	POST("free", GRBfreemodel(master));
	GRBfreeenv(env_master);
}

/**
* @brief Solve LP relaxation (master) of k-rank relaxations
*/
void solveLP_clique_cycle(BitGraph& G, bool odd_wheel = false) {
	// Random generator
	std::mt19937 gen(13);  // Seed = 13
	std::uniform_real_distribution<> perturba(1e-07, 1e-05); // white noise

	// Copy graph
	int n = G.num_vertices();

	int      status;
	int      n_iter = 0;
	double   alpha_lp;

	// Pointers for Gurobi
	GRBenv* env_master;
	GRBmodel* master;
	// Arrays for passing data to Gurobi
	char* le = (char*)malloc(sizeof(char) * (n + 1));
	int* ind = (int*)malloc(sizeof(int) * (n + 1));
	double* val = (double*)malloc(sizeof(double) * (n + 1));
	double* ub = (double*)malloc(sizeof(double) * n);
	double* xbar = (double*)malloc(sizeof(double) * n);
	double* oj = (double*)malloc(sizeof(double) * n);
	double* pert = (double*)malloc(sizeof(double) * n);

	for (int i = 0; i < n; ++i) {
		ub[i] = 1.0;
		pert[i] = 1;
		oj[i] = 1.0;
	}

	// Build master problem
	GRBloadenv(&env_master, NULL);
	GRBsetintparam(env_master, GRB_INT_PAR_OUTPUTFLAG, 0);
	GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
	GRBsetintparam(env_master, GRB_INT_PAR_THREADS, 1);

	GRBnewmodel(env_master, &master, "mss", n, oj, NULL, ub, NULL, NULL);

	GRBsetintattr(master, "ModelSense", GRB_MAXIMIZE);

	// Solve the problem
	int ccs = 0;
	int U_card = 0;
	double rhs = 0.0;
	double violation = 0.0;
	auto start = std::chrono::steady_clock::now();
	auto end = std::chrono::steady_clock::now();
	double elapsed;

	/* Edge constraints */
	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			if (G.is_edge(i, j)) {
				ind[0] = i;
				val[0] = 1.0;
				ind[1] = j;
				val[1] = 1.0;

				POST("add edge constraint", GRBaddconstr(master, 2, ind, val, GRB_LESS_EQUAL, 1, ""));
			}
		}
	}

	size_t idx = 0;
	elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	bool isInt = false;
	double nl = 0;

	//G.print();
	SeparatorOddCycle cutter(G);
	SeparatorOddWheel cutter_wheel(G);
	vector<double> cutvs = { 0.5, 0.25, 0.1, 0.05, 0.01 };

	while (elapsed < TIMEOUT) {
		GRBoptimize(master);

		GRBgetintattr(master, "Status", &status);
		if (status == GRB_UNBOUNDED || status == GRB_INFEASIBLE) {
			fprintf(stdout, "Unbounded or Infeasible\n");
			goto QUIT;
		}
		// Take the current LP decision vector
		POST("get X", GRBgetdblattrarray(master, "X", 0, n, xbar));
		POST("get obj", GRBgetdblattr(master, "ObjVal", &alpha_lp));

		// Check if the solution is integer
		if (isInteger(xbar, n)) {
			isInt = true;
			break;
		}

		end = std::chrono::steady_clock::now();
		elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

		std::string msg = "LogLocal: it " + std::to_string(n_iter++) + " Bound " + std::to_string(alpha_lp);

		// Incremental timeout
		double t0 = TIMEOUT - elapsed;
		if (t0 < 1e-03) break;

		for (int i = 0; i < n; ++i)
			xbar[i] += perturba(gen);  // Diversity a random in [1e-05 ... 1e0-7]

		int res_code = 0;
		res_code = separatorBnB(&G, xbar, master, 1, CUT_VIOL, start, msg);

		if (res_code == 0) {
			idx = 0;
			res_code = cutter.separate(G, xbar, master, msg);
			if (res_code == 0 && odd_wheel == false)
				break;

			if (res_code == 0) {
				res_code = cutter_wheel.separate(G, xbar, master, msg);
				if (res_code == 0)
					break;
			}
		}
	}

	end = std::chrono::steady_clock::now();
	elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	printf("LogGlobal: it %d Bound %.4f Time %.3f integer: %d\n",
		n_iter++, alpha_lp, elapsed, isInt);

QUIT:
	free(ind);
	free(val);
	free(xbar);
	free(oj);
	free(le);
	free(ub);
	free(pert);

	POST("free", GRBfreemodel(master));
	GRBfreeenv(env_master);
}


/*******************************************************************************************/
int main(int argc, char** argv) {
	if (argc < 3) {
		fprintf(stdout, "\nusage:  $ ./mssRank <filename> <alpha> <complement>\n\n");
		exit(EXIT_SUCCESS);
	}

	auto start = std::chrono::steady_clock::now();

	BitGraph G;
	read_graph(argv[1], G);

	BitGraph H;
	auto n = G.num_vertices();
	H.init(n);
	for (int i = 0; i < n - 1; i++)
		for (int j = i + 1; j < n; j++) {
			if (argc == 4) { // COMPLEMENT argv[3]
				if (!G.is_edge(i, j))
					H.add_edge(i, j);
			}
			else { // Copy graph, no complement
				if (G.is_edge(i, j))
					H.add_edge(i, j);
			}
		}

	// TODO: un codice unico per gli algoritmi con un flag vero per il complemento

	// list di cose da testare:

	double algo = atof(argv[2]);
	if (algo < 1) {
		fprintf(stdout, "\nusage:  $ ./mssRank <filename> <alpha> <complement>  ===> with alpha>1\n\n");
		exit(EXIT_SUCCESS);
	}

	if (algo <= 5) 
		solveLP(H, algo, false); // False ==> BaB
	else {
		if (algo <= 10)
			solveLP(H, algo - 5, true); // True ==> BaC 
		else {
			if (algo == 11)
				solveLP_oddcycle(H, false); // only odd-hole 
			if (algo == 12)
				solveLP_oddcycle(H, true); // only odd-hole + odd-wheel
			if (algo == 13)
				solveLP_zero_half(H); // zero-one half aggressive Gurobi
			if (algo == 14)
				solveLP_clique_cycle(H, false); // only odd-hole + odd-wheel
		}
	}

	//solveBranchAndCut(H, atof(argv[2]));

	auto end = std::chrono::steady_clock::now();
	auto elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

	return EXIT_SUCCESS;
}