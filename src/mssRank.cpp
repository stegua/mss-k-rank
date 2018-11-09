/*
*  Main authors:
*     Stefano Gualandi <stefano.gualandi@gmail.com>
*
*  Copyright:
*     Stefano Gualandi, 2017
*
*  Last update: November, 2018
*/

#include "Utils.h"

#include "separators.h"

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
      if (x_bar[i] > INT_TOL && x_bar[i] < 1.0 - INT_TOL)
         return false;
   return true;
}

/**
* @brief Solve LP relaxation (master)
*/
void solveLP(BitGraph & G, double alpha, bool BnC) {
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
   GRBenv   *env_master;
   GRBmodel *master;
   // Arrays for passing data to Gurobi
   char*     le = (char*)malloc(sizeof(char)     * (n + 1));
   int*     ind = (int*)malloc(sizeof(int)       * (n + 1));
   double*  val = (double*)malloc(sizeof(double) * (n + 1));
   double*   ub = (double*)malloc(sizeof(double) * n);
   double* xbar = (double*)malloc(sizeof(double) * n);
   double*   oj = (double*)malloc(sizeof(double) * n);
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
   auto end   = std::chrono::steady_clock::now();
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
            } else
               break;
         }
      }
   }

   end = std::chrono::steady_clock::now();
   elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
   printf("LogGlobal: it %d Bound %.4f Time %.3f alpha: %.f integer: %d\n",
          n_iter++, alpha_lp, elapsed, std::min(alpha0,alpha), isInt);

   for (int i = 0; i < n; ++i)
      POST("set to binary", GRBsetcharattrelement(master, GRB_CHAR_ATTR_VTYPE, i, GRB_BINARY));
   GRBupdatemodel(master);

   GRBoptimize(master);

   end = std::chrono::steady_clock::now();
   elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   POST("get obj", GRBgetdblattr(master, "ObjVal", &alpha_lp));
   POST("get obj p", GRBgetdblattr(master, "NodeCount", &nl));
   fprintf(stdout, "ott: %f time: %f nodes: %.f\n", alpha_lp, elapsed, nl);
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
int main(int argc, char **argv) {
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
         } else { // Copy graph, no complement
            if (G.is_edge(i, j))
               H.add_edge(i, j);
         }
      }

   solveLP(H, atof(argv[2]), true);

   //solveBranchAndCut(H, atof(argv[2]));

   auto end = std::chrono::steady_clock::now();
   auto elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   return EXIT_SUCCESS;
}