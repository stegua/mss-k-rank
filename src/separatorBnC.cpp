/*
*  Main authors:
*     Stefano Gualandi <stefano.gualandi@gmail.com>
*
*  Copyright:
*     Stefano Gualandi, 2017
*
*  Last update: November, 2018
*/

#include "separators.h"

// Use Cliquer to find all maximal clique in the conflict graph
extern "C" {
#include "cliquer.h"
#include "set.h"
#include "graph.h"
#include "reorder.h"
   /* Records a clique into the clique list using dynamic allocation.
    * Used as opts->user_function.  */
   static int clique_count = 0;
   static set_t *clique_list;
   static int clique_list_size = 0;
   boolean record_clique_func(set_t s, graph_t *g, clique_options * opts) {
      if (clique_count >= clique_list_size) {
         clique_list = (set_t*)realloc(clique_list, (clique_list_size + 512) * sizeof(set_t));
         clique_list_size += 512;
      }
      clique_list[clique_count] = set_duplicate(s);
      clique_count++;
      return 1;
   }
}

struct callback_data {
   graph_t* g;
   int n;
   double alpha;
	double runtime;
};


/**
* @brief Separation lazycallback for the separation procedure
*/
static set_t
find_alpha_cut(
   int       n0, // Number of vertices in the original graph
   int       n,  // Number of vertices in vector V
   graph_t*  h,  // Complement of the input graph (where clique are stable set)
   int*      V,  // Subset of vertices that gives the subgraph
//      set_t     s,  /// New violated stable set
   double alpha,
	double runtime) {
   int i, j;
   graph_t* F = graph_new(n);
   set_t    r = NULL;
   set_t    s = NULL;
   clique_options* opts;
   // Set the options for using Cliquer
   opts = (clique_options*)malloc(sizeof(clique_options));
   // Set incremental timeout
   double* timeout;
   timeout = (double*)malloc(sizeof(double));
   timeout[0] = runtime;
   opts->user_data = timeout;
   opts->time_function = clique_time_out;
   opts->reorder_function = reorder_by_greedy_coloring;
   opts->reorder_map = NULL;
   // Build Support Graph
   for (i = 0; i < n - 1; ++i)
      for (j = i + 1; j < n; ++j)
         if (GRAPH_IS_EDGE_FAST(h, V[i], V[j]))
            GRAPH_ADD_EDGE(F, i, j);
   // Find violated cut
   r = clique_find_single(F, alpha + 1, 0, TRUE, opts);
   // Map Clique to original graph
   if (r != NULL && set_size(r) > 0) {
      s = set_new(n0);
      for (i = 0; i < n; ++i)
         if (SET_CONTAINS_FAST(r, i))
            SET_ADD_ELEMENT(s, V[i]);
      maximalize_clique_random(s, h);
   }

   // Free memory
   if (r != NULL)
      set_free(r);
   graph_free(F);
   free(opts);

   return s;
}


/**
* @brief Lazycut callback for generating violated inequalities
*        in the separator: looks on the support graph
*        stable set of cardinalities \alpha+1
*/
int __stdcall
subsetelim(GRBmodel *model,
           void     *cbdata,
           int       where,
           void     *usrdata) {
   struct  callback_data *mydata = (struct callback_data *) usrdata;
   double  alpha = mydata->alpha;
   double  runtime = mydata->runtime;
   int     n = mydata->n;
   int*    V;
   double* xbar;
   int i, len, nz;
   int error = 0;

   if (where == GRB_CB_MIPSOL) {
      xbar = (double *)malloc((n + 1) * sizeof(double));
      V = (int*)malloc((n + 1) * sizeof(int));
      if (xbar == NULL) {
         fprintf(stderr, "Out of memory\n");
         exit(1);
      }

      GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, xbar);

      // The solution must be integer
      nz = 0;
      for (i = 0; i < n; ++i) {
         if (xbar[i] > 0.98)
            V[nz++] = i;
      }

      if (nz > 0) {
         set_t s = find_alpha_cut(n, nz, mydata->g, V, alpha, runtime);
         if (s != NULL && set_size(s) > alpha) {
            //set_print(s);
            int    *ind = NULL;
            double *val = NULL;
            len = set_size(s);
            ind = (int *)malloc(len * sizeof(int));
            val = (double *)malloc(len * sizeof(double));
            if (ind == NULL || val == NULL) {
               fprintf(stderr, "Out of memory\n");
               exit(1);
            }
            // Add violated ineq.
            nz = 0;
            for (i = 0; i < n; ++i)
               if (SET_CONTAINS_FAST(s, i)) {
                  val[nz] = 1.0;
                  ind[nz] = i;
                  nz++;
               }
            POST("add lazy", GRBcblazy(cbdata, nz, ind, val, GRB_LESS_EQUAL, alpha));

            free(ind);
            free(val);
         }
         if (s != NULL)
            set_free(s);
      }
      free(xbar);
      free(V);
   }
   return error;
}

/**
* @brief
*/
int
separatorBnC(graph_t* g, graph_t* h, double* x_bar, GRBmodel* master, double timeout, int iter,
             double alpha, double cut_viol,
             tpoint start_time, const std::string& msg) {
   int found = 1;

   // Auxiliary data structure
   int      n = g->n;
   int      status;
   struct   callback_data mydata;
   // Solve the problem
   double violation = 0.0;

   GRBenv   *env_separa;
   GRBmodel *separa;

   // Arrays for passing data to Gurobi
   char*     le = (char*)malloc(sizeof(char)     * (n + 1));
   int*     ind = (int*)malloc(sizeof(int)       * (n + 1));
   double*  val = (double*)malloc(sizeof(double) * (n + 1));
   double*   ub = (double*)malloc(sizeof(double) * n);
   double*   lb = (double*)malloc(sizeof(double) * n);
   double*   oj = (double*)malloc(sizeof(double) * n);

   // Build pricing subproblem enviroment
   GRBloadenv(&env_separa, NULL);
   GRBsetintparam(env_separa, GRB_INT_PAR_OUTPUTFLAG, 0);
   GRBsetintparam(env_separa, GRB_INT_PAR_METHOD, 1);   // Force dual simplex
   GRBsetintparam(env_separa, GRB_INT_PAR_LAZYCONSTRAINTS, 1);
   GRBsetintparam(env_separa, GRB_INT_PAR_SOLUTIONLIMIT,   1);
   GRBsetdblparam(env_separa, GRB_DBL_PAR_TIMELIMIT, timeout);

   int zeros = 0;
   for (int i = 0; i < n; ++i) {
      le[i] = GRB_BINARY;
      val[i] = 1;
      ind[i] = i;
      lb[i] = 0.0;
      ub[i] = 1.0;
      oj[i] = x_bar[i];
      if (oj[i] < INT_TOL || oj[i] > 1.0 - INT_TOL) {
         ub[i] = 0.0;
         oj[i] = 0.0;
         zeros++;
      }
   }

   GRBnewmodel(env_separa, &separa, "mss_sep", n, oj, lb, ub, le, NULL);
   GRBsetintattr(separa, "ModelSense", GRB_MAXIMIZE);
   GRBupdatemodel(separa);
   // lazy constraint
   mydata.n = n;
   mydata.g = h;   // stable set are computed as clique on the complement graph
   mydata.alpha = alpha;
   mydata.runtime = timeout;
   POST("add callback", GRBsetcallbackfunc(separa, subsetelim, (void *)&mydata));
   GRBupdatemodel(separa);

   // Solve the separation problem
   int solCount = 2;
   double runtime = 0;
   while (true) {
      GRBoptimize(separa);
      POST("get status", GRBgetintattr(separa, "Status", &status));
      if (status == GRB_OPTIMAL)
         break;
      double rnt;
      if (status == GRB_TIME_LIMIT) {
         POST("get time", GRBgetdblattr(separa, "Runtime", &rnt));
         runtime += rnt;
         if (runtime >= timeout)
            break;
         mydata.runtime = timeout - runtime;
         POST("set time", GRBsetdblattr(separa, "TimeLimit", timeout - runtime));
      } else {
         if (status == GRB_SOLUTION_LIMIT) {
            POST("get obj p", GRBgetdblattr(separa, "ObjVal", &violation));
            if (violation - alpha > cut_viol) {
               // SET THE NODE LIMIT 5000 FROM NOW
               double nl = 0;
               POST("get obj p", GRBgetdblattr(separa, "NodeCount", &nl));
               POST("set nodelimit", GRBsetdblparam(GRBgetenv(separa), GRB_DBL_PAR_NODELIMIT, nl+5000));
            }
            POST("set timelimit", GRBsetintparam(GRBgetenv(separa), GRB_INT_PAR_SOLUTIONLIMIT, solCount++));
         } else {
            if (status == GRB_NODE_LIMIT) {
               break;
            }
            fprintf(stdout, "FATAL ERROR - STATUS %d\n", status);
            exit(0);
         }
      }
      runtime += rnt;
      POST("get time", GRBgetdblattr(separa, "Runtime", &rnt));
      mydata.runtime = timeout - runtime;
	fprintf(stdout, "%f %f\n", timeout, runtime);
	fflush(stdout);
      //POST("set time", GRBsetdblattr(separa, "TimeLimit", timeout - runtime));
   }
   solCount = 0;

   if (status == GRB_OPTIMAL || status == GRB_SOLUTION_LIMIT
         || status == GRB_TIME_LIMIT || status == GRB_NODE_LIMIT) {
      POST("solcount", GRBgetintattr(separa, "SolCount", &solCount));
      if (solCount == 0) {
         found = 0;
         fprintf(stdout, "%s Violation: %.3f |U| = %d alpha=%.0f zeros: %d cut_viol: %.3f ",
                 msg.c_str(), 0.0, 0, alpha, n - zeros, cut_viol);
         goto QUIT;
      }

      POST("get obj p", GRBgetdblattr(separa, "ObjVal", &violation));
      if (violation - alpha < CUT_VIOL) {
         found = 0;
         fprintf(stdout, "%s Violation: %.3f |U| = %d alpha=%.0f zeros: %d cut_viol: %.3f ",
                 msg.c_str(), 0.0, 0, alpha, n - zeros, cut_viol);
         goto QUIT;
      }
      if (violation - alpha >= CUT_VIOL && violation - alpha <= cut_viol)
         found = 2;

      // Add the new cut
      int nz = 0;
      double u_i = 0.0;
      for (int i = 0; i < n; ++i) {
         POST("get X", GRBgetdblattrelement(separa, "X", i, &u_i));
         if (u_i > 0.99) {  // is u_i == 1.0 ?
            val[nz] = 1.0;
            ind[nz] = i;
            nz++;
         }
      }

      POST("add cut", GRBaddconstr(master, nz, ind, val, GRB_LESS_EQUAL, alpha, (const char*)NULL));
      fprintf(stdout, "%s Violation: %.3f |U| = %d alpha=%.0f zeros: %d cut_viol: %.3f ",
              msg.c_str(), violation, nz, alpha, n-zeros, cut_viol);
   } else {
      found = 0;
      fprintf(stdout, "%s Violation: %.3f |U| = %d alpha=%.0f zeros: %d cut_viol: %.3f ",
              msg.c_str(), 0.0, 0, alpha, n-zeros, cut_viol);
   }

QUIT:
   free(ind);
   free(val);
   free(oj);
   free(le);
   free(ub);
   free(lb);

   POST("free", GRBfreemodel(separa));
   GRBfreeenv(env_separa);

   auto end = std::chrono::steady_clock::now();
   auto elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start_time).count()) / 1000;
   printf("alpha=%.0f Time %.3f\n", alpha, elapsed);

   return found;
}
