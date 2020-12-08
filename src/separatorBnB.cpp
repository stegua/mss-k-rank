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

#include "Find_1_rank.h"
#include "Find_2_rank.h"
#include "Find_3_rank.h"
#include "Find_4_rank.h"
#include "Find_5_rank.h"

/**
* @brief Separate rank ineq. by branch-and-bound
*/
int separatorBnB(BitGraph * G, double * x_bar, GRBmodel * master,
                 double alpha, double cut_viol, tpoint start_time,
                 const std::string& msg) {
   int n = G->num_vertices();
   int zeros = 0;
   for (int i = 0; i < n; ++i) {
      if (x_bar[i] < INT_TOL || x_bar[i] > 1.0 - INT_TOL) {
         G->set_node_weight(i, 0.0);
         zeros++;
      } else
         G->set_node_weight(i, x_bar[i]);
   }

   typedef std::pair<double, vertex_set_t> mypair;

   mypair result;

   if (alpha == 1)
      result = best_1_rank(G, cut_viol, start_time);
   if (alpha == 2)
      result = best_2_rank(G, cut_viol, start_time);
   if (alpha == 3)
      result = best_3_rank(G, cut_viol, start_time);
   if (alpha == 4)
      result = best_4_rank(G, cut_viol, start_time);
   if (alpha == 5)
      result = best_5_rank(G, cut_viol, start_time);

   int ret_code = 0;
   if (result.first >= alpha + CUT_VIOL) {
      fprintf(stdout, "%s Violation: %.3f |U| = %d alpha=%.0f zeros: %d cut_viol: %.3f ",
              msg.c_str(), result.first, (int)result.second.size(), alpha, n-zeros, cut_viol);
      auto end = std::chrono::steady_clock::now();
      auto elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start_time).count()) / 1000;
      printf("Time %.3f\n", elapsed);

      int* ind = (int*)malloc(n * sizeof(int));
      double* val = (double*)malloc(n * sizeof(double));
      int nz = 0;
      for (int i : result.second) {
         val[nz] = 1.0;
         ind[nz] = i;
         nz++;
      }
      POST("add cut", GRBaddconstr(master, nz, ind, val, GRB_LESS_EQUAL, alpha, (const char*)NULL));

      free(ind);
      free(val);
      if (result.first <= alpha + cut_viol)
         ret_code = 2;
      else
         ret_code = 1;
   } //else
     // fprintf(stdout, "%s Violation: %.3f |U| = %d alpha=%.0f zeros: %d cut_viol: %.3f ",
     //         msg.c_str(), 0.0, 0, alpha, n-zeros, cut_viol);


   return ret_code;
}



/* Define structure to pass data to the callback function */
struct callback_data_cg {
   int n;
   BitGraph* G;

   double alpha;
   tpoint start_time;
};

/* Rank separation subroutine */
int __stdcall
rank_separation(GRBmodel *model,
                void     *cbdata,
                int       where,
                void     *usrdata) {
   int error = 0;

   if (where == GRB_CB_MIPNODE) {
      struct callback_data_cg *mydata = (struct callback_data_cg *) usrdata;
      int n = mydata->n;
      double *x_bar = NULL;

      double nodecount;
      GRBcbget(cbdata, where, GRB_CB_MIPNODE_NODCNT, (void *)&nodecount);
      if (nodecount < 100) {
         double cut_viol = CUT_VIOL;
         if (nodecount < 2)
            cut_viol = 0.1;

         x_bar = (double *)malloc(n * sizeof(double));
         if (x_bar == NULL) {
            fprintf(stderr, "Out of memory\n");
            exit(1);
         }

         double tmp_lp_obj = 0;
         GRBcbget(cbdata, where, GRB_CB_MIPNODE_OBJBND, &tmp_lp_obj);
         GRBcbget(cbdata, where, GRB_CB_MIPNODE_REL, x_bar);
         int zeros = 0;
         for (int i = 0; i < n; ++i) {
            if (x_bar[i] < INT_TOL || x_bar[i] > 1.0 - INT_TOL) {
               mydata->G->set_node_weight(i, 0.0);
               zeros++;
            } else
               mydata->G->set_node_weight(i, x_bar[i]);
         }

         std::pair<double, vertex_set_t> best_rank;

         double alpha = mydata->alpha;
         auto start_time = std::chrono::steady_clock::now();
         if (alpha == 1)
            best_rank = best_1_rank(mydata->G, cut_viol, start_time);
         if (alpha == 2)
            best_rank = best_2_rank(mydata->G, cut_viol, start_time);
         if (alpha == 3)
            best_rank = best_3_rank(mydata->G, cut_viol, start_time);
         if (alpha == 4)
            best_rank = best_4_rank(mydata->G, cut_viol, start_time);
         if (alpha == 5)
            best_rank = best_5_rank(mydata->G, cut_viol, start_time);


         if (best_rank.first > alpha + CUT_VIOL) {
            printf("LB: %.3f Violation: %.3f  |U| = %d  alpha=%.0f  zeros: %d\n",
                   tmp_lp_obj, best_rank.first, (int)best_rank.second.size(), alpha, zeros);

            int* ind = (int*)malloc(n * sizeof(int));
            double* val = (double*)malloc(n * sizeof(double));
            int nz = 0;

            for (int i : best_rank.second) {
               val[nz] = 1.0;
               ind[nz] = i;
               nz++;
            }

            POST("add cut", GRBcbcut(cbdata, nz, ind, val, GRB_LESS_EQUAL, alpha));

            free(ind);
            free(val);
         }

         free(x_bar);
      }
   }

   return error;
}


/*
* @brief Solve max stable set by branch and cut using our separation algorithm
*/
int solveBranchAndCut(BitGraph & G, double alpha) {
   auto start = std::chrono::steady_clock::now();
   double elapsed;

   int n = G.num_vertices();
   int      status;
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
   char	 *vtype = (char*)malloc(sizeof(char)     * n);

   // Build master problem
   GRBloadenv(&env_master, NULL);
   GRBsetintparam(env_master, GRB_INT_PAR_OUTPUTFLAG, 0);
   GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
   GRBsetintparam(env_master, GRB_INT_PAR_PRECRUSH, 1);   // TO BE SET TO 1 SINCE WE USE USER CUTS
   GRBsetintparam(env_master, GRB_INT_PAR_THREADS, 1);
   GRBsetdblparam(env_master, GRB_DBL_PAR_TIMELIMIT, TIMEOUT);
   GRBsetdblparam(env_master, GRB_DBL_PAR_NODELIMIT, 1);

   //GRBsetintparam(env_master, GRB_INT_PAR_CUTS, 0);
   //GRBsetintparam(env_master, GRB_INT_PAR_CLIQUECUTS, 0);

   for (int i = 0; i < n; ++i) {
      ub[i] = 1.0;
      oj[i] = 1.0;
      vtype[i] = GRB_BINARY;
   }
   GRBnewmodel(env_master, &master, "mss", n, oj, NULL, ub, vtype, NULL);

   GRBsetintattr(master, "ModelSense", GRB_MAXIMIZE);

   auto end = std::chrono::steady_clock::now();
   elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

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

   /* Set callback function */
   if (alpha > 0) {
      struct callback_data_cg mydata;
      mydata.n = n;
      mydata.G = &G;
      mydata.alpha = alpha;
      mydata.start_time = std::chrono::steady_clock::now();

      int error = 0;
      error = GRBsetcallbackfunc(master, rank_separation, (void *)&mydata);
      if (error) goto QUIT;
   }

   // Solve the problem
   start = std::chrono::steady_clock::now();
   elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   GRBoptimize(master);

   GRBgetintattr(master, "Status", &status);
   if (status == GRB_UNBOUNDED || status == GRB_INFEASIBLE) {
      fprintf(stdout, "Unbounded or Infeasible\n");
      goto QUIT;
   }
   // Take the current LP decision vector
   POST("get obj", GRBgetdblattr(master, "ObjVal", &alpha_lp));

   end = std::chrono::steady_clock::now();
   elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
   printf("Time %.3f - Value:%.3f\n", elapsed, alpha_lp);

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
