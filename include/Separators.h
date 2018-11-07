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

#include "BitGraph.h"

// Use Cliquer to find all maximal clique in the conflict graph
extern "C" {
#include "graph.h"
#include "set.h"
#include "graph.h"
}

/*
* @brief Solve separation of rank ineq. by branch-and-bound
*/
int separatorBnB(BitGraph * G, double * x_bar, GRBmodel * master,
                 double alpha, double cut_viol,
                 tpoint start_time, const std::string& msg);

/*
* @brief Solve separation of rank ineq. by branch-and-cut
*/
int separatorBnC(graph_t* g, graph_t* h, double* x_bar, GRBmodel* master, double timeout, int iter,
                 double alpha, double cut_viol,
                 tpoint start_time, const std::string& msg);

int solveBranchAndCut(BitGraph & G, double alpha);