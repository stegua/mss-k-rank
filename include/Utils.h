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

#include <chrono>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <random>
#include <future>

const double INT_TOL = 0.0001;

const uint64_t NODE_LIMIT = std::numeric_limits<uint64_t>::max();

const double CUT_VIOL = 0.01;
const double TIMEOUT = 7200;
const int    PERTURB = 2;        // {0: none, 1: -eps, 2: +eps, 3: coord}


extern "C" {
#include "gurobi_c.h"
#include "math.h"
}

#define POST(s, x) if (x) { fprintf(stdout,"%s: GRB error %d\n", s, x); exit(0); }

typedef std::chrono::steady_clock::time_point tpoint;

#include <cassert>
#include <vector>
using std::vector;

typedef uint16_t          vertex_t;
typedef vector<vertex_t>  vertex_set_t;
