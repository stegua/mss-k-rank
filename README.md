# Maximum Stable Set by k-Rank Cut Inequalities
This repository contains the source code used for the computational experiments in the following paper:

* S. Congiglio, S. Gualandi. *Optimizing over the Closure of Rank Inequalities with a Small Right-Hand Side for the Maximum Stable Set problem via Bilevel Programming*. Accepted for publication at **INFORMS Journal of Computing** (to appear).

## Requiriments and Installation
The source code relies on the three following external libraries:

1. [Gurobi](https://www.gurobi.com/), a very fast commercial mathematical optimization solver that can be downloaded and used with an academic free license.
2. [Cliquer 1.21](https://users.aalto.fi/~pat/cliquer.html), a maximum weight clique solver, developed S. Niskanen and Patric R. J. Östergård. We have included a slightly modified version of *cliquer-1.21* under the directory `externs`.
3. [OpenMP](https://www.openmp.org/): this is optional, but if you have it installed, the code should run faster. If you don't have OpenMP on your computer, you should be able to use the code anyway.

### Compiling under Windows 10
We have developed the code under Windows 10, using MS Visual Studio Community 2019. 
The MSVC project files are available under the directory `msvc`.
Please, in the project file, check the file path where you have installed Gurobi on your computer.
The code compiles in *Release* mode at 64 bits.

### Compiling under Linux
We run all our tests on a Linux server running CentOS.
The repository contains a `Makefile`, where you have to set the path for the Gurobi `include` and `lib` directories correctly.
Note that we are compiling using the C++17 standard (`-std=c++17`), we are setting some optimization flags (`-march=native -mavx2 -mfma`), and we use the OpenMP library (`-fopenmp`).

## Separation routines
You may be interested in the following files, which contains the cut separator we have implemented in our solver:

1. `Find_n_rank.h`: they contain the exact separation algorithm for *k-rank inequalities*, with n=2,3,4,5.
2. `DigraphSpp.h`: it contains the exact separation algorithms (nested loops) for small cliques (with cardinalities smaller or equal than 4);
the separation algorithm *OddCycle* and *OddWheel* rely on a custom implementation of our shortest path algorithm using binary heaps.

The implementation of the separation algorithms is self-contained (it does not depend on external libraries), and it would be possible to reuse the code in other contexts.


### Datasets
In the experiments, we use the following set of graphs, which are present in the subfolder `data`:

1. `data/dimacs`: A subset of the graphs from the DIMACS challenge on the maximum clique problem, which are complemented to obtain maximum stable set instances. We include only graphs of moderate size.
2. `data/sparse`: Very sparse graph with up to 400 vertices.
3. `data/rnd`: Erdos/Reny random graphs with 50, 75, and 100 vertices and an edge density ranging from 10% up to 90%.

All the graphs are stored in the ascii text [DIMACS format](http://www.diag.uniroma1.it/challenge9/format.shtml).

## Reproducing the results
In order to facilitate the reproduction of our results, we include three python script that we use to automatize our tests:

1. `run_test_dimacs.py`
2. `run_test_sparse.py`
3. `run_test_ijoc.py`

Before running the scripts, you should check the variables pointing to the directory where the instances are saved.

Please, contact us by email if you encounter any issues.

### Citing the paper
If you ever find our paper useful for your research, or if you reuse any snippet from our code, please cite our work as follows:

```
@article{SCSG2021,
  author  = {Stefano Coniglio and Stefano Gualandi}, 
  title   = {Optimizing over the Closure of Rank Inequalities with a Small Right-Hand Side for the Maximum Stable Set problem via Bilevel Programming},
  journal = {{INFORMS} Journal on Computing},
  year    = 2021,
  note    = {Accepted for publication, to apper}
}
```