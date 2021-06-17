# Maximum Stable Set by k-Rank Cut Inequalities
This repositories contains the source code used for the computational experiments in the following paper:

* S. Congiglio, S. Gualandi. *Optimizing over the Closure of Rank Inequalities with a Small Right-Hand Side for the Maximum Stable Set problem via Bilevel Programming*. Accepted for publication at INFORMS Journal of Computing (to appear).

## Requiriments
The source code relies on the three following external libraries:

1. [Gurobi](https://www.gurobi.com/), a very fast commercial mathematical optimization solver that can be downloaded and used with an academic free license.
2. [Cliquer 1.21](https://users.aalto.fi/~pat/cliquer.html), a maximum weight clique solver, developed S. Niskanen and Patric R. J. Östergård. We have included a slightly modified version of cliquer under the directory `externs`.
3. [OpenMP](https://www.openmp.org/): this is optional, but if you have it installed, the code should run faster. If you don't have openmp on your library, you should be able to use the code anyway.

### Windows
We have developed the code under Windows 10, using MS Visual Studio Community 2019. 
The MSVC project files are available under the directory `msvc`.
Please, in the project file check the file path where you have installed Gurobi on your computer.
The code compile in *Release* mode at 64 bits.

### Linux
We run all our test on a linux server running CentOS.
The repository contains a `Makefile`, where you have to set correctly the path for the Gurobi `include` and `lib` directories.
Note that we are compiling using the C++17 standard (`-std=c++17`), we are setting some optimization flags (`-march=native -mavx2 -mfma`), and we use the OpenMP library (`-fopenmp`).

## Datasets
In the experiments, we use the following set of graphs, which are present in the subfolder `data`:

1. A subset of the graphs from the DIMACS challenge on the maximum clique problem, which are complemented to obtain maximum stable set instances. We include only graph of moderate size.
2. Very sparse graph with up to 400 vertices.
3. Erdos/Reny random graphs with 50, 75 and 100 vertices, and an edge density ranging from 10% up to 90%.

All the graphs are stored in the ascii text DIMACS format.

## Separation routines
You may be interested in the following files, which containt the cut separator we have implemented in our solver:

1. `Find_n_rank.h`: they contain the exact separation algorithm for *k-rank inequalities*, with n=2,3,4,5.
2. `DigraphSpp.h`: it contains the exact separation algorithms (nested loops) for small cliques (with cardinalities smaller or equal than 4);
the separation algorithm *OddCycle* and *OddWheel*, which are based on a custom implementation of a shortest path algorithm using binary heaps.

The implementation of the separation algorithms is self-contained (it does not depend on external libraries) and it would be possible to reuse the code in other context.

## Citing the paper
If you ever find our paper useful for your research, or if you reuse any snippet from our code, please cite our work as follows:

```
@article{SCSG2021,
  author  = {Stefano Conilgio and Stefano Gualandi}, 
  title   = {Optimizing over the Closure of Rank Inequalities with a Small Right-Hand Side for the Maximum Stable Set problem via Bilevel Programming},
  journal = {{INFORMS} Journal on Computing},
  year    = 2021,
  note    = {Accepted for publication, to apper}
}

```