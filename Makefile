#
# @fileoverview Copyright (c) 2019, by Stefano Gualandi, UniPv,
#               via Ferrata, 1, Pavia, Italy, 27100
#
# @author stefano.gualandi@gmail.com (Stefano Gualandi)
#
#

COMPILER    = g++ ${OPTFLAG}
LINKER      = g++ ${LDFLAGS}

# Directory for my files
MYHOME          = ${PWD}
BIN             = ${MYHOME}/bin
INCLUDE         = ${MYHOME}/include/
LIB             = ${MYHOME}/lib
SRC             = ${MYHOME}/src

OPTFLAG = -O3 -ffast-math -march=native -mavx2 -mfma -DNDEBUG -Wall -std=c++17 -fopenmp -DLINUX
LDFLAGS = -O3 -DNDEBUG -lm -pthread -std=c++17 -fopenmp

# Gurobi
GUROBI_INC = -I/cm/shared/apps/gurobi910/linux64/include/
GUROBI_LIB = -L/cm/shared/apps/gurobi910/linux64/lib/ -lgurobi91

# Cliquer
CLIQUER_INC = -I${MYHOME}/externs/cliquer-1.21/
CLIQUER_LIB = ${MYHOME}/externs/cliquer-1.21/cliquer.o ${MYHOME}/externs/cliquer-1.21/graph.o ${MYHOME}/externs/cliquer-1.21/reorder.o

# Directory for output files
OUT_DIR=bin lib

bitgraph: ${OUT_DIR} ${SRC}/bitGraph.cpp
	${COMPILER} -c -g -pg ${SRC}/bitGraph.cpp -o ${LIB}/bitGraph.o -I${INCLUDE} ${CLIQUER_INC} ${GUROBI_INC}

separatorBnB: ${OUT_DIR} ${SRC}/separatorBnB.cpp
	${COMPILER} -c -g -pg ${SRC}/separatorBnB.cpp -o ${LIB}/separatorBnB.o -I${INCLUDE} ${GUROBI_INC} ${CLIQUER_INC}

separatorBnC: ${OUT_DIR} ${SRC}/separatorBnC.cpp
	${COMPILER} -c -g -pg ${SRC}/separatorBnC.cpp -o ${LIB}/separatorBnC.o -I${INCLUDE} ${GUROBI_INC} ${CLIQUER_INC}

mss: ${OUT_DIR} ${SRC}/mssRank.cpp
	${COMPILER} -c -g -pg ${SRC}/mssRank.cpp -o ${LIB}/mssRank.o -I${INCLUDE} ${GUROBI_INC} ${CLIQUER_INC}
	${LINKER} -o ${BIN}/mssRank ${CLIQUER_LIB} ${LIB}/bitGraph.o ${LIB}/separatorBnB.o ${LIB}/separatorBnC.o ${LIB}/mssRank.o ${GUROBI_LIB}

# Create subdirectory for output files (bin,lib)
MKDIR_P = mkdir -p

lib:
	${MKDIR_P} lib

bin:
	${MKDIR_P} bin

# Be careful!
clean::
	rm -f *.o
	rm -f ${LIB}/*.o
	rm -f *~
	rm -f ${SRC}/*~ ${INCLUDE}/*~
