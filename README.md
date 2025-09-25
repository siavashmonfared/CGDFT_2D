![screenshot](2D_condensation.png)

# CGDFT_2D: 2D Coarse-Grained Lattice Gas Density Functional Theory

[![GitHub](https://img.shields.io/badge/GitHub-CGDFT__2D-blue?logo=github)](https://github.com/siavashmonfared/CGDFT_2D)
[![C++](https://img.shields.io/badge/C%2B%2B-11-blue?logo=cplusplus)](https://en.cppreference.com/w/cpp/11)
[![MPI](https://img.shields.io/badge/MPI-Parallelized-green?logo=openmpi)](https://www.open-mpi.org/)
[![Python](https://img.shields.io/badge/Python-Visualization-orange?logo=python)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A 2D implementation of Coarse-Grained Lattice Gas Density Functional Theory with MPI parallelization. To run, ensure that proper number of processors are assigned in the input file. 

## Building 

```
mpic++ -std=c++0x code.cpp -o code
mpiexec -np NumProc code 
```


## Visualizing 

To visualize, modify the exec.sh file included based on the number of time steps and frequency of data output. Make sure the dimensions are correct in the python file. Then: 

```
./exec.sh 

```

