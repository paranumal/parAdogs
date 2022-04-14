## ParADoGs
ParADoGs (**Par**allel **A**ccelerated **D**istribution of **G**raph**s**) is an experiemental distributed and fine-grain parallel mesh partitioner. It is designed for, and uses libraries from, [libParanumal](https://github.com/paranumal/libparanumal/), a finite element testbed funded in part by the US Department of Energy as part of the activities of the [Center for Efficient Exscale Discretizations](http://ceed.exascaleproject.org).

---
### 1. Overview

ParADoGs implements some classic graph partitioining heuristics, favoring algorithms that have large amounts of fine-grain parallism. Graphs are generated from the connectivity of an input mesh. ParADoGs can either generate simple box meshes, or use a GMSH `.msh` file.

A. Supported elements:
  - Triangles, quadrilaterals, tetrahedra, hexahedra.

B. Parititoners:
  - Resursive Inhertial Partitioning.
  - Resursive Multilevel Spectral Partitioning.

C. Local Ordering:
  - Cuthill-Mckee.

---
### 2. Dependencies
- Message Passing Interface (MPI v3.0 or higher).
  * The ParADoGs makefiles assume that mpic++ is installed and visible in your path.
- OpenBLAS library.
  * By default, the build system will look for `libopenblas` in your default library search paths.
  * The library paths can also be manually specified in `make.top` with the `LIBP_BLAS_DIR` variable.
  * Some Linux distributions will package the OpenBLAS library. For example, on Ubuntu systems the library can be installed via `sudo apt install libopenblas-dev`
- Open Concurrent Compute Abstraction (OCCA)
  * OCCA is tracked in ParADoGs as a git submodule.
  * OCCA will try to detect if any of these execution models are installed: OpenMP, CUDA, HIP, OpenCL, and/or SYCL.
  * By default, if OCCA does not detect a chosen mode of execution it will default to Serial execution.
  * You will need to adjust the libParanumal setup input files to choose the execution model and compute device appropriate for your system.
  * The OCCA github repo is [here](https://github.com/libocca/occa)
  * The OCCA webpage is [here](http://libocca.org)

Optional:
- Paraview
- Gmsh

---
### 3. Building
```
git clone --recursive https://github.com/paranumal/parAdogs
cd parAdogs
make -j `nproc` 
```

---
### 4. Running Simple Tests
A helper run script for running and immediately viewing the resulting mesh partitioning is provided as
```
runParAdogs.sh <np> <setup.rc>
```

For example, a test which partitions a simple box mesh of Hexahedral elements into 16 paritions can be run as
```
runParAdogs.sh 16 setups/setupHex3D.rc
```

To run directly with MPI, the syntax is
```
mpirun --np <np> paradogsMain <setup.rc>
```

---
### 5. Stress tests
Some larger tests are provided in the `tests/` directory. The `.msh` files must first be generated with Gmsh,
```
gmsh tests/FenceTet3D.geo -3 -format msh22 -o tests/FenceTet3D.msh
```
and the tests can be run as usual,
```
runParAdogs.sh <np> tests/FenceTet3D.rc
```

---

### 6. License

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
