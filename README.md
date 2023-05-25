# dwarf-p-cloudsc2-tl-ad

`dwarf-p-cloudsc2-tl-ad` is intended to test tangent-linear and adjoint
version of the CLOUDSC2 cloud microphysics scheme of the IFS.

*This package is made available to support research collaborations and is not
officially supported by ECMWF*

## Contact

Michael Lange (michael.lange@ecmwf.int),
Willem Deconinck (willem.deconinck@ecmwf.int)
Balthasar Reuter (balthasar.reuter@ecmwf.int),

## Licence

`dwarf-p-cloudsc2-tl-ad` is distributed under the Apache Licence Version 2.0.
See [LICENSE](LICENSE) file for details.

## Prototypes available

- **dwarf-cloudsc2-nl**: The nonlinear forward run only. This has been created
  from the original CLOUDSC dwarf. _The reference results have not been updated,
  so errors are expected._
- **dwarf-cloudsc2-tl**: Tangent linear version of CLOUDSC2 that performs a Taylor
  test to validate the TL code.
- **dwarf-cloudsc2-ad**: Adjoint test of CLOUDSC2 that validates adjoint symmetry.
- **dwarf-cloudsc2-nl-loki**: Experimental version of Loki port of CLOUDSC2 NL

## Download and Installation

The preferred method to install the CLOUDSC dwarf uses the bundle
definition shipped in the main repository. For this please
install the bundle via:
```
./cloudsc-bundle create  # Checks out dependency packages
./cloudsc-bundle build [--build-type=debug|bit|release] [--arch=$PWD/arch/ecmwf/machine/compiler/version/env.sh]
```
## Loki variant build
At the moment, Loki variant cannot be compiled with GNU 11 compiler. On ATOS HPC, Loki variant may be built using
```
./cloudsc-bundle build --clean --with-loki --loki-frontend=fp --arch=./arch/ecmwf/hpc2020/gnu/9.3.0/
./cloudsc-bundle build --clean --with-loki --loki-frontend=fp --arch=./arch/ecmwf/hpc2020/intel/2021.4.0
```
Targetting GPU, one should use:
```
 ./cloudsc-bundle build --with-gpu --with-loki --loki-frontend=fp --arch=./arch/ecmwf/hpc2020/nvhpc/22.1 
```
Running on a GPU node requires specificaton of the CUDA heapsize, e.g. 
```
# Go into build directory and enable environment
cd build && . env.sh

# Run Loki-SCC with vector-level kernels
NV_ACC_CUDA_HEAPSIZE=9G ./bin/dwarf-cloudsc2-nl-loki-scc 1 256000 128

# Run Loki-SCC-hoist with hoisted temporaries
./bin/dwarf-cloudsc2-nl-loki-scc-hoist 1 256000 128
```
NVHPC will not build Loki variant with --build-type=debug due to the compiler error. Possibly, look into the compiler flags is needed here.

## Example usage and verification

Following the build, please run the following to set up the environment:
```
cd build
source env.sh
```

Example run of the non-linear baseline run on CPU (4 threads):
```
./bin/dwarf-cloudsc2-nl 4 160000 32
```

Verify the correctness of the tangent linear version by running a
Taylor test for each column with:
```
./bin/dwarf-cloudsc2-tl 1 100 1
```
_Note that this is not yet ready for performance evaluation._

Verify the correctness of the adjoint, please run:
```
./bin/dwarf-cloudsc2-ad 1 100 100
```
_Note that this is not yet ready for performance evaluation._

## Performance
ATOS node was allocated with:
```sh
   export OMP_NUM_THREADS=64
   OMP_PLACES="{$(seq -s '},{' 0 $(($OMP_NUM_THREADS-1)) )}" srun -q np --ntasks=1 --hint=nomultithread --cpus-per-task=$OMP_NUM_THREADS --pty /bin/bash
```
Command:
```./bin/dwarf-cloudsc2-nl $OMP_NUM_THREADS 163840 32```
Current performance on gcc/11.2.0:  76 MFlops/s
Current performance on intel/2021.4.0: 94k MFlops/s
 
