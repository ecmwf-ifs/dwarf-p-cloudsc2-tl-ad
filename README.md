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
- **dwarf-cloudsc2-nl-loki-scc**: Experimental version of Loki port of CLOUDSC2 NL with Loki-SCC transformation
- **dwarf-cloudsc2-nl-loki-scc-hoist**: Experimental version of Loki port of CLOUDSC2 NL with Loki-SCC-H transformation

## Download and Installation

The preferred method to install the CLOUDSC2-TL/AD dwarf uses the bundle
definition shipped in the main repository. For this please
install the bundle via:
```
./cloudsc-bundle create  # Checks out dependency packages
./cloudsc-bundle build [--build-type=debug|bit|release] [--arch=$PWD/arch/ecmwf/machine/compiler/version/env.sh]
```

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

## Targetting GPU accelerators via Loki

The CLOUDSC2-TL/AD dwarf currently has limited GPU capabilities that
are restricted to the non-linear variant only, via the Loki
source-to-source transformation package. The respective variants can
be enabled, using the Nvidia compiler toolchain, via
```
./cloudsc-bundle build --with-loki --arch=./arch/ecmwf/hpc2020/nvhpc/22.1
```

When the Loki variant is enabled, two different transformations are run and built:
* **Loki-SCC** is the most conservative supported Loki GPU
  transformation, which maps the outer block loop to SMs via `!$acc
  loop gang` and marks kernel routines as `!$acc routine vector`.
* **Loki-SCC-H**, or Loki-SCC-hoist, is similar to Loki-SCC, but
  hoists temporary array allocations, as well as the innter vector
  loop to the driver layer. These optimisations improve on-device
  memory movement and yield higher overall compute throughput.

Running on a GPU node requires specificaton of the CUDA heapsize for the Loki-SCC variant, e.g.
```
# Go into build directory and enable environment
cd build && . env.sh

# Run Loki-SCC with vector-level kernels
NV_ACC_CUDA_HEAPSIZE=9G ./bin/dwarf-cloudsc2-nl-loki-scc 1 256000 128

# Run Loki-SCC-hoist with hoisted temporaries
./bin/dwarf-cloudsc2-nl-loki-scc-hoist 1 256000 128
```

To debug the Loki-generated code, one may also suppresse OpenACC during the build via:
```
# Using NVHPC-22.11
./cloudsc-bundle build --clean --with-loki --cmake=ENABLE_ACC=off --arch=./arch/ecmwf/hpc2020/nvhpc/22.1/

# Using GNU-11
./cloudsc-bundle build --clean --with-loki --cmake=ENABLE_ACC=off --arch=./arch/ecmwf/hpc2020/gnu/11.2.0/
```

## Performance

_Note that this is not yet ready for performance evaluation._
