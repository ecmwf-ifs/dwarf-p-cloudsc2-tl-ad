Dwarf-P-cloudsc2-tl-ad
-----------------------------------
Contact: Michael Lange (michael.lange@ecmwf.int),
Filip Vana (filip.vana@ecmwf.int),
Willem Deconinck (willem.deconinck@ecmwf.int)

Dwarf-P-cloudsc2-tl-ad is intended to test tangent-linear and adjoint
version of the CLOUDSC2 cloud microphysics scheme.

Prototypes available
--------------------
- **dwarf-cloudsc2-nl**: The nonlinear forward run only. This has been created
  from the original CLOUDSC dwarf. _The reference results have not been updated,
  so errors are expected._
- **dwarf-cloudsc2-tl**: Tangent linear version of CLOUDSC2 that performs a Taylor
  test to validate the TL code.
- **dwarf-cloudsc2-ad**: Adjoint test of CLOUDSC2 that validates adjoint symmetry.

Download and Installation
-------------------------
The preferred method to install the CLOUDSC dwarf uses the bundle
definition shipped in the main repository. For this please
```
git clone ssh://git@git.ecmwf.int/escape/dwarf-p-cloudsc2-tl-ad.git cloudsc2-tl-ad
cd cloudsc2-tl-ad
```
Then simply install the bundle via:
```
./cloudsc-bundle create  # Checks out dependency packages
./cloudsc-bundle build [--build-type=debug|bit|release] [--arch=$PWD/arch/ecmwf/machine/compiler/version/env.sh]
```
