Current implementation provides working idem and scc variants, the latter with 3.5 GFlops/s. scc-hoist currently requires manual edit of 
cloudsc_driver_loki_mod.scc_hoist.F90 and remove empty entry enter data create and exit data delete. It also get only 37 Mflops/s. 
Workaround for *scc* variants: static allocation of the CETA array, otherwise runtime errors.

