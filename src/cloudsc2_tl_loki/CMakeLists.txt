# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC2_TL_LOKI DEFAULT ON
	DESCRIPTION "Build the Fortran CLOUDSC2 tangent-linear dwarf Loki"
)

if( HAVE_CLOUDSC2_TL_LOKI )

    set( COMMON_MODULE "${CMAKE_CURRENT_SOURCE_DIR}/../common/module" )
    set( COMMON_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../common/include" )

#   # Define the binary build target for this variant

    set( LOKI_FRONTEND "fp" CACHE STRING "Frontend parser for Loki transforms" )

    # OFP frontend cannot deal with statement functions, so we toggle them here
    set( CLOUDSC_DEFINE_STMT_FUNC "" )

    ####################################################
    ##  Idempotence mode:                             ##
    ##   * Internal "do-nothing" mode for Loki debug  ##
    ####################################################


    loki_transform_convert(
        MODE idem FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS ${COMMON_MODULE}/yomphyder.F90
        INCLUDES ${COMMON_INCLUDE}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
	OUTPUT 
	       loki-idem/satur.idem.F90 
               loki-idem/error.idem.F90
	       loki-idem/cloudsc2.idem.F90 
	       loki-idem/cloudsc2tl.idem.F90 
	       loki-idem/cloudsc_driver_tl_loki_mod.idem.F90
	       DEPENDS satur.F90 error.F90 cloudsc2.F90 cloudsc2tl.F90 cloudsc_driver_tl_loki_mod.F90 
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc2-tl-loki-idem
        SOURCES
            dwarf_cloudsc.F90
            loki-idem/cloudsc_driver_tl_loki_mod.idem.F90
            loki-idem/cloudsc2tl.idem.F90 
            loki-idem/cloudsc2.idem.F90
            loki-idem/satur.idem.F90
            loki-idem/error.idem.F90
        LIBS
            cloudsc2-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc2-tl-loki-idem-serial
        COMMAND bin/dwarf-cloudsc2-tl-loki-idem
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

#    ####################################################
#    ##  "Single Column Coalesced" (SCC) mode          ##
#    ##   * Removes horizontal vector loops            ##
#    ##   * Invokes compute kernel as `!$acc vector`   ##
#    ####################################################
#
#
#    loki_transform_convert(
#        MODE scc FRONTEND ${LOKI_FRONTEND} CPP
#        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
#        PATH ${CMAKE_CURRENT_SOURCE_DIR}
#        HEADERS ${COMMON_MODULE}/yomphyder.F90 ${COMMON_MODULE}/yoecldp.F90
#        INCLUDES ${COMMON_INCLUDE}
#        DEFINITIONS CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
#        DATA_OFFLOAD REMOVE_OPENMP
#        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc
#        OUTPUT  loki-scc/satur.scc.F90 loki-scc/cloudsc2.scc.F90 loki-scc/cloudsc2tl.scc.F90 loki-scc/cloudsc_driver_tl_loki_mod.scc.F90
#        DEPENDS satur.F90 cloudsc2.F90 cloudsc2tl.F90 cloudsc_driver_tl_loki_mod.F90
#    )
#
#    ecbuild_add_executable( TARGET dwarf-cloudsc2-tl-loki-scc
#        SOURCES
#            dwarf_cloudsc.F90
#            loki-scc/cloudsc_driver_tl_loki_mod.scc.F90
#            loki-scc/cloudsc2tl.scc.F90
#            loki-scc/cloudsc2.scc.F90
#            loki-scc/satur.scc.F90
##           loki-scc/error.scc.F90
#        LIBS
#            cloudsc2-common-lib
#        DEFINITIONS ${CLOUDSC_DEFINITIONS}
#    )
#
#    ecbuild_add_test(
#        TARGET dwarf-cloudsc2-tl--loki-scc-serial
#        COMMAND bin/dwarf-cloudsc2-tl-loki-scc
#        ARGS 1 1280 128
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
#        OMP 1
#        ENVIRONMENT NV_ACC_CUDA_HEAPSIZE=19G
#    ) 
#
#    ####################################################
#    ##  SCC-hoist mode                                ##
#    ##   * SCC with vector loop hoisted               ##
#    ##   * Kernel is "seq, but args are full blocks   ##
#    ##   * Temporary arrays hoisted to driver         ##
#    ####################################################
#
#
#    loki_transform_convert(
#        MODE scc-hoist FRONTEND ${LOKI_FRONTEND} CPP
#        DIRECTIVE none
#        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
#        PATH ${CMAKE_CURRENT_SOURCE_DIR}
#        HEADERS ${COMMON_MODULE}/yomphyder.F90 ${COMMON_MODULE}/yoecldp.F90
#        INCLUDES ${COMMON_INCLUDE}
#        DEFINITIONS CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
#        DATA_OFFLOAD REMOVE_OPENMP
#        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-hoist
#        OUTPUT  loki-scc-hoist/satur.scc_hoist.F90 
#                loki-scc-hoist/cloudsc2.scc_hoist.F90
#                loki-scc-hoist/cloudsc2tl.scc_hoist.F90
#                loki-scc-hoist/cloudsc_driver_tl_loki_mod.scc_hoist.F90
#        DEPENDS satur.F90 cloudsc2.F90 cloudsc2tl.F90 cloudsc_driver_tl_loki_mod.F90
#    )
#
#
#    ecbuild_add_executable( TARGET dwarf-cloudsc2-tl-loki-scc-hoist
#        SOURCES
#            dwarf_cloudsc.F90
#            loki-scc-hoist/cloudsc_driver_tl_loki_mod.scc_hoist.F90
#            loki-scc-hoist/cloudsc2tl.scc_hoist.F90
#            loki-scc-hoist/cloudsc2.scc_hoist.F90
#            loki-scc-hoist/satur.scc_hoist.F90
##           loki-scc-hoist/error.scc_hoist.F90
#        LIBS
#            cloudsc2-common-lib
#        DEFINITIONS ${CLOUDSC_DEFINITIONS}
#    )
#
#    ecbuild_add_test(
#        TARGET dwarf-cloudsc2-tl-loki-scc-hoist-serial
#        COMMAND bin/dwarf-cloudsc2-tl-loki-scc-hoist
#        ARGS 1 1280 128
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
#        OMP 1
#        ENVIRONMENT NV_ACC_CUDA_HEAPSIZE=19G
#    ) 

endif()
