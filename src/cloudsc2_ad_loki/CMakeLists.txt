# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC2_AD_LOKI DEFAULT OFF
    DESCRIPTION "Build the Fortran CLOUDSC2 adjoint dwarf LOKI"
)

if( HAVE_CLOUDSC2_AD_LOKI )

    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )

    set( COMMON_MODULE "${CMAKE_CURRENT_SOURCE_DIR}/../common/module" )
    set( COMMON_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../common/include" )

    set( LOKI_FRONTEND "fp" CACHE STRING "Frontend parser for Loki transforms" )


    ####################################################
    ##  Idempotence mode:                             ##
    ##   * Internal "do-nothing" mode for Loki debug  ##
    ####################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-idem/satur.idem.F90
            loki-idem/cloudsc2ad.idem.F90
            loki-idem/cloudsc2tl.idem.F90
            loki-idem/cloudsc_driver_ad_loki_mod.idem.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        DEPENDS
            satur.F90
            cloudsc2ad.F90
            cloudsc2tl.F90
            cloudsc_driver_ad_loki_mod.F90
        MODE idem
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc2-ad-loki-idem
        SOURCES
            dwarf_cloudsc.F90
            loki-idem/cloudsc_driver_ad_loki_mod.idem.F90
            loki-idem/cloudsc2tl.idem.F90
            loki-idem/cloudsc2ad.idem.F90
            loki-idem/satur.idem.F90
        LIBS
            cloudsc2-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc2-ad-loki-idem-serial
        COMMAND bin/dwarf-cloudsc2-ad-loki-idem
        ARGS
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

     ####################################################
     ##  "Single Column Coalesced" (SCC) mode          ##
     ##   * Removes horizontal vector loops            ##
     ##   * Invokes compute kernel as `!$acc vector`   ##
     ####################################################

     loki_transform(
        COMMAND convert
        OUTPUT
            loki-scc/satur.scc.F90
            loki-scc/cloudsc2ad.scc.F90
            loki-scc/cloudsc2tl.scc.F90
            loki-scc/cloudsc_driver_ad_loki_mod.scc.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc
        DEPENDS
            satur.F90
            cloudsc2ad.F90
            cloudsc2tl.F90
            cloudsc_driver_ad_loki_mod.F90
        MODE scc
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        DATA_OFFLOAD
        REMOVE_OPENMP
        CPP
        DEFINITIONS
            ${CLOUDSC_DEFINITIONS}
            CLOUDSC_GPU_TIMING
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
    )

     ecbuild_add_executable( TARGET dwarf-cloudsc2-ad-loki-scc
         SOURCES
             dwarf_cloudsc.F90
             loki-scc/cloudsc_driver_ad_loki_mod.scc.F90
             loki-scc/cloudsc2tl.scc.F90
             loki-scc/cloudsc2ad.scc.F90
             loki-scc/satur.scc.F90
         LIBS
             cloudsc2-common-lib
         DEFINITIONS ${CLOUDSC_DEFINITIONS}
     )

     ecbuild_add_test(
         TARGET dwarf-cloudsc2-ad--loki-scc-serial
         COMMAND bin/dwarf-cloudsc2-ad-loki-scc
         ARGS
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
         OMP 1
         ENVIRONMENT NV_ACC_CUDA_HEAPSIZE=9G
     )
    ####################################################
    ##  SCC-hoist mode                                ##
    ##   * SCC with vector loop hoisted               ##
    ##   * Kernel is "seq, but args are full blocks   ##
    ##   * Temporary arrays hoisted to driver         ##
    ####################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-scc-hoist/satur.scc_hoist.F90
            loki-scc-hoist/cloudsc2ad.scc_hoist.F90
            loki-scc-hoist/cloudsc2tl.scc_hoist.F90
            loki-scc-hoist/cloudsc_driver_ad_loki_mod.scc_hoist.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-hoist
        DEPENDS
            satur.F90
            cloudsc2ad.F90
            cloudsc2tl.F90
            cloudsc_driver_ad_loki_mod.F90
        MODE scc-hoist
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        DATA_OFFLOAD
        REMOVE_OPENMP
        CPP
        DEFINITIONS
            ${CLOUDSC_DEFINITIONS}
            CLOUDSC_GPU_TIMING
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
    )

     ecbuild_add_executable( TARGET dwarf-cloudsc2-ad-loki-scc-hoist
         SOURCES
             dwarf_cloudsc.F90
             loki-scc-hoist/cloudsc_driver_ad_loki_mod.scc_hoist.F90
             loki-scc-hoist/cloudsc2tl.scc_hoist.F90
             loki-scc-hoist/cloudsc2ad.scc_hoist.F90
             loki-scc-hoist/satur.scc_hoist.F90
         LIBS
             cloudsc2-common-lib
         DEFINITIONS ${CLOUDSC_DEFINITIONS}
     )

     ecbuild_add_test(
         TARGET dwarf-cloudsc2-ad--loki-scc-hoist-serial
         COMMAND bin/dwarf-cloudsc2-ad-loki-scc-hoist
         ARGS
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
         OMP 1
         ENVIRONMENT NV_ACC_CUDA_HEAPSIZE=9G
     )

endif()
