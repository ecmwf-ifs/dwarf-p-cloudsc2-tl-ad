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


    loki_transform_convert(
        MODE idem FRONTEND ${LOKI_FRONTEND}
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS ${COMMON_MODULE}/yomphyder.F90
        INCLUDES ${COMMON_INCLUDE}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        OUTPUT
               loki-idem/satur.idem.F90
               loki-idem/cloudsc2ad.idem.F90
               loki-idem/cloudsc2tl.idem.F90
               loki-idem/cloudsc_driver_ad_loki_mod.idem.F90
               DEPENDS satur.F90 cloudsc2ad.F90 cloudsc2tl.F90 cloudsc_driver_ad_loki_mod.F90
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
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
endif()
