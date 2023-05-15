# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC2_AD_LOKI DEFAULT ON
    DESCRIPTION "Build the Fortran CLOUDSC2 adjoint dwarf LOKI"
)

if( HAVE_CLOUDSC2_AD_LOKI )

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
               loki-idem/cloudsc2ad.idem.F90
               loki-idem/cloudsc2tl.idem.F90
               loki-idem/cloudsc_driver_ad_loki_mod.idem.F90
               DEPENDS satur.F90 cloudsc2.F90 cloudsc2tl.F90 cloudsc_driver_ad_loki_mod.F90
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc2-tl-loki-idem
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
#   # Define the binary build target for this variant
#   ecbuild_add_executable( TARGET dwarf-cloudsc2-ad
#       SOURCES
#           dwarf_cloudsc.F90
#           cloudsc_driver_ad_mod.F90
#           satur.F90
#           cuadjtqs.F90
#           cuadjtqstl.F90
#           cuadjtqsad.F90
#           cloudsc2ad.F90
#           cloudsc2tl.F90
#       DEFINITIONS ${CLOUDSC_DEFINITIONS}
#   )
#   target_link_libraries( dwarf-cloudsc2-ad PRIVATE cloudsc2-common-lib )

#   if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
#       target_link_libraries( dwarf-cloudsc2-ad PRIVATE OpenMP::OpenMP_Fortran )
#   endif()

#   if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
#       target_link_libraries( dwarf-cloudsc2-ad PRIVATE MPI::MPI_Fortran )
#   endif()

#   execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
#     ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
#   execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
#     ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )

#   ecbuild_add_test(
#       TARGET dwarf-cloudsc2-ad
#       COMMAND bin/dwarf-cloudsc2-ad
#       ARGS
#       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
#       OMP 1
#   )
endif()
