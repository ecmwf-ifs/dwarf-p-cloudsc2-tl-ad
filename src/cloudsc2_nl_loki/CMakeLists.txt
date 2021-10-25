# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC2_NL_LOKI DEFAULT ON
    DESCRIPTION "Build the Fortran CLOUDSC2 non-linear dwarf with Loki transformations"
)

if( HAVE_CLOUDSC2_NL_LOKI )

    ####################################################
    ##  Define various pre-processing modes via Loki  ##
    ####################################################

    set( COMMON_MODULE "${CMAKE_CURRENT_SOURCE_DIR}/../common/module" )
    set( COMMON_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../common/include" )
    set( XMOD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/xmod" )

    set( LOKI_FRONTEND "fp" CACHE STRING "Frontend parser for Loki transforms" )

    # Ensure xmod directory for OMNI frontend
    file(MAKE_DIRECTORY ${XMOD_DIR})

    ####################################################
    ##  Idempotence mode:                             ##
    ##   * Internal "do-nothing" mode for Loki debug  ##
    ####################################################
    add_custom_command(
        OUTPUT
            loki-idem/cloudsc_driver_mod.idem.F90 loki-idem/cloudsc2.idem.F90
            loki-idem/satur.idem.F90 loki-idem/cuadjtqs.idem.F90
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        COMMAND loki-transform.py convert --mode idem --frontend ${LOKI_FRONTEND}
                      --config ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
                      --path ${CMAKE_CURRENT_SOURCE_DIR}
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --cpp --include ${COMMON_INCLUDE}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        DEPENDS cloudsc2.F90 cloudsc_driver_mod.F90
	COMMENT "[Loki] Pre-processing: mode=idem frontend=${LOKI_FRONTEND}"
    )

    # Define the binary build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc2-nl-loki-idem
        SOURCES
            dwarf_cloudsc.F90
            loki-idem/cloudsc_driver_mod.idem.F90
	    loki-idem/cloudsc2.idem.F90
	    loki-idem/satur.idem.F90
	    loki-idem/cuadjtqs.idem.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc2-nl-loki-idem PRIVATE cloudsc2-common-lib )

    if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc2-nl-loki-idem PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )

endif()
