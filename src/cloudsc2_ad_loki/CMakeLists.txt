# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC2_AD DEFAULT ON
    DESCRIPTION "Build the Fortran CLOUDSC2 adjoint dwarf"
)

if( HAVE_CLOUDSC2_AD )

    # Define the binary build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc2-ad
        SOURCES
            dwarf_cloudsc.F90
            cloudsc_driver_ad_mod.F90
            satur.F90
            cuadjtqs.F90
            cuadjtqstl.F90
            cuadjtqsad.F90
            cloudsc2ad.F90
            cloudsc2tl.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc2-ad PRIVATE cloudsc2-common-lib )

    if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc2-ad PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
        target_link_libraries( dwarf-cloudsc2-ad PRIVATE MPI::MPI_Fortran )
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )

    ecbuild_add_test(
        TARGET dwarf-cloudsc2-ad
        COMMAND bin/dwarf-cloudsc2-ad
        ARGS
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
endif()
