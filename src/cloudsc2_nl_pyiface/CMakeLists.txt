# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC2_NL DEFAULT ON
    DESCRIPTION "Build the Fortran CLOUDSC2 non-linear dwarf"
)

if( HAVE_CLOUDSC2_NL )

    # Define the binary build target for this variant
    ecbuild_add_library( TARGET dwarf-cloudsc2-nl-lib
        SOURCES
            ./fortransrc/dwarf_cloudsc.F90
            ./fortransrc/cloudsc_driver_mod.F90
            ./fortransrc/satur.F90
            ./fortransrc/cuadjtqs.F90
            ./fortransrc/cloudsc2.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    target_link_libraries( dwarf-cloudsc2-nl-lib PRIVATE cloudsc-common-lib HDF5::HDF5 )

    ecbuild_add_executable( TARGET dwarf-cloudsc2-nl
        SOURCES
            dwarf_cloudsc.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    target_link_libraries( dwarf-cloudsc2-nl PRIVATE cloudsc-common-lib dwarf-cloudsc2-nl-lib)

    if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc2-nl PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
        target_link_libraries( dwarf-cloudsc2-nl PRIVATE MPI::MPI_Fortran )
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )
    
endif()
