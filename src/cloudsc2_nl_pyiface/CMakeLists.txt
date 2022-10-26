# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC2_NL DEFAULT ON
    DESCRIPTION "Build the Fortran CLOUDSC2 non-linear dwarf"
)

if( HAVE_CLOUDSC2_NL )
   message (STATUS "Setting up module directory")
   ecbuild_enable_fortran(MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../module)
   message (STATUS  ${CMAKE_CURRENT_BINARY_DIR}/../../module)
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

    ecbuild_add_executable( TARGET dwarf-cloudsc2-nl-piref
        SOURCES
            ./fortransrc/dwarf_cloudsc.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    target_link_libraries( dwarf-cloudsc2-nl-piref PRIVATE cloudsc-common-lib dwarf-cloudsc2-nl-lib)

    if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc2-nl-piref PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
        target_link_libraries( dwarf-cloudsc2-nl-piref PRIVATE MPI::MPI_Fortran )
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )

# F90wrap section   
     set (PYTHON_MODN cloudsc2)
     set(F90WRAP_SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/../common/module/yomphyder.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/../common/module/parkind1.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/fortransrc/cloudsc_driver_mod.F90
        )
      add_custom_command(TARGET dwarf-cloudsc2-nl-piref POST_BUILD
         COMMAND echo "Current directory:"
         COMMAND pwd
         COMMAND echo ${CMAKE_CURRENT_SOURCE_DIR}  
         COMMAND echo ${CMAKE_CURRENT_BINARY_DIR}  
         COMMAND echo ${CMAKE_CURRENT_BINARY_DIR}/../../module
         COMMAND  f90wrap -m${PYTHON_MODN} ${F90WRAP_SOURCES} -k ${CMAKE_CURRENT_SOURCE_DIR}/kind_map -P
#How to access ecbuild MODULE_DIRECTORY ?
         COMMAND  f2py-f90wrap -c -m _${PYTHON_MODN} --include-paths ${CMAKE_CURRENT_BINARY_DIR}/../../module -I${CMAKE_CURRENT_BINARY_DIR}/../../module -I${CMAKE_CURRENT_BINARY_DIR}/../common/module -L${CMAKE_CURRENT_BINARY_DIR}/../../lib -lcloudsc-common-lib  -ldwarf-cloudsc2-nl-lib f90wrap_cloudsc_driver_mod.f90 f90wrap_yomphyder.f90
         COMMENT "Running F90wrap"
         DEPENDS dwarf-cloudsc2-nl-lib
         VERBATIM
      )
endif()
