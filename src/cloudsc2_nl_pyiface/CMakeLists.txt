# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC2_NL_PYIFACE DEFAULT ON
    DESCRIPTION "Build the Fortran CLOUDSC2 non-linear dwarf with an interface to Python"
)

if( HAVE_CLOUDSC2_NL_PYIFACE )
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

     set (PYTHON_MODN cloudsc2)

     set(FORTRAN_PYTHON_COMMON_SOURCES
            yomphyder
            yoecld
            yoecldp
            yoethf
            yomcst
            yoephli
            yophnc
        )
     set(FORTRAN_PYTHON_DWARF_SOURCES
            cloudsc_driver_mod
        )
     set (F90WRAP_COMMON_SOURCES ${FORTRAN_PYTHON_COMMON_SOURCES})
     set (F90WRAP_DWARF_SOURCES  ${FORTRAN_PYTHON_DWARF_SOURCES})
     set (F2PY_COMMON_SOURCES ${FORTRAN_PYTHON_COMMON_SOURCES})
     set (F2PY_DWARF_SOURCES  ${FORTRAN_PYTHON_DWARF_SOURCES})
     list(TRANSFORM F90WRAP_COMMON_SOURCES PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/../common/module/)
     list(TRANSFORM F90WRAP_COMMON_SOURCES APPEND .F90 )
     list(TRANSFORM F90WRAP_DWARF_SOURCES PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/fortransrc/)
     list(TRANSFORM F90WRAP_DWARF_SOURCES APPEND .F90 )
     list(TRANSFORM F2PY_COMMON_SOURCES PREPEND f90wrap_)
     list(TRANSFORM F2PY_COMMON_SOURCES APPEND .f90 )
     list(TRANSFORM F2PY_DWARF_SOURCES PREPEND f90wrap_)
     list(TRANSFORM F2PY_DWARF_SOURCES APPEND .f90 )
     message("${F90WRAP_COMMON_SOURCES}")
     message("${F90WRAP_PYTHON_SOURCES}")
     message("${F2PY_COMMON_SOURCES}")
     message("${F2PY_DWARF_SOURCES}")


     set(PYTHON_SOURCES
            ./pythonsrc/dwarfdriver.py
            ./pythonsrc/cloudsc2_data.py
        )

     add_custom_target( pythonexecs ALL SOURCES ${PYTHON_SOURCES}
         COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pythonsrc  ${CMAKE_CURRENT_BINARY_DIR}/pythonsrc
         COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/pythonsrc  ../../bin/pythonexec 
         COMMAND sh -c "rm -f ../../bin/pythonexec/*.f90" 
         COMMENT "Copy python scripts to executable directory"
        )


    target_link_libraries( dwarf-cloudsc2-nl-lib PRIVATE cloudsc-common-lib HDF5::HDF5 )
    ecbuild_add_executable( TARGET dwarf-cloudsc2-nl-piref
        SOURCES
            ./fortransrc/dwarf_cloudsc.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    add_dependencies(dwarf-cloudsc2-nl-piref pythonexecs)

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

# F90wrap execution section  
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pythonsrc  ${CMAKE_CURRENT_BINARY_DIR}/pythonsrc)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/pythonsrc/input.h5 )
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/pythonsrc/reference.h5 )


      add_custom_command(TARGET dwarf-cloudsc2-nl-piref POST_BUILD 
         WORKING_DIRECTORY ./pythonsrc
         COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pythonsrc  ${CMAKE_CURRENT_BINARY_DIR}/pythonsrc
         COMMAND echo "Executing f90wrap"
         COMMAND  f90wrap -m${PYTHON_MODN} ${F90WRAP_COMMON_SOURCES} ${F90WRAP_DWARF_SOURCES} -k ${CMAKE_CURRENT_SOURCE_DIR}/pythonsrc/kind_map -P > f90wrap_log.txt 
         COMMAND echo "Executing f2py-f90wrap"
         COMMAND pwd
         COMMAND  f2py-f90wrap -c -m _${PYTHON_MODN} --include-paths ${CMAKE_CURRENT_BINARY_DIR}/../../module
                                                                   -I${CMAKE_CURRENT_BINARY_DIR}/../../module
                                                                   -I${CMAKE_CURRENT_BINARY_DIR}/../common/module
                                                                   -L${CMAKE_CURRENT_BINARY_DIR}/../../lib 
                                                                   -lcloudsc-common-lib  -ldwarf-cloudsc2-nl-lib ${F2PY_COMMON_SOURCES} ${F2PY_DWARF_SOURCES} > f2py_log.txt
         COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/pythonsrc/  ../../bin/pythonexec/ 
         COMMAND echo ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/../../lib/  ${CMAKE_CURRENT_BINARY_DIR}/../../bin/pythonexec 
         COMMAND sh -c "rm -f ../../bin/pythonexec/*.f90" 
         COMMENT "Running F90wrap"
         DEPENDS dwarf-cloudsc2-nl-lib
         VERBATIM
      )

endif()
