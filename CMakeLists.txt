###############################################################################
# (C) Copyright 1996-2015 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
###############################################################################

# define minimum version of cmake required
cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

find_package( ecbuild REQUIRED )

# define the project
project(dwarf-p-cloudsc2-tl-ad VERSION 1.0 LANGUAGES C Fortran )

include( cmake/compat.cmake )

### OpenMP
if( NOT DEFINED ENABLE_OMP OR ENABLE_OMP )
    find_package( OpenMP COMPONENTS Fortran C )
    # Ensure targets are available
    if( OpenMP_Fortran_FOUND AND NOT TARGET OpenMP::OpenMP_Fortran )
        ecbuild_critical( "OpenMP::OpenMP_Fortran should have been defined as a CMake target" )
    endif()
    if( OpenMP_C_FOUND AND NOT TARGET OpenMP::OpenMP_C )
        ecbuild_critical( "OpenMP::OpenMP_C should have been defined as a CMake target" )
    endif()
endif()
ecbuild_add_option( FEATURE OMP
    DESCRIPTION "OpenMP" DEFAULT ON
    CONDITION OpenMP_C_FOUND AND OpenMP_Fortran_FOUND )

# HDF5: This dwarf fully depends on HDF5 for input and reference data
ecbuild_find_package( NAME HDF5 REQUIRED COMPONENTS Fortran )
list(APPEND CLOUDSC_DEFINITIONS HAVE_HDF5 ${HDF5_Fortran_DEFINITIONS})

# build executables
add_subdirectory(src)

# finalize
ecbuild_install_project(NAME dwarf-P-cloudsc2-tl-ad)

# add symlinks to input/validation data to the bin directory 
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ../../config-files/input.h5 ../build/bin/input.h5 )
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
     ../../config-files/reference.h5 ../build/bin/reference.h5 )

# print summary
ecbuild_print_summary()
