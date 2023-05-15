# (C) Copyright 1988- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# define minimum version of cmake required
cmake_minimum_required( VERSION 3.17 FATAL_ERROR )

# Disable warnings about setting `ENABLE_ACC` variable for ecbuild_add_option
cmake_policy( SET CMP0077 NEW )

find_package( ecbuild REQUIRED )

# define the project
project(dwarf-p-cloudsc2-tl-ad VERSION 1.0 LANGUAGES C Fortran )

include( cmake/compat.cmake )
if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  ecbuild_add_fortran_flags("-ffree-line-length-none")
  if( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0")
    ecbuild_add_fortran_flags("-fallow-argument-mismatch")
  endif()
  if( CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "11.0")
    set( ENABLE_ACC OFF )
  endif()
endif()

include( cloudsc_compile_options )

if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(DIRECTIVE none)
endif()
### OpenACC
if( NOT DEFINED ENABLE_ACC OR ENABLE_ACC )
    # Incredibly inconvenient: FindOpenACC does _not_ set OpenACC_FOUND, only
    # the language-specific components OpenACC_Fortran_FOUND and OpenACC_C_FOUND.
    # This means, even internally CMake considers OpenACC as not found.
    # (See eg get_property(... GLOBAL PROPERTY PACKAGES_NOT_FOUND))
    # Therefore, we search for OpenACC, set OpenACC_FOUND ourselves according to
    # the result, and then, trigger a second find_package via ecbuild_add_option.
    # This then conveniently takes the previously set OpenACC_FOUND into account
    # and rectifies CMake's internal bookkeeping in the process.
    find_package( OpenACC )
    if( OpenACC_Fortran_FOUND AND OpenACC_C_FOUND )
        set( OpenACC_FOUND ON )
    endif()
endif()
ecbuild_add_option( FEATURE ACC
    DESCRIPTION "OpenACC" DEFAULT ON
    REQUIRED_PACKAGES "OpenACC" )


### OpenMP
ecbuild_add_option( FEATURE OMP
    DESCRIPTION "OpenMP" DEFAULT ON
    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran C" )
include(features/OMP)

### MPI
ecbuild_add_option( FEATURE MPI
    DESCRIPTION "MPI" DEFAULT OFF
    REQUIRED_PACKAGES "MPI COMPONENTS Fortran" )
if( HAVE_MPI )
    list(APPEND CLOUDSC_DEFINITIONS HAVE_MPI )
endif()

### HDF5
ecbuild_add_option( FEATURE HDF5
    DESCRIPTION "Use HDF5 to read input and reference data"
    REQUIRED_PACKAGES "HDF5 COMPONENTS Fortran"
    DEFAULT ON )
if( HAVE_HDF5 )
    list(APPEND CLOUDSC_DEFINITIONS HAVE_HDF5 )
endif()

ecbuild_find_package( NAME loki REQUIRED )

# build executables
add_subdirectory(src)

# finalize
ecbuild_install_project(NAME dwarf-P-cloudsc2-tl-ad)

# print summary
ecbuild_print_summary()
