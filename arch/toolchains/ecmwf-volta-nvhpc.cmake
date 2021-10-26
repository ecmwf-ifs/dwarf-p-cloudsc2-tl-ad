####################################################################
# COMPILER
####################################################################

set( ECBUILD_FIND_MPI ON )

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS             "-mp -mp=bind,allcores,numa" )
set( OpenMP_CXX_FLAGS           "-mp -mp=bind,allcores,numa" )
set( OpenMP_Fortran_FLAGS       "-mp -mp=bind,allcores,numa" )

####################################################################
# OpenAcc FLAGS
####################################################################

set( OpenACC_Fortran_FLAGS "-acc -ta=tesla:lineinfo,deepcopy,fastmath" )
set( OpenACC_Fortran_FLAGS "${OpenACC_Fortran_FLAGS} -Minfo=accel" )
set( OpenACC_Fortran_FLAGS "${OpenACC_Fortran_FLAGS} -Mcuda" )

set( ECBUILD_Fortran_LINK_FLAGS "-acc -ta=tesla:pinned,lineinfo,deepcopy,fastmath" )

####################################################################
# COMMON FLAGS
####################################################################

set(ECBUILD_Fortran_FLAGS "-fpic")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mframe")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mbyteswapio")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mstack_arrays")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mrecursive")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Ktrap=fp")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Kieee")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mdaz")


set( ECBUILD_Fortran_FLAGS_BIT "-O2 -gopt" )

set( ECBUILD_C_FLAGS "-O2 -gopt -traceback" )

set( ECBUILD_CXX_FLAGS "-O2 -gopt" )

# Fix for C++ template headers needed for Serialbox
set( GNU_HEADER_INCLUDE "-I/usr/local/apps/gcc/7.3.0/lib/gcc/x86_64-linux-gnu/7.3.0/include-fixed" )
set( ECBUILD_CXX_FLAGS "${ECBUILD_CXX_FLAGS} ${GNU_HEADER_INCLUDE}" )
