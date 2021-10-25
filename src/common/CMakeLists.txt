ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)

list(APPEND CLOUDSC2_COMMON_SOURCES
  module/parkind1.F90
  module/yoecldp.F90
  module/yomcst.F90
  module/yoethf.F90
  module/yoephli.F90
  module/yoecld.F90
  module/yoecldp.F90
  module/yophnc.F90
  module/yomncl.F90
  module/yomphyder.F90
  module/fcttre_mod.F90
  module/fccld_mod.F90
  module/routines.F90
  module/abor1.F90
  module/timer_mod.F90
  module/mycpu.c
  module/ec_pmon_mod.F90
  module/expand_mod.F90
  module/file_io_mod.F90
  module/hdf5_file_mod.F90
  module/validate_mod.F90
  module/cloudsc2_array_state_mod.F90
  module/cloudsc_mpi_mod.F90
)


ecbuild_add_library( TARGET cloudsc2-common-lib
    SOURCES ${CLOUDSC2_COMMON_SOURCES}
    DEFINITIONS ${CLOUDSC_DEFINITIONS}
)
target_include_directories( cloudsc2-common-lib PUBLIC $<INSTALL_INTERFACE:include> )
target_include_directories( cloudsc2-common-lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_include_directories( cloudsc2-common-lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/module>)

# We do all I/O via HDF% in this dwarf
target_link_libraries( cloudsc2-common-lib PUBLIC HDF5::HDF5 )

if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
    target_link_libraries( cloudsc2-common-lib PRIVATE OpenMP::OpenMP_Fortran )
endif()

if( HAVE_ACC AND TARGET OpenACC::OpenACC_Fortran )
  target_link_libraries( cloudsc2-common-lib PRIVATE OpenACC::OpenACC_Fortran )
endif()

if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
    target_sources( cloudsc2-common-lib PRIVATE module/cloudsc_mpif.F90 )
    target_link_libraries( cloudsc2-common-lib PRIVATE MPI::MPI_Fortran )
endif()
