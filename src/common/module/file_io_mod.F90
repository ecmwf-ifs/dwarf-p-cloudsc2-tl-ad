! Copyright (C) 2003- ECMWF
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
module file_io_mod
  USE PARKIND1 , ONLY : JPIM, JPRB, JPRD

#ifdef HAVE_SERIALBOX
  USE m_serialize, ONLY: &
   fs_create_savepoint, &
   fs_add_serializer_metainfo, &
   fs_get_serializer_metainfo, &
   fs_read_field, &
   fs_write_field
  USE utils_ppser, ONLY:  &
   ppser_initialize, &
   ppser_finalize, &
   ppser_serializer, &
   ppser_serializer_ref, &
   ppser_set_mode, &
   ppser_savepoint
#endif

#ifdef HAVE_HDF5
  USE hdf5_file_mod, only: hdf5_file
#endif

  implicit none

#ifdef HAVE_HDF5
  TYPE(hdf5_file) :: data_file
#endif

  interface load_scalar
    procedure load_scalar_real, load_scalar_int, load_scalar_log
  end interface load_scalar

  interface load_array
    procedure load_array_l1, load_array_i1, load_array_r1, load_array_r2, load_array_r3
  end interface load_array

  interface write_scalar
    procedure write_scalar_real, write_scalar_int
  end interface write_scalar

  interface write_array
    procedure write_array_r1, write_array_r2, write_array_r3
  end interface write_array

contains

  subroutine input_initialize(name, read_only)
    ! Initialize input file for a given data set name
    character(len=*), intent(in) :: name
    logical, intent(in) :: read_only
#ifdef HAVE_SERIALBOX
    ! Get dimensions information from stored metadata
    call ppser_initialize(directory='data', prefix='dummy', prefix_ref=NAME)
    call fs_create_savepoint(NAME, ppser_savepoint)
    call ppser_set_mode(1)
#elif defined(HAVE_HDF5)
    call data_file%open_file(NAME//'.h5', rdwr=.not.read_only)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine input_initialize

  subroutine input_finalize()
#ifdef HAVE_SERIALBOX
    call ppser_finalize()
#elif defined(HAVE_HDF5)
    call data_file%close_file()
#else
    call abor1('ERROR: Serialbox not found.')
#endif
  end subroutine input_finalize

    
  subroutine load_scalar_real(name, variable)
    character(len=*), intent(in) :: name
    real(kind=JPRB), intent(inout) :: variable
    real(kind=JPRD) :: buffer

#ifdef HAVE_SERIALBOX
    call fs_get_serializer_metainfo(ppser_serializer_ref, name, buffer)
    variable = buffer
#elif defined(HAVE_HDF5)
    call data_file%load(name, buffer)
    variable = buffer
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_scalar_real

  subroutine load_scalar_int(name, variable)
    character(len=*), intent(in) :: name
    integer(kind=JPIM), intent(inout) :: variable

#ifdef HAVE_SERIALBOX
    call fs_get_serializer_metainfo(ppser_serializer_ref, name, variable)
#elif defined(HAVE_HDF5)
    call data_file%load(name, variable)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_scalar_int

  subroutine load_scalar_log(name, variable)
    character(len=*), intent(in) :: name
    logical, intent(inout) :: variable
    integer(kind=4) :: tmp

#ifdef HAVE_SERIALBOX
    call fs_get_serializer_metainfo(ppser_serializer_ref, name, variable)
#elif defined(HAVE_HDF5)
    call data_file%load(name, tmp)
    variable = merge(.TRUE., .FALSE., tmp > 0)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_scalar_log

  
  subroutine load_array_i1(name, start, end, size, nlon, buffer)
    ! Load data from file into the local memory buffer
    character(len=*), intent(in) :: name
    integer(kind=jpim), intent(in) :: start, end, size, nlon
    integer(kind=jpim), intent(out) :: buffer(size)
    integer(kind=jpim), allocatable :: rbuf(:)

#ifdef HAVE_SERIALBOX
    allocate(rbuf(nlon))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, rbuf)
    buffer(:) = rbuf(start:end)
    deallocate(rbuf)
#elif defined(HAVE_HDF5)
    call data_file%load(name, buffer, start, size)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_array_i1

  subroutine load_array_l1(name, start, end, size, nlon, buffer)
    ! Load data from file into the local memory buffer
    character(len=*), intent(in) :: name
    integer(kind=jpim), intent(in) :: start, end, size, nlon
    logical, intent(out) :: buffer(size)
    logical, allocatable :: rbuf(:)

#ifdef HAVE_SERIALBOX
    allocate(rbuf(nlon))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, rbuf)
    buffer(:) = rbuf(start:end)
    deallocate(rbuf)
#elif defined(HAVE_HDF5)
    call data_file%load(name, buffer, start, size)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_array_l1

  subroutine load_array_r1(name, start, end, size, nlon, buffer)
    ! Load data from file into the local memory buffer
    character(len=*), intent(in) :: name
    integer(kind=jpim), intent(in) :: start, end, size, nlon
    real(kind=jprb), intent(out) :: buffer(size)
    real(kind=jprd), allocatable :: rbuf(:)

#ifdef HAVE_SERIALBOX
    allocate(rbuf(nlon))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, rbuf)
    buffer(:) = rbuf(start:end)
    deallocate(rbuf)
#elif defined(HAVE_HDF5)
    allocate(rbuf(size))
    call data_file%load(name, rbuf, start, size)
    buffer(:) = rbuf(:)
    deallocate(rbuf)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_array_r1

  subroutine load_array_r2(name, start, end, size, nlon, nlev, buffer)
    ! Load data from file into the local memory buffer
    character(len=*), intent(in) :: name
    integer(kind=jpim), intent(in) :: start, end, size, nlon, nlev
    real(kind=jprb), intent(out) :: buffer(size,nlev)
    integer(kind=jpim) :: istart(2), isize(2)
    real(kind=jprd), allocatable :: rbuf(:,:)

#ifdef HAVE_SERIALBOX
    allocate(rbuf(nlon,nlev))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, rbuf)
    buffer(:,:) = rbuf(start:end,:)
    deallocate(rbuf)
#elif defined(HAVE_HDF5)
    istart(1) = start
    istart(2) = 1
    isize(1) = size
    isize(2) = nlev
    allocate(rbuf(size,nlev))
    call data_file%load(name, rbuf, istart, isize)
    buffer(:,:) = rbuf(:,:)
    deallocate(rbuf)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_array_r2

  subroutine load_array_r3(name, start, end, size, nlon, nlev, ndim, buffer)
    ! Load data from file into the local memory buffer
    character(len=*), intent(in) :: name
    integer(kind=jpim), intent(in) :: start, end, size, nlon, nlev, ndim
    real(kind=jprb), intent(out) :: buffer(size,nlev,ndim)
    integer(kind=jpim) :: istart(3), isize(3)
    real(kind=jprd), allocatable :: rbuf(:,:,:)

#ifdef HAVE_SERIALBOX
    allocate(rbuf(nlon,nlev,ndim))
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, name, rbuf)
    buffer(:,:,:) = rbuf(start:end,:,:)
    deallocate(rbuf)
#elif defined(HAVE_HDF5)
    istart(1) = start
    istart(2) = 1
    istart(3) = 1
    isize(1) = size
    isize(2) = nlev
    isize(3) = ndim
    allocate(rbuf(size,nlev,ndim))
    call data_file%load(name, rbuf, istart, isize)
    buffer(:,:,:) = rbuf(:,:,:)
    deallocate(rbuf)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine load_array_r3


  subroutine write_scalar_real(name, variable)
    character(len=*), intent(in) :: name
    real(kind=JPRB), intent(inout) :: variable
    real(kind=JPRD) :: buffer

#ifdef HAVE_SERIALBOX
    buffer = variable
    call fs_add_serializer_metainfo(ppser_serializer, name, buffer)
#elif defined(HAVE_HDF5)
    call abor1('ERROR: HDF5 writes not yet supported.')
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine write_scalar_real

  subroutine write_scalar_int(name, variable)
    character(len=*), intent(in) :: name
    integer(kind=JPIM), intent(inout) :: variable
#ifdef HAVE_SERIALBOX
    call fs_add_serializer_metainfo(ppser_serializer, name, buffer)
#elif defined(HAVE_HDF5)
    call data_file%write(name, variable)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine write_scalar_int


  subroutine write_array_r1(name, nlon, buffer)
    ! Load data from file into the local memory buffer
    character(len=*), intent(in) :: name
    integer(kind=jpim), intent(in) :: nlon
    real(kind=jprb), intent(out) :: buffer(nlon)
    real(kind=jprd), allocatable :: wbuf(:)

#ifdef HAVE_SERIALBOX
    allocate(wbuf(nlon))
    wbuf(:) = buffer(:nlon)
    call fs_write_field(ppser_serializer_ref, name, wbuf)
    deallocate(wbuf)
#elif defined(HAVE_HDF5)
    allocate(wbuf(nlon))
    wbuf(:) = buffer(:nlon)
    call data_file%write(name, wbuf)
    deallocate(wbuf)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine write_array_r1

  subroutine write_array_r2(name, nlon, nlev, buffer)
    ! Load data from file into the local memory buffer
    character(len=*), intent(in) :: name
    integer(kind=jpim), intent(in) :: nlon, nlev
    real(kind=jprb), intent(out) :: buffer(nlon, nlev)
    real(kind=jprd), allocatable :: wbuf(:,:)

#ifdef HAVE_SERIALBOX
    allocate(wbuf(nlon,nlev))
    wbuf(:,:) = buffer(:nlon,:nlev)
    call fs_write_field(ppser_serializer_ref, name, wbuf)
    deallocate(wbuf)
#elif defined(HAVE_HDF5)
    allocate(wbuf(nlon,nlev))
    wbuf(:,:) = buffer(:nlon,:nlev)
    call data_file%write(name, wbuf)
    deallocate(wbuf)
#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine write_array_r2

  subroutine write_array_r3(name, nlon, nlev, ndim, buffer)
    ! Load data from file into the local memory buffer
    character(len=*), intent(in) :: name
    integer(kind=jpim), intent(in) :: nlon, nlev, ndim
    real(kind=jprb), intent(out) :: buffer(nlon, nlev, ndim)
    real(kind=jprd), allocatable :: wbuf(:,:,:)

#ifdef HAVE_SERIALBOX
    allocate(wbuf(nlon,nlev,ndim))
    wbuf(:) = buffer(:nlon,:nlev,:ndim)
    call fs_write_field(ppser_serializer_ref, name, wbuf)
    deallocate(wbuf)
#elif defined(HAVE_HDF5)
    allocate(wbuf(nlon,nlev,ndim))
    wbuf(:,:,:) = buffer(:nlon,:nlev,:ndim)
    call data_file%write(name, wbuf)
    deallocate(wbuf)

#else
    call abor1('ERROR: Serialbox and HDF5 not found.')
#endif
  end subroutine write_array_r3

end module file_io_mod
