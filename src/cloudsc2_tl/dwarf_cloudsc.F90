! Copyright (C) 2003- ECMWF
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
PROGRAM DWARF_CLOUDSC

USE PARKIND1, ONLY: JPIM, JPIB, JPRB
USE CLOUDSC_MPI_MOD, ONLY: CLOUDSC_MPI_INIT, CLOUDSC_MPI_END, NUMPROC, IRANK
USE CLOUDSC2_ARRAY_STATE_MOD, ONLY: CLOUDSC2_ARRAY_STATE
USE CLOUDSC_DRIVER_TL_MOD, ONLY: CLOUDSC_DRIVER_TL
USE EC_PMON_MOD, ONLY: EC_PMON
USE YOECLD   , ONLY : YRECLD
USE YOPHNC   , ONLY : YRPHNC
USE YOEPHLI  , ONLY : YREPHLI
USE YOMNCL   , ONLY : YRNCL


IMPLICIT NONE

CHARACTER(LEN=20) :: CLARG
INTEGER(KIND=JPIM) :: IARGS, LENARG, JARG, I, JK

INTEGER(KIND=JPIM) :: NUMOMP   = 1     ! Number of OpenMP threads for this run
INTEGER(KIND=JPIM) :: NGPTOTG  = 16384 ! Number of grid points (as read from command line)
INTEGER(KIND=JPIM) :: NPROMA   = 32    ! NPROMA blocking factor (currently active)
INTEGER(KIND=JPIM) :: NGPTOT           ! Local number of grid points

REAL(KIND=JPRB) :: ZPRES(0:200),ZPRESF(200) ! KLEV should be bellow 200
REAL(KIND=JPRB) :: VP00

TYPE(CLOUDSC2_ARRAY_STATE) :: GLOBAL_STATE

INTEGER(KIND=JPIB) :: ENERGY, POWER
CHARACTER(LEN=1)   :: CLEC_PMON

CALL GET_ENVIRONMENT_VARIABLE('EC_PMON', CLEC_PMON)
IF (CLEC_PMON == '1') THEN
  CALL EC_PMON(ENERGY, POWER)
  print *, "EC_PMON:: Initial (idle) power: ", POWER
END IF

IARGS = COMMAND_ARGUMENT_COUNT()

! Get the number of OpenMP threads to use for the benchmark
if (IARGS >= 1) then
   CALL GET_COMMAND_ARGUMENT(1, CLARG, LENARG)
   READ(CLARG(1:LENARG),*) NUMOMP
end if

! Initialize MPI environment
CALL CLOUDSC_MPI_INIT(NUMOMP)

! Get total number of grid points (NGPTOT) with which to run the benchmark
IF (IARGS >= 2) THEN
  CALL GET_COMMAND_ARGUMENT(2, CLARG, LENARG)
  READ(CLARG(1:LENARG),*) NGPTOTG
END IF

! Determine local number of grid points
NGPTOT = (NGPTOTG - 1) / NUMPROC + 1
if (IRANK == NUMPROC - 1) then
  NGPTOT = NGPTOTG - (NUMPROC - 1) * NGPTOT
end if

! Get the block size (NPROMA) for which to run the benchmark  
IF (IARGS >= 3) THEN
  CALL GET_COMMAND_ARGUMENT(3, CLARG, LENARG)
  READ(CLARG(1:LENARG),*) NPROMA
ENDIF

! TODO: Create a global global memory state from serialized input data
CALL GLOBAL_STATE%LOAD(NPROMA, NGPTOT, NGPTOTG)

!cloudsc2
!set up other modules not initialized in cloudsc
!YRECLD
allocate(YRECLD)
allocate(YRECLD%CETA(GLOBAL_STATE%KLEV))
! security
if (GLOBAL_STATE%KLEV>200) then
  print *, ' Dimension of ZPRES/ZPRESF is too short. '
  stop
endif
! Here we do bit shortcut wrt what it should be when following the model
!VP00 = 101325._JPRB
!ZPRES(GLOBAL_STATE%KLEV)=VP00
!CALL GPHPRE(1,GLOBAL_STATE%KLEV,1,1,YDVAB,ZPRES,PRESF=ZPRESF)
!DO JK=1,GLOBAL_STATE%KLEV
!  YRECLD%CETA(JK)= ZPRESF(JK)/ZPRES(NFLEVG)
!ENDDO
! Instead we do this similar computation involving values already available
DO JK=1,GLOBAL_STATE%KLEV
  YRECLD%CETA(JK)= GLOBAL_STATE%PAP(1,JK,1)/GLOBAL_STATE%PAPH(1,GLOBAL_STATE%KLEV+1,1)
ENDDO
! YRPHNC
allocate(YRPHNC)
YRPHNC%LEVAPLS2=.false.
!YRNCL
allocate(YRNCL)
YRNCL%LREGCL=.false.  ! no regularization in TL test
! overload LPHYLIN
YREPHLI%LPHYLIN=.true.

! Call the driver to perform the parallel loop over our kernel
CALL CLOUDSC_DRIVER_TL(NUMOMP, NPROMA, GLOBAL_STATE%KLEV, NGPTOT, NGPTOTG, GLOBAL_STATE%PTSPHY, &
     & GLOBAL_STATE%PT, GLOBAL_STATE%PQ, &
     & GLOBAL_STATE%TENDENCY_CML, GLOBAL_STATE%TENDENCY_LOC, &
     & GLOBAL_STATE%PAP,      GLOBAL_STATE%PAPH, &
     & GLOBAL_STATE%PLU,      GLOBAL_STATE%PLUDE, &
     & GLOBAL_STATE%PMFU,     GLOBAL_STATE%PMFD, &
     & GLOBAL_STATE%PA,       GLOBAL_STATE%PCLV,  GLOBAL_STATE%PSUPSAT,&
     & GLOBAL_STATE%PCOVPTOT, &
     & GLOBAL_STATE%PFPLSL,   GLOBAL_STATE%PFPLSN,   GLOBAL_STATE%PFHPSL,   GLOBAL_STATE%PFHPSN &
     & )

!! Validate the output against serialized reference data
!CALL GLOBAL_STATE%VALIDATE(NPROMA, NGPTOT, NGPTOTG)

! Tear down MPI environment
CALL CLOUDSC_MPI_END()

END PROGRAM DWARF_CLOUDSC
