! Copyright (C) 2003- ECMWF
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CLOUDSC_DRIVER_PYIFACE_MOD
  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM
  USE EC_PMON_MOD, ONLY: EC_PMON
  USE CLOUDSC2_ARRAY_STATE_MOD, ONLY: CLOUDSC2_ARRAY_STATE
  USE YOPHNC   , ONLY : TPHNC
  USE YOECLD   , ONLY : TECLD
  USE YOECLDP  , ONLY : TECLDP
  USE YOEPHLI  , ONLY : TEPHLI
  USE YOMCST   , ONLY : TOMCST
  USE YOETHF   , ONLY : TOETHF
  IMPLICIT NONE
  SAVE
  TYPE(CLOUDSC2_ARRAY_STATE) :: GLOBAL_STATE


CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_NO_DERV_TPES( &
     & NUMOMP, NPROMA, NLEV, NCLV, NCLDQL, NCLDQI, & 
     & NGPTOT, NGPBLKS,  NGPTOTG, PTSPHY, &
     & PT, PQ, &
     & BUFFER_CML, BUFFER_LOC, &
     & PAP,      PAPH, &
     & PLU,      PLUDE,    PMFU,     PMFD, &
     & PA,       PCLV,     PSUPSAT,&
     & PCOVPTOT, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN , & 
     & YDOMCST, YDOETHF, YDECLD, & 
     & YDECLDP, YDEPHLI, YDPHNC )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC2 kernel

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NCLV, NCLDQL, NCLDQI, NGPTOT, NGPBLKS, NGPTOTG
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: PQ(:,:,:)    ! Q at start of callpar
    REAL(KIND=JPRB), INTENT(INOUT)    :: BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB), INTENT(INOUT)    :: BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB),    INTENT(IN)    :: PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(INOUT) :: PLUDE(:,:,:) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFD(:,:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(INOUT) :: PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(IN)    :: PCLV(:,:,:,:) 
    REAL(KIND=JPRB),    INTENT(IN)    :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB),    INTENT(INOUT) :: PCOVPTOT(:,:,:) ! Precip fraction
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSL(:,:,:) ! liq+rain sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSN(:,:,:) ! ice+snow sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSL(:,:,:) ! Enthalpy flux for liq
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSN(:,:,:) ! Enthalp flux for ice
    TYPE(TOMCST)   ,    INTENT(INOUT) :: YDOMCST 
    TYPE(TOETHF)   ,    INTENT(INOUT) :: YDOETHF 
    TYPE(TECLD)    ,    INTENT(INOUT) :: YDECLD 
    TYPE(TECLDP)   ,    INTENT(INOUT) :: YDECLDP
    TYPE(TEPHLI)   ,    INTENT(INOUT) :: YDEPHLI
    TYPE(TPHNC)    ,    INTENT(INOUT) :: YDPHNC

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND

    TYPE(PERFORMANCE_TIMER) :: TIMER
    REAL(KIND=JPRD), PARAMETER :: ZHPM = 3996006.0_JPRD  ! The nominal number of flops per 100 columns

    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    LOGICAL            :: LDRAIN1D = .FALSE.
    REAL(KIND=JPRB)    :: ZQSAT(NPROMA,NLEV) ! local array

1003 format(5x,'NUMPROC=',i0', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
   if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
   end if

    ! Global timer for the parallel region
   CALL TIMER%START(NUMOMP)

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID) &
    !$omp& private(ZQSAT) &
    !$omp& num_threads(NUMOMP)

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)
    CALL INITIALIZE_FCTTRE_PARAMETERS(YDOMCST, YDOETHF, YDECLD, YDECLDP, YDEPHLI, YDPHNC )

    !$omp do schedule(runtime)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         PCOVPTOT(:,:,IBL) = 0.0_JPRB

         ! Fill in ZQSAT
         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:), 2 )
 
         CALL CLOUDSC2 ( &
              &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
              & PTSPHY,&
              & PAPH(:,:,IBL),  PAP(:,:,IBL), &
              & PQ(:,:,IBL), ZQSAT(:,:), PT(:,:,IBL), &
              & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
              & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
              &  BUFFER_LOC(:,:,1,IBL), BUFFER_CML(:,:,1,IBL), &
              &  BUFFER_LOC(:,:,3,IBL), BUFFER_CML(:,:,3,IBL), &
              &  BUFFER_LOC(:,:,3+NCLDQL,IBL), BUFFER_CML(:,:,3+NCLDQL,IBL),  &
              &  BUFFER_LOC(:,:,3+NCLDQI,IBL), BUFFER_CML(:,:,3+NCLDQI,IBL),  &
              &  PSUPSAT(:,:,IBL), &
              &  PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
              &  PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL), &
              &  YDOMCST, YDOETHF, YDECLD, YDECLDP, YDEPHLI, YDPHNC )

         ! Log number of columns processed by this thread
         CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
      ENDDO

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc
      !$omp end do nowait

       CALL TIMER%THREAD_END(TID)

      !$omp end parallel

      CALL TIMER%END()

       CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, ZHPM, NGPTOT)
    
  CONTAINS
SUBROUTINE INITIALIZE_FCTTRE_PARAMETERS(YDOMCST, YDOETHF, YDECLD, YDECLDP, YDEPHLI, YDPHNC )
USE FCTTRE_MOD, ONLY : FCTTRE_CONSTANTS_SET

TYPE(TOMCST)      ,INTENT(IN) :: YDOMCST
TYPE(TOETHF)      ,INTENT(IN) :: YDOETHF
TYPE(TECLD)       ,INTENT(IN) :: YDECLD
TYPE(TECLDP)      ,INTENT(IN) :: YDECLDP
TYPE(TEPHLI)      ,INTENT(IN) :: YDEPHLI
TYPE(TPHNC)       ,INTENT(IN) :: YDPHNC

CALL FCTTRE_CONSTANTS_SET(  &
 &                 RG_IN=YDOMCST%RG, RD_IN=YDOMCST%RD, RCPD_IN=YDOMCST%RCPD, &
 &                 RETV_IN=YDOMCST%RETV, RLVTT_IN=YDOMCST%RLVTT, &
 &                 RLSTT_IN=YDOMCST%RLSTT, RLMLT_IN=YDOMCST%RLMLT, RTT_IN=YDOMCST%RTT, & 
 &                 RV_IN=YDOMCST%RV, RA_IN=YDOMCST%RA, RPI_IN=YDOMCST%RPI, &
 &                 R2ES_IN=YDOETHF%R2ES, R3LES_IN=YDOETHF%R3LES, R3IES_IN=YDOETHF%R3IES, &
 &                 R4LES_IN=YDOETHF%R4LES, R4IES_IN=YDOETHF%R4IES, &
 &                 R5LES_IN=YDOETHF%R5LES, R5IES_IN=YDOETHF%R5IES, &
 &                 R5ALVCP_IN=YDOETHF%R5ALVCP, R5ALSCP_IN=YDOETHF%R5ALSCP,&
 &                 RALVDCP_IN=YDOETHF%RALVDCP, RALSDCP_IN=YDOETHF%RALSDCP, &
 &                 RALFDCP_IN=YDOETHF%RALFDCP, RTWAT_IN=YDOETHF%RTWAT,&
 &                 RTICE_IN=YDOETHF%RTICE, RTICECU_IN=YDOETHF%RTICECU, &
 &                 RTWAT_RTICE_R_IN=YDOETHF%RTWAT_RTICE_R,&
 &                 RTWAT_RTICECU_R_IN=YDOETHF%RTWAT_RTICECU_R, &
 &                 RKOOP1_IN=YDOETHF%RKOOP1, RKOOP2_IN=YDOETHF%RKOOP2)

END SUBROUTINE INITIALIZE_FCTTRE_PARAMETERS
  END SUBROUTINE CLOUDSC_DRIVER_NO_DERV_TPES

END MODULE CLOUDSC_DRIVER_PYIFACE_MOD
