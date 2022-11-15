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
! USE YOECLD   , ONLY : YRECLD
! USE YOEPHLI  , ONLY : YREPHLI
  IMPLICIT NONE
  SAVE
  TYPE(CLOUDSC2_ARRAY_STATE) :: GLOBAL_STATE


CONTAINS
! In this file the legacy driver with derived type argument removed (CLOUDSC_DRIVER_NO_DERV_TYPES) is
! supplemented with the experiment init (CLOUDSC_DRIVER_INIT) and validation(CLOUDSC_DRIVER_VALIDATION)
! extracted from the topmost dwarf_cloudsc.F90 experiment driver

  SUBROUTINE CLOUDSC_DRIVER_INIT( &
     & NUMOMP, NPROMA, NLEV, NCLV,  &
     & NGPTOT, NGPBLKS, NGPTOTG, &
       PTSPHY, &
     & OUTPUT_PT, OUTPUT_PQ, &
     & OUTPUT_B_CML, OUTPUT_B_LOC, &
     & OUTPUT_PAP,      OUTPUT_PAPH, &
     & OUTPUT_PLU,      OUTPUT_PLUDE,    OUTPUT_PMFU,     OUTPUT_PMFD, &
     & OUTPUT_PA,       OUTPUT_PCLV,     OUTPUT_PSUPSAT,&
     & OUTPUT_PCOVPTOT , &
     & OUTPUT_YDMCST, OUTPUT_YDOETHF, OUTPUT_YDECLD,  &
     & OUTPUT_YDECLDP, OUTPUT_YDEPHLI, OUTPUT_YDPHNC)  

   USE YOMCST   , ONLY : RETV, RG  ,RCPD     ,&
    & RLVTT    ,RLSTT    ,RLMLT    ,RTT     ,RD, &
      RA, RPI, RV
   USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
    & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
    & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU   ,&
    & RTWAT_RTICE_R      ,RTWAT_RTICECU_R    ,RVTMP2, &
    & RALFDCP, RKOOP1, RKOOP2
  USE YOECLD   , ONLY : YRECLD
  USE YOECLDP  , ONLY : YRECLDP
  USE YOEPHLI  , ONLY : YREPHLI
  USE YOPHNC   , ONLY : YRPHNC

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NCLV, NGPTOT,  NGPTOTG, NGPBLKS
    REAL(KIND=JPRB)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PQ(:,:,:)    ! Q at start of callpar
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_B_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_B_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PLUDE(:,:,:) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PMFD(:,:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PCLV(:,:,:,:) 
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PSUPSAT(:,:,:)
    REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PCOVPTOT(:,:,:) ! Precip fraction
    INTEGER(KIND=JPIM) ::  JK
    TYPE(TOMCST)   , INTENT(INOUT) :: OUTPUT_YDMCST 
    TYPE(TOETHF)   , INTENT(INOUT) :: OUTPUT_YDOETHF 
    TYPE(TECLD)    , INTENT(INOUT) :: OUTPUT_YDECLD 
    TYPE(TECLDP)   , INTENT(INOUT) :: OUTPUT_YDECLDP
    TYPE(TEPHLI)   , INTENT(INOUT) :: OUTPUT_YDEPHLI
    TYPE(TPHNC)    , INTENT(INOUT) :: OUTPUT_YDPHNC

     print *, 'TEST 0'
    CALL GLOBAL_STATE%LOAD(NPROMA, NGPTOT, NGPTOTG)

    !cloudsc2
    !set up other modules not initialized in cloudsc
    !YRECLD
    !print *,'Test i3'
    allocate(YRECLD)
!   !print *,'Test i4'
    allocate(YRECLD%CETA(GLOBAL_STATE%KLEV))
!   ! security
!   if (GLOBAL_STATE%KLEV>200) then
!     print *, ' Dimension of ZPRES/ZPRESF is too short. '
!     stop
!   endif
    DO JK=1,GLOBAL_STATE%KLEV
      YRECLD%CETA(JK)= GLOBAL_STATE%PAP(1,JK,1)/GLOBAL_STATE%PAPH(1,GLOBAL_STATE%KLEV+1,1)
    ENDDO
!   ! YRPHNC
    allocate(YRPHNC)
    YRPHNC%LEVAPLS2=.false.
!   ! overload LPHYLIN
    YREPHLI%LPHYLIN=.true.
     print *, 'TEST 1'
    OUTPUT_PT      (:,:,:) = GLOBAL_STATE%PT      (:,:,:)
    OUTPUT_PQ      (:,:,:) = GLOBAL_STATE%PQ      (:,:,:)
    OUTPUT_B_CML (:,:,:,:) = GLOBAL_STATE%B_CML   (:,:,:,:)
    OUTPUT_B_LOC (:,:,:,:) = GLOBAL_STATE%B_LOC   (:,:,:,:)
    OUTPUT_PAP     (:,:,:) = GLOBAL_STATE%PAP     (:,:,:)
    OUTPUT_PAPH    (:,:,:) = GLOBAL_STATE%PAPH    (:,:,:)
    OUTPUT_PLU     (:,:,:) = GLOBAL_STATE%PLU     (:,:,:)
    OUTPUT_PLUDE   (:,:,:) = GLOBAL_STATE%PLUDE   (:,:,:)
    OUTPUT_PMFU    (:,:,:) = GLOBAL_STATE%PMFU    (:,:,:)
    OUTPUT_PMFD    (:,:,:) = GLOBAL_STATE%PMFD    (:,:,:)
    OUTPUT_PA      (:,:,:) = GLOBAL_STATE%PA      (:,:,:)
    OUTPUT_PCLV  (:,:,:,:) = GLOBAL_STATE%PCLV    (:,:,:,:)
    OUTPUT_PSUPSAT (:,:,:) = GLOBAL_STATE%PSUPSAT (:,:,:)
    OUTPUT_PCOVPTOT(:,:,:) = GLOBAL_STATE%PCOVPTOT(:,:,:)

!    allocate(OUTPUT_YDECLD%CETA(GLOBAL_STATE%KLEV))
     print *, 'TEST 1'
    OUTPUT_YDECLDP%RCLCRIT =YRECLDP%RCLCRIT
     print *, 'TEST 2'
    OUTPUT_YDECLDP%RKCONV  =YRECLDP%RKCONV
    OUTPUT_YDECLDP%RLMIN  =YRECLDP%RLMIN
    OUTPUT_YDECLDP%RPECONS=YRECLDP%RPECONS
    OUTPUT_YDEPHLI%LPHYLIN=YREPHLI%LPHYLIN
    OUTPUT_YDEPHLI%RLPTRC =YREPHLI%RLPTRC
    OUTPUT_YDPHNC%LEVAPLS2=YRPHNC%LEVAPLS2
    OUTPUT_YDOETHF%RVTMP2 = RVTMP2
    OUTPUT_YDEPHLI%LPHYLIN=.true.
    OUTPUT_YDPHNC%LEVAPLS2=.false.
     print *, 'TEST 3'
    OUTPUT_YDECLD%CETA(1:GLOBAL_STATE%KLEV)=YRECLD%CETA(1:GLOBAL_STATE%KLEV)
     print *, 'TEST 4'
    OUTPUT_YDMCST%RG   =RG
    OUTPUT_YDMCST%RG   =RD
    OUTPUT_YDMCST%RCPD =RCPD
    OUTPUT_YDMCST%RETV =RETV
    OUTPUT_YDMCST%RLVTT=RLVTT
    OUTPUT_YDMCST%RLSTT=RLSTT
    OUTPUT_YDMCST%RLMLT=RLMLT 
    OUTPUT_YDMCST%RTT  =RTT
    OUTPUT_YDMCST%RV   =RV 
    OUTPUT_YDMCST%RA   =RA
    OUTPUT_YDMCST%RPI  =RPI

    OUTPUT_YDOETHF%R2ES=R2ES
    OUTPUT_YDOETHF%R3LES=R3LES
    OUTPUT_YDOETHF%R3IES=R3IES
    OUTPUT_YDOETHF%R4LES=R4LES
    OUTPUT_YDOETHF%R4IES=R4IES
    OUTPUT_YDOETHF%R5LES=R5LES
    OUTPUT_YDOETHF%R5IES=R5IES
    OUTPUT_YDOETHF%R5ALVCP=R5ALVCP
    OUTPUT_YDOETHF%R5ALSCP=R5ALSCP
    OUTPUT_YDOETHF%RALVDCP=RALVDCP
    OUTPUT_YDOETHF%RALSDCP=RALSDCP
    OUTPUT_YDOETHF%RALFDCP=RALFDCP
    OUTPUT_YDOETHF%RTWAT=RTWAT
    OUTPUT_YDOETHF%RTICE=RTICE
    OUTPUT_YDOETHF%RTICECU=RTICECU
    OUTPUT_YDOETHF%RTWAT_RTICE_R=RTWAT_RTICE_R
    OUTPUT_YDOETHF%RTWAT_RTICECU_R=RTWAT_RTICECU_R
    OUTPUT_YDOETHF%RKOOP1=RKOOP1
    OUTPUT_YDOETHF%RKOOP2=RKOOP2

     print *, 'Finishing driver init '
  END SUBROUTINE CLOUDSC_DRIVER_INIT

  SUBROUTINE CLOUDSC_DRIVER_NO_DERV_TPES( &
     & YDOMCST, YDOETHF, YDECLD, & 
     & YDECLDP, YDEPHLI, YDPHNC, & 
     & NUMOMP, NPROMA, NLEV,   & 
     & NCLV, NCLDQL, NCLDQI, NGPTOT, NGPBLKS,  NGPTOTG, PTSPHY, &
     & PT, PQ, &
     & BUFFER_CML, BUFFER_LOC, &
     & PAP,      PAPH, &
     & PLU,      PLUDE,    PMFU,     PMFD, &
     & PA,       PCLV,     PSUPSAT,&
     & PCOVPTOT, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN )
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
    REAL(KIND=JPRB),    INTENT(INOUT)    :: PA(:,:,:)    ! Original Cloud fraction (t)
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
!   if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
!   end if

    ! Global timer for the parallel region
!   CALL TIMER%START(NUMOMP)

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID) &
    !$omp& private(ZQSAT) &
    !$omp& num_threads(NUMOMP)

    ! Local timer for each thread
!   TID = GET_THREAD_NUM()
!   CALL TIMER%THREAD_START(TID)
    print *, 'Cloudsc compute: before calling fcttre  set '
    CALL INITIALIZE_FCTTRE_PARAMETERS(YDOMCST, YDOETHF, YDECLD, YDECLDP, YDEPHLI, YDPHNC )
    print *, 'Cloudsc compute: after calling fcttre  set '

    !$omp do schedule(runtime)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         PCOVPTOT(:,:,IBL) = 0.0_JPRB
!         TENDENCY_LOC(IBL)%cld(:,:,NCLV) = 0.0_JPRB

         ! Fill in ZQSAT
         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:), 2) 

    print *, 'Cloudsc compute: befoe calling cloudsc2 '
         CALL CLOUDSC2 ( &
              &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
              & PTSPHY,&
              & PAPH(:,:,IBL),  PAP(:,:,IBL), &
              & PQ(:,:,IBL), ZQSAT(:,:), PT(:,:,IBL), &
              & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
              & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
!             &  TENDENCY_LOC(IBL)%T, TENDENCY_CML(IBL)%T, &
!             &  TENDENCY_LOC(IBL)%Q, TENDENCY_CML(IBL)%Q, &
!             &  TENDENCY_LOC(IBL)%CLD(:,:,NCLDQL), TENDENCY_CML(IBL)%CLD(:,:,NCLDQL), &
!             &  TENDENCY_LOC(IBL)%CLD(:,:,NCLDQI), TENDENCY_CML(IBL)%CLD(:,:,NCLDQI), &
              &  BUFFER_LOC(:,:,1,IBL), BUFFER_CML(:,:,1,IBL), &
              &  BUFFER_LOC(:,:,3,IBL), BUFFER_CML(:,:,3,IBL), &
              &  BUFFER_LOC(:,:,3+NCLDQL,IBL), BUFFER_CML(:,:,3+NCLDQL,IBL),  &
              &  BUFFER_LOC(:,:,3+NCLDQI,IBL), BUFFER_CML(:,:,3+NCLDQI,IBL),  &
              &  PSUPSAT(:,:,IBL), &
              &  PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
              &  PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL), &
              &  YDOMCST, YDOETHF, YDECLD, YDECLDP, YDEPHLI, YDPHNC )
    print *, 'After compute: befoe calling cloudsc2 '

         ! Log number of columns processed by this thread
!         CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
      ENDDO

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc
      !$omp end do nowait

!      CALL TIMER%THREAD_END(TID)

      !$omp end parallel

!     CALL TIMER%END()

!      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, ZHPM, NGPTOT)
    
  CONTAINS
SUBROUTINE INITIALIZE_FCTTRE_PARAMETERS(YDOMCST, YDOETHF, YDECLD, YDECLDP, YDEPHLI, YDPHNC )
USE FCTTRE_MOD, ONLY : FCTTRE_CONSTANTS_SET
!USE YOMCST   , ONLY : RETV, RG  ,RCPD     ,&
! & RLVTT    ,RLSTT    ,RLMLT    ,RTT     ,RD, &
!   RA, RPI, RV
!USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
! & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
! & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU   ,&
! & RTWAT_RTICE_R      ,RTWAT_RTICECU_R    ,RVTMP2, &
! & RALFDCP, RKOOP1, RKOOP2
!USE YOECLD   , ONLY : YRECLD
!USE YOECLDP  , ONLY : YRECLDP
!USE YOEPHLI  , ONLY : YREPHLI
!USE YOPHNC   , ONLY : YRPHNC
!USE YOETHF   , ONLY : RVTMP2
TYPE(TOMCST)      ,INTENT(IN) :: YDOMCST
TYPE(TOETHF)      ,INTENT(IN) :: YDOETHF
TYPE(TECLD)       ,INTENT(IN) :: YDECLD
TYPE(TECLDP)      ,INTENT(IN) :: YDECLDP
TYPE(TEPHLI)      ,INTENT(IN) :: YDEPHLI
TYPE(TPHNC)       ,INTENT(IN) :: YDPHNC
!YDECLDP%RCLCRIT =YRECLDP%RCLCRIT
!YDECLDP%RKCONV  =YRECLDP%RKCONV
!YDECLDP%RLMIN  =YRECLDP%RLMIN
!YDECLDP%RPECONS=YRECLDP%RPECONS
!YDEPHLI%LPHYLIN=YREPHLI%LPHYLIN
!YDEPHLI%RLPTRC =YREPHLI%RLPTRC
!YDPHNC%LEVAPLS2=YRPHNC%LEVAPLS2
!YDOETHF%RVTMP2 = RVTMP2
!YDEPHLI%LPHYLIN=.true.
!YDPHNC%LEVAPLS2=.false.
!YDECLD%CETA(:)=YRECLD%CETA(:)
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

  SUBROUTINE CLOUDSC_DRIVER_VALIDATE( &
     & NUMOMP, NPROMA, NLEV, NCLV, NGPTOT, NGPBLKS,  NGPTOTG, PTSPHY, &
     & INPUT_PT, INPUT_PQ, &
     & INPUT_BUFFER_CML, INPUT_BUFFER_LOC, &
     & INPUT_PAP,      INPUT_PAPH, &
     & INPUT_PLU,      INPUT_PLUDE,    INPUT_PMFU,     INPUT_PMFD, &
     & INPUT_PA,       INPUT_PCLV,     INPUT_PSUPSAT,&
     & INPUT_PCOVPTOT, &
     & INPUT_PFPLSL,   INPUT_PFPLSN,   INPUT_PFHPSL,   INPUT_PFHPSN &
     & )

    INTEGER(KIND=JPIM), INTENT(IN) :: NUMOMP, NPROMA, NLEV, NCLV, NGPTOT, NGPBLKS, NGPTOTG
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PQ(:,:,:)    ! Q at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PLUDE(:,:,:) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PMFD(:,:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PCLV(:,:,:,:) 
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PSUPSAT(:,:,:)
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PCOVPTOT(:,:,:) ! Precip fraction
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PFPLSL(:,:,:)
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PFPLSN(:,:,:)
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PFHPSL(:,:,:)
    REAL(KIND=JPRB),    INTENT(IN)    :: INPUT_PFHPSN(:,:,:)

    GLOBAL_STATE%PT      (:,:,:) = INPUT_PT        (:,:,:)
    GLOBAL_STATE%PQ      (:,:,:) = INPUT_PQ        (:,:,:)
    GLOBAL_STATE%B_CML (:,:,:,:) = INPUT_BUFFER_CML(:,:,:,:)
    GLOBAL_STATE%B_LOC (:,:,:,:) = INPUT_BUFFER_LOC(:,:,:,:)
    GLOBAL_STATE%PAP     (:,:,:) = INPUT_PAP       (:,:,:)
    GLOBAL_STATE%PAPH    (:,:,:) = INPUT_PAPH      (:,:,:)
    GLOBAL_STATE%PLU     (:,:,:) = INPUT_PLU       (:,:,:)
    GLOBAL_STATE%PLUDE   (:,:,:) = INPUT_PLUDE     (:,:,:)
    GLOBAL_STATE%PMFU    (:,:,:) = INPUT_PMFU      (:,:,:)
    GLOBAL_STATE%PMFD    (:,:,:) = INPUT_PMFD      (:,:,:)
    GLOBAL_STATE%PA      (:,:,:) = INPUT_PA        (:,:,:)
    GLOBAL_STATE%PCLV  (:,:,:,:) = INPUT_PCLV      (:,:,:,:)
    GLOBAL_STATE%PSUPSAT (:,:,:) = INPUT_PSUPSAT   (:,:,:)
    GLOBAL_STATE%PCOVPTOT(:,:,:) = INPUT_PCOVPTOT  (:,:,:)
    GLOBAL_STATE%PFPLSL  (:,:,:) = INPUT_PFPLSL    (:,:,:)
    GLOBAL_STATE%PFPLSN  (:,:,:) = INPUT_PFPLSN    (:,:,:)
    GLOBAL_STATE%PFHPSL  (:,:,:) = INPUT_PFHPSL    (:,:,:)
    GLOBAL_STATE%PFHPSN  (:,:,:) = INPUT_PFHPSN    (:,:,:)

    CALL GLOBAL_STATE%VALIDATE(NPROMA, NGPTOT, NGPTOTG)

  END SUBROUTINE CLOUDSC_DRIVER_VALIDATE

! SUBROUTINE CLOUDSC_DRIVER_INIT_COMPUTE( &
!    & NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS,  NGPTOTG, PTSPHY, &
!    & OUTPUT_PT      , OUTPUT_PQ, &
!    & OUTPUT_B_CML   , OUTPUT_B_LOC, &
!    & OUTPUT_PAP     , OUTPUT_PAPH, &
!    & OUTPUT_PLU     , OUTPUT_PLUDE,    OUTPUT_PMFU,     OUTPUT_PMFD, &
!    & OUTPUT_PA      , OUTPUT_PCLV,     OUTPUT_PSUPSAT,&
!    & OUTPUT_PCOVPTOT, & 
!    & OUTPUT_PFPLSL  , OUTPUT_PFPLSN,   OUTPUT_PFHPSL,   OUTPUT_PFHPSN &
!    & )
!   ! Driver routine that performans the parallel NPROMA-blocking and
!   ! invokes the CLOUDSC2 kernel

!   INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, NGPTOTG
!   REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PT(:,:,:)    ! T at start of callpar
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PQ(:,:,:)    ! Q at start of callpar
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_B_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_CML
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_B_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! Storage buffer for TENDENCY_LOC
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PAP(:,:,:)   ! Pressure on full levels
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PAPH(:,:,:)  ! Pressure on half levels
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PLU(:,:,:)   ! Conv. condensate
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PLUDE(:,:,:) ! Conv. detrained water
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PMFU(:,:,:)  ! Conv. mass flux up
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PMFD(:,:,:)  ! Conv. mass flux down
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PA(:,:,:)    ! Original Cloud fraction (t)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PCLV(:,:,:,:) 
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PSUPSAT(:,:,:)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PCOVPTOT(:,:,:) ! Precip fraction
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFPLSL(:,:,:)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFPLSN(:,:,:)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFHPSL(:,:,:)
!   REAL(KIND=JPRB),    INTENT(OUT)    :: OUTPUT_PFHPSN(:,:,:)

!   INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND
!   INTEGER(KIND=JPIM) ::  JK

!   TYPE(PERFORMANCE_TIMER) :: TIMER
!   REAL(KIND=JPRD), PARAMETER :: ZHPM = 3996006.0_JPRD  ! The nominal number of flops per 100 columns

!   INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
!   LOGICAL            :: LDRAIN1D = .FALSE.
!   REAL(KIND=JPRB)    :: ZQSAT(NPROMA,NLEV) ! local array

!    CALL GLOBAL_STATE%LOAD(NPROMA, NGPTOT, NGPTOTG)

!    !cloudsc2
!    !set up other modules not initialized in cloudsc
!    !YRECLD
!    allocate(YRECLD)
!    allocate(YRECLD%CETA(GLOBAL_STATE%KLEV))
!    ! security
!    if (GLOBAL_STATE%KLEV>200) then
!     print *, ' Dimension of ZPRES/ZPRESF is too short. '
!     stop
!    endif
!    DO JK=1,GLOBAL_STATE%KLEV
!     YRECLD%CETA(JK)= GLOBAL_STATE%PAP(1,JK,1)/GLOBAL_STATE%PAPH(1,GLOBAL_STATE%KLEV+1,1)
!    ENDDO
!    ! YRPHNC
!    allocate(YRPHNC)
!    YRPHNC%LEVAPLS2=.false.
!    ! overload LPHYLIN
!    YREPHLI%LPHYLIN=.true.
!    CALL CLOUDSC_DRIVER_NO_DERV_TPES(NUMOMP, NPROMA, GLOBAL_STATE%KLEV, NGPTOT, &
!       & GLOBAL_STATE%NBLOCKS , NGPTOTG, GLOBAL_STATE%PTSPHY,                   &
!       & GLOBAL_STATE%PT      , GLOBAL_STATE%PQ,                                &
!       & GLOBAL_STATE%B_CML   , GLOBAL_STATE%B_LOC,                             &
!       & GLOBAL_STATE%PAP     , GLOBAL_STATE%PAPH,                              &
!       & GLOBAL_STATE%PLU     , GLOBAL_STATE%PLUDE,                             &
!       & GLOBAL_STATE%PMFU    , GLOBAL_STATE%PMFD,                              &
!       & GLOBAL_STATE%PA      , GLOBAL_STATE%PCLV,  GLOBAL_STATE%PSUPSAT,       &
!       & GLOBAL_STATE%PCOVPTOT,                                                 &
!       & GLOBAL_STATE%PFPLSL  ,   GLOBAL_STATE%PFPLSN,   GLOBAL_STATE%PFHPSL,   &
!       & GLOBAL_STATE%PFHPSN &
!       & )
!    OUTPUT_PT      (:,:,:)   =  GLOBAL_STATE%PT      (:,:,:)
!    OUTPUT_PQ      (:,:,:)   =  GLOBAL_STATE%PQ      (:,:,:)
!    OUTPUT_B_CML   (:,:,:,:) =  GLOBAL_STATE%B_CML   (:,:,:,:)
!    OUTPUT_B_LOC   (:,:,:,:) =  GLOBAL_STATE%B_LOC   (:,:,:,:)
!    OUTPUT_PAP     (:,:,:)   =  GLOBAL_STATE%PAP     (:,:,:)
!    OUTPUT_PAPH    (:,:,:)   =  GLOBAL_STATE%PAPH    (:,:,:)
!    OUTPUT_PLU     (:,:,:)   =  GLOBAL_STATE%PLU     (:,:,:)
!    OUTPUT_PLUDE   (:,:,:)   =  GLOBAL_STATE%PLUDE   (:,:,:)
!    OUTPUT_PMFU    (:,:,:)   =  GLOBAL_STATE%PMFU    (:,:,:)
!    OUTPUT_PMFD    (:,:,:)   =  GLOBAL_STATE%PMFD    (:,:,:)
!    OUTPUT_PA      (:,:,:)   =  GLOBAL_STATE%PA      (:,:,:)
!    OUTPUT_PCLV    (:,:,:,:) =  GLOBAL_STATE%PCLV    (:,:,:,:)
!    OUTPUT_PSUPSAT (:,:,:)   =  GLOBAL_STATE%PSUPSAT (:,:,:)
!    OUTPUT_PCOVPTOT(:,:,:)   =  GLOBAL_STATE%PCOVPTOT(:,:,:)
!    OUTPUT_PFPLSL  (:,:,:)   =  GLOBAL_STATE%PFPLSL  (:,:,:)
!    OUTPUT_PFPLSN  (:,:,:)   =  GLOBAL_STATE%PFPLSN  (:,:,:)
!    OUTPUT_PFHPSL  (:,:,:)   =  GLOBAL_STATE%PFHPSL  (:,:,:)
!    OUTPUT_PFHPSN  (:,:,:)   =  GLOBAL_STATE%PFHPSN  (:,:,:)
!      
!    CALL GLOBAL_STATE%VALIDATE(NPROMA, NGPTOT, NGPTOTG)
! END SUBROUTINE CLOUDSC_DRIVER_INIT_COMPUTE

!  SUBROUTINE CLOUDSC_DRIVER_INIT_COMPUTE_VALIDATE( &
!             PTSPHY, &
!             NUMOMP, NPROMA, NLEV, NGPTOT, NGPBLKS, & 
!             NGPTOTG)
!     USE PARKIND1, ONLY: JPIM
!     USE YOECLD   , ONLY : YRECLD
!     USE YOPHNC   , ONLY : YRPHNC
!     USE YOEPHLI  , ONLY : YREPHLI
!         INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT,  NGPTOTG, NGPBLKS
!         REAL(KIND=JPRB)    :: PTSPHY       ! Physics timestep
!         INTEGER(KIND=JPIM) ::  JK
!         INTEGER(KIND=JPIM) ::  OUTPUT_NLEV,OUTPUT_NGPBLKS 
!     !Compute the number of blocks based on NGPTOT and NPROMA
!     ! TODO: Create a global global memory state from serialized input data
!     CALL GLOBAL_STATE%LOAD(NPROMA, NGPTOT, NGPTOTG)
!
!     !cloudsc2
!     !set up other modules not initialized in cloudsc
!     !YRECLD
!     allocate(YRECLD)
!     allocate(YRECLD%CETA(GLOBAL_STATE%KLEV))
!     ! security
!     if (GLOBAL_STATE%KLEV>200) then
!       print *, ' Dimension of ZPRES/ZPRESF is too short. '
!       stop
!     endif
!     DO JK=1,GLOBAL_STATE%KLEV
!       YRECLD%CETA(JK)= GLOBAL_STATE%PAP(1,JK,1)/GLOBAL_STATE%PAPH(1,GLOBAL_STATE%KLEV+1,1)
!     ENDDO
!     ! YRPHNC
!     allocate(YRPHNC)
!     YRPHNC%LEVAPLS2=.false.
!     ! overload LPHYLIN
!     YREPHLI%LPHYLIN=.true.
! 
!     !Computation
!     CALL CLOUDSC_DRIVER_NO_DERV_TPES(NUMOMP, NPROMA, &
!          & GLOBAL_STATE%KLEV, NGPTOT, GLOBAL_STATE%NBLOCKS, NGPTOTG, &
!          & GLOBAL_STATE%PTSPHY, &
!          & GLOBAL_STATE%PT, GLOBAL_STATE%PQ, &
!          & GLOBAL_STATE%B_CML,    GLOBAL_STATE%B_LOC, &
!          & GLOBAL_STATE%PAP,      GLOBAL_STATE%PAPH, &
!          & GLOBAL_STATE%PLU,      GLOBAL_STATE%PLUDE, &
!          & GLOBAL_STATE%PMFU,     GLOBAL_STATE%PMFD, &
!          & GLOBAL_STATE%PA,       GLOBAL_STATE%PCLV,  GLOBAL_STATE%PSUPSAT,&
!          & GLOBAL_STATE%PCOVPTOT, &
!          & GLOBAL_STATE%PFPLSL,   GLOBAL_STATE%PFPLSN,   GLOBAL_STATE%PFHPSL, &
!          & GLOBAL_STATE%PFHPSN, &
!          & YOMCST, YOETHF, YOECLD, YOECLDP, YEOPHLI, YEOPHNC )
!     !Validation
!       CALL GLOBAL_STATE%VALIDATE(NPROMA, NGPTOT, NGPTOTG)
!  END SUBROUTINE CLOUDSC_DRIVER_INIT_COMPUTE_VALIDATE
END MODULE CLOUDSC_DRIVER_PYIFACE_MOD
