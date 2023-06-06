! Copyright (C) 2003- ECMWF
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CLOUDSC_DRIVER_AD_MOD
  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, NCLDQL, NCLDQI 
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM
  USE EC_PMON_MOD, ONLY: EC_PMON

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_AD( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NGPBLKS, PTSPHY, &
     & PT, PQ, BUFFER_CML, BUFFER_LOC, &
     & PAP,      PAPH, &
     & PLU,      PLUDE,    PMFU,     PMFD, &
     & PA,       PCLV,     PSUPSAT,&
     & PCOVPTOT, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN, &
     & YDCST, YDTHF, YHNC, YPHLI, YCLD, YCLDP, YNCL, LCETA )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC2 kernel
    USE YOMNCL   , ONLY : TNCL
    USE YOMCST   , ONLY : TOMCST
    USE YOETHF   , ONLY : TOETHF
    USE YOPHNC   , ONLY : TPHNC
    USE YOEPHLI  , ONLY : TEPHLI
    USE YOECLD   , ONLY : TECLD
    USE YOECLDP  , ONLY : TECLDP
    USE CLOUDSC2AD_MOD, ONLY : CLOUDSC2AD
    USE CLOUDSC2TL_MOD, ONLY : CLOUDSC2TL
    USE SATUR_MOD, ONLY : SATUR 

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: LCETA(NLEV)
    REAL(KIND=JPRB),    INTENT(IN)    :: PT(NPROMA,NLEV,NGPBLKS)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: PQ(NPROMA,NLEV,NGPBLKS)    ! Q at start of callpar
    REAL(KIND=JPRB),    INTENT(INOUT) :: BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! TENDENCY_CML storage buffer
    REAL(KIND=JPRB),    INTENT(INOUT) :: BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! TENDENCY_LOC storage buffer
    REAL(KIND=JPRB),    INTENT(IN)    :: PAP(NPROMA,NLEV,NGPBLKS)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PAPH(NPROMA,NLEV+1,NGPBLKS)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PLU(NPROMA,NLEV,NGPBLKS)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(INOUT) :: PLUDE(NPROMA,NLEV,NGPBLKS) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFU(NPROMA,NLEV,NGPBLKS)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFD(NPROMA,NLEV,NGPBLKS)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(INOUT)  :: PA(NPROMA,NLEV,NGPBLKS)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(IN)    :: PCLV(NPROMA,NLEV,NCLV,NGPBLKS)
    REAL(KIND=JPRB),    INTENT(IN)    :: PSUPSAT(NPROMA,NLEV,NGPBLKS)
    REAL(KIND=JPRB),    INTENT(INOUT) :: PCOVPTOT(NPROMA,NLEV,NGPBLKS) ! Precip fraction
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSL(NPROMA,NLEV+1,NGPBLKS) ! liq+rain sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSN(NPROMA,NLEV+1,NGPBLKS) ! ice+snow sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSL(NPROMA,NLEV+1,NGPBLKS) ! Enthalpy flux for liq
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSN(NPROMA,NLEV+1,NGPBLKS) ! Enthalp flux for ice

    INTEGER(KIND=JPIM) :: JKGLO,IBL,JROF,ICEND,NGPBLKS

    TYPE(PERFORMANCE_TIMER) :: TIMER
    REAL(KIND=JPRD), PARAMETER :: ZHPM = 3996006.0_JPRD  ! The nominal number of flops per 100 columns

    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    LOGICAL            :: LDRAIN1D = .FALSE.
    REAL(KIND=JPRB)    :: ZQSAT(NPROMA,NLEV,NGPBLKS) ! local array

    REAL(KIND=JPRB)    :: ZNORMG
    REAL(KIND=JPRB)    :: ZNORM1(NPROMA), ZNORM2(NPROMA), ZNORM3(NPROMA)
    REAL(KIND=JPRB)    ::   ZAPH(NPROMA,NLEV+1,NGPBLKS), &
     &                       ZAP(NPROMA,NLEV  ,NGPBLKS), &
     &                        ZQ(NPROMA,NLEV  ,NGPBLKS), &
     &                    ZZQSAT(NPROMA,NLEV  ,NGPBLKS), &
     &                        ZT(NPROMA,NLEV  ,NGPBLKS), &
     &                       ZLU(NPROMA,NLEV  ,NGPBLKS), &
     &                      ZMFU(NPROMA,NLEV  ,NGPBLKS), &
     &                      ZMFD(NPROMA,NLEV  ,NGPBLKS), &
     &                   ZTENI_T(NPROMA,NLEV  ,NGPBLKS), &
     &                   ZTENI_Q(NPROMA,NLEV  ,NGPBLKS), &
     &                   ZTENI_L(NPROMA,NLEV  ,NGPBLKS), &
     &                   ZTENI_I(NPROMA,NLEV  ,NGPBLKS), &
     &                   ZSUPSAT(NPROMA,NLEV  ,NGPBLKS)
    REAL(KIND=JPRB)    ::ZTENO_T(NPROMA,NLEV  ,NGPBLKS), &
     &                   ZTENO_Q(NPROMA,NLEV  ,NGPBLKS), &
     &                   ZTENO_L(NPROMA,NLEV  ,NGPBLKS), & 
     &                   ZTENO_I(NPROMA,NLEV  ,NGPBLKS), &
     &                      ZCLC(NPROMA,NLEV  ,NGPBLKS), &
     &                    ZFPLSL(NPROMA,NLEV+1,NGPBLKS), &
     &                    ZFPLSN(NPROMA,NLEV+1,NGPBLKS), &
     &                    ZFHPSL(NPROMA,NLEV+1,NGPBLKS), &
     &                    ZFHPSN(NPROMA,NLEV+1,NGPBLKS), &
     &                  ZCOVPTOT(NPROMA,NLEV  ,NGPBLKS)
    REAL(KIND=JPRB)    ::  ZAPH0(NPROMA,NLEV+1,NGPBLKS), &
     &                      ZAP0(NPROMA,NLEV  ,NGPBLKS), &
     &                       ZQ0(NPROMA,NLEV  ,NGPBLKS), &
     &                   ZZQSAT0(NPROMA,NLEV  ,NGPBLKS), &
     &                       ZT0(NPROMA,NLEV  ,NGPBLKS), &
     &                      ZLU0(NPROMA,NLEV  ,NGPBLKS), &
     &                     ZMFU0(NPROMA,NLEV  ,NGPBLKS), &
     &                     ZMFD0(NPROMA,NLEV  ,NGPBLKS), &
     &                  ZTENI_T0(NPROMA,NLEV  ,NGPBLKS), &
     &                  ZTENI_Q0(NPROMA,NLEV  ,NGPBLKS), &
     &                  ZTENI_L0(NPROMA,NLEV  ,NGPBLKS), &
     &                  ZTENI_I0(NPROMA,NLEV  ,NGPBLKS), &
     &                  ZSUPSAT0(NPROMA,NLEV  ,NGPBLKS)
    REAL(KIND=JPRB)    ::     ZDRVL(NPROMA,NLEV  ,NGPBLKS), &
     &                       ZL0(NPROMA,NLEV  ,NGPBLKS), &
     &                        ZDRVI(NPROMA,NLEV  ,NGPBLKS), &
     &                       ZI0(NPROMA,NLEV  ,NGPBLKS), &
     &                     ZDRVLUDE(NPROMA,NLEV  ,NGPBLKS), &
     &                    ZLUDE0(NPROMA,NLEV  ,NGPBLKS)
    TYPE(TNCL)      :: YNCL
    TYPE(TOMCST)    :: YDCST
    TYPE(TOETHF)    :: YDTHF
    TYPE(TPHNC)     :: YHNC
    TYPE(TEPHLI)    :: YPHLI
    TYPE(TECLD)     :: YCLD
    TYPE(TECLDP)    :: YCLDP

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)
    !$loki data 
    ZNORMG=0._JPRB

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         PCOVPTOT(:,:,IBL) = 0.0_JPRB

         ! Fill in ZQSAT
         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:,IBL), 2, YDCST, YDTHF) 


         ! Preparation for TL

         ! Increments (IN)
         ZAPH(:,:,IBL)    =PAPH(:,:,IBL)*0.01_JPRB
         ZAP(:,:,IBL)     =PAP(:,:,IBL)*0.01_JPRB
         ZQ(:,:,IBL)      =PQ(:,:,IBL)*0.01_JPRB
         ZZQSAT(:,:,IBL)  =ZQSAT(:,:,IBL)     *0.01_JPRB
         ZT(:,:,IBL)      = PT(:,:,IBL)*0.01_JPRB
         ZDRVL(:,:,IBL)      = PCLV(:,:,NCLDQL,IBL)*0.01_JPRB
         ZDRVI(:,:,IBL)      = PCLV(:,:,NCLDQI,IBL)*0.01_JPRB
         ZDRVLUDE(:,:,IBL)   = PLUDE(:,:,IBL)*0.01_JPRB
         ZLU(:,:,IBL)     = PLU(:,:,IBL)*0.01_JPRB
         ZMFU(:,:,IBL)    = PMFU(:,:,IBL)*0.01_JPRB
         ZMFD(:,:,IBL)    = PMFD(:,:,IBL)*0.01_JPRB
         ZTENI_T(:,:,IBL) = BUFFER_CML(:,:,1,IBL)*0.01_JPRB
         ZTENI_Q(:,:,IBL) = BUFFER_CML(:,:,3,IBL)*0.01_JPRB
         ZTENI_L(:,:,IBL) = BUFFER_CML(:,:,3+NCLDQL,IBL)*0.01_JPRB
         ZTENI_I(:,:,IBL) = BUFFER_CML(:,:,3+NCLDQI,IBL)*0.01_JPRB
         ZSUPSAT(:,:,IBL) = 0.00_JPRB  ! obsolette, beter not use

         ! Storage of initial increments
            ZAPH0(:,:,IBL) = ZAPH(:,:,IBL)
             ZAP0(:,:,IBL) = ZAP(:,:,IBL)
              ZQ0(:,:,IBL) = ZQ(:,:,IBL)
          ZZQSAT0(:,:,IBL) = ZZQSAT(:,:,IBL)
              ZT0(:,:,IBL) = ZT(:,:,IBL)
              ZL0(:,:,IBL) = ZDRVL(:,:,IBL)
              ZI0(:,:,IBL) = ZDRVI(:,:,IBL)
           ZLUDE0(:,:,IBL) = ZDRVLUDE(:,:,IBL)
             ZLU0(:,:,IBL) = ZLU(:,:,IBL)
            ZMFU0(:,:,IBL) = ZMFU(:,:,IBL)
            ZMFD0(:,:,IBL) = ZMFD(:,:,IBL)
         ZTENI_T0(:,:,IBL) = ZTENI_T(:,:,IBL)
         ZTENI_Q0(:,:,IBL) = ZTENI_Q(:,:,IBL)
         ZTENI_L0(:,:,IBL) = ZTENI_L(:,:,IBL)
         ZTENI_I0(:,:,IBL) = ZTENI_I(:,:,IBL)
         ZSUPSAT0(:,:,IBL) = ZSUPSAT(:,:,IBL)

         ! Tangent linear integration
         CALL  CLOUDSC2TL ( &
            &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
            & PTSPHY,LCETA, &
            ! trajectory
            & PAPH(:,:,IBL),  PAP(:,:,IBL), &
            & PQ(:,:,IBL), ZQSAT(:,:,IBL), PT(:,:,IBL), &
            & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
            & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
            & BUFFER_LOC(:,:,1,IBL), BUFFER_CML(:,:,1,IBL), &
            & BUFFER_LOC(:,:,3,IBL), BUFFER_CML(:,:,3,IBL), &
            & BUFFER_LOC(:,:,3+NCLDQL,IBL), BUFFER_CML(:,:,3+NCLDQL,IBL), &
            & BUFFER_LOC(:,:,3+NCLDQI,IBL), BUFFER_CML(:,:,3+NCLDQI,IBL), &
            & PSUPSAT(:,:,IBL), &
            & PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
            & PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL), &
            ! increments
            & ZAPH(:,:,IBL), ZAP(:,:,IBL), ZQ(:,:,IBL), ZZQSAT(:,:,IBL), ZT(:,:,IBL), ZDRVL(:,:,IBL), ZDRVI(:,:,IBL), &
            & ZDRVLUDE(:,:,IBL), ZLU(:,:,IBL), ZMFU(:,:,IBL), ZMFD(:,:,IBL), &
            & ZTENO_T(:,:,IBL), ZTENI_T(:,:,IBL), ZTENO_Q(:,:,IBL), ZTENI_Q(:,:,IBL), &   ! o,i,o,i
            & ZTENO_L(:,:,IBL), ZTENI_L(:,:,IBL), ZTENO_I(:,:,IBL), ZTENI_I(:,:,IBL), ZSUPSAT(:,:,IBL), &  ! o,i,o,i
            & ZCLC(:,:,IBL)   , ZFPLSL(:,:,IBL)   , ZFPLSN(:,:,IBL) ,&        ! o
            & ZFHPSL(:,:,IBL) , ZFHPSN(:,:,IBL)   , ZCOVPTOT(:,:,IBL), &
            & YDCST, YDTHF, YHNC, YPHLI, YCLD, YCLDP, YNCL )       ! o

         ! First norm
         DO JROF=1,ICEND
           ZNORM1(JROF)=SUM( ZTENO_T(JROF,1:NLEV  ,IBL)* ZTENO_T(JROF,1:NLEV  ,IBL)) &
           &          + SUM( ZTENO_Q(JROF,1:NLEV  ,IBL)* ZTENO_Q(JROF,1:NLEV  ,IBL)) &
           &          + SUM( ZTENO_L(JROF,1:NLEV  ,IBL)* ZTENO_L(JROF,1:NLEV  ,IBL)) &
           &          + SUM( ZTENO_I(JROF,1:NLEV  ,IBL)* ZTENO_I(JROF,1:NLEV  ,IBL)) &
           &          + SUM(    ZCLC(JROF,1:NLEV  ,IBL)*    ZCLC(JROF,1:NLEV  ,IBL)) &
           &          + SUM(  ZFPLSL(JROF,1:NLEV+1,IBL)*  ZFPLSL(JROF,1:NLEV+1,IBL)) &
           &          + SUM(  ZFPLSN(JROF,1:NLEV+1,IBL)*  ZFPLSN(JROF,1:NLEV+1,IBL)) &
           &          + SUM(  ZFHPSL(JROF,1:NLEV+1,IBL)*  ZFHPSL(JROF,1:NLEV+1,IBL)) &
           &          + SUM(  ZFHPSN(JROF,1:NLEV+1,IBL)*  ZFHPSN(JROF,1:NLEV+1,IBL)) &
           &          + SUM(ZCOVPTOT(JROF,1:NLEV  ,IBL)*ZCOVPTOT(JROF,1:NLEV  ,IBL))
         ENDDO

         ! Initiaslization of output variables
             ZAPH(:,:,IBL) = 0.0_JPRB
              ZAP(:,:,IBL) = 0.0_JPRB
               ZQ(:,:,IBL) = 0.0_JPRB
           ZZQSAT(:,:,IBL) = 0.0_JPRB
               ZT(:,:,IBL) = 0.0_JPRB
               ZDRVL(:,:,IBL) = 0.0_JPRB
               ZDRVI(:,:,IBL) = 0.0_JPRB
            ZDRVLUDE(:,:,IBL) = 0.0_JPRB
              ZLU(:,:,IBL) = 0.0_JPRB
             ZMFU(:,:,IBL) = 0.0_JPRB
             ZMFD(:,:,IBL) = 0.0_JPRB
         ZTENI_T (:,:,IBL) = 0.0_JPRB
         ZTENI_Q (:,:,IBL) = 0.0_JPRB
         ZTENI_L (:,:,IBL) = 0.0_JPRB
         ZTENI_I (:,:,IBL) = 0.0_JPRB
         ZSUPSAT (:,:,IBL) = 0.0_JPRB

         ! Adjoint integration
         CALL  CLOUDSC2AD ( &
            &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
            & PTSPHY, LCETA,&
            ! trajectory
            & PAPH(:,:,IBL),  PAP(:,:,IBL), &
            & PQ(:,:,IBL), ZQSAT(:,:,IBL), PT(:,:,IBL), &
            & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
            & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
            & BUFFER_LOC(:,:,1,IBL), BUFFER_CML(:,:,1,IBL), &
            & BUFFER_LOC(:,:,3,IBL), BUFFER_CML(:,:,3,IBL), &
            & BUFFER_LOC(:,:,3+NCLDQL,IBL), BUFFER_CML(:,:,3+NCLDQL,IBL), &
            & BUFFER_LOC(:,:,3+NCLDQI,IBL), BUFFER_CML(:,:,3+NCLDQI,IBL), &
            & PSUPSAT(:,:,IBL), &
            & PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
            & PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL), &
            ! increments
            & ZAPH(:,:,IBL), ZAP(:,:,IBL), ZQ(:,:,IBL), ZZQSAT(:,:,IBL), ZT(:,:,IBL), ZDRVL(:,:,IBL), ZDRVI(:,:,IBL), &
            & ZDRVLUDE(:,:,IBL), ZLU(:,:,IBL), ZMFU(:,:,IBL), ZMFD(:,:,IBL), &
            & ZTENO_T(:,:,IBL), ZTENI_T(:,:,IBL), ZTENO_Q(:,:,IBL), ZTENI_Q(:,:,IBL), &   ! o,i,o,i
            & ZTENO_L(:,:,IBL), ZTENI_L(:,:,IBL), ZTENO_I(:,:,IBL), ZTENI_I(:,:,IBL), ZSUPSAT(:,:,IBL), &  ! o,i,o,i
            & ZCLC(:,:,IBL)   , ZFPLSL(:,:,IBL)   , ZFPLSN(:,:,IBL) ,&        ! o
            & ZFHPSL(:,:,IBL) , ZFHPSN(:,:,IBL)   , ZCOVPTOT(:,:,IBL), &
            & YDCST, YDTHF, YHNC, YPHLI, YCLD, YCLDP, YNCL)       ! o

         ! Second norm
         DO JROF=1,ICEND
           ZNORM2(JROF)=SUM(   ZAPH0(JROF,1:NLEV+1,IBL) *    ZAPH(JROF,1:NLEV+1,IBL)) &
           &          + SUM(    ZAP0(JROF,1:NLEV  ,IBL) *     ZAP(JROF,1:NLEV  ,IBL)) &
           &          + SUM(     ZQ0(JROF,1:NLEV  ,IBL) *      ZQ(JROF,1:NLEV  ,IBL)) &
           &          + SUM( ZZQSAT0(JROF,1:NLEV  ,IBL) *  ZZQSAT(JROF,1:NLEV  ,IBL)) &
           &          + SUM(     ZT0(JROF,1:NLEV  ,IBL) *      ZT(JROF,1:NLEV  ,IBL)) &
           &          + SUM(     ZL0(JROF,1:NLEV  ,IBL) *   ZDRVL(JROF,1:NLEV  ,IBL)) &
           &          + SUM(     ZI0(JROF,1:NLEV  ,IBL) *   ZDRVI(JROF,1:NLEV  ,IBL)) &
           &          + SUM(  ZLUDE0(JROF,1:NLEV  ,IBL) *ZDRVLUDE(JROF,1:NLEV  ,IBL)) &
           &          + SUM(    ZLU0(JROF,1:NLEV  ,IBL) *     ZLU(JROF,1:NLEV  ,IBL)) &
           &          + SUM(   ZMFU0(JROF,1:NLEV  ,IBL) *    ZMFU(JROF,1:NLEV  ,IBL)) &
           &          + SUM(   ZMFD0(JROF,1:NLEV  ,IBL) *    ZMFD(JROF,1:NLEV  ,IBL)) &
           &          + SUM(ZTENI_T0(JROF,1:NLEV  ,IBL) * ZTENI_T(JROF,1:NLEV  ,IBL)) &
           &          + SUM(ZTENI_Q0(JROF,1:NLEV  ,IBL) * ZTENI_Q(JROF,1:NLEV  ,IBL)) &
           &          + SUM(ZTENI_L0(JROF,1:NLEV  ,IBL) * ZTENI_L(JROF,1:NLEV  ,IBL)) &
           &          + SUM(ZTENI_I0(JROF,1:NLEV  ,IBL) * ZTENI_I(JROF,1:NLEV  ,IBL)) &
           &          + SUM(ZSUPSAT0(JROF,1:NLEV  ,IBL) * ZSUPSAT(JROF,1:NLEV  ,IBL))
           ! Third norm
           ! Note the machine precision is defined here as strictly 64bits
           ! as we assume at worst 12 digits agreements in norms.
           IF (ZNORM2(JROF) == 0._JPRB ) THEN
             ZNORM3(JROF)=ABS(ZNORM1(JROF)-ZNORM2(JROF))/EPSILON(1._8)
           ELSE
             ZNORM3(JROF)=ABS(ZNORM1(JROF)-ZNORM2(JROF))/EPSILON(1._8)/ZNORM2(JROF)
           ENDIF
         ENDDO

         ZNORMG=MAX(ZNORMG,MAXVAL(ZNORM3(1:ICEND)))

#ifndef CLOUDSC_GPU_TIMING
         ! Log number of columns processed by this thread
         CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
#endif
      ENDDO

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc

      CALL TIMER%THREAD_END(TID)

      !$loki end data 

      CALL TIMER%END()

#ifdef CLOUDSC_GPU_TIMING
      ! On GPUs, adding block-level column totals is cumbersome and
      ! error prone, and of little value due to the large number of
      ! processing "thread teams". Instead we register the total here.
      CALL TIMER % THREAD_LOG(TID=TID, IGPC=NGPTOT)
#endif

      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, ZHPM, NGPTOT)

      ! Print final test results
      print *, ' AD TEST '
      print *, ' The maximum error is ',ZNORMG,' times the zero of the machine. '
      print *, '   =============================  '
      IF (ZNORMG < 10000._JPRB) THEN
        print *, '   =           TEST OK         = '
      ELSE
        print *, '   =        TEST FAILED        = '
      ENDIF
      print *, '   =============================  '
      
    
  END SUBROUTINE CLOUDSC_DRIVER_AD

END MODULE CLOUDSC_DRIVER_AD_MOD
