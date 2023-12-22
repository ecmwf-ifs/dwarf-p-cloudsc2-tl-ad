! Copyright (C) 2003- ECMWF
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CLOUDSC_DRIVER_TL_MOD
  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, NCLDQL, NCLDQI 
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_TL( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NGPBLKS, PTSPHY, &
     & PT, PQ, BUFFER_CML, BUFFER_LOC, &
     & PAP,      PAPH, &
     & PLU,      PLUDE,    PMFU,     PMFD, &
     & PA,       PCLV,     PSUPSAT,&
     & PCOVPTOT, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN,  &
     & YDCST, YDTHF, YHNC, YPHLI, YCLD, YCLDP, YNCL, LCETA)
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC2 kernel
    USE CLOUDSC2_MOD, ONLY: CLOUDSC2
    USE CLOUDSC2TL_MOD, ONLY: CLOUDSC2TL
    USE SATUR_MOD, ONLY: SATUR 
    USE ERROR_MOD, ONLY: VALIDATE_TAYLOR_TEST
    USE YOMCST   , ONLY : TOMCST
    USE YOETHF   , ONLY : TOETHF
    USE YOPHNC   , ONLY : TPHNC
    USE YOEPHLI  , ONLY : TEPHLI
    USE YOECLD   , ONLY : TECLD
    USE YOECLDP  , ONLY : TECLDP
    USE YOMNCL   , ONLY : TNCL

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NGPBLKS
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: LCETA(NLEV)
    REAL(KIND=JPRB),    INTENT(IN)    :: PT(NPROMA,NLEV,NGPBLKS)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: PQ(NPROMA,NLEV,NGPBLKS)    ! Q at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: BUFFER_CML(NPROMA,NLEV,3+NCLV,NGPBLKS) ! TENDENCY_CML storage buffer
    REAL(KIND=JPRB),    INTENT(OUT)   :: BUFFER_LOC(NPROMA,NLEV,3+NCLV,NGPBLKS) ! TENDENCY_LOC storage buffer
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
    INTEGER(KIND=JPIM) :: JKGLO,IBL,IBLT,ICEND,ISTART,ITEST,INEGAT,ITEMPNEGAT

    TYPE(PERFORMANCE_TIMER) :: TIMER
    REAL(KIND=JPRD), PARAMETER :: ZHPM = 3996006.0_JPRD  ! The nominal number of flops per 100 columns

    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    LOGICAL            :: LDRAIN1D = .FALSE.
    REAL(KIND=JPRB)    :: ZQSAT(NPROMA,NLEV,NGPBLKS) ! local array

    INTEGER(KIND=JPIB) :: ILAM
    REAL(KIND=JPRB)    :: ZLAMBDA, ZCOUNT, ZNORM, ZNORMG(10)
    REAL(KIND=JPRB)    :: ZAPH(NPROMA,NLEV+1,NGPBLKS), ZAP(NPROMA,NLEV,NGPBLKS), ZQ(NPROMA,NLEV,NGPBLKS), &
     & ZZQSAT(NPROMA,NLEV,NGPBLKS), ZT(NPROMA,NLEV,NGPBLKS), &
     & ZLU(NPROMA,NLEV,NGPBLKS), ZMFU(NPROMA,NLEV,NGPBLKS), ZMFD(NPROMA,NLEV,NGPBLKS), &
     & ZTENI_T(NPROMA,NLEV,NGPBLKS), ZTENI_Q(NPROMA,NLEV,NGPBLKS), ZTENI_L(NPROMA,NLEV,NGPBLKS), &
     & ZTENI_I(NPROMA,NLEV,NGPBLKS), ZSUPSAT(NPROMA,NLEV,NGPBLKS)
    REAL(KIND=JPRB)    :: ZTENO_T(NPROMA,NLEV,NGPBLKS), ZTENO_Q(NPROMA,NLEV,NGPBLKS), &
     & ZTENO_L(NPROMA,NLEV,NGPBLKS), ZTENO_I(NPROMA,NLEV,NGPBLKS), ZCLC(NPROMA,NLEV,NGPBLKS), &
     & ZFPLSL(NPROMA,NLEV+1,NGPBLKS), ZFPLSN(NPROMA,NLEV+1,NGPBLKS), ZFHPSL(NPROMA,NLEV+1,NGPBLKS), &
     & ZFHPSN(NPROMA,NLEV+1,NGPBLKS), ZCOVPTOT(NPROMA,NLEV,NGPBLKS)
    REAL(KIND=JPRB)    :: PAPH5(NPROMA,NLEV+1,NGPBLKS), PAP5(NPROMA,NLEV,NGPBLKS), &
     & PQ5(NPROMA,NLEV,NGPBLKS), ZQSAT5(NPROMA,NLEV,NGPBLKS), PT5(NPROMA,NLEV,NGPBLKS), &
     & PCLVL5(NPROMA,NLEV,NGPBLKS), PCLVI5(NPROMA,NLEV,NGPBLKS), PLUDE5(NPROMA,NLEV,NGPBLKS), &
     & PLU5(NPROMA,NLEV,NGPBLKS), PMFU5(NPROMA,NLEV,NGPBLKS), PMFD5(NPROMA,NLEV,NGPBLKS), &
     & ZTENI_T5(NPROMA,NLEV,NGPBLKS), ZTENI_Q5(NPROMA,NLEV,NGPBLKS), ZTENI_L5(NPROMA,NLEV,NGPBLKS), &
     & ZTENI_I5(NPROMA,NLEV,NGPBLKS), PSUPSAT5(NPROMA,NLEV,NGPBLKS)
    INTEGER(KIND=JPIM), PARAMETER :: NLAM=10
    REAL(KIND=JPRB)    :: ZTENO_T5(NPROMA,NLEV,NLAM,NGPBLKS), &
     &                    ZTENO_Q5(NPROMA,NLEV,NLAM,NGPBLKS), &
     &                    ZTENO_L5(NPROMA,NLEV,NLAM,NGPBLKS), &
     &                    ZTENO_I5(NPROMA,NLEV,NLAM,NGPBLKS), &
     &                         PA5(NPROMA,NLEV,NLAM,NGPBLKS), &
     &                    PFPLSL5(NPROMA,NLEV+1,NLAM,NGPBLKS), &
     &                    PFPLSN5(NPROMA,NLEV+1,NLAM,NGPBLKS), &
     &                    PFHPSL5(NPROMA,NLEV+1,NLAM,NGPBLKS), &
     &                    PFHPSN5(NPROMA,NLEV+1,NLAM,NGPBLKS), &
     &                    PCOVPTOT5(NPROMA,NLEV,NLAM,NGPBLKS)
    REAL(KIND=JPRB)    :: ZDRVL(NPROMA,NLEV,NGPBLKS)
    REAL(KIND=JPRB)    :: ZDRVI(NPROMA,NLEV,NGPBLKS)
    REAL(KIND=JPRB)    :: ZDRVLUDE(NPROMA,NLEV,NGPBLKS)
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
    
    ZNORMG(:)=0.

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         PCOVPTOT(:,:,IBL) = 0.0_JPRB
         BUFFER_LOC(:,:,2,IBL) = 0.0_JPRB
         BUFFER_LOC(:,:,4:3+NCLV,IBL) = 0.0_JPRB

         ! Fill in ZQSAT
         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:,IBL), 2, YDCST, YDTHF) 

         CALL CLOUDSC2 ( &
              &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
              & PTSPHY, LCETA, &
              & PAPH(:,:,IBL),  PAP(:,:,IBL), &
              & PQ(:,:,IBL), ZQSAT(:,:,IBL), PT(:,:,IBL), &
              & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
              & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
              & BUFFER_LOC(:,:,1,IBL), BUFFER_CML(:,:,1,IBL), &
              & BUFFER_LOC(:,:,3,IBL), BUFFER_CML(:,:,3,IBL), &
              & BUFFER_LOC(:,:,3+NCLDQL,IBL), BUFFER_CML(:,:,3+NCLDQL,IBL), &
              & BUFFER_LOC(:,:,3+NCLDQI,IBL), BUFFER_CML(:,:,3+NCLDQI,IBL), &
              &  PSUPSAT(:,:,IBL), &
              &  PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
              &  PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL), & 
              &  YDCST, YDTHF, YHNC, YPHLI, YCLD, YCLDP)

         ! Preparation for TL

         ! Increments (IN)
         ZAPH(:,:,IBL)    =PAPH(:,:,IBL)*0.01_JPRB
         ZAP(:,:,IBL)     =PAP(:,:,IBL)*0.01_JPRB
         ZQ(:,:,IBL)      =PQ(:,:,IBL)*0.01_JPRB
         ZZQSAT(:,:,IBL)  =ZQSAT(:,:,IBL)     *0.01_JPRB
         ZT(:,:,IBL)      = PT(:,:,IBL)*0.01_JPRB
         ZDRVL(:,:,IBL)   = PCLV(:,:,NCLDQL,IBL)*0.01_JPRB
         ZDRVI(:,:,IBL)   = PCLV(:,:,NCLDQI,IBL)*0.01_JPRB
         ZDRVLUDE(:,:,IBL)= PLUDE(:,:,IBL)*0.01_JPRB
         ZLU(:,:,IBL)     = PLU(:,:,IBL)*0.01_JPRB
         ZMFU(:,:,IBL)    = PMFU(:,:,IBL)*0.01_JPRB
         ZMFD(:,:,IBL)    = PMFD(:,:,IBL)*0.01_JPRB
         ZTENI_T(:,:,IBL) = BUFFER_CML(:,:,1,IBL)*0.01_JPRB
         ZTENI_Q(:,:,IBL) = BUFFER_CML(:,:,3,IBL)*0.01_JPRB
         ZTENI_L(:,:,IBL) = BUFFER_CML(:,:,3+NCLDQL,IBL)*0.01_JPRB
         ZTENI_I(:,:,IBL) = BUFFER_CML(:,:,3+NCLDQI,IBL)*0.01_JPRB
         ZSUPSAT(:,:,IBL) = PSUPSAT(:,:,IBL)*0.01_JPRB
         ! Call TL
         CALL  CLOUDSC2TL ( &
            &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
            & PTSPHY,LCETA,&
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
            & ZFHPSL(:,:,IBL) , ZFHPSN(:,:,IBL)   , ZCOVPTOT(:,:,IBL),&
            & YDCST, YDTHF, YHNC, YPHLI, YCLD, YCLDP, YNCL )       ! o

         ! Loop over incrementing states
         DO ILAM=1,10
           ZLAMBDA=10._JPRB**(-REAL(ILAM,JPRB))
           ! Define perturbed NL state
           PAPH5(:,:,IBL) = PAPH(:,:,IBL) + ZLAMBDA*ZAPH(:,:,IBL)
           PAP5(:,:,IBL)  = PAP(:,:,IBL)  + ZLAMBDA*ZAP(:,:,IBL)
           PQ5(:,:,IBL)   = PQ(:,:,IBL)   + ZLAMBDA*ZQ(:,:,IBL)
           ZQSAT5(:,:,IBL)= ZQSAT(:,:,IBL)    + ZLAMBDA*ZZQSAT(:,:,IBL)
           PT5(:,:,IBL)   = PT(:,:,IBL)   + ZLAMBDA*ZT(:,:,IBL)
           PCLVL5(:,:,IBL)= PCLV(:,:,NCLDQL,IBL) + ZLAMBDA*ZDRVL(:,:,IBL)
           PCLVI5(:,:,IBL)= PCLV(:,:,NCLDQI,IBL) + ZLAMBDA*ZDRVI(:,:,IBL)
           PLUDE5(:,:,IBL)= PLUDE(:,:,IBL)+ ZLAMBDA*ZDRVLUDE(:,:,IBL)
           PLU5(:,:,IBL)  = PLU(:,:,IBL)  + ZLAMBDA*ZLU(:,:,IBL)
           PMFU5(:,:,IBL) = PMFU(:,:,IBL) + ZLAMBDA*ZMFU(:,:,IBL)
           PMFD5(:,:,IBL) = PMFD(:,:,IBL) + ZLAMBDA*ZMFD(:,:,IBL)
           ZTENI_T5(:,:,IBL)=BUFFER_CML(:,:,1,IBL) + ZLAMBDA*ZTENI_T(:,:,IBL)
           ZTENI_Q5(:,:,IBL)=BUFFER_CML(:,:,3,IBL) + ZLAMBDA*ZTENI_Q(:,:,IBL)
           ZTENI_L5(:,:,IBL)=BUFFER_CML(:,:,3+NCLDQL,IBL) + ZLAMBDA*ZTENI_L(:,:,IBL)
           ZTENI_I5(:,:,IBL)=BUFFER_CML(:,:,3+NCLDQI,IBL) + ZLAMBDA*ZTENI_I(:,:,IBL)
           PSUPSAT5(:,:,IBL)=PSUPSAT(:,:,IBL) + ZLAMBDA*ZSUPSAT(:,:,IBL)
           ! Call the NL code
           CALL CLOUDSC2 ( &
              &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
              & PTSPHY, LCETA,&
              & PAPH5(:,:,IBL),  PAP5(:,:,IBL), &
              & PQ5(:,:,IBL), ZQSAT5(:,:,IBL), PT5(:,:,IBL), &
              & PCLVL5(:,:,IBL), PCLVI5(:,:,IBL), &
              & PLUDE5(:,:,IBL), PLU5(:,:,IBL), PMFU5(:,:,IBL), PMFD5(:,:,IBL),&
              & ZTENO_T5(:,:,ILAM,IBL), ZTENI_T5(:,:,IBL), &
              & ZTENO_Q5(:,:,ILAM,IBL), ZTENI_Q5(:,:,IBL), &
              & ZTENO_L5(:,:,ILAM,IBL), ZTENI_L5(:,:,IBL), &
              & ZTENO_I5(:,:,ILAM,IBL), ZTENI_I5(:,:,IBL), &
              &                         PSUPSAT5(:,:,IBL), &
              &      PA5(:,:,ILAM,IBL), &
              &  PFPLSL5(:,:,ILAM,IBL), &
              &  PFPLSN5(:,:,ILAM,IBL), &
              &  PFHPSL5(:,:,ILAM,IBL), &
              &  PFHPSN5(:,:,ILAM,IBL), &
              &PCOVPTOT5(:,:,ILAM,IBL), &
              & YDCST, YDTHF, YHNC, YPHLI, YCLD, YCLDP)


         ENDDO  ! end of lambda loops

#ifndef CLOUDSC_GPU_TIMING
         ! Log number of columns processed by this thread
         CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
#endif
      ENDDO

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

      CALL VALIDATE_TAYLOR_TEST(NPROMA, NLEV, NLAM, NGPTOT, &
       & BUFFER_LOC(:,:,1,:)       ,  ZTENO_T5(:,:,:,:),  ZTENO_T(:,:,:), &
       & BUFFER_LOC(:,:,3,:)       ,  ZTENO_Q5(:,:,:,:),  ZTENO_Q(:,:,:), &
       & BUFFER_LOC(:,:,3+NCLDQL,:),  ZTENO_L5(:,:,:,:),  ZTENO_L(:,:,:), &
       & BUFFER_LOC(:,:,3+NCLDQI,:),  ZTENO_I5(:,:,:,:),  ZTENO_I(:,:,:), &
       & PA(:,:,:)                 ,       PA5(:,:,:,:),     ZCLC(:,:,:), &
       & PFPLSL(:,:,:)             ,   PFPLSL5(:,:,:,:),   ZFPLSL(:,:,:), &
       & PFPLSN(:,:,:)             ,   PFPLSN5(:,:,:,:),   ZFPLSN(:,:,:), &
       & PFHPSL(:,:,:)             ,   PFHPSL5(:,:,:,:),   ZFHPSL(:,:,:), &
       & PFHPSN(:,:,:)             ,   PFHPSN5(:,:,:,:),   ZFHPSN(:,:,:), &
       & PCOVPTOT(:,:,:)           , PCOVPTOT5(:,:,:,:), ZCOVPTOT(:,:,:)  &
       & )
    
  END SUBROUTINE CLOUDSC_DRIVER_TL

END MODULE CLOUDSC_DRIVER_TL_MOD
