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
    USE ERROR_MOD, ONLY: ERROR_NORM 
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
    INTEGER(KIND=JPIB), PARAMETER :: NLAM=10
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
    REAL(KIND=JPRB)    :: ZL(NPROMA,NLEV,NGPBLKS)
    REAL(KIND=JPRB)    :: ZI(NPROMA,NLEV,NGPBLKS)
    REAL(KIND=JPRB)    :: ZLUDE(NPROMA,NLEV,NGPBLKS)
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
         ZL(:,:,IBL)      = PCLV(:,:,NCLDQL,IBL)*0.01_JPRB
         ZI(:,:,IBL)      = PCLV(:,:,NCLDQI,IBL)*0.01_JPRB
         ZLUDE(:,:,IBL)   = PLUDE(:,:,IBL)*0.01_JPRB
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
            & ZAPH(:,:,IBL), ZAP(:,:,IBL), ZQ(:,:,IBL), ZZQSAT(:,:,IBL), ZT(:,:,IBL), ZL(:,:,IBL), ZI(:,:,IBL), &
            & ZLUDE(:,:,IBL), ZLU(:,:,IBL), ZMFU(:,:,IBL), ZMFD(:,:,IBL), &
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
           PCLVL5(:,:,IBL)= PCLV(:,:,NCLDQL,IBL) + ZLAMBDA*ZL(:,:,IBL)
           PCLVI5(:,:,IBL)= PCLV(:,:,NCLDQI,IBL) + ZLAMBDA*ZI(:,:,IBL)
           PLUDE5(:,:,IBL)= PLUDE(:,:,IBL)+ ZLAMBDA*ZLUDE(:,:,IBL)
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

         ! Log number of columns processed by this thread
 !       CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
      ENDDO

      CALL TIMER%THREAD_END(TID)
      !$loki end data

      CALL TIMER%END()

      CALL TIMER%THREAD_LOG(TID, IGPC=NGPTOT)
      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, ZHPM, NGPTOT)

    DO JKGLO=1,NGPTOT,NPROMA
       IBLT=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
         ! Loop over incrementing states
         DO ILAM=1,10
           ZLAMBDA=10._JPRB**(-REAL(ILAM,JPRB))
           ! Compute final test norm
           ZCOUNT=0._JPRB
           ZNORM= 0._JPRB 
         
           CALL ERROR_NORM(1,ICEND,ICEND, BUFFER_LOC(:,:,1,IBLT)       ,  ZTENO_T5(:,:,ILAM,IBLT),  ZTENO_T(:,:,IBLT), ZNORM, ZCOUNT, ZLAMBDA)
           CALL ERROR_NORM(1,ICEND,ICEND, BUFFER_LOC(:,:,3,IBLT)       ,  ZTENO_Q5(:,:,ILAM,IBLT),  ZTENO_Q(:,:,IBLT), ZNORM, ZCOUNT, ZLAMBDA)
           CALL ERROR_NORM(1,ICEND,ICEND, BUFFER_LOC(:,:,3+NCLDQL,IBLT),  ZTENO_L5(:,:,ILAM,IBLT),  ZTENO_L(:,:,IBLT), ZNORM, ZCOUNT, ZLAMBDA)
           CALL ERROR_NORM(1,ICEND,ICEND, BUFFER_LOC(:,:,3+NCLDQI,IBLT),  ZTENO_I5(:,:,ILAM,IBLT),  ZTENO_I(:,:,IBLT), ZNORM, ZCOUNT, ZLAMBDA)
           CALL ERROR_NORM(1,ICEND,ICEND, PA(:,:,IBLT)                 ,       PA5(:,:,ILAM,IBLT),     ZCLC(:,:,IBLT), ZNORM, ZCOUNT, ZLAMBDA)
           CALL ERROR_NORM(1,ICEND,ICEND, PFPLSL(:,:,IBLT)             ,   PFPLSL5(:,:,ILAM,IBLT),   ZFPLSL(:,:,IBLT), ZNORM, ZCOUNT, ZLAMBDA)
           CALL ERROR_NORM(1,ICEND,ICEND, PFPLSN(:,:,IBLT)             ,   PFPLSN5(:,:,ILAM,IBLT),   ZFPLSN(:,:,IBLT), ZNORM, ZCOUNT, ZLAMBDA)
           CALL ERROR_NORM(1,ICEND,ICEND, PFHPSL(:,:,IBLT)             ,   PFHPSL5(:,:,ILAM,IBLT),   ZFHPSL(:,:,IBLT), ZNORM, ZCOUNT, ZLAMBDA)
           CALL ERROR_NORM(1,ICEND,ICEND, PFHPSN(:,:,IBLT)             ,   PFHPSN5(:,:,ILAM,IBLT),   ZFHPSN(:,:,IBLT), ZNORM, ZCOUNT, ZLAMBDA)
           CALL ERROR_NORM(1,ICEND,ICEND, PCOVPTOT(:,:,IBLT)           , PCOVPTOT5(:,:,ILAM,IBLT), ZCOVPTOT(:,:,IBLT), ZNORM, ZCOUNT, ZLAMBDA)

           ! Global norm (normalize by number of active statistics)
           IF (ZNORM == 0._JPRB .OR. ZCOUNT == 0._JPRB) THEN
             print *, ' TL is totally wrong !!! ',ZNORM,ZCOUNT
             stop
           ELSE
             ZNORMG(ILAM)=MAX(ZNORMG(ILAM),ZNORM/ZCOUNT)
           ENDIF
         ENDDO  ! end of lambda loops
      ENDDO
        
      ! Evaluate the test and print the otput
      print *, ' TL Taylor test '
      print *, '                Lambda   Result'
      istart=0
      DO ILAM=1,10
         print *, ILAM, ZNORMG(ILAM)
         ! Redefine ZNORMG
         ZNORMG(ILAM)=ABS(1._JPRB - ZNORMG(ILAM))
         ! filter out first members with strong NL departures
         if (istart == 0 .AND.  ZNORMG(ILAM) < 0.5_JPRB ) istart=ILAM 
      ENDDO

      print *, '   ==============================================   '
      IF (ISTART == 0 .OR. ISTART > 4 ) THEN
        print *, '       TEST FAILLED, err 13 '
      ELSE
        ! V-shape test
        ITEST=-10
        INEGAT=1
        DO ILAM=ISTART,10-1
          IF (ZNORMG(ILAM+1)/ZNORMG(ILAM) < 1._JPRB ) THEN
            ITEMPNEGAT = 1
          ELSE
            ITEMPNEGAT = 0
          ENDIF
          IF (INEGAT > ITEMPNEGAT) ITEST=ITEST+10
          INEGAT=ITEMPNEGAT
        ENDDO
        IF (ITEST == -10) ITEST = 11 ! no change of sign at all
        ! Accuracy test
        IF (MINVAL(ZNORMG(ISTART:10)) > 0.00001_JPRB) ITEST=ITEST+7  ! Hard limit
        IF (MINVAL(ZNORMG(ISTART:10)) > 0.000001_JPRB) ITEST=ITEST+5  ! Soft limit
        ! Final prints
        IF (ITEST > 5) THEN
          print *, '       TEST FAILLED, err ',ITEST
        ELSE
          print *, '       TEST PASSED, penalty ',ITEST
        ENDIF 
      ENDIF
      print *, '   ==============================================   '
    
  END SUBROUTINE CLOUDSC_DRIVER_TL

END MODULE CLOUDSC_DRIVER_TL_MOD
