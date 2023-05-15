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
    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND,ISTART,ITEST,INEGAT,ITEMPNEGAT

    TYPE(PERFORMANCE_TIMER) :: TIMER
    REAL(KIND=JPRD), PARAMETER :: ZHPM = 3996006.0_JPRD  ! The nominal number of flops per 100 columns

    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    LOGICAL            :: LDRAIN1D = .FALSE.
    REAL(KIND=JPRB)    :: ZQSAT(NPROMA,NLEV) ! local array

    INTEGER(KIND=JPIB) :: ILAM
    REAL(KIND=JPRB)    :: ZLAMBDA, ZCOUNT, ZNORM, ZNORMG(10)
    REAL(KIND=JPRB)    :: ZAPH(NPROMA,NLEV+1), ZAP(NPROMA,NLEV), ZQ(NPROMA,NLEV), &
     & ZZQSAT(NPROMA,NLEV), ZT(NPROMA,NLEV), ZL(NPROMA,NLEV), ZI(NPROMA,NLEV), &
     & ZLUDE(NPROMA,NLEV), ZLU(NPROMA,NLEV), ZMFU(NPROMA,NLEV), ZMFD(NPROMA,NLEV), &
     & ZTENI_T(NPROMA,NLEV), ZTENI_Q(NPROMA,NLEV), ZTENI_L(NPROMA,NLEV), &
     & ZTENI_I(NPROMA,NLEV), ZSUPSAT(NPROMA,NLEV)
    REAL(KIND=JPRB)    :: ZTENO_T(NPROMA,NLEV), ZTENO_Q(NPROMA,NLEV), &
     & ZTENO_L(NPROMA,NLEV), ZTENO_I(NPROMA,NLEV), ZCLC(NPROMA,NLEV), &
     & ZFPLSL(NPROMA,NLEV+1), ZFPLSN(NPROMA,NLEV+1), ZFHPSL(NPROMA,NLEV+1), &
     & ZFHPSN(NPROMA,NLEV+1), ZCOVPTOT(NPROMA,NLEV)
    REAL(KIND=JPRB)    :: PAPH5(NPROMA,NLEV+1), PAP5(NPROMA,NLEV), &
     & PQ5(NPROMA,NLEV), ZQSAT5(NPROMA,NLEV), PT5(NPROMA,NLEV), &
     & PCLVL5(NPROMA,NLEV), PCLVI5(NPROMA,NLEV), PLUDE5(NPROMA,NLEV), &
     & PLU5(NPROMA,NLEV), PMFU5(NPROMA,NLEV), PMFD5(NPROMA,NLEV), &
     & ZTENI_T5(NPROMA,NLEV), ZTENI_Q5(NPROMA,NLEV), ZTENI_L5(NPROMA,NLEV), &
     & ZTENI_I5(NPROMA,NLEV), PSUPSAT5(NPROMA,NLEV)
    REAL(KIND=JPRB)    :: ZTENO_T5(NPROMA,NLEV), ZTENO_Q5(NPROMA,NLEV), &
     & ZTENO_L5(NPROMA,NLEV), ZTENO_I5(NPROMA,NLEV), PA5(NPROMA,NLEV), &
     & PFPLSL5(NPROMA,NLEV+1), PFPLSN5(NPROMA,NLEV+1), PFHPSL5(NPROMA,NLEV+1), &
     & PFHPSN5(NPROMA,NLEV+1), PCOVPTOT5(NPROMA,NLEV)

    TYPE(TOMCST)    :: YDCST, LOCAL_YDCST
    TYPE(TOETHF)    :: YDTHF, LOCAL_YDTHF
    TYPE(TPHNC)     :: YHNC, LOCAL_YHNC
    TYPE(TEPHLI)    :: YPHLI, LOCAL_YPHLI
    TYPE(TECLD)     :: YCLD, LOCAL_YCLD
    TYPE(TECLDP)    :: YCLDP, LOCAL_YCLDP
    TYPE(TNCL)      :: YNCL, LOCAL_YNCL

    LOCAL_YDCST=YDCST
    LOCAL_YDTHF=YDTHF
    LOCAL_YHNC=YHNC
    LOCAL_YPHLI=YPHLI
    LOCAL_YCLD=YCLD
    LOCAL_YCLDP=YCLDP
    LOCAL_YNCL=YNCL

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)
    
    ZNORMG(:)=0.

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         PCOVPTOT(:,:,IBL) = 0.0_JPRB
!        TENDENCY_LOC(IBL)%cld(:,:,NCLV) = 0.0_JPRB

         ! Fill in ZQSAT
         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:), 2, LOCAL_YDCST, LOCAL_YDTHF) 

         CALL CLOUDSC2 ( &
              &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
              & PTSPHY, LCETA, &
              & PAPH(:,:,IBL),  PAP(:,:,IBL), &
              & PQ(:,:,IBL), ZQSAT(:,:), PT(:,:,IBL), &
              & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
              & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
              & BUFFER_LOC(:,:,1,IBL), BUFFER_CML(:,:,1,IBL), &
              & BUFFER_LOC(:,:,3,IBL), BUFFER_CML(:,:,3,IBL), &
              & BUFFER_LOC(:,:,3+NCLDQL,IBL), BUFFER_CML(:,:,3+NCLDQL,IBL), &
              & BUFFER_LOC(:,:,3+NCLDQI,IBL), BUFFER_CML(:,:,3+NCLDQI,IBL), &
              &  PSUPSAT(:,:,IBL), &
              &  PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
              &  PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL), & 
              &  LOCAL_YDCST, LOCAL_YDTHF, LOCAL_YHNC, LOCAL_YPHLI, LOCAL_YCLD, LOCAL_YCLDP)

         ! Preparation for TL

         ! Increments (IN)
         ZAPH    =PAPH(:,:,IBL)*0.01_JPRB
         ZAP     =PAP(:,:,IBL)*0.01_JPRB
         ZQ      =PQ(:,:,IBL)*0.01_JPRB
         ZZQSAT  =ZQSAT     *0.01_JPRB
         ZT      = PT(:,:,IBL)*0.01_JPRB
         ZL      = PCLV(:,:,NCLDQL,IBL)*0.01_JPRB
         ZI      = PCLV(:,:,NCLDQI,IBL)*0.01_JPRB
         ZLUDE   = PLUDE(:,:,IBL)*0.01_JPRB
         ZLU     = PLU(:,:,IBL)*0.01_JPRB
         ZMFU    = PMFU(:,:,IBL)*0.01_JPRB
         ZMFD    = PMFD(:,:,IBL)*0.01_JPRB
         ZTENI_T = BUFFER_CML(:,:,1,IBL)*0.01_JPRB
         ZTENI_Q = BUFFER_CML(:,:,3,IBL)*0.01_JPRB
         ZTENI_L = BUFFER_CML(:,:,3+NCLDQL,IBL)*0.01_JPRB
         ZTENI_I = BUFFER_CML(:,:,3+NCLDQI,IBL)*0.01_JPRB
         ZSUPSAT = PSUPSAT(:,:,IBL)*0.01_JPRB
         ! Call TL
         CALL  CLOUDSC2TL ( &
            &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
            & PTSPHY,LCETA,&
            ! trajectory
            & PAPH(:,:,IBL),  PAP(:,:,IBL), &
            & PQ(:,:,IBL), ZQSAT(:,:), PT(:,:,IBL), &
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
            & ZAPH, ZAP, ZQ, ZZQSAT, ZT, ZL, ZI, &
            & ZLUDE, ZLU, ZMFU, ZMFD, &
            & ZTENO_T, ZTENI_T, ZTENO_Q, ZTENI_Q, &   ! o,i,o,i
            & ZTENO_L, ZTENI_L, ZTENO_I, ZTENI_I, ZSUPSAT, &  ! o,i,o,i
            & ZCLC   , ZFPLSL   , ZFPLSN ,&        ! o
            & ZFHPSL , ZFHPSN   , ZCOVPTOT,&
            & LOCAL_YDCST, LOCAL_YDTHF, LOCAL_YHNC, LOCAL_YPHLI, LOCAL_YCLD, LOCAL_YCLDP, LOCAL_YNCL )       ! o

         ! Loop over incrementing states
         DO ILAM=1,10
           ZLAMBDA=10._JPRB**(-REAL(ILAM,JPRB))
           ! Define perturbed NL state
           PAPH5(:,:) = PAPH(:,:,IBL) + ZLAMBDA*ZAPH(:,:)
           PAP5(:,:)  = PAP(:,:,IBL)  + ZLAMBDA*ZAP(:,:)
           PQ5(:,:)   = PQ(:,:,IBL)   + ZLAMBDA*ZQ(:,:)
           ZQSAT5(:,:)= ZQSAT(:,:)    + ZLAMBDA*ZZQSAT(:,:)
           PT5(:,:)   = PT(:,:,IBL)   + ZLAMBDA*ZT(:,:)
           PCLVL5(:,:)= PCLV(:,:,NCLDQL,IBL) + ZLAMBDA*ZL(:,:)
           PCLVI5(:,:)= PCLV(:,:,NCLDQI,IBL) + ZLAMBDA*ZI(:,:)
           PLUDE5(:,:)= PLUDE(:,:,IBL)+ ZLAMBDA*ZLUDE(:,:)
           PLU5(:,:)  = PLU(:,:,IBL)  + ZLAMBDA*ZLU(:,:)
           PMFU5(:,:) = PMFU(:,:,IBL) + ZLAMBDA*ZMFU(:,:)
           PMFD5(:,:) = PMFD(:,:,IBL) + ZLAMBDA*ZMFD(:,:)
           ZTENI_T5(:,:)=BUFFER_CML(:,:,1,IBL) + ZLAMBDA*ZTENI_T(:,:)
           ZTENI_Q5(:,:)=BUFFER_CML(:,:,3,IBL) + ZLAMBDA*ZTENI_Q(:,:)
           ZTENI_L5(:,:)=BUFFER_CML(:,:,3+NCLDQL,IBL) + ZLAMBDA*ZTENI_L(:,:)
           ZTENI_I5(:,:)=BUFFER_CML(:,:,3+NCLDQI,IBL) + ZLAMBDA*ZTENI_I(:,:)
           PSUPSAT5(:,:)=PSUPSAT(:,:,IBL) + ZLAMBDA*ZSUPSAT(:,:)
           ! Call the NL code
           CALL CLOUDSC2 ( &
              &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
              & PTSPHY, LCETA,&
              & PAPH5,  PAP5, &
              & PQ5, ZQSAT5, PT5, &
              & PCLVL5, PCLVI5, &
              & PLUDE5, PLU5, PMFU5, PMFD5,&
              & ZTENO_T5, ZTENI_T5, &
              & ZTENO_Q5, ZTENI_Q5, &
              & ZTENO_L5, ZTENI_L5, &
              & ZTENO_I5, ZTENI_I5, &
              & PSUPSAT5, &
              & PA5(:,:), PFPLSL5(:,:),   PFPLSN5(:,:), &
              & PFHPSL5(:,:),   PFHPSN5(:,:), PCOVPTOT5(:,:), &
              & LOCAL_YDCST, LOCAL_YDTHF, LOCAL_YHNC, LOCAL_YPHLI, LOCAL_YCLD, LOCAL_YCLDP)

           ! Compute final test norm
           ZCOUNT=0._JPRB
           ZNORM= 0._JPRB 

           ! Global norm (normalize by number of active statistics)
           IF (ZNORM == 0._JPRB .OR. ZCOUNT == 0._JPRB) THEN
             print *, ' TL is totally wrong !!! ',ZNORM,ZCOUNT
             stop
           ELSE
             ZNORMG(ILAM)=MAX(ZNORMG(ILAM),ZNORM/ZCOUNT)
           ENDIF

         ENDDO  ! end of lambda loops

         ! Log number of columns processed by this thread
         CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
      ENDDO

      CALL TIMER%THREAD_END(TID)

      CALL TIMER%END()

      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, ZHPM, NGPTOT)

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
