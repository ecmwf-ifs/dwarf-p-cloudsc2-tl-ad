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
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, PTSPHY, &
     & PT, PQ, TENDENCY_CML, TENDENCY_LOC, &
     & PAP,      PAPH, &
     & PLU,      PLUDE,    PMFU,     PMFD, &
     & PA,       PCLV,     PSUPSAT,&
     & PCOVPTOT, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN &
     & )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC2 kernel

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: PQ(:,:,:)    ! Q at start of callpar
    TYPE(STATE_TYPE),   INTENT(IN)    :: TENDENCY_CML(:) ! cumulative tendency used for final output
    TYPE(STATE_TYPE),   INTENT(OUT)   :: TENDENCY_LOC(:) ! local tendency from cloud scheme
    REAL(KIND=JPRB),    INTENT(IN)    :: PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(INOUT) :: PLUDE(:,:,:) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFD(:,:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(IN)    :: PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(IN)    :: PCLV(:,:,:,:) 
    REAL(KIND=JPRB),    INTENT(IN)    :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB),    INTENT(INOUT) :: PCOVPTOT(:,:,:) ! Precip fraction
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSL(:,:,:) ! liq+rain sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSN(:,:,:) ! ice+snow sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSL(:,:,:) ! Enthalpy flux for liq
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSN(:,:,:) ! Enthalp flux for ice

    INTEGER(KIND=JPIM) :: JKGLO,IBL,JROF,ICEND,NGPBLKS

    TYPE(PERFORMANCE_TIMER) :: TIMER
    INTEGER(KIND=JPIM) :: TID ! thread id from 0 .. NUMOMP - 1
    LOGICAL            :: LDRAIN1D = .FALSE.
    REAL(KIND=JPRB)    :: ZQSAT(NPROMA,NLEV) ! local array

    REAL(KIND=JPRB)    :: ZNORMG
    REAL(KIND=JPRB)    :: ZNORM1(NPROMA), ZNORM2(NPROMA), ZNORM3(NPROMA)
    REAL(KIND=JPRB)    :: ZAPH(NPROMA,NLEV+1), ZAP(NPROMA,NLEV), ZQ(NPROMA,NLEV), &
     & ZZQSAT(NPROMA,NLEV), ZT(NPROMA,NLEV), ZL(NPROMA,NLEV), ZI(NPROMA,NLEV), &
     & ZLUDE(NPROMA,NLEV), ZLU(NPROMA,NLEV), ZMFU(NPROMA,NLEV), ZMFD(NPROMA,NLEV), &
     & ZTENI_T(NPROMA,NLEV), ZTENI_Q(NPROMA,NLEV), ZTENI_L(NPROMA,NLEV), &
     & ZTENI_I(NPROMA,NLEV), ZSUPSAT(NPROMA,NLEV)
    REAL(KIND=JPRB)    :: ZTENO_T(NPROMA,NLEV), ZTENO_Q(NPROMA,NLEV), &
     & ZTENO_L(NPROMA,NLEV), ZTENO_I(NPROMA,NLEV), ZCLC(NPROMA,NLEV), &
     & ZFPLSL(NPROMA,NLEV+1), ZFPLSN(NPROMA,NLEV+1), ZFHPSL(NPROMA,NLEV+1), &
     & ZFHPSN(NPROMA,NLEV+1), ZCOVPTOT(NPROMA,NLEV)
    REAL(KIND=JPRB)    :: ZAPH0(NPROMA,NLEV+1), ZAP0(NPROMA,NLEV), ZQ0(NPROMA,NLEV), &
     & ZZQSAT0(NPROMA,NLEV), ZT0(NPROMA,NLEV), ZL0(NPROMA,NLEV), ZI0(NPROMA,NLEV), &
     & ZLUDE0(NPROMA,NLEV), ZLU0(NPROMA,NLEV), ZMFU0(NPROMA,NLEV), ZMFD0(NPROMA,NLEV), &
     & ZTENI_T0(NPROMA,NLEV), ZTENI_Q0(NPROMA,NLEV), ZTENI_L0(NPROMA,NLEV), &
     & ZTENI_I0(NPROMA,NLEV), ZSUPSAT0(NPROMA,NLEV)

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMPROC=',i0', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)
    
    ZNORMG=0._JPRB

    !$omp parallel default(shared) private(JKGLO,IBL,JROF,ICEND,TID) &
    !$omp& private(ZQSAT) &
    !$omp& private(ZNORM1,ZNORM2,ZNORM3) &
    !$omp& private(ZAPH,ZAP,ZQ,ZZQSAT,ZT,ZL,ZI,ZLUDE,ZLU,ZMFU,ZMFD) &
    !$omp& private(ZTENI_T,ZTENI_Q,ZTENI_L,ZTENI_I,ZSUPSAT) &
    !$omp& private(ZTENO_T,ZTENO_Q,ZTENO_L,ZTENO_I) &
    !$omp& private(ZCLC,ZFPLSL,ZFPLSN,ZFHPSL,ZFHPSN,ZCOVPTOT) &
    !$omp& private(ZAPH0,ZAP0,ZQ0,ZZQSAT0,ZT0,ZL0,ZI0,ZLUDE0,ZLU0,ZMFU0,ZMFD0) &
    !$omp& private(ZTENI_T0,ZTENI_Q0,ZTENI_L0,ZTENI_I0,ZSUPSAT0) &
    !$omp& num_threads(NUMOMP)

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    !$omp do schedule(runtime) reduction(max:znormg)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         PCOVPTOT(:,:,IBL) = 0.0_JPRB
         TENDENCY_LOC(IBL)%cld(:,:,NCLV) = 0.0_JPRB

         ! Fill in ZQSAT
         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:), 2) 


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
         ZTENI_T = TENDENCY_CML(IBL)%T*0.01_JPRB
         ZTENI_Q = TENDENCY_CML(IBL)%Q*0.01_JPRB
         ZTENI_L = TENDENCY_CML(IBL)%CLD(:,:,NCLDQL)*0.01_JPRB
         ZTENI_I = TENDENCY_CML(IBL)%CLD(:,:,NCLDQI)*0.01_JPRB
         ZSUPSAT = 0.00_JPRB  ! obsolette, beter not use

         ! Storage of initial increments
         ZAPH0   = ZAPH
         ZAP0    = ZAP
         ZQ0     = ZQ
         ZZQSAT0 = ZZQSAT
         ZT0     = ZT
         ZL0     = ZL
         ZI0     = ZI
         ZLUDE0  = ZLUDE
         ZLU0    = ZLU
         ZMFU0   = ZMFU
         ZMFD0   = ZMFD
         ZTENI_T0= ZTENI_T
         ZTENI_Q0= ZTENI_Q
         ZTENI_L0= ZTENI_L
         ZTENI_I0= ZTENI_I
         ZSUPSAT0= ZSUPSAT

         ! Tangent linear integration
         CALL  CLOUDSC2TL ( &
            &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
            & PTSPHY,&
            ! trajectory
            & PAPH(:,:,IBL),  PAP(:,:,IBL), &
            & PQ(:,:,IBL), ZQSAT(:,:), PT(:,:,IBL), &
            & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
            & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
            & TENDENCY_LOC(IBL)%T, TENDENCY_CML(IBL)%T, &
            & TENDENCY_LOC(IBL)%Q, TENDENCY_CML(IBL)%Q, &
            & TENDENCY_LOC(IBL)%CLD(:,:,NCLDQL), TENDENCY_CML(IBL)%CLD(:,:,NCLDQL), &
            & TENDENCY_LOC(IBL)%CLD(:,:,NCLDQI), TENDENCY_CML(IBL)%CLD(:,:,NCLDQI), &
            & PSUPSAT(:,:,IBL), &
            & PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
            & PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL), &
            ! increments
            & ZAPH, ZAP, ZQ, ZZQSAT, ZT, ZL, ZI, &
            & ZLUDE, ZLU, ZMFU, ZMFD, &
            & ZTENO_T, ZTENI_T, ZTENO_Q, ZTENI_Q, &   ! o,i,o,i
            & ZTENO_L, ZTENI_L, ZTENO_I, ZTENI_I, ZSUPSAT, &  ! o,i,o,i
            & ZCLC   , ZFPLSL   , ZFPLSN ,&        ! o
            & ZFHPSL , ZFHPSN   , ZCOVPTOT )       ! o

         ! First norm
         DO JROF=1,ICEND
           ZNORM1(JROF)=SUM(ZTENO_T(JROF,1:NLEV)*ZTENO_T(JROF,1:NLEV)) &
           &  + SUM(ZTENO_Q(JROF,1:NLEV)*ZTENO_Q(JROF,1:NLEV)) &
           &  + SUM(ZTENO_L(JROF,1:NLEV)*ZTENO_L(JROF,1:NLEV)) &
           &  + SUM(ZTENO_I(JROF,1:NLEV)*ZTENO_I(JROF,1:NLEV)) &
           &  + SUM(ZCLC(JROF,1:NLEV)*ZCLC(JROF,1:NLEV)) &
           &  + SUM(ZFPLSL(JROF,1:NLEV+1)*ZFPLSL(JROF,1:NLEV+1)) &
           &  + SUM(ZFPLSN(JROF,1:NLEV+1)*ZFPLSN(JROF,1:NLEV+1)) &
           &  + SUM(ZFHPSL(JROF,1:NLEV+1)*ZFHPSL(JROF,1:NLEV+1)) &
           &  + SUM(ZFHPSN(JROF,1:NLEV+1)*ZFHPSN(JROF,1:NLEV+1)) &
           &  + SUM(ZCOVPTOT(JROF,1:NLEV)*ZCOVPTOT(JROF,1:NLEV))
         ENDDO

         ! Initiaslization of output variables
         ZAPH    = 0.0_JPRB
         ZAP     = 0.0_JPRB
         ZQ      = 0.0_JPRB
         ZZQSAT  = 0.0_JPRB
         ZT      = 0.0_JPRB
         ZL      = 0.0_JPRB
         ZI      = 0.0_JPRB
         ZLUDE   = 0.0_JPRB
         ZLU     = 0.0_JPRB
         ZMFU    = 0.0_JPRB
         ZMFD    = 0.0_JPRB
         ZTENI_T = 0.0_JPRB
         ZTENI_Q = 0.0_JPRB
         ZTENI_L = 0.0_JPRB
         ZTENI_I = 0.0_JPRB
         ZSUPSAT = 0.0_JPRB

         ! Adjoint integration
         CALL  CLOUDSC2AD ( &
            &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
            & PTSPHY,&
            ! trajectory
            & PAPH(:,:,IBL),  PAP(:,:,IBL), &
            & PQ(:,:,IBL), ZQSAT(:,:), PT(:,:,IBL), &
            & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
            & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
            & TENDENCY_LOC(IBL)%T, TENDENCY_CML(IBL)%T, &
            & TENDENCY_LOC(IBL)%Q, TENDENCY_CML(IBL)%Q, &
            & TENDENCY_LOC(IBL)%CLD(:,:,NCLDQL), TENDENCY_CML(IBL)%CLD(:,:,NCLDQL), &
            & TENDENCY_LOC(IBL)%CLD(:,:,NCLDQI), TENDENCY_CML(IBL)%CLD(:,:,NCLDQI), &
            & PSUPSAT(:,:,IBL), &
            & PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
            & PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL), &
            ! increments
            & ZAPH, ZAP, ZQ, ZZQSAT, ZT, ZL, ZI, &
            & ZLUDE, ZLU, ZMFU, ZMFD, &
            & ZTENO_T, ZTENI_T, ZTENO_Q, ZTENI_Q, &   ! o,i,o,i
            & ZTENO_L, ZTENI_L, ZTENO_I, ZTENI_I, ZSUPSAT, &  ! o,i,o,i
            & ZCLC   , ZFPLSL   , ZFPLSN ,&        ! o
            & ZFHPSL , ZFHPSN   , ZCOVPTOT )       ! o

         ! Second norm
         DO JROF=1,ICEND
           ZNORM2(JROF)=SUM(ZAPH0(JROF,1:NLEV+1)*ZAPH(JROF,1:NLEV+1)) &
           &  + SUM(ZAP0(JROF,1:NLEV)*ZAP(JROF,1:NLEV)) &
           &  + SUM(ZQ0(JROF,1:NLEV)*ZQ(JROF,1:NLEV)) &
           &  + SUM(ZZQSAT0(JROF,1:NLEV)*ZZQSAT(JROF,1:NLEV)) &
           &  + SUM(ZT0(JROF,1:NLEV)*ZT(JROF,1:NLEV)) &
           &  + SUM(ZL0(JROF,1:NLEV)*ZL(JROF,1:NLEV)) &
           &  + SUM(ZI0(JROF,1:NLEV)*ZI(JROF,1:NLEV)) &
           &  + SUM(ZLUDE0(JROF,1:NLEV)*ZLUDE(JROF,1:NLEV)) &
           &  + SUM(ZLU0(JROF,1:NLEV)*ZLU(JROF,1:NLEV)) &
           &  + SUM(ZMFU0(JROF,1:NLEV)*ZMFU(JROF,1:NLEV)) &
           &  + SUM(ZMFD0(JROF,1:NLEV)*ZMFD(JROF,1:NLEV)) &
           &  + SUM(ZTENI_T0(JROF,1:NLEV)*ZTENI_T(JROF,1:NLEV)) &
           &  + SUM(ZTENI_Q0(JROF,1:NLEV)*ZTENI_Q(JROF,1:NLEV)) &
           &  + SUM(ZTENI_L0(JROF,1:NLEV)*ZTENI_L(JROF,1:NLEV)) &
           &  + SUM(ZTENI_I0(JROF,1:NLEV)*ZTENI_I(JROF,1:NLEV)) &
           &  + SUM(ZSUPSAT0(JROF,1:NLEV)*ZSUPSAT(JROF,1:NLEV))
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

         ! Log number of columns processed by this thread
         CALL TIMER%THREAD_LOG(TID, IGPC=ICEND)
      ENDDO

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc
      !$omp end do nowait

      CALL TIMER%THREAD_END(TID)

      !$omp end parallel

      CALL TIMER%END()

      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

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
