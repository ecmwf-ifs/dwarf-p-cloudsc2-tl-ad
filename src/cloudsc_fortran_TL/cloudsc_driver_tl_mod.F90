MODULE CLOUDSC_DRIVER_TL_MOD
  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMPHYDER, ONLY: STATE_TYPE
  USE YOECLDP, ONLY : NCLV, NCLDQL, NCLDQI 
  USE CLOUDSC_MPI_MOD, ONLY: NUMPROC, IRANK
  USE TIMER_MOD, ONLY : PERFORMANCE_TIMER, GET_THREAD_NUM
  USE EC_PMON_MOD, ONLY: EC_PMON

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CLOUDSC_DRIVER_TL( &
     & NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, KFLDX, PTSPHY, &
     & PT, PQ, TENDENCY_CML, TENDENCY_TMP, TENDENCY_LOC, &
     & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
     & PHRSW,    PHRLW, &
     & PVERVEL,  PAP,      PAPH, &
     & PLSM,     LDCUM,    KTYPE, &
     & PLU,      PLUDE,    PSNDE,    PMFU,     PMFD, &
     & PA,       PCLV,     PSUPSAT,&
     & PLCRIT_AER,PICRIT_AER, PRE_ICE, &
     & PCCN,     PNICE,&
     & PCOVPTOT, PRAINFRAC_TOPRFZ, &
     & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG, &
     & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG, &
     & PFSQLTUR, PFSQITUR, &
     & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN &
     & )
    ! Driver routine that performans the parallel NPROMA-blocking and
    ! invokes the CLOUDSC2 kernel

    INTEGER(KIND=JPIM), INTENT(IN)    :: NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG
    INTEGER(KIND=JPIM), INTENT(IN)    :: KFLDX      ! NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY       ! Physics timestep
    REAL(KIND=JPRB),    INTENT(IN)    :: PT(:,:,:)    ! T at start of callpar
    REAL(KIND=JPRB),    INTENT(IN)    :: PQ(:,:,:)    ! Q at start of callpar
    TYPE(STATE_TYPE),   INTENT(IN)    :: TENDENCY_CML(:) ! cumulative tendency used for final output
    TYPE(STATE_TYPE),   INTENT(IN)    :: TENDENCY_TMP(:) !  NOT_USED
    TYPE(STATE_TYPE),   INTENT(OUT)   :: TENDENCY_LOC(:) ! local tendency from cloud scheme
    REAL(KIND=JPRB),    INTENT(IN)    :: PVFA(:,:,:)  ! CC from VDF scheme NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PVFL(:,:,:)  ! Liq from VDF scheme NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PVFI(:,:,:)  ! Ice from VDF scheme NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PDYNA(:,:,:) ! CC from Dynamics  NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PDYNL(:,:,:) ! Liq from Dynamics  NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PDYNI(:,:,:) ! Liq from Dynamics  NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PHRSW(:,:,:) ! Short-wave heating rate NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PHRLW(:,:,:) ! Long-wave heating rate NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PVERVEL(:,:,:) !Vertical velocity NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PAP(:,:,:)   ! Pressure on full levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PAPH(:,:,:)  ! Pressure on half levels
    REAL(KIND=JPRB),    INTENT(IN)    :: PLSM(:,:)    ! Land fraction (0-1) NOT_USED
    LOGICAL        ,    INTENT(IN)    :: LDCUM(:,:)   ! Convection active  NOT_USED
    INTEGER(KIND=JPIM), INTENT(IN)    :: KTYPE(:,:)   ! Convection type 0,1,2  NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PLU(:,:,:)   ! Conv. condensate
    REAL(KIND=JPRB),    INTENT(INOUT) :: PLUDE(:,:,:) ! Conv. detrained water
    REAL(KIND=JPRB),    INTENT(IN)    :: PSNDE(:,:,:) ! Conv. detrained snow NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFU(:,:,:)  ! Conv. mass flux up
    REAL(KIND=JPRB),    INTENT(IN)    :: PMFD(:,:,:)  ! Conv. mass flux down
    REAL(KIND=JPRB),    INTENT(IN)    :: PA(:,:,:)    ! Original Cloud fraction (t)
    REAL(KIND=JPRB),    INTENT(IN)    :: PCLV(:,:,:,:) 
    REAL(KIND=JPRB),    INTENT(IN)    :: PSUPSAT(:,:,:)
    REAL(KIND=JPRB),    INTENT(IN)    :: PLCRIT_AER(:,:,:) ! NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PICRIT_AER(:,:,:) ! NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PRE_ICE(:,:,:)    ! NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PCCN(:,:,:)     ! liquid cloud condensation nuclei NOT_USED
    REAL(KIND=JPRB),    INTENT(IN)    :: PNICE(:,:,:)    ! ice number concentration (cf. CCN) NOT_USED

    REAL(KIND=JPRB),    INTENT(INOUT) :: PCOVPTOT(:,:,:) ! Precip fraction
    REAL(KIND=JPRB),    INTENT(OUT)   :: PRAINFRAC_TOPRFZ(:,:) ! NOT_USED
    ! Flux diagnostics for DDH budget
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQLF(:,:,:)  ! Flux of liquid NOT_USED
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQIF(:,:,:)  ! Flux of ice  NOT_USED
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQLNG(:,:,:) ! -ve corr for liq  NOT_USED
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQNNG(:,:,:) ! -ve corr for ice  NOT_USED
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQRF(:,:,:)  ! Flux diagnostics  NOT_USED
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQSF(:,:,:)  !    for DDH, generic  NOT_USED
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQRNG(:,:,:) ! rain  NOT_USED
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFCQSNG(:,:,:) ! snow  NOT_USED
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQLTUR(:,:,:) ! liquid flux due to VDF NOT_USED
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFSQITUR(:,:,:) ! ice flux due to VDF NOT_USED
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSL(:,:,:) ! liq+rain sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFPLSN(:,:,:) ! ice+snow sedim flux
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSL(:,:,:) ! Enthalpy flux for liq
    REAL(KIND=JPRB),    INTENT(OUT)   :: PFHPSN(:,:,:) ! Enthalp flux for ice

    INTEGER(KIND=JPIM) :: JKGLO,IBL,ICEND,NGPBLKS,ISTART,ITEST,INEGAT,ITEMPNEGAT

    TYPE(PERFORMANCE_TIMER) :: TIMER
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

    NGPBLKS = (NGPTOT / NPROMA) + MIN(MOD(NGPTOT,NPROMA), 1)
1003 format(5x,'NUMPROC=',i0', NUMOMP=',i0,', NGPTOTG=',i0,', NPROMA=',i0,', NGPBLKS=',i0)
    if (irank == 0) then
      write(0,1003) NUMPROC,NUMOMP,NGPTOTG,NPROMA,NGPBLKS
    end if

    ! Global timer for the parallel region
    CALL TIMER%START(NUMOMP)
    
    ZNORMG(:)=0.

    !$omp parallel default(shared) private(JKGLO,IBL,ICEND,TID,energy,power) &
    !$omp& private(ZQSAT) &
    !$omp& private(ILAM,ZLAMBDA,ZCOUNT,ZNORM) &
    !$omp& private(ZAPH,ZAP,ZQ,ZZQSAT,ZT,ZL,ZI,ZLUDE,ZLU,ZMFU,ZMFD) &
    !$omp& private(ZTENI_T,ZTENI_Q,ZTENI_L,ZTENI_I, ZSUPSAT) &
    !$omp& private(ZTENO_T,ZTENO_Q,ZTENO_L,ZTENO_I) &
    !$omp& private(ZCLC,ZFPLSL,ZFPLSN,ZFHPSL,ZFHPSN,ZCOVPTOT) &
    !$omp& private(PAPH5,PAP5,PQ5,ZQSAT5,PT5,PCLVL5,PCLVI5,PLUDE5,PLU5) &
    !$omp& private(PMFU5,PMFD5,ZTENI_T5,ZTENI_Q5,ZTENI_L5,ZTENI_I5,PSUPSAT5) &
    !$omp& private(ZTENO_T5,ZTENO_Q5,ZTENO_L5,ZTENO_I5,PA5) &
    !$omp& private(PFPLSL5, PFPLSN5, PFHPSL5,PFHPSN5, PCOVPTOT5) &
    !$omp& num_threads(NUMOMP)

    ! Local timer for each thread
    TID = GET_THREAD_NUM()
    CALL TIMER%THREAD_START(TID)

    !$omp do schedule(runtime) reduction(+:power_total,power_count) reduction(max:power_max) &
    !$OMP& REDUCTION(MAX:ZNORMG)
    DO JKGLO=1,NGPTOT,NPROMA
       IBL=(JKGLO-1)/NPROMA+1
       ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)

         !-- These were uninitialized : meaningful only when we compare error differences
         PCOVPTOT(:,:,IBL) = 0.0_JPRB
         TENDENCY_LOC(IBL)%cld(:,:,NCLV) = 0.0_JPRB

         ! Fill in ZQSAT
         CALL SATUR (1, ICEND, NPROMA, 1, NLEV, .TRUE., &
              & PAP(:,:,IBL), PT(:,:,IBL), ZQSAT(:,:), 2) 

         CALL CLOUDSC2 ( &
              &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
              & PTSPHY,&
              & PAPH(:,:,IBL),  PAP(:,:,IBL), &
              & PQ(:,:,IBL), ZQSAT(:,:), PT(:,:,IBL), &
              & PCLV(:,:,NCLDQL,IBL), PCLV(:,:,NCLDQI,IBL), &
              & PLUDE(:,:,IBL), PLU(:,:,IBL), PMFU(:,:,IBL), PMFD(:,:,IBL),&
              &  TENDENCY_LOC(IBL)%T, TENDENCY_CML(IBL)%T, &
              &  TENDENCY_LOC(IBL)%Q, TENDENCY_CML(IBL)%Q, &
              &  TENDENCY_LOC(IBL)%CLD(:,:,NCLDQL), TENDENCY_CML(IBL)%CLD(:,:,NCLDQL), &
              &  TENDENCY_LOC(IBL)%CLD(:,:,NCLDQI), TENDENCY_CML(IBL)%CLD(:,:,NCLDQI), &
              &  PSUPSAT(:,:,IBL), &
              &  PA(:,:,IBL), PFPLSL(:,:,IBL),   PFPLSN(:,:,IBL), &
              &  PFHPSL(:,:,IBL),   PFHPSN(:,:,IBL), PCOVPTOT(:,:,IBL))

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
         ZSUPSAT = PSUPSAT(:,:,IBL)*0.01_JPRB
         ! Call TL
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
           ZTENI_T5(:,:)=TENDENCY_CML(IBL)%T(:,:) + ZLAMBDA*ZTENI_T(:,:)
           ZTENI_Q5(:,:)=TENDENCY_CML(IBL)%Q(:,:) + ZLAMBDA*ZTENI_Q(:,:)
           ZTENI_L5(:,:)=TENDENCY_CML(IBL)%CLD(:,:,NCLDQL) + ZLAMBDA*ZTENI_L(:,:)
           ZTENI_I5(:,:)=TENDENCY_CML(IBL)%CLD(:,:,NCLDQI) + ZLAMBDA*ZTENI_I(:,:)
           PSUPSAT5(:,:)=PSUPSAT(:,:,IBL) + ZLAMBDA*ZSUPSAT(:,:)
           ! Call the NL code
           CALL CLOUDSC2 ( &
              &  1, ICEND, NPROMA, 1, NLEV, LDRAIN1D, &
              & PTSPHY,&
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
              & PFHPSL5(:,:),   PFHPSN5(:,:), PCOVPTOT5(:,:))

           ! Compute final test norm
           ZCOUNT=0._JPRB
           ZNORM= 0._JPRB 
           IF (SUM(ZTENO_T(1:ICEND,1:NLEV)*ZLAMBDA) /= 0._JPRB) THEN
             ZCOUNT=ZCOUNT+1._JPRB
             ZNORM= ZNORM &
             & +ABS(SUM(TENDENCY_LOC(IBL)%T(1:ICEND,1:NLEV)-ZTENO_T5(1:ICEND,1:NLEV)) &
             &  /SUM(ZTENO_T(1:ICEND,1:NLEV)*ZLAMBDA))
           ENDIF
           IF (SUM(ZTENO_Q(1:ICEND,1:NLEV)*ZLAMBDA) /= 0._JPRB) THEN
             ZCOUNT=ZCOUNT+1._JPRB
             ZNORM= ZNORM &
             & +ABS(SUM(TENDENCY_LOC(IBL)%Q(1:ICEND,1:NLEV)-ZTENO_Q5(1:ICEND,1:NLEV)) &
             &  /SUM(ZTENO_Q(1:ICEND,1:NLEV)*ZLAMBDA))
           ENDIF
           IF (SUM(ZTENO_L(1:ICEND,1:NLEV)*ZLAMBDA) /= 0._JPRB) THEN
             ZCOUNT=ZCOUNT+1._JPRB
             ZNORM= ZNORM &
             & +ABS(SUM(TENDENCY_LOC(IBL)%CLD(1:ICEND,1:NLEV,NCLDQL) &
             &            -ZTENO_L5(1:ICEND,1:NLEV)) &
             &  /SUM(ZTENO_L(1:ICEND,1:NLEV)*ZLAMBDA))
           ENDIF
           IF (SUM(ZCLC(1:ICEND,1:NLEV)*ZLAMBDA) /= 0._JPRB ) THEN
             ZCOUNT=ZCOUNT+1._JPRB
             ZNORM= ZNORM &
             & +ABS(SUM(TENDENCY_LOC(IBL)%CLD(1:ICEND,1:NLEV,NCLDQI) &
             &            -ZTENO_I5(1:ICEND,1:NLEV)) &
             &  /SUM(ZTENO_I(1:ICEND,1:NLEV)*ZLAMBDA))
           ENDIF
           IF (SUM(ZCLC(1:ICEND,1:NLEV)*ZLAMBDA) /= 0._JPRB ) THEN
             ZCOUNT=ZCOUNT+1._JPRB
             ZNORM= ZNORM &
             & +ABS(SUM(PA(1:ICEND,1:NLEV,IBL)-PA5(1:ICEND,1:NLEV)) &
             &  /SUM(ZCLC(1:ICEND,1:NLEV)*ZLAMBDA))
           ENDIF
           IF (SUM(ZFPLSL(1:ICEND,1:NLEV+1)*ZLAMBDA) /= 0._JPRB ) THEN
             ZCOUNT=ZCOUNT+1._JPRB
             ZNORM= ZNORM &
             & +ABS(SUM(PFPLSL(1:ICEND,1:NLEV+1,IBL)-PFPLSL5(1:ICEND,1:NLEV+1)) &
             &  /SUM(ZFPLSL(1:ICEND,1:NLEV+1)*ZLAMBDA))
           ENDIF
           IF (SUM(ZFPLSN(1:ICEND,1:NLEV+1)*ZLAMBDA) /= 0._JPRB ) THEN
             ZCOUNT=ZCOUNT+1._JPRB
             ZNORM= ZNORM &
             & +ABS(SUM(PFPLSN(1:ICEND,1:NLEV+1,IBL)-PFPLSN5(1:ICEND,1:NLEV+1)) &
             &  /SUM(ZFPLSN(1:ICEND,1:NLEV+1)*ZLAMBDA))
           ENDIF
           IF (SUM(ZFHPSL(1:ICEND,1:NLEV+1)*ZLAMBDA) /= 0._JPRB ) THEN
             ZCOUNT=ZCOUNT+1._JPRB
             ZNORM= ZNORM &
             & +ABS(SUM(PFHPSL(1:ICEND,1:NLEV+1,IBL)-PFHPSL5(1:ICEND,1:NLEV+1)) &
             &  /SUM(ZFHPSL(1:ICEND,1:NLEV+1)*ZLAMBDA))
           ENDIF
           IF (SUM(ZFHPSN(1:ICEND,1:NLEV+1)*ZLAMBDA) /= 0._JPRB ) THEN
             ZCOUNT=ZCOUNT+1._JPRB
             ZNORM= ZNORM &
             & +ABS(SUM(PFHPSN(1:ICEND,1:NLEV+1,IBL)-PFHPSN5(1:ICEND,1:NLEV+1)) &
             &  /SUM(ZFHPSN(1:ICEND,1:NLEV+1)*ZLAMBDA))
           ENDIF
           IF (SUM(ZCOVPTOT(1:ICEND,1:NLEV)*ZLAMBDA) /= 0._JPRB ) THEN
             ZCOUNT=ZCOUNT+1._JPRB
             ZNORM= ZNORM &
             & +ABS(SUM(PCOVPTOT(1:ICEND,1:NLEV,IBL)-PCOVPTOT5(1:ICEND,1:NLEV)) &
             &  /SUM(ZCOVPTOT(1:ICEND,1:NLEV)*ZLAMBDA))
           ENDIF

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

      !-- The "nowait" is here to get correct local timings (tloc) per thread
      !   i.e. we should not wait for slowest thread to finish before measuring tloc
      !$omp end do nowait

      CALL TIMER%THREAD_END(TID)

      !$omp end parallel

      CALL TIMER%END()

      CALL TIMER%PRINT_PERFORMANCE(NPROMA, NGPBLKS, NGPTOT)

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
