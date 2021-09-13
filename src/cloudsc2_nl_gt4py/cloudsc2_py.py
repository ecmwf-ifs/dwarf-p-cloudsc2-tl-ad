import numpy as np
import math

def foealfa(ptare, yrethf):
  return min(1.0,((max(yrethf.rtice,min(yrethf.rtwat,ptare))-yrethf.rtice)*yrethf.rtwat_rtice_r)**2)

def foeewm(ptare, yrethf, yrmcst):
  return yrethf.r2es*(foealfa(ptare, yrethf)*np.exp(yrethf.r3les*(ptare-yrmcst.rtt)/(ptare-yrethf.r4les)) + \
                      (1.0-foealfa(ptare, yrethf))*np.exp(yrethf.r3ies*(ptare-yrmcst.rtt)/(ptare-yrethf.r4ies)))


def satur(kidia, kfdia, klon, ktdia, klev, ldphylin, paprsf, pt, pqsat, kflag, yrethf, yrmcst):

    #----------------------------------------------------------------------
    #*    1.           DEFINE CONSTANTS
    zqmax = 0.5

    #     *
    #----------------------------------------------------------------------
    #     *    2.           CALCULATE SATURATION SPECIFIC HUMIDITY
    #                       --------------------------------------

    if ldphylin:
        for jk in range(klev):
            for jl in range(kfdia):
                ztarg = pt[jk,jl]
                zalfa = foealfa(ztarg, yrethf)

                zfoeewl = yrethf.r2es*np.exp(yrethf.r3les*(ztarg-yrmcst.rtt)/(ztarg-yrethf.r4les))
                zfoeewi = yrethf.r2es*np.exp(yrethf.r3ies*(ztarg-yrmcst.rtt)/(ztarg-yrethf.r4ies))
                zfoeew = zalfa*zfoeewl+(1.0-zalfa)*zfoeewi

                zqs    = zfoeew/paprsf[jk,jl]
                if zqs > zqmax:
                    zqs=zqmax

                zcor = 1.0/(1.0-yrmcst.retv*zqs)
                pqsat[jk,jl]=zqs*zcor

    else:
        for jk in range(klev):
            for jl in range(kfdia):
                if(kflag == 1):
                    zew  = foeewmcu(pt[jk,jl])
                else:
                    zew  = foeewm(pt[jk,jl])

                zqs  = zew/paprsf[jk,jl]
                zqs  = min(zqmax,zqs)
                zcor = 1.0/(1.0-yrmcst.retv*zqs)
                pqsat[jk,jl]=zqs*zcor


def cloudsc2_py(kidia: None, kfdia: None, klon: None, ktdia: None, klev: None, ldrain1d: None, \
  ptsphy: None, paphp1: None, papp1: None, pqm1: None, pqs: None, ptm1: None, pl: None, pi: None, \
  plude: None, plu: None, pmfu: None, pmfd: None, ptent: None, pgtent: None, ptenq: None, \
  pgtenq: None, ptenl: None, pgtenl: None, pteni: None, pgteni: None, psupsat: None, pclc: None, \
  pfplsl: None, pfplsn: None, pfhpsl: None, pfhpsn: None, pcovptot: None, \
  yrecldp: None, yrecld: None, yrmcst: None, yrethf:None, yrephli: None):
  
  #     -----------------------------------------------------------------
  
  #*       0.1   ARGUMENTS.
  #              ----------
  
  #     -----------------------------------------------------------------
  
  #*       0.2   LOCAL ARRAYS.
  #              -------------
  
  
  #=======================================================================
  #     TUNABLE CONSTANTS (to be moved to include files later)
  #=======================================================================
  
  # zscal is a scale factor that linearly reduces the variance between
  # qv=qv-crit and qv-qsat
  # 0 = No scaling
  # 1 = full scaling (i.e. variance=0 when qv=qsat)
  
  zscal = 0.9
  
  #=======================================================================
  
  zrfl = np.ndarray(order="C", shape=(klon,))
  zsfl = np.ndarray(order="C", shape=(klon,))
  zrfln = np.ndarray(order="C", shape=(klon,))
  zsfln = np.ndarray(order="C", shape=(klon,))
  zgdp = np.ndarray(order="C", shape=(klon,))
  zdqdt = np.ndarray(order="C", shape=(klon,))
  zdtdt = np.ndarray(order="C", shape=(klon,))
  zdldt = np.ndarray(order="C", shape=(klon,))
  zdidt = np.ndarray(order="C", shape=(klon,))
  zqcrit = np.ndarray(order="C", shape=(klon,))
  zcovpclr = np.ndarray(order="C", shape=(klon,))
  zcovptot = np.ndarray(order="C", shape=(klon,))
  zdqsdtemp = np.ndarray(order="C", shape=(klon,))
  zcorqs = np.ndarray(order="C", shape=(klon,))
  zdtgdp = np.ndarray(order="C", shape=(klon,))
  zqold = np.ndarray(order="C", shape=(klon,))
  zpp = np.ndarray(order="C", shape=(klon,))
  zdq = np.ndarray(order="C", shape=(klon,))
  zqlim = np.ndarray(order="C", shape=(klon,))
  zscalm = np.ndarray(order="C", shape=(klev,))
  zqsat = np.ndarray(order="C", shape=(klon,))
  zfoeew = np.ndarray(order="C", shape=(klon,))
  zfwat = np.ndarray(order="C", shape=(klon,))
  
  
  ztrpaus = np.ndarray(order="C", shape=(klon,))
  
  
  llflag = np.ndarray(order="C", shape=(klon,))
  #REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
  
  #     ------------------------------------------------------------------
  ztp1 = np.ndarray(order="F", shape=(klev, klon,))
  zqp1 = np.ndarray(order="F", shape=(klev, klon,))
  zl = np.ndarray(order="F", shape=(klev, klon,))
  zi = np.ndarray(order="F", shape=(klev, klon,))
  zlude = np.ndarray(order="F", shape=(klev, klon,))
  zqc = np.ndarray(order="F", shape=(klev, klon,))
  zqlwc = np.ndarray(order="F", shape=(klev, klon,))
  zqiwc = np.ndarray(order="F", shape=(klev, klon,))
  zdp = np.ndarray(order="F", shape=(klev, klon,))
  zlsdcp = np.ndarray(order="F", shape=(klev, klon,))
  zlfdcp = np.ndarray(order="F", shape=(klev, klon,))
  zlvdcp = np.ndarray(order="F", shape=(klev, klon,))
  zrfreeze = np.ndarray(order="F", shape=(klev, klon,))
  zcondl = np.ndarray(order="F", shape=(klev, klon,))
  zcondi = np.ndarray(order="F", shape=(klev, klon,))
  zevapr = np.ndarray(order="F", shape=(klev, klon,))
  zevaps = np.ndarray(order="F", shape=(klev, klon,))
  #     ------------------------------------------------------------------
  
  #IF (LHOOK) CALL DR_HOOK('CLOUDSC2',0,ZHOOK_HANDLE)
  
  #*         1.     SET-UP INPUT QUANTITIES
  #                 -----------------------
  
  #*         1.1    Set-up tunning parameters
  
  # set up constants required
  
  zckcodtl = 2.0*yrecldp.rkconv*ptsphy
  zckcodti = 5.0*yrecldp.rkconv*ptsphy
  zcons2 = 1.0 / ((ptsphy*yrmcst.rg))
  zcons3 = yrmcst.rlvtt / yrmcst.rcpd
  zmeltp2 = yrmcst.rtt + 2.0
  zqtmst = 1.0 / ptsphy
  
  zqmax = 0.5
  zeps1 = 1.E-12
  zeps2 = 1.E-10
  
  #     --------------------------------------------------------------------
  
  #*         2.1    COMPUTE CRITICAL RELATIVE HUMIDITY AND RELATIVE HUMIDITY
  #                 --------------------------------------------------------
  
  # first guess values for T, q, ql and qi
  
  for jk in range(klev):
    for jl in range(kfdia):
      ztp1[jk, jl] = ptm1[jk, jl] + ptsphy*pgtent[jk, jl]
      zqp1[jk, jl] = \
        pqm1[jk, jl] + ptsphy*pgtenq[jk, jl] + psupsat[jk, jl]
      zl[jk, jl] = pl[jk, jl] + ptsphy*pgtenl[jk, jl]
      zi[jk, jl] = pi[jk, jl] + ptsphy*pgteni[jk, jl]
  
  for jk in range(klev):
    
    # Parameter for cloud formation

    zscalm[jk] = zscal*max((yrecld.ceta[jk] - 0.2), zeps1)**0.2
    
    for jl in range(kfdia):
      
      # thermodynamic constants
      
      zdp[jk, jl] = paphp1[jk+1, jl] - paphp1[jk, jl]
      zzz = 1.0 / (yrmcst.rcpd + yrmcst.rcpd*yrethf.rvtmp2*zqp1[jk, jl])
      zlfdcp[jk, jl] = yrmcst.rlmlt*zzz
      zlsdcp[jk, jl] = yrmcst.rlstt*zzz
      zlvdcp[jk, jl] = yrmcst.rlvtt*zzz
      llflag[jl] = True
  
  #     ------------------------------------------------------------------
  
  #*         2.2    INITIALIZATION OF CLOUD AND PRECIPITATION ARRAYS
  #                 ------------------------------------------------
  
  #       Clear cloud and freezing arrays
  
  for jk in range(klev):
    for jl in range(kfdia):
      pclc[jk, jl] = 0.0
      zqc[jk, jl] = 0.0
      zqlwc[jk, jl] = 0.0
      zqiwc[jk, jl] = 0.0
      zrfreeze[jk, jl] = 0.0
      zcondl[jk, jl] = 0.0
      zcondi[jk, jl] = 0.0
      zevapr[jk, jl] = 0.0
      zevaps[jk, jl] = 0.0
      pcovptot[jk, jl] = 0.0
  
  #       Set to zero precipitation fluxes at the top
  
  for jl in range(kfdia):
    zrfl[jl] = 0.0
    zsfl[jl] = 0.0
    pfplsl[0, jl] = 0.0
    pfplsn[0, jl] = 0.0
    zcovptot[jl] = 0.0
    zcovpclr[jl] = 0.0
  
  # Eta value at tropopause
  for jl in range(kfdia):
    ztrpaus[jl] = 0.1
  for jk in range(klev - 1):
    for jl in range(kfdia):
      llo1 = \
        yrecld.ceta[jk] > 0.1 and yrecld.ceta[jk] < 0.4 and ztp1[jk, jl] > ztp1[jk+1, jl]
      if llo1:
        ztrpaus[jl] = yrecld.ceta[jk]
  
  #     ------------------------------------------------------------------
  
  #*        3. COMPUTE LAYER CLOUD AMOUNTS
  #            ---------------------------
  
  # Large loop over KLEV
  # Calculates
  #   1. diagnostic CC and QL
  #   2. Convective CC and QL
  #   3. Rainfall
  
  for jk in range(klev):
    
    #       3.1   INITIALIZATION
    
    for jl in range(kfdia):
      
      #-----------------------------------
      # calculate dqs/dT correction factor
      #-----------------------------------
      
      if yrephli.lphylin or ldrain1d:
        zoealfaw = 0.545*(np.tanh(0.17*(ztp1[jk, jl] - yrephli.rlptrc)) + 1.0)
        if ztp1[jk, jl] < yrmcst.rtt:
          zfwat[jl] = zoealfaw
          z3es = yrethf.r3ies
          z4es = yrethf.r4ies
        else:
          zfwat[jl] = 1.0
          z3es = yrethf.r3les
          z4es = yrethf.r4les
        zfoeew[jl] = \
          yrethf.r2es*np.exp((z3es*(ztp1[jk, jl] - yrmcst.rtt)) / (ztp1[jk, jl] - z4es))
        zesdp = zfoeew[jl] / papp1[jk, jl]
        if zesdp > zqmax:
          zesdp = zqmax
      else:
        zfwat[jl] = foealfa(ztp1[jk,jl], yrethf)
        zfoeew[jl] = foeewm(ztp1[jk,jl], yrethf, yrmcst)
        zesdp = zfoeew[jl] / papp1[jk, jl]
      zfacw = yrethf.r5les / ((ztp1[jk, jl] - yrethf.r4les)**2)
      zfaci = yrethf.r5ies / ((ztp1[jk, jl] - yrethf.r4ies)**2)
      zfac = zfwat[jl]*zfacw + (1.0 - zfwat[jl])*zfaci
      zcor = 1.0 / (1.0 - yrmcst.retv*zesdp)
      zdqsdtemp[jl] = zfac*zcor*pqs[jk, jl]
      zcorqs[jl] = 1.0 + zcons3*zdqsdtemp[jl]
      
      # use clipped state
      
      zqlim[jl] = zqp1[jk, jl]
      if zqp1[jk, jl] > pqs[jk, jl]:
        zqlim[jl] = pqs[jk, jl]
      
      # set up critical value of humidity
      
      zeta3 = ztrpaus[jl]
      zrh1 = 1.0
      zrh2 = 0.35 + 0.14*((zeta3 - 0.25) / 0.15)**2 + (0.04*min(zeta3 - 0.25, 0.0)) / 0.15
      zrh3 = 1.0
      zdeta2 = 0.3
      zdeta1 = 0.09 + (0.16*(0.4 - zeta3)) / 0.3
      if yrecld.ceta[jk] < zeta3:
        zcrh2 = zrh3
      elif yrecld.ceta[jk] >= zeta3 and yrecld.ceta[jk] < (zeta3 + zdeta2):
        zcrh2 = zrh3 + (zrh2 - zrh3)*((yrecld.ceta[jk] - zeta3) / zdeta2)
      elif yrecld.ceta[jk] >= (zeta3 + zdeta2) and yrecld.ceta[jk] < (1.0 - zdeta1):
        zcrh2 = zrh2
      elif yrecld.ceta[jk] >= (1.0 - zdeta1):
        zcrh2 = zrh1 + (zrh2 - zrh1)*((1.0 - yrecld.ceta[jk]) / zdeta1)**0.5
      # Allow ice supersaturation at cold temperatures
      if ztp1[jk, jl] < yrethf.rtice:
        zsupsat = 1.8 - 3.E-03*ztp1[jk, jl]
      else:
        zsupsat = 1.0
      zqsat[jl] = pqs[jk, jl]*zsupsat
      zqcrit[jl] = zcrh2*zqsat[jl]
    
    # Simple UNIFORM distribution of total water from Letreut & Li (90)
    
    for jl in range(kfdia):
      zqt = zqp1[jk, jl] + zl[jk, jl] + zi[jk, jl]
      if zqt <= zqcrit[jl]:
        pclc[jk, jl] = 0.0
        zqc[jk, jl] = 0.0
      elif zqt >= zqsat[jl]:
        pclc[jk, jl] = 1.0
        zqc[jk, jl] = (1.0 - zscalm[jk])*(zqsat[jl] - zqcrit[jl])
      else:
        zqpd = zqsat[jl] - zqt
        zqcd = zqsat[jl] - zqcrit[jl]
        pclc[jk, jl] = \
          1.0 - np.sqrt(zqpd / (zqcd - zscalm[jk]*(zqt - zqcrit[jl])))
        zqc[jk, jl] = \
          (zscalm[jk]*zqpd + (1.0 - zscalm[jk])*zqcd)*pclc[jk, jl]**2
    
    # Add convective component
    
    for jl in range(kfdia):
      zgdp[jl] = yrmcst.rg / (paphp1[jk+1, jl] - paphp1[jk, jl])
      zlude[jk, jl] = plude[jk, jl]*ptsphy*zgdp[jl]
      if jk < klev:
        llo1 = zlude[jk, jl] >= yrecldp.rlmin and plu[jk+1, jl] >= zeps2
      else:
        llo1 = False
      if llo1:
        pclc[jk, jl] = pclc[jk, jl] + (1.0 - pclc[jk, jl])*(1.0 - \
          np.exp(-zlude[jk, jl] / plu[jk+1, jl]))
        zqc[jk, jl] = zqc[jk, jl] + zlude[jk, jl]
    
    # Add compensating subsidence component
    
    for jl in range(kfdia):
      zfac1 = 1.0 / ((yrmcst.rd*ztp1[jk, jl]))
      zrho = papp1[jk, jl]*zfac1
      zfac2 = 1.0 / (papp1[jk, jl] - yrmcst.retv*zfoeew[jl])
      zrodqsdp = -zrho*pqs[jk, jl]*zfac2
      zldcp = \
        zfwat[jl]*zlvdcp[jk, jl] + (1.0 - zfwat[jl])*zlsdcp[jk, jl]
      zfac3 = 1.0 / (1.0 + zldcp*zdqsdtemp[jl])
      dtdzmo = yrmcst.rg*(1.0 / yrmcst.rcpd - zldcp*zrodqsdp)*zfac3
      zdqsdz = zdqsdtemp[jl]*dtdzmo - yrmcst.rg*zrodqsdp
      zfac4 = 1.0 / zrho
      zdqc = min(zdqsdz*(pmfu[jk, jl] + pmfd[jk, jl])*ptsphy*zfac4, zqc[jk, jl])
      zqc[jk, jl] = zqc[jk, jl] - zdqc
    
    # New cloud liquid/ice contents and condensation rates (liquid/ice)
    
    for jl in range(kfdia):
      zqlwc[jk, jl] = zqc[jk, jl]*zfwat[jl]
      zqiwc[jk, jl] = zqc[jk, jl]*(1.0 - zfwat[jl])
      zcondl[jk, jl] = (zqlwc[jk, jl] - zl[jk, jl])*zqtmst
      zcondi[jk, jl] = (zqiwc[jk, jl] - zi[jk, jl])*zqtmst
    
    
    # Calculate precipitation overlap.
    # Simple form based on Maximum Overlap.
    
    for jl in range(kfdia):
      if pclc[jk, jl] > zcovptot[jl]:
        # total rain frac
        zcovptot[jl] = pclc[jk, jl]
      zcovpclr[jl] = zcovptot[jl] - pclc[jk, jl]        # clear sky frac
      zcovpclr[jl] = max(zcovpclr[jl], 0.0)
    
    #*         3.3    CALCULATE PRECIPITATION
    
    # Melting of incoming snow
    
    for jl in range(kfdia):
      if zsfl[jl] != 0.0:
        zcons = (zcons2*zdp[jk, jl]) / zlfdcp[jk, jl]
        zsnmlt = min(zsfl[jl], zcons*max(0.0, (ztp1[jk, jl] - zmeltp2)))
        zrfln[jl] = zrfl[jl] + zsnmlt
        zsfln[jl] = zsfl[jl] - zsnmlt
        ztp1[jk, jl] = ztp1[jk, jl] - zsnmlt / zcons
      else:
        zrfln[jl] = zrfl[jl]
        zsfln[jl] = zsfl[jl]
    
    for jl in range(kfdia):
      
      #   Diagnostic calculation of rain production from cloud liquid water
      
      if pclc[jk, jl] > zeps2:
        # if yrphnc.levapls2 or ldrain1d:
        if False or ldrain1d:
          zlcrit = 1.9*yrecldp.rclcrit
        else:
          zlcrit = yrecldp.rclcrit*2.
        zcldl = zqlwc[jk, jl] / pclc[jk, jl]          # in-cloud liquid
        zd = zckcodtl*(1.0 - np.exp(-(zcldl / zlcrit)**2))
        zlnew = pclc[jk, jl]*zcldl*np.exp(-zd)
        zprr = zqlwc[jk, jl] - zlnew
        zqlwc[jk, jl] = zqlwc[jk, jl] - zprr
      else:
        zprr = 0.0
      
      #   Diagnostic calculation of snow production from cloud ice
      
      if pclc[jk, jl] > zeps2:
        # if yrphnc.levapls2 or ldrain1d:
        if False or ldrain1d:
          zlcrit = 1.E-04
        else:
          zlcrit = yrecldp.rclcrit*2.
        zcldi = zqiwc[jk, jl] / pclc[jk, jl]          # in-cloud ice
        zd = zckcodti*np.exp(0.025*(ztp1[jk, jl] - yrmcst.rtt))*(1.0 - np.exp(-(zcldi / \
          zlcrit)**2))
        zinew = pclc[jk, jl]*zcldi*np.exp(-zd)
        zprs = zqiwc[jk, jl] - zinew
        zqiwc[jk, jl] = zqiwc[jk, jl] - zprs
      else:
        zprs = 0.0
      
      #   New precipitation (rain + snow)
      
      zdr = zcons2*zdp[jk, jl]*(zprr + zprs)
      
      #   Rain fraction (different from cloud liquid water fraction!)
      
      if ztp1[jk, jl] < yrmcst.rtt:
        zrfreeze[jk, jl] = zcons2*zdp[jk, jl]*zprr
        zfwatr = 0.0
      else:
        zfwatr = 1.0
      
      zrn = zfwatr*zdr
      zsn = (1.0 - zfwatr)*zdr
      zrfln[jl] = zrfln[jl] + zrn
      zsfln[jl] = zsfln[jl] + zsn
      
      #   Precip evaporation
      
      zprtot = zrfln[jl] + zsfln[jl]
      llo2 = zprtot > zeps2 and zcovpclr[jl] > zeps2 and ldrain1d
      if llo2:
        
        zpreclr = (zprtot*zcovpclr[jl]) / zcovptot[jl]
        
        #     This is the humidity in the moistest zcovpclr region
        
        zqe = pqs[jk, jl] - ((pqs[jk, jl] - zqlim[jl])*zcovpclr[jl]) \
          / (1.0 - pclc[jk, jl])**2
        zbeta = yrmcst.rg*yrecldp.rpecons*(((np.sqrt(papp1[jk, jl] / paphp1[klev, jl]) / \
          5.09E-3)*zpreclr) / zcovpclr[jl])**0.5777
        
        #     implicit solution:
        zb = (ptsphy*zbeta*(pqs[jk, jl] - zqe)) / (1.0 + zbeta*ptsphy*zcorqs[jl])
        
        #     exact solution:
        #     ZB=(PQS(JL,JK)-ZQE)*(_ONE_-EXP(-ZBETA*ZCORQS(JL)*PTSPHY))/ZCORQS(JL)
        
        zdtgdp[jl] = (ptsphy*yrmcst.rg) / (paphp1[jk+1, jl] - paphp1[jk, jl])
        
        zdpr = (zcovpclr[jl]*zb) / zdtgdp[jl]
        zdpr = min(zdpr, zpreclr)
        zpreclr = zpreclr - zdpr          # take away from clr sky flux
        if zpreclr <= 0.0:
          zcovptot[jl] = pclc[jk, jl]
        #reset
        pcovptot[jk, jl] = zcovptot[jl]
        
        # warm proportion
        zevapr[jk, jl] = (zdpr*zrfln[jl]) / zprtot
        zrfln[jl] = zrfln[jl] - zevapr[jk, jl]
        
        # ice proportion
        zevaps[jk, jl] = (zdpr*zsfln[jl]) / zprtot
        zsfln[jl] = zsfln[jl] - zevaps[jk, jl]
      
    
    # Update of T and Q tendencies due to:
    #  - condensation/evaporation of cloud liquid water/ice
    #  - detrainment of convective cloud condensate
    #  - evaporation of precipitation
    #  - freezing of rain (impact on T only).
    
    for jl in range(kfdia):
      zdqdt[jl] = -(zcondl[jk, jl] + zcondi[jk, jl]) + (plude[jk, \
        jl] + zevapr[jk, jl] + zevaps[jk, jl])*zgdp[jl]
      
      zdtdt[jl] = zlvdcp[jk, jl]*zcondl[jk, jl] + zlsdcp[jk, jl] \
        *zcondi[jk, jl] - (zlvdcp[jk, jl]*zevapr[jk, jl] + \
        zlsdcp[jk, jl]*zevaps[jk, jl] + plude[jk, jl]*(zfwat[jl] \
        *zlvdcp[jk, jl] + (1.0 - zfwat[jl])*zlsdcp[jk, jl]) - \
        (zlsdcp[jk, jl] - zlvdcp[jk, jl])*zrfreeze[jk, jl])*zgdp[jl]
      
      # first guess T and Q
      ztp1[jk, jl] = ztp1[jk, jl] + ptsphy*zdtdt[jl]
      zqp1[jk, jl] = zqp1[jk, jl] + ptsphy*zdqdt[jl]
      
      zpp[jl] = papp1[jk, jl]
      zqold[jl] = zqp1[jk, jl]
    
    # clipping of final qv
    
    # -----------------------------------
    # IK=JK
    # ICALL=0
    # CALL CUADJTQS ( KIDIA, KFDIA, KLON, KLEV, IK,&
    #   & ZPP  , ZTP1  , ZQP1 , LLFLAG, ICALL  )
    # -----------------------------------
    # Manually inlined CUADJTQS
    # -----------------------------------
    zqmax = 0.5
    for jl in range(kfdia):
      if ztp1[jk, jl] > yrmcst.rtt:
        z3es = yrethf.r3les
        z4es = yrethf.r4les
        z5alcp = yrethf.r5alvcp
        zaldcp = yrethf.ralvdcp
      else:
        z3es = yrethf.r3ies
        z4es = yrethf.r4ies
        z5alcp = yrethf.r5alscp
        zaldcp = yrethf.ralsdcp
      
      zqp = 1.0 / zpp[jl]
      ztarg = ztp1[jk, jl]
      zfoeew[jl] = yrethf.r2es*np.exp((z3es*(ztarg - yrmcst.rtt)) / (ztarg - z4es))
      zqsat[jl] = zqp*zfoeew[jl]
      if zqsat[jl] > zqmax:
        zqsat[jl] = zqmax
      zcor = 1.0 / (1.0 - yrmcst.retv*zqsat[jl])
      zqsat[jl] = zqsat[jl]*zcor
      z2s = z5alcp / (ztarg - z4es)**2
      zcond1 = (zqp1[jk, jl] - zqsat[jl]) / (1.0 + zqsat[jl]*zcor*z2s)
      ztp1[jk, jl] = ztp1[jk, jl] + zaldcp*zcond1
      zqp1[jk, jl] = zqp1[jk, jl] - zcond1
      ztarg = ztp1[jk, jl]
      zfoeew[jl] = yrethf.r2es*np.exp((z3es*(ztarg - yrmcst.rtt)) / (ztarg - z4es))
      zqsat[jl] = zqp*zfoeew[jl]
      if zqsat[jl] > zqmax:
        zqsat[jl] = zqmax
      zcor = 1.0 / (1.0 - yrmcst.retv*zqsat[jl])
      zqsat[jl] = zqsat[jl]*zcor
      z2s = z5alcp / (ztarg - z4es)**2
      zcond1 = (zqp1[jk, jl] - zqsat[jl]) / (1.0 + zqsat[jl]*zcor*z2s)
      ztp1[jk, jl] = ztp1[jk, jl] + zaldcp*zcond1
      zqp1[jk, jl] = zqp1[jk, jl] - zcond1
    # -----------------------------------
    
    for jl in range(kfdia):
      zdq[jl] = max(0.0, zqold[jl] - zqp1[jk, jl])
      zdr2 = zcons2*zdp[jk, jl]*zdq[jl]
      # Update rain fraction and freezing.
      # Note: impact of new temperature ZTP1 on ZFWAT is neglected here.
      if ztp1[jk, jl] < yrmcst.rtt:
        zrfreeze2 = zfwat[jl]*zdr2
        zfwatr = 0.0
      else:
        zrfreeze2 = 0.0
        zfwatr = 1.0
      zrn = zfwatr*zdr2
      zsn = (1.0 - zfwatr)*zdr2
      # Note: The extra condensation due to the adjustment goes directly to precipitation
      zcondl[jk, jl] = zcondl[jk, jl] + zfwatr*zdq[jl]*zqtmst
      zcondi[jk, jl] = zcondi[jk, jl] + (1.0 - zfwatr)*zdq[jl]*zqtmst
      zrfln[jl] = zrfln[jl] + zrn
      zsfln[jl] = zsfln[jl] + zsn
      zrfreeze[jk, jl] = zrfreeze[jk, jl] + zrfreeze2
    
    for jl in range(kfdia):
      zdqdt[jl] = -(zcondl[jk, jl] + zcondi[jk, jl]) + (plude[jk, \
        jl] + zevapr[jk, jl] + zevaps[jk, jl])*zgdp[jl]
      
      zdtdt[jl] = zlvdcp[jk, jl]*zcondl[jk, jl] + zlsdcp[jk,jl] * \
        zcondi[jk, jl] - (zlvdcp[jk, jl]*zevapr[jk, jl] + \
        zlsdcp[jk, jl]*zevaps[jk, jl] + plude[jk, jl]*(zfwat[jl] * \
        zlvdcp[jk, jl] + (1.0 - zfwat[jl])*zlsdcp[jk, jl]) - \
        (zlsdcp[jk, jl] - zlvdcp[jk, jl])*zrfreeze[jk, jl])*zgdp[jl]
      
      zdldt[jl] = (zqlwc[jk, jl] - zl[jk, jl])*zqtmst
      
      zdidt[jl] = (zqiwc[jk, jl] - zi[jk, jl])*zqtmst
      
      ptenq[jk, jl] = zdqdt[jl]
      ptent[jk, jl] = zdtdt[jl]
      ptenl[jk, jl] = zdldt[jl]
      pteni[jk, jl] = zdidt[jl]
      
      pfplsl[jk+1, jl] = zrfln[jl]
      pfplsn[jk+1, jl] = zsfln[jl]
    
    # record rain flux for next level
    
    for jl in range(kfdia):
      zrfl[jl] = zrfln[jl]
      zsfl[jl] = zsfln[jl]
    
  #jk
  
  #*     ENTHALPY FLUXES DUE TO PRECIPITATION
  #      ------------------------------------
  
  for jk in range(klev + 1):
    for jl in range(kfdia):
      pfhpsl[jk, jl] = -pfplsl[jk, jl]*yrmcst.rlvtt
      pfhpsn[jk, jl] = -pfplsn[jk, jl]*yrmcst.rlstt
  
  #     ------------------------------------------------------------------
  
  #IF (LHOOK) CALL DR_HOOK('CLOUDSC2',1,ZHOOK_HANDLE)
  return 
