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
        for jk in range(1, klev + 1):
            for jl in range(kidia, kfdia + 1):
                ztarg = pt[jk-1,jl-1]
                zalfa = foealfa(ztarg, yrethf)

                zfoeewl = yrethf.r2es*np.exp(yrethf.r3les*(ztarg-yrmcst.rtt)/(ztarg-yrethf.r4les))
                zfoeewi = yrethf.r2es*np.exp(yrethf.r3ies*(ztarg-yrmcst.rtt)/(ztarg-yrethf.r4ies))
                zfoeew = zalfa*zfoeewl+(1.0-zalfa)*zfoeewi

                zqs    = zfoeew/paprsf[jk-1,jl-1]
                if zqs > zqmax:
                    zqs=zqmax

                zcor = 1.0/(1.0-yrmcst.retv*zqs)
                pqsat[jk-1,jl-1]=zqs*zcor

    else:
        for jk in range(1, klev + 1):
            for jl in range(kidia, kfdia + 1):
                if(kflag == 1):
                    zew  = foeewmcu(pt[jk-1,jl-1])
                else:
                    zew  = foeewm(pt[jk-1,jl-1])

                zqs  = zew/paprsf[jk-1,jl-1]
                zqs  = min(zqmax,zqs)
                zcor = 1.0/(1.0-yrmcst.retv*zqs)
                pqsat[jk-1,jl-1]=zqs*zcor


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
  
  for jk in range(1, klev + 1):
    for jl in range(kidia, kfdia + 1):
      ztp1[jk-1, jl-1] = ptm1[jk-1, jl-1] + ptsphy*pgtent[jk-1, jl-1]
      zqp1[jk-1, jl-1] = \
        pqm1[jk-1, jl-1] + ptsphy*pgtenq[jk-1, jl-1] + psupsat[jk-1, jl-1]
      zl[jk-1, jl-1] = pl[jk-1, jl-1] + ptsphy*pgtenl[jk-1, jl-1]
      zi[jk-1, jl-1] = pi[jk-1, jl-1] + ptsphy*pgteni[jk-1, jl-1]
  
  for jk in range(1, klev + 1):
    
    # Parameter for cloud formation

    zscalm[jk-1] = zscal*max((yrecld.ceta[jk-1] - 0.2), zeps1)**0.2
    
    for jl in range(kidia, kfdia + 1):
      
      # thermodynamic constants
      
      zdp[jk-1, jl-1] = paphp1[jk, jl-1] - paphp1[jk-1, jl-1]
      zzz = 1.0 / (yrmcst.rcpd + yrmcst.rcpd*yrethf.rvtmp2*zqp1[jk-1, jl-1])
      zlfdcp[jk-1, jl-1] = yrmcst.rlmlt*zzz
      zlsdcp[jk-1, jl-1] = yrmcst.rlstt*zzz
      zlvdcp[jk-1, jl-1] = yrmcst.rlvtt*zzz
      llflag[jl-1] = True
  
  #     ------------------------------------------------------------------
  
  #*         2.2    INITIALIZATION OF CLOUD AND PRECIPITATION ARRAYS
  #                 ------------------------------------------------
  
  #       Clear cloud and freezing arrays
  
  for jk in range(1, klev + 1):
    for jl in range(kidia, kfdia + 1):
      pclc[jk-1, jl-1] = 0.0
      zqc[jk-1, jl-1] = 0.0
      zqlwc[jk-1, jl-1] = 0.0
      zqiwc[jk-1, jl-1] = 0.0
      zrfreeze[jk-1, jl-1] = 0.0
      zcondl[jk-1, jl-1] = 0.0
      zcondi[jk-1, jl-1] = 0.0
      zevapr[jk-1, jl-1] = 0.0
      zevaps[jk-1, jl-1] = 0.0
      pcovptot[jk-1, jl-1] = 0.0
  
  #       Set to zero precipitation fluxes at the top
  
  for jl in range(kidia, kfdia + 1):
    zrfl[jl-1] = 0.0
    zsfl[jl-1] = 0.0
    pfplsl[0, jl-1] = 0.0
    pfplsn[0, jl-1] = 0.0
    zcovptot[jl-1] = 0.0
    zcovpclr[jl-1] = 0.0
  
  # Eta value at tropopause
  for jl in range(kidia, kfdia + 1):
    ztrpaus[jl-1] = 0.1
  for jk in range(1, klev - 1 + 1):
    for jl in range(kidia, kfdia + 1):
      llo1 = \
        yrecld.ceta[jk-1] > 0.1 and yrecld.ceta[jk-1] < 0.4 and ztp1[jk-1, jl-1] > ztp1[jk, jl-1]
      if llo1:
        ztrpaus[jl-1] = yrecld.ceta[jk-1]
  
  #     ------------------------------------------------------------------
  
  #*        3. COMPUTE LAYER CLOUD AMOUNTS
  #            ---------------------------
  
  # Large loop over KLEV
  # Calculates
  #   1. diagnostic CC and QL
  #   2. Convective CC and QL
  #   3. Rainfall
  
  for jk in range(ktdia, klev + 1):
    
    #       3.1   INITIALIZATION
    
    for jl in range(kidia, kfdia + 1):
      
      #-----------------------------------
      # calculate dqs/dT correction factor
      #-----------------------------------
      
      if yrephli.lphylin or ldrain1d:
        zoealfaw = 0.545*(np.tanh(0.17*(ztp1[jk-1, jl-1] - yrephli.rlptrc)) + 1.0)
        if ztp1[jk-1, jl-1] < yrmcst.rtt:
          zfwat[jl-1] = zoealfaw
          z3es = yrethf.r3ies
          z4es = yrethf.r4ies
        else:
          zfwat[jl-1] = 1.0
          z3es = yrethf.r3les
          z4es = yrethf.r4les
        zfoeew[jl-1] = \
          yrethf.r2es*np.exp((z3es*(ztp1[jk-1, jl-1] - yrmcst.rtt)) / (ztp1[jk-1, jl-1] - z4es))
        zesdp = zfoeew[jl-1] / papp1[jk-1, jl-1]
        if zesdp > zqmax:
          zesdp = zqmax
      else:
        zfwat[jl-1] = foealfa(ztp1[jk-1,jl-1], yrethf)
        zfoeew[jl-1] = foeewm(ztp1[jk-1,jl-1], yrethf, yrmcst)
        zesdp = zfoeew[jl-1] / papp1[jk-1, jl-1]
      zfacw = yrethf.r5les / ((ztp1[jk-1, jl-1] - yrethf.r4les)**2)
      zfaci = yrethf.r5ies / ((ztp1[jk-1, jl-1] - yrethf.r4ies)**2)
      zfac = zfwat[jl-1]*zfacw + (1.0 - zfwat[jl-1])*zfaci
      zcor = 1.0 / (1.0 - yrmcst.retv*zesdp)
      zdqsdtemp[jl-1] = zfac*zcor*pqs[jk-1, jl-1]
      zcorqs[jl-1] = 1.0 + zcons3*zdqsdtemp[jl-1]
      
      # use clipped state
      
      zqlim[jl-1] = zqp1[jk-1, jl-1]
      if zqp1[jk-1, jl-1] > pqs[jk-1, jl-1]:
        zqlim[jl-1] = pqs[jk-1, jl-1]
      
      # set up critical value of humidity
      
      zeta3 = ztrpaus[jl-1]
      zrh1 = 1.0
      zrh2 = 0.35 + 0.14*((zeta3 - 0.25) / 0.15)**2 + (0.04*min(zeta3 - 0.25, 0.0)) / 0.15
      zrh3 = 1.0
      zdeta2 = 0.3
      zdeta1 = 0.09 + (0.16*(0.4 - zeta3)) / 0.3
      if yrecld.ceta[jk-1] < zeta3:
        zcrh2 = zrh3
      elif yrecld.ceta[jk-1] >= zeta3 and yrecld.ceta[jk-1] < (zeta3 + zdeta2):
        zcrh2 = zrh3 + (zrh2 - zrh3)*((yrecld.ceta[jk-1] - zeta3) / zdeta2)
      elif yrecld.ceta[jk-1] >= (zeta3 + zdeta2) and yrecld.ceta[jk-1] < (1.0 - zdeta1):
        zcrh2 = zrh2
      elif yrecld.ceta[jk-1] >= (1.0 - zdeta1):
        zcrh2 = zrh1 + (zrh2 - zrh1)*((1.0 - yrecld.ceta[jk-1]) / zdeta1)**0.5
      # Allow ice supersaturation at cold temperatures
      if ztp1[jk-1, jl-1] < yrethf.rtice:
        zsupsat = 1.8 - 3.E-03*ztp1[jk-1, jl-1]
      else:
        zsupsat = 1.0
      zqsat[jl-1] = pqs[jk-1, jl-1]*zsupsat
      zqcrit[jl-1] = zcrh2*zqsat[jl-1]
    
    # Simple UNIFORM distribution of total water from Letreut & Li (90)
    
    for jl in range(kidia, kfdia + 1):
      zqt = zqp1[jk-1, jl-1] + zl[jk-1, jl-1] + zi[jk-1, jl-1]
      if zqt <= zqcrit[jl-1]:
        pclc[jk-1, jl-1] = 0.0
        zqc[jk-1, jl-1] = 0.0
      elif zqt >= zqsat[jl-1]:
        pclc[jk-1, jl-1] = 1.0
        zqc[jk-1, jl-1] = (1.0 - zscalm[jk-1])*(zqsat[jl-1] - zqcrit[jl-1])
      else:
        zqpd = zqsat[jl-1] - zqt
        zqcd = zqsat[jl-1] - zqcrit[jl-1]
        pclc[jk-1, jl-1] = \
          1.0 - np.sqrt(zqpd / (zqcd - zscalm[jk-1]*(zqt - zqcrit[jl-1])))
        zqc[jk-1, jl-1] = \
          (zscalm[jk-1]*zqpd + (1.0 - zscalm[jk-1])*zqcd)*pclc[jk-1, jl-1]**2
    
    # Add convective component
    
    for jl in range(kidia, kfdia + 1):
      zgdp[jl-1] = yrmcst.rg / (paphp1[jk, jl-1] - paphp1[jk-1, jl-1])
      zlude[jk-1, jl-1] = plude[jk-1, jl-1]*ptsphy*zgdp[jl-1]
      if jk < klev:
        llo1 = zlude[jk-1, jl-1] >= yrecldp.rlmin and plu[jk, jl-1] >= zeps2
      else:
        llo1 = False
      if llo1:
        pclc[jk-1, jl-1] = pclc[jk-1, jl-1] + (1.0 - pclc[jk-1, jl-1])*(1.0 - \
          np.exp(-zlude[jk-1, jl-1] / plu[jk, jl-1]))
        zqc[jk-1, jl-1] = zqc[jk-1, jl-1] + zlude[jk-1, jl-1]
    
    # Add compensating subsidence component
    
    for jl in range(kidia, kfdia + 1):
      zfac1 = 1.0 / ((yrmcst.rd*ztp1[jk-1, jl-1]))
      zrho = papp1[jk-1, jl-1]*zfac1
      zfac2 = 1.0 / (papp1[jk-1, jl-1] - yrmcst.retv*zfoeew[jl-1])
      zrodqsdp = -zrho*pqs[jk-1, jl-1]*zfac2
      zldcp = \
        zfwat[jl-1]*zlvdcp[jk-1, jl-1] + (1.0 - zfwat[jl-1])*zlsdcp[jk-1, jl-1]
      zfac3 = 1.0 / (1.0 + zldcp*zdqsdtemp[jl-1])
      dtdzmo = yrmcst.rg*(1.0 / yrmcst.rcpd - zldcp*zrodqsdp)*zfac3
      zdqsdz = zdqsdtemp[jl-1]*dtdzmo - yrmcst.rg*zrodqsdp
      zfac4 = 1.0 / zrho
      zdqc = min(zdqsdz*(pmfu[jk-1, jl-1] + pmfd[jk-1, jl-1])*ptsphy*zfac4, zqc[-1 + \
        jk, jl-1])
      zqc[jk-1, jl-1] = zqc[jk-1, jl-1] - zdqc
    
    # New cloud liquid/ice contents and condensation rates (liquid/ice)
    
    for jl in range(kidia, kfdia + 1):
      zqlwc[jk-1, jl-1] = zqc[jk-1, jl-1]*zfwat[jl-1]
      zqiwc[jk-1, jl-1] = zqc[jk-1, jl-1]*(1.0 - zfwat[jl-1])
      zcondl[jk-1, jl-1] = (zqlwc[jk-1, jl-1] - zl[jk-1, jl-1])*zqtmst
      zcondi[jk-1, jl-1] = (zqiwc[jk-1, jl-1] - zi[jk-1, jl-1])*zqtmst
    
    
    # Calculate precipitation overlap.
    # Simple form based on Maximum Overlap.
    
    for jl in range(kidia, kfdia + 1):
      if pclc[jk-1, jl-1] > zcovptot[jl-1]:
        # total rain frac
        zcovptot[jl-1] = pclc[jk-1, jl-1]
      zcovpclr[jl-1] = zcovptot[jl-1] - pclc[jk-1, jl-1]        # clear sky frac
      zcovpclr[jl-1] = max(zcovpclr[jl-1], 0.0)
    
    #*         3.3    CALCULATE PRECIPITATION
    
    # Melting of incoming snow
    
    for jl in range(kidia, kfdia + 1):
      if zsfl[jl-1] != 0.0:
        zcons = (zcons2*zdp[jk-1, jl-1]) / zlfdcp[jk-1, jl-1]
        zsnmlt = min(zsfl[jl-1], zcons*max(0.0, (ztp1[jk-1, jl-1] - zmeltp2)))
        zrfln[jl-1] = zrfl[jl-1] + zsnmlt
        zsfln[jl-1] = zsfl[jl-1] - zsnmlt
        ztp1[jk-1, jl-1] = ztp1[jk-1, jl-1] - zsnmlt / zcons
      else:
        zrfln[jl-1] = zrfl[jl-1]
        zsfln[jl-1] = zsfl[jl-1]
    
    for jl in range(kidia, kfdia + 1):
      
      #   Diagnostic calculation of rain production from cloud liquid water
      
      if pclc[jk-1, jl-1] > zeps2:
        # if yrphnc.levapls2 or ldrain1d:
        if False or ldrain1d:
          zlcrit = 1.9*yrecldp.rclcrit
        else:
          zlcrit = yrecldp.rclcrit*2.
        zcldl = zqlwc[jk-1, jl-1] / pclc[jk-1, jl-1]          # in-cloud liquid
        zd = zckcodtl*(1.0 - np.exp(-(zcldl / zlcrit)**2))
        zlnew = pclc[jk-1, jl-1]*zcldl*np.exp(-zd)
        zprr = zqlwc[jk-1, jl-1] - zlnew
        zqlwc[jk-1, jl-1] = zqlwc[jk-1, jl-1] - zprr
      else:
        zprr = 0.0
      
      #   Diagnostic calculation of snow production from cloud ice
      
      if pclc[jk-1, jl-1] > zeps2:
        # if yrphnc.levapls2 or ldrain1d:
        if False or ldrain1d:
          zlcrit = 1.E-04
        else:
          zlcrit = yrecldp.rclcrit*2.
        zcldi = zqiwc[jk-1, jl-1] / pclc[jk-1, jl-1]          # in-cloud ice
        zd = zckcodti*np.exp(0.025*(ztp1[jk-1, jl-1] - yrmcst.rtt))*(1.0 - np.exp(-(zcldi / \
          zlcrit)**2))
        zinew = pclc[jk-1, jl-1]*zcldi*np.exp(-zd)
        zprs = zqiwc[jk-1, jl-1] - zinew
        zqiwc[jk-1, jl-1] = zqiwc[jk-1, jl-1] - zprs
      else:
        zprs = 0.0
      
      #   New precipitation (rain + snow)
      
      zdr = zcons2*zdp[jk-1, jl-1]*(zprr + zprs)
      
      #   Rain fraction (different from cloud liquid water fraction!)
      
      if ztp1[jk-1, jl-1] < yrmcst.rtt:
        zrfreeze[jk-1, jl-1] = zcons2*zdp[jk-1, jl-1]*zprr
        zfwatr = 0.0
      else:
        zfwatr = 1.0
      
      zrn = zfwatr*zdr
      zsn = (1.0 - zfwatr)*zdr
      zrfln[jl-1] = zrfln[jl-1] + zrn
      zsfln[jl-1] = zsfln[jl-1] + zsn
      
      #   Precip evaporation
      
      zprtot = zrfln[jl-1] + zsfln[jl-1]
      llo2 = zprtot > zeps2 and zcovpclr[jl-1] > zeps2 and ldrain1d
      if llo2:
        
        zpreclr = (zprtot*zcovpclr[jl-1]) / zcovptot[jl-1]
        
        #     This is the humidity in the moistest zcovpclr region
        
        zqe = pqs[jk-1, jl-1] - ((pqs[jk-1, jl-1] - zqlim[jl-1])*zcovpclr[-1 + \
          jl]) / (1.0 - pclc[jk-1, jl-1])**2
        zbeta = yrmcst.rg*yrecldp.rpecons*(((np.sqrt(papp1[jk-1, jl-1] / paphp1[klev, jl-1]) / \
          5.09E-3)*zpreclr) / zcovpclr[jl-1])**0.5777
        
        #     implicit solution:
        zb = (ptsphy*zbeta*(pqs[jk-1, jl-1] - zqe)) / (1.0 + zbeta*ptsphy*zcorqs[jl-1])
        
        #     exact solution:
        #     ZB=(PQS(JL,JK)-ZQE)*(_ONE_-EXP(-ZBETA*ZCORQS(JL)*PTSPHY))/ZCORQS(JL)
        
        zdtgdp[jl-1] = (ptsphy*yrmcst.rg) / (paphp1[jk, jl-1] - paphp1[jk-1, jl-1])
        
        zdpr = (zcovpclr[jl-1]*zb) / zdtgdp[jl-1]
        zdpr = min(zdpr, zpreclr)
        zpreclr = zpreclr - zdpr          # take away from clr sky flux
        if zpreclr <= 0.0:
          zcovptot[jl-1] = pclc[jk-1, jl-1]
        #reset
        pcovptot[jk-1, jl-1] = zcovptot[jl-1]
        
        # warm proportion
        zevapr[jk-1, jl-1] = (zdpr*zrfln[jl-1]) / zprtot
        zrfln[jl-1] = zrfln[jl-1] - zevapr[jk-1, jl-1]
        
        # ice proportion
        zevaps[jk-1, jl-1] = (zdpr*zsfln[jl-1]) / zprtot
        zsfln[jl-1] = zsfln[jl-1] - zevaps[jk-1, jl-1]
      
    
    # Update of T and Q tendencies due to:
    #  - condensation/evaporation of cloud liquid water/ice
    #  - detrainment of convective cloud condensate
    #  - evaporation of precipitation
    #  - freezing of rain (impact on T only).
    
    for jl in range(kidia, kfdia + 1):
      zdqdt[jl-1] = -(zcondl[jk-1, jl-1] + zcondi[jk-1, jl-1]) + (plude[jk-1, \
        jl-1] + zevapr[jk-1, jl-1] + zevaps[jk-1, jl-1])*zgdp[jl-1]
      
      zdtdt[jl-1] = zlvdcp[jk-1, jl-1]*zcondl[jk-1, jl-1] + zlsdcp[jk-1, -1 + \
        jl]*zcondi[jk-1, jl-1] - (zlvdcp[jk-1, jl-1]*zevapr[jk-1, jl-1] + \
        zlsdcp[jk-1, jl-1]*zevaps[jk-1, jl-1] + plude[jk-1, jl-1]*(zfwat[-1 + \
        jl]*zlvdcp[jk-1, jl-1] + (1.0 - zfwat[jl-1])*zlsdcp[jk-1, jl-1]) - \
        (zlsdcp[jk-1, jl-1] - zlvdcp[jk-1, jl-1])*zrfreeze[jk-1, jl-1])*zgdp[-1 \
        + jl]
      
      # first guess T and Q
      ztp1[jk-1, jl-1] = ztp1[jk-1, jl-1] + ptsphy*zdtdt[jl-1]
      zqp1[jk-1, jl-1] = zqp1[jk-1, jl-1] + ptsphy*zdqdt[jl-1]
      
      zpp[jl-1] = papp1[jk-1, jl-1]
      zqold[jl-1] = zqp1[jk-1, jl-1]
    
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
    for jl in range(kidia, kfdia + 1):
      if ztp1[jk-1, jl-1] > yrmcst.rtt:
        z3es = yrethf.r3les
        z4es = yrethf.r4les
        z5alcp = yrethf.r5alvcp
        zaldcp = yrethf.ralvdcp
      else:
        z3es = yrethf.r3ies
        z4es = yrethf.r4ies
        z5alcp = yrethf.r5alscp
        zaldcp = yrethf.ralsdcp
      
      zqp = 1.0 / zpp[jl-1]
      ztarg = ztp1[jk-1, jl-1]
      zfoeew[jl-1] = yrethf.r2es*np.exp((z3es*(ztarg - yrmcst.rtt)) / (ztarg - z4es))
      zqsat[jl-1] = zqp*zfoeew[jl-1]
      if zqsat[jl-1] > zqmax:
        zqsat[jl-1] = zqmax
      zcor = 1.0 / (1.0 - yrmcst.retv*zqsat[jl-1])
      zqsat[jl-1] = zqsat[jl-1]*zcor
      z2s = z5alcp / (ztarg - z4es)**2
      zcond1 = (zqp1[jk-1, jl-1] - zqsat[jl-1]) / (1.0 + zqsat[jl-1]*zcor*z2s)
      ztp1[jk-1, jl-1] = ztp1[jk-1, jl-1] + zaldcp*zcond1
      zqp1[jk-1, jl-1] = zqp1[jk-1, jl-1] - zcond1
      ztarg = ztp1[jk-1, jl-1]
      zfoeew[jl-1] = yrethf.r2es*np.exp((z3es*(ztarg - yrmcst.rtt)) / (ztarg - z4es))
      zqsat[jl-1] = zqp*zfoeew[jl-1]
      if zqsat[jl-1] > zqmax:
        zqsat[jl-1] = zqmax
      zcor = 1.0 / (1.0 - yrmcst.retv*zqsat[jl-1])
      zqsat[jl-1] = zqsat[jl-1]*zcor
      z2s = z5alcp / (ztarg - z4es)**2
      zcond1 = (zqp1[jk-1, jl-1] - zqsat[jl-1]) / (1.0 + zqsat[jl-1]*zcor*z2s)
      ztp1[jk-1, jl-1] = ztp1[jk-1, jl-1] + zaldcp*zcond1
      zqp1[jk-1, jl-1] = zqp1[jk-1, jl-1] - zcond1
    # -----------------------------------
    
    for jl in range(kidia, kfdia + 1):
      zdq[jl-1] = max(0.0, zqold[jl-1] - zqp1[jk-1, jl-1])
      zdr2 = zcons2*zdp[jk-1, jl-1]*zdq[jl-1]
      # Update rain fraction and freezing.
      # Note: impact of new temperature ZTP1 on ZFWAT is neglected here.
      if ztp1[jk-1, jl-1] < yrmcst.rtt:
        zrfreeze2 = zfwat[jl-1]*zdr2
        zfwatr = 0.0
      else:
        zrfreeze2 = 0.0
        zfwatr = 1.0
      zrn = zfwatr*zdr2
      zsn = (1.0 - zfwatr)*zdr2
      # Note: The extra condensation due to the adjustment goes directly to precipitation
      zcondl[jk-1, jl-1] = zcondl[jk-1, jl-1] + zfwatr*zdq[jl-1]*zqtmst
      zcondi[jk-1, jl-1] = zcondi[jk-1, jl-1] + (1.0 - zfwatr)*zdq[jl-1]*zqtmst
      zrfln[jl-1] = zrfln[jl-1] + zrn
      zsfln[jl-1] = zsfln[jl-1] + zsn
      zrfreeze[jk-1, jl-1] = zrfreeze[jk-1, jl-1] + zrfreeze2
    
    for jl in range(kidia, kfdia + 1):
      zdqdt[jl-1] = -(zcondl[jk-1, jl-1] + zcondi[jk-1, jl-1]) + (plude[jk-1, \
        jl-1] + zevapr[jk-1, jl-1] + zevaps[jk-1, jl-1])*zgdp[jl-1]
      
      zdtdt[jl-1] = zlvdcp[jk-1, jl-1]*zcondl[jk-1, jl-1] + zlsdcp[jk-1, -1 + \
        jl]*zcondi[jk-1, jl-1] - (zlvdcp[jk-1, jl-1]*zevapr[jk-1, jl-1] + \
        zlsdcp[jk-1, jl-1]*zevaps[jk-1, jl-1] + plude[jk-1, jl-1]*(zfwat[-1 + \
        jl]*zlvdcp[jk-1, jl-1] + (1.0 - zfwat[jl-1])*zlsdcp[jk-1, jl-1]) - \
        (zlsdcp[jk-1, jl-1] - zlvdcp[jk-1, jl-1])*zrfreeze[jk-1, jl-1])*zgdp[-1 \
        + jl]
      
      zdldt[jl-1] = (zqlwc[jk-1, jl-1] - zl[jk-1, jl-1])*zqtmst
      
      zdidt[jl-1] = (zqiwc[jk-1, jl-1] - zi[jk-1, jl-1])*zqtmst
      
      ptenq[jk-1, jl-1] = zdqdt[jl-1]
      ptent[jk-1, jl-1] = zdtdt[jl-1]
      ptenl[jk-1, jl-1] = zdldt[jl-1]
      pteni[jk-1, jl-1] = zdidt[jl-1]
      
      pfplsl[jk, jl-1] = zrfln[jl-1]
      pfplsn[jk, jl-1] = zsfln[jl-1]
    
    # record rain flux for next level
    
    for jl in range(kidia, kfdia + 1):
      zrfl[jl-1] = zrfln[jl-1]
      zsfl[jl-1] = zsfln[jl-1]
    
  #jk
  
  #*     ENTHALPY FLUXES DUE TO PRECIPITATION
  #      ------------------------------------
  
  for jk in range(1, klev + 1 + 1):
    for jl in range(kidia, kfdia + 1):
      pfhpsl[jk-1, jl-1] = -pfplsl[jk-1, jl-1]*yrmcst.rlvtt
      pfhpsn[jk-1, jl-1] = -pfplsn[jk-1, jl-1]*yrmcst.rlstt
  
  #     ------------------------------------------------------------------
  
  #IF (LHOOK) CALL DR_HOOK('CLOUDSC2',1,ZHOOK_HANDLE)
  return 
