import numpy as np
from pathlib import Path
from collections import OrderedDict
from itertools import product

import gt4py
import gt4py.gtscript as gtscript
import gt4py.storage as gt_storage

from cloudsc2_inputs import load_input_parameters


backend = "numpy"
dtype = np.float64
origin = (0, 0, 0)

# Get the constant parameters from file
rootpath = Path(__file__).resolve().parents[2]
input_path = rootpath/'config-files/input.h5'
yrecldp, yrmcst, yrethf, yrephli, yrecld = load_input_parameters(path=input_path)


def wrap_input_arrays(satur_args, cloudsc_args):
    wrapped_fields = OrderedDict()

    def wrap(data):
        # shape = data.shape
        # ijk_data = np.reshape(data, (10, 10, shape[0]))
        # return gt_storage.from_array(ijk_data, backend, default_origin=origin, dtype=dtype)

        shape = data.shape
        ijk_data = np.reshape(data, (1, shape[0], shape[1]))
        return gt_storage.from_array(ijk_data, backend, default_origin=origin, dtype=dtype)

        # return gt_storage.from_array(data, backend, default_origin=origin, dtype=dtype)

    cloudsc_args['paphp1'] = wrap(cloudsc_args['paphp1'])
    cloudsc_args['papp1'] = wrap(cloudsc_args['papp1'])
    cloudsc_args['pqm1'] = wrap(cloudsc_args['pqm1'])
    cloudsc_args['pqs'] = wrap(cloudsc_args['pqs'])
    cloudsc_args['ptm1'] = wrap(cloudsc_args['ptm1'])
    cloudsc_args['pl'] = wrap(cloudsc_args['pl'])
    cloudsc_args['pi'] = wrap(cloudsc_args['pi'])
    cloudsc_args['plude'] = wrap(cloudsc_args['plude'])
    cloudsc_args['plu'] = wrap(cloudsc_args['plu'])
    cloudsc_args['pmfu'] = wrap(cloudsc_args['pmfu'])
    cloudsc_args['pmfd'] = wrap(cloudsc_args['pmfd'])
    cloudsc_args['ptent'] = wrap(cloudsc_args['ptent'])
    cloudsc_args['pgtent'] = wrap(cloudsc_args['pgtent'])
    cloudsc_args['ptenq'] = wrap(cloudsc_args['ptenq'])
    cloudsc_args['pgtenq'] = wrap(cloudsc_args['pgtenq'])
    cloudsc_args['ptenl'] = wrap(cloudsc_args['ptenl'])
    cloudsc_args['pgtenl'] = wrap(cloudsc_args['pgtenl'])
    cloudsc_args['pteni'] = wrap(cloudsc_args['pteni'])
    cloudsc_args['pgteni'] = wrap(cloudsc_args['pgteni'])
    cloudsc_args['psupsat'] = wrap(cloudsc_args['psupsat'])
    cloudsc_args['pclc'] = wrap(cloudsc_args['pclc'])
    cloudsc_args['pfplsl'] = wrap(cloudsc_args['pfplsl'])
    cloudsc_args['pfplsn'] = wrap(cloudsc_args['pfplsn'])
    cloudsc_args['pfhpsl'] = wrap(cloudsc_args['pfhpsl'])
    cloudsc_args['pfhpsn'] = wrap(cloudsc_args['pfhpsn'])
    cloudsc_args['pcovptot'] = wrap(cloudsc_args['pcovptot'])

    satur_args['paprsf'] = cloudsc_args['papp1']
    satur_args['pt'] = cloudsc_args['ptm1']
    satur_args['pqsat'] = cloudsc_args['pqs']

    return satur_args, cloudsc_args


def foealfa(ptare, yrethf):
  return min(1.0,((max(yrethf.rtice,min(yrethf.rtwat,ptare))-yrethf.rtice)*yrethf.rtwat_rtice_r)**2)


def foeewm(ptare, yrethf, yrmcst):
  return yrethf.r2es*(foealfa(ptare, yrethf)*np.exp(yrethf.r3les*(ptare-yrmcst.rtt)/(ptare-yrethf.r4les)) + \
                      (1.0-foealfa(ptare, yrethf))*np.exp(yrethf.r3ies*(ptare-yrmcst.rtt)/(ptare-yrethf.r4ies)))


def satur_py_gt4py(kidia, kfdia, klon, ktdia, klev, ldphylin, paprsf, pt, pqsat, kflag):

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
                ztarg = pt[0,jk,jl]
                zalfa = foealfa(ztarg, yrethf)

                zfoeewl = yrethf.r2es*np.exp(yrethf.r3les*(ztarg-yrmcst.rtt)/(ztarg-yrethf.r4les))
                zfoeewi = yrethf.r2es*np.exp(yrethf.r3ies*(ztarg-yrmcst.rtt)/(ztarg-yrethf.r4ies))
                zfoeew = zalfa*zfoeewl+(1.0-zalfa)*zfoeewi

                zqs    = zfoeew/paprsf[0,jk,jl]
                if zqs > zqmax:
                    zqs=zqmax

                zcor = 1.0/(1.0-yrmcst.retv*zqs)
                pqsat[0,jk,jl]=zqs*zcor

    else:
        for jk in range(klev):
            for jl in range(kfdia):
                if(kflag == 1):
                    zew  = foeewmcu(pt[0,jk,jl])
                else:
                    zew  = foeewm(pt[0,jk,jl])

                zqs  = zew/paprsf[0,jk,jl]
                zqs  = min(zqmax,zqs)
                zcor = 1.0/(1.0-yrmcst.retv*zqs)
                pqsat[0,jk,jl]=zqs*zcor


                
def cloudsc2_py_gt4py(
        kidia, kfdia, klon, ktdia, klev, ldrain1d, ptsphy, paphp1,
        papp1, pqm1, pqs, ptm1, pl, pi, plude, plu, pmfu, pmfd, ptent,
        pgtent, ptenq, pgtenq, ptenl, pgtenl, pteni, pgteni, psupsat,
        pclc, pfplsl, pfplsn, pfhpsl, pfhpsn, pcovptot):
  
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
  
  zrfl = np.ndarray(order="C", shape=(1,1,klon))
  zsfl = np.ndarray(order="C", shape=(1,1,klon))
  zrfln = np.ndarray(order="C", shape=(1,1,klon))
  zsfln = np.ndarray(order="C", shape=(1,1,klon))
  zgdp = np.ndarray(order="C", shape=(1,1,klon))
  zdqdt = np.ndarray(order="C", shape=(1,1,klon))
  zdtdt = np.ndarray(order="C", shape=(1,1,klon))
  zdldt = np.ndarray(order="C", shape=(1,1,klon))
  zdidt = np.ndarray(order="C", shape=(1,1,klon))
  zqcrit = np.ndarray(order="C", shape=(1,1,klon))
  zcovpclr = np.ndarray(order="C", shape=(1,1,klon))
  zcovptot = np.ndarray(order="C", shape=(1,1,klon))
  zdqsdtemp = np.ndarray(order="C", shape=(1,1,klon))
  zcorqs = np.ndarray(order="C", shape=(1,1,klon))
  zdtgdp = np.ndarray(order="C", shape=(1,1,klon))
  zqold = np.ndarray(order="C", shape=(1,1,klon))
  zpp = np.ndarray(order="C", shape=(1,1,klon))
  zdq = np.ndarray(order="C", shape=(1,1,klon))
  zqlim = np.ndarray(order="C", shape=(1,1,klon))
  zscalm = np.ndarray(order="C", shape=(klev,))
  zqsat = np.ndarray(order="C", shape=(1,1,klon))
  zfoeew = np.ndarray(order="C", shape=(1,1,klon))
  zfwat = np.ndarray(order="C", shape=(1,1,klon))
  
  
  ztrpaus = np.ndarray(order="C", shape=(1,1,klon))
  
  
  llflag = np.ndarray(order="C", shape=(1,1,klon))
  #REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
  
  #     ------------------------------------------------------------------
  ztp1 = np.ndarray(order="F", shape=(1,klev,klon))
  zqp1 = np.ndarray(order="F", shape=(1,klev,klon))
  zl = np.ndarray(order="F", shape=(1,klev,klon))
  zi = np.ndarray(order="F", shape=(1,klev,klon))
  zlude = np.ndarray(order="F", shape=(1,klev,klon))
  zqc = np.ndarray(order="F", shape=(1,klev,klon))
  zqlwc = np.ndarray(order="F", shape=(1,klev,klon))
  zqiwc = np.ndarray(order="F", shape=(1,klev,klon))
  zdp = np.ndarray(order="F", shape=(1,klev,klon))
  zlsdcp = np.ndarray(order="F", shape=(1,klev,klon))
  zlfdcp = np.ndarray(order="F", shape=(1,klev,klon))
  zlvdcp = np.ndarray(order="F", shape=(1,klev,klon))
  zrfreeze = np.ndarray(order="F", shape=(1,klev,klon))
  zcondl = np.ndarray(order="F", shape=(1,klev,klon))
  zcondi = np.ndarray(order="F", shape=(1,klev,klon))
  zevapr = np.ndarray(order="F", shape=(1,klev,klon))
  zevaps = np.ndarray(order="F", shape=(1,klev,klon))
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
      ztp1[0,jk,jl] = ptm1[0,jk,jl] + ptsphy*pgtent[0,jk,jl]
      zqp1[0,jk,jl] = \
        pqm1[0,jk,jl] + ptsphy*pgtenq[0,jk,jl] + psupsat[0,jk,jl]
      zl[0,jk,jl] = pl[0,jk,jl] + ptsphy*pgtenl[0,jk,jl]
      zi[0,jk,jl] = pi[0,jk,jl] + ptsphy*pgteni[0,jk,jl]
  
  for jk in range(klev):
    
    # Parameter for cloud formation

    zscalm[jk] = zscal*max((yrecld.ceta[jk] - 0.2), zeps1)**0.2
    
    for jl in range(kfdia):
      
      # thermodynamic constants
      zdp[0,jk,jl] = paphp1[0,jk+1,jl] - paphp1[0,jk,jl]
      zzz = 1.0 / (yrmcst.rcpd + yrmcst.rcpd*yrethf.rvtmp2*zqp1[0,jk,jl])
      zlfdcp[0,jk,jl] = yrmcst.rlmlt*zzz
      zlsdcp[0,jk,jl] = yrmcst.rlstt*zzz
      zlvdcp[0,jk,jl] = yrmcst.rlvtt*zzz
      llflag[0,0,jl] = True
  
  #     ------------------------------------------------------------------
  
  #*         2.2    INITIALIZATION OF CLOUD AND PRECIPITATION ARRAYS
  #                 ------------------------------------------------
  
  #       Clear cloud and freezing arrays
  
  for jk in range(klev):
    for jl in range(kfdia):
      pclc[0,jk,jl] = 0.0
      zqc[0,jk,jl] = 0.0
      zqlwc[0,jk,jl] = 0.0
      zqiwc[0,jk,jl] = 0.0
      zrfreeze[0,jk,jl] = 0.0
      zcondl[0,jk,jl] = 0.0
      zcondi[0,jk,jl] = 0.0
      zevapr[0,jk,jl] = 0.0
      zevaps[0,jk,jl] = 0.0
      pcovptot[0,jk,jl] = 0.0
  
  #       Set to zero precipitation fluxes at the top
  
  for jl in range(kfdia):
    zrfl[0,0,jl] = 0.0
    zsfl[0,0,jl] = 0.0
    pfplsl[0,0,jl] = 0.0
    pfplsn[0,0,jl] = 0.0
    zcovptot[0,0,jl] = 0.0
    zcovpclr[0,0,jl] = 0.0
  
  # Eta value at tropopause
  for jl in range(kfdia):
    ztrpaus[0,0,jl] = 0.1
  for jk in range(klev - 1):
    for jl in range(kfdia):
      llo1 = \
        yrecld.ceta[jk] > 0.1 and yrecld.ceta[jk] < 0.4 and ztp1[0,jk,jl] > ztp1[0,jk+1,jl]
      if llo1:
        ztrpaus[0,0,jl] = yrecld.ceta[jk]
  
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
        zoealfaw = 0.545*(np.tanh(0.17*(ztp1[0,jk,jl] - yrephli.rlptrc)) + 1.0)
        if ztp1[0,jk,jl] < yrmcst.rtt:
          zfwat[0,0,jl] = zoealfaw
          z3es = yrethf.r3ies
          z4es = yrethf.r4ies
        else:
          zfwat[0,0,jl] = 1.0
          z3es = yrethf.r3les
          z4es = yrethf.r4les
        zfoeew[0,0,jl] = \
          yrethf.r2es*np.exp((z3es*(ztp1[0,jk,jl] - yrmcst.rtt)) / (ztp1[0,jk,jl] - z4es))
        zesdp = zfoeew[0,0,jl] / papp1[0,jk,jl]
        if zesdp > zqmax:
          zesdp = zqmax
      else:
        zfwat[0,0,jl] = foealfa(ztp1[0,jk,jl], yrethf)
        zfoeew[0,0,jl] = foeewm(ztp1[0,jk,jl], yrethf, yrmcst)
        zesdp = zfoeew[0,0,jl] / papp1[0,jk,jl]
      zfacw = yrethf.r5les / ((ztp1[0,jk,jl] - yrethf.r4les)**2)
      zfaci = yrethf.r5ies / ((ztp1[0,jk,jl] - yrethf.r4ies)**2)
      zfac = zfwat[0,0,jl]*zfacw + (1.0 - zfwat[0,0,jl])*zfaci
      zcor = 1.0 / (1.0 - yrmcst.retv*zesdp)
      zdqsdtemp[0,0,jl] = zfac*zcor*pqs[0,jk,jl]
      zcorqs[0,0,jl] = 1.0 + zcons3*zdqsdtemp[0,0,jl]
      
      # use clipped state
      
      zqlim[0,0,jl] = zqp1[0,jk,jl]
      if zqp1[0,jk,jl] > pqs[0,jk,jl]:
        zqlim[0,0,jl] = pqs[0,jk,jl]
      
      # set up critical value of humidity
      
      zeta3 = ztrpaus[0,0,jl]
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
      if ztp1[0,jk,jl] < yrethf.rtice:
        zsupsat = 1.8 - 3.E-03*ztp1[0,jk,jl]
      else:
        zsupsat = 1.0
      zqsat[0,0,jl] = pqs[0,jk,jl]*zsupsat
      zqcrit[0,0,jl] = zcrh2*zqsat[0,0,jl]
    
    # Simple UNIFORM distribution of total water from Letreut & Li (90)
    
    for jl in range(kfdia):
      zqt = zqp1[0,jk,jl] + zl[0,jk,jl] + zi[0,jk,jl]
      if zqt <= zqcrit[0,0,jl]:
        pclc[0,jk,jl] = 0.0
        zqc[0,jk,jl] = 0.0
      elif zqt >= zqsat[0,0,jl]:
        pclc[0,jk,jl] = 1.0
        zqc[0,jk,jl] = (1.0 - zscalm[jk])*(zqsat[0,0,jl] - zqcrit[0,0,jl])
      else:
        zqpd = zqsat[0,0,jl] - zqt
        zqcd = zqsat[0,0,jl] - zqcrit[0,0,jl]
        pclc[0,jk,jl] = \
          1.0 - np.sqrt(zqpd / (zqcd - zscalm[jk]*(zqt - zqcrit[0,0,jl])))
        zqc[0,jk,jl] = \
          (zscalm[jk]*zqpd + (1.0 - zscalm[jk])*zqcd)*pclc[0,jk,jl]**2
    
    # Add convective component
    
    for jl in range(kfdia):
      zgdp[0,0,jl] = yrmcst.rg / (paphp1[0,jk+1,jl] - paphp1[0,jk,jl])
      zlude[0,jk,jl] = plude[0,jk,jl]*ptsphy*zgdp[0,0,jl]
      if jk < klev:
        llo1 = zlude[0,jk,jl] >= yrecldp.rlmin and plu[0,jk+1,jl] >= zeps2
      else:
        llo1 = False
      if llo1:
        pclc[0,jk,jl] = pclc[0,jk,jl] + (1.0 - pclc[0,jk,jl])*(1.0 - \
          np.exp(-zlude[0,jk,jl] / plu[0,jk+1,jl]))
        zqc[0,jk,jl] = zqc[0,jk,jl] + zlude[0,jk,jl]
    
    # Add compensating subsidence component
    
    for jl in range(kfdia):
      zfac1 = 1.0 / ((yrmcst.rd*ztp1[0,jk,jl]))
      zrho = papp1[0,jk,jl]*zfac1
      zfac2 = 1.0 / (papp1[0,jk,jl] - yrmcst.retv*zfoeew[0,0,jl])
      zrodqsdp = -zrho*pqs[0,jk,jl]*zfac2
      zldcp = \
        zfwat[0,0,jl]*zlvdcp[0,jk,jl] + (1.0 - zfwat[0,0,jl])*zlsdcp[0,jk,jl]
      zfac3 = 1.0 / (1.0 + zldcp*zdqsdtemp[0,0,jl])
      dtdzmo = yrmcst.rg*(1.0 / yrmcst.rcpd - zldcp*zrodqsdp)*zfac3
      zdqsdz = zdqsdtemp[0,0,jl]*dtdzmo - yrmcst.rg*zrodqsdp
      zfac4 = 1.0 / zrho
      zdqc = min(zdqsdz*(pmfu[0,jk,jl] + pmfd[0,jk,jl])*ptsphy*zfac4, zqc[0,jk,jl])
      zqc[0,jk,jl] = zqc[0,jk,jl] - zdqc
    
    # New cloud liquid/ice contents and condensation rates (liquid/ice)
    
    for jl in range(kfdia):
      zqlwc[0,jk,jl] = zqc[0,jk,jl]*zfwat[0,0,jl]
      zqiwc[0,jk,jl] = zqc[0,jk,jl]*(1.0 - zfwat[0,0,jl])
      zcondl[0,jk,jl] = (zqlwc[0,jk,jl] - zl[0,jk,jl])*zqtmst
      zcondi[0,jk,jl] = (zqiwc[0,jk,jl] - zi[0,jk,jl])*zqtmst
    
    
    # Calculate precipitation overlap.
    # Simple form based on Maximum Overlap.
    
    for jl in range(kfdia):
      if pclc[0,jk,jl] > zcovptot[0,0,jl]:
        # total rain frac
        zcovptot[0,0,jl] = pclc[0,jk,jl]
      zcovpclr[0,0,jl] = zcovptot[0,0,jl] - pclc[0,jk,jl]        # clear sky frac
      zcovpclr[0,0,jl] = max(zcovpclr[0,0,jl], 0.0)
    
    #*         3.3    CALCULATE PRECIPITATION
    
    # Melting of incoming snow
    
    for jl in range(kfdia):
      if zsfl[0,0,jl] != 0.0:
        zcons = (zcons2*zdp[0,jk,jl]) / zlfdcp[0,jk,jl]
        zsnmlt = min(zsfl[0,0,jl], zcons*max(0.0, (ztp1[0,jk,jl] - zmeltp2)))
        zrfln[0,0,jl] = zrfl[0,0,jl] + zsnmlt
        zsfln[0,0,jl] = zsfl[0,0,jl] - zsnmlt
        ztp1[0,jk,jl] = ztp1[0,jk,jl] - zsnmlt / zcons
      else:
        zrfln[0,0,jl] = zrfl[0,0,jl]
        zsfln[0,0,jl] = zsfl[0,0,jl]
    
    for jl in range(kfdia):
      
      #   Diagnostic calculation of rain production from cloud liquid water
      
      if pclc[0,jk,jl] > zeps2:
        # if yrphnc.levapls2 or ldrain1d:
        if False or ldrain1d:
          zlcrit = 1.9*yrecldp.rclcrit
        else:
          zlcrit = yrecldp.rclcrit*2.
        zcldl = zqlwc[0,jk,jl] / pclc[0,jk,jl]          # in-cloud liquid
        zd = zckcodtl*(1.0 - np.exp(-(zcldl / zlcrit)**2))
        zlnew = pclc[0,jk,jl]*zcldl*np.exp(-zd)
        zprr = zqlwc[0,jk,jl] - zlnew
        zqlwc[0,jk,jl] = zqlwc[0,jk,jl] - zprr
      else:
        zprr = 0.0
      
      #   Diagnostic calculation of snow production from cloud ice
      
      if pclc[0,jk,jl] > zeps2:
        # if yrphnc.levapls2 or ldrain1d:
        if False or ldrain1d:
          zlcrit = 1.E-04
        else:
          zlcrit = yrecldp.rclcrit*2.
        zcldi = zqiwc[0,jk,jl] / pclc[0,jk,jl]          # in-cloud ice
        zd = zckcodti*np.exp(0.025*(ztp1[0,jk,jl] - yrmcst.rtt))*(1.0 - np.exp(-(zcldi / \
          zlcrit)**2))
        zinew = pclc[0,jk,jl]*zcldi*np.exp(-zd)
        zprs = zqiwc[0,jk,jl] - zinew
        zqiwc[0,jk,jl] = zqiwc[0,jk,jl] - zprs
      else:
        zprs = 0.0
      
      #   New precipitation (rain + snow)
      
      zdr = zcons2*zdp[0,jk,jl]*(zprr + zprs)
      
      #   Rain fraction (different from cloud liquid water fraction!)
      
      if ztp1[0,jk,jl] < yrmcst.rtt:
        zrfreeze[0,jk,jl] = zcons2*zdp[0,jk,jl]*zprr
        zfwatr = 0.0
      else:
        zfwatr = 1.0
      
      zrn = zfwatr*zdr
      zsn = (1.0 - zfwatr)*zdr
      zrfln[0,0,jl] = zrfln[0,0,jl] + zrn
      zsfln[0,0,jl] = zsfln[0,0,jl] + zsn
      
      #   Precip evaporation
      
      zprtot = zrfln[0,0,jl] + zsfln[0,0,jl]
      llo2 = zprtot > zeps2 and zcovpclr[0,0,jl] > zeps2 and ldrain1d
      if llo2:
        
        zpreclr = (zprtot*zcovpclr[0,0,jl]) / zcovptot[0,0,jl]
        
        #     This is the humidity in the moistest zcovpclr region
        
        zqe = pqs[0,jk,jl] - ((pqs[0,jk,jl] - zqlim[0,0,jl])*zcovpclr[0,0,jl]) \
          / (1.0 - pclc[0,jk,jl])**2
        zbeta = yrmcst.rg*yrecldp.rpecons*(((np.sqrt(papp1[0,jk,jl] / paphp1[klev, jl]) / \
          5.09E-3)*zpreclr) / zcovpclr[0,0,jl])**0.5777
        
        #     implicit solution:
        zb = (ptsphy*zbeta*(pqs[0,jk,jl] - zqe)) / (1.0 + zbeta*ptsphy*zcorqs[0,0,jl])
        
        #     exact solution:
        #     ZB=(PQS(JL,JK)-ZQE)*(_ONE_-EXP(-ZBETA*ZCORQS(JL)*PTSPHY))/ZCORQS(JL)
        
        zdtgdp[0,0,jl] = (ptsphy*yrmcst.rg) / (paphp1[0,jk+1,jl] - paphp1[0,jk,jl])
        
        zdpr = (zcovpclr[0,0,jl]*zb) / zdtgdp[0,0,jl]
        zdpr = min(zdpr, zpreclr)
        zpreclr = zpreclr - zdpr          # take away from clr sky flux
        if zpreclr <= 0.0:
          zcovptot[0,0,jl] = pclc[0,jk,jl]
        #reset
        pcovptot[0,jk,jl] = zcovptot[0,0,jl]
        
        # warm proportion
        zevapr[0,jk,jl] = (zdpr*zrfln[0,0,jl]) / zprtot
        zrfln[0,0,jl] = zrfln[0,0,jl] - zevapr[0,jk,jl]
        
        # ice proportion
        zevaps[0,jk,jl] = (zdpr*zsfln[0,0,jl]) / zprtot
        zsfln[0,0,jl] = zsfln[0,0,jl] - zevaps[0,jk,jl]
      
    
    # Update of T and Q tendencies due to:
    #  - condensation/evaporation of cloud liquid water/ice
    #  - detrainment of convective cloud condensate
    #  - evaporation of precipitation
    #  - freezing of rain (impact on T only).
    
    for jl in range(kfdia):
      zdqdt[0,0,jl] = -(zcondl[0,jk,jl] + zcondi[0,jk,jl]) + (plude[0,jk,jl] \
        + zevapr[0,jk,jl] + zevaps[0,jk,jl])*zgdp[0,0,jl]
      
      zdtdt[0,0,jl] = zlvdcp[0,jk,jl]*zcondl[0,jk,jl] + zlsdcp[0,jk,jl] \
        *zcondi[0,jk,jl] - (zlvdcp[0,jk,jl]*zevapr[0,jk,jl] + \
        zlsdcp[0,jk,jl]*zevaps[0,jk,jl] + plude[0,jk,jl]*(zfwat[0,0,jl] \
        *zlvdcp[0,jk,jl] + (1.0 - zfwat[0,0,jl])*zlsdcp[0,jk,jl]) - \
        (zlsdcp[0,jk,jl] - zlvdcp[0,jk,jl])*zrfreeze[0,jk,jl])*zgdp[0,0,jl]
      
      # first guess T and Q
      ztp1[0,jk,jl] = ztp1[0,jk,jl] + ptsphy*zdtdt[0,0,jl]
      zqp1[0,jk,jl] = zqp1[0,jk,jl] + ptsphy*zdqdt[0,0,jl]
      
      zpp[0,0,jl] = papp1[0,jk,jl]
      zqold[0,0,jl] = zqp1[0,jk,jl]
    
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
      if ztp1[0,jk,jl] > yrmcst.rtt:
        z3es = yrethf.r3les
        z4es = yrethf.r4les
        z5alcp = yrethf.r5alvcp
        zaldcp = yrethf.ralvdcp
      else:
        z3es = yrethf.r3ies
        z4es = yrethf.r4ies
        z5alcp = yrethf.r5alscp
        zaldcp = yrethf.ralsdcp
      
      zqp = 1.0 / zpp[0,0,jl]
      ztarg = ztp1[0,jk,jl]
      zfoeew[0,0,jl] = yrethf.r2es*np.exp((z3es*(ztarg - yrmcst.rtt)) / (ztarg - z4es))
      zqsat[0,0,jl] = zqp*zfoeew[0,0,jl]
      if zqsat[0,0,jl] > zqmax:
        zqsat[0,0,jl] = zqmax
      zcor = 1.0 / (1.0 - yrmcst.retv*zqsat[0,0,jl])
      zqsat[0,0,jl] = zqsat[0,0,jl]*zcor
      z2s = z5alcp / (ztarg - z4es)**2
      zcond1 = (zqp1[0,jk,jl] - zqsat[0,0,jl]) / (1.0 + zqsat[0,0,jl]*zcor*z2s)
      ztp1[0,jk,jl] = ztp1[0,jk,jl] + zaldcp*zcond1
      zqp1[0,jk,jl] = zqp1[0,jk,jl] - zcond1
      ztarg = ztp1[0,jk,jl]
      zfoeew[0,0,jl] = yrethf.r2es*np.exp((z3es*(ztarg - yrmcst.rtt)) / (ztarg - z4es))
      zqsat[0,0,jl] = zqp*zfoeew[0,0,jl]
      if zqsat[0,0,jl] > zqmax:
        zqsat[0,0,jl] = zqmax
      zcor = 1.0 / (1.0 - yrmcst.retv*zqsat[0,0,jl])
      zqsat[0,0,jl] = zqsat[0,0,jl]*zcor
      z2s = z5alcp / (ztarg - z4es)**2
      zcond1 = (zqp1[0,jk,jl] - zqsat[0,0,jl]) / (1.0 + zqsat[0,0,jl]*zcor*z2s)
      ztp1[0,jk,jl] = ztp1[0,jk,jl] + zaldcp*zcond1
      zqp1[0,jk,jl] = zqp1[0,jk,jl] - zcond1
    # -----------------------------------
    
    for jl in range(kfdia):
      zdq[0,0,jl] = max(0.0, zqold[0,0,jl] - zqp1[0,jk,jl])
      zdr2 = zcons2*zdp[0,jk,jl]*zdq[0,0,jl]
      # Update rain fraction and freezing.
      # Note: impact of new temperature ZTP1 on ZFWAT is neglected here.
      if ztp1[0,jk,jl] < yrmcst.rtt:
        zrfreeze2 = zfwat[0,0,jl]*zdr2
        zfwatr = 0.0
      else:
        zrfreeze2 = 0.0
        zfwatr = 1.0
      zrn = zfwatr*zdr2
      zsn = (1.0 - zfwatr)*zdr2
      # Note: The extra condensation due to the adjustment goes directly to precipitation
      zcondl[0,jk,jl] = zcondl[0,jk,jl] + zfwatr*zdq[0,0,jl]*zqtmst
      zcondi[0,jk,jl] = zcondi[0,jk,jl] + (1.0 - zfwatr)*zdq[0,0,jl]*zqtmst
      zrfln[0,0,jl] = zrfln[0,0,jl] + zrn
      zsfln[0,0,jl] = zsfln[0,0,jl] + zsn
      zrfreeze[0,jk,jl] = zrfreeze[0,jk,jl] + zrfreeze2
    
    for jl in range(kfdia):
      zdqdt[0,0,jl] = -(zcondl[0,jk,jl] + zcondi[0,jk,jl]) + (plude[0,jk,jl] \
        + zevapr[0,jk,jl] + zevaps[0,jk,jl])*zgdp[0,0,jl]
      
      zdtdt[0,0,jl] = zlvdcp[0,jk,jl]*zcondl[0,jk,jl] + zlsdcp[0,jk,jl] * \
        zcondi[0,jk,jl] - (zlvdcp[0,jk,jl]*zevapr[0,jk,jl] + \
        zlsdcp[0,jk,jl]*zevaps[0,jk,jl] + plude[0,jk,jl]*(zfwat[0,0,jl] * \
        zlvdcp[0,jk,jl] + (1.0 - zfwat[0,0,jl])*zlsdcp[0,jk,jl]) - \
        (zlsdcp[0,jk,jl] - zlvdcp[0,jk,jl])*zrfreeze[0,jk,jl])*zgdp[0,0,jl]
      
      zdldt[0,0,jl] = (zqlwc[0,jk,jl] - zl[0,jk,jl])*zqtmst
      
      zdidt[0,0,jl] = (zqiwc[0,jk,jl] - zi[0,jk,jl])*zqtmst
      
      ptenq[0,jk,jl] = zdqdt[0,0,jl]
      ptent[0,jk,jl] = zdtdt[0,0,jl]
      ptenl[0,jk,jl] = zdldt[0,0,jl]
      pteni[0,jk,jl] = zdidt[0,0,jl]
      
      pfplsl[0,jk+1,jl] = zrfln[0,0,jl]
      pfplsn[0,jk+1,jl] = zsfln[0,0,jl]
    
    # record rain flux for next level
    
    for jl in range(kfdia):
      zrfl[0,0,jl] = zrfln[0,0,jl]
      zsfl[0,0,jl] = zsfln[0,0,jl]
    
  #jk
  
  #*     ENTHALPY FLUXES DUE TO PRECIPITATION
  #      ------------------------------------
  
  for jk in range(klev + 1):
    for jl in range(kfdia):
      pfhpsl[0,jk,jl] = -pfplsl[0,jk,jl]*yrmcst.rlvtt
      pfhpsn[0,jk,jl] = -pfplsn[0,jk,jl]*yrmcst.rlstt
  
  #     ------------------------------------------------------------------
  
  #IF (LHOOK) CALL DR_HOOK('CLOUDSC2',1,ZHOOK_HANDLE)
  return 
