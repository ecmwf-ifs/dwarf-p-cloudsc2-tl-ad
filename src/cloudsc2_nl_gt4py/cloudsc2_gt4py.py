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
        shape = data.shape
        # ijk_data = np.reshape(data, (1, shape[0], shape[1]))
        ijk_data = np.reshape(np.transpose(data), (10, 10, shape[0]))
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


def satur_py_gt4py(iend, jend, klon, ktdia, klev, ldphylin, paprsf, pt, pqsat, kflag):

    #----------------------------------------------------------------------
    #*    1.           DEFINE CONSTANTS
    zqmax = 0.5

    #     *
    #----------------------------------------------------------------------
    #     *    2.           CALCULATE SATURATION SPECIFIC HUMIDITY
    #                       --------------------------------------

    if ldphylin:
        for jk in range(klev):
            for i, j in product(range(iend), range(jend)):
                ztarg = pt[i,j,jk]
                zalfa = foealfa(ztarg, yrethf)

                zfoeewl = yrethf.r2es*np.exp(yrethf.r3les*(ztarg-yrmcst.rtt)/(ztarg-yrethf.r4les))
                zfoeewi = yrethf.r2es*np.exp(yrethf.r3ies*(ztarg-yrmcst.rtt)/(ztarg-yrethf.r4ies))
                zfoeew = zalfa*zfoeewl+(1.0-zalfa)*zfoeewi

                zqs    = zfoeew/paprsf[i,j,jk]
                if zqs > zqmax:
                    zqs=zqmax

                zcor = 1.0/(1.0-yrmcst.retv*zqs)
                pqsat[i,j,jk]=zqs*zcor

    else:
        for jk in range(klev):
            for i, j in product(range(iend), range(jend)):
                if(kflag == 1):
                    zew  = foeewmcu(pt[i,j,jk])
                else:
                    zew  = foeewm(pt[i,j,jk])

                zqs  = zew/paprsf[i,j,jk]
                zqs  = min(zqmax,zqs)
                zcor = 1.0/(1.0-yrmcst.retv*zqs)
                pqsat[i,j,jk]=zqs*zcor


                
def cloudsc2_py_gt4py(
        iend, jend, klon, ktdia, klev, ldrain1d, ptsphy, paphp1,
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
  
  zrfl = np.ndarray(order="C", shape=(iend,jend,1))
  zsfl = np.ndarray(order="C", shape=(iend,jend,1))
  zrfln = np.ndarray(order="C", shape=(iend,jend,1))
  zsfln = np.ndarray(order="C", shape=(iend,jend,1))
  zgdp = np.ndarray(order="C", shape=(iend,jend,1))
  zdqdt = np.ndarray(order="C", shape=(iend,jend,1))
  zdtdt = np.ndarray(order="C", shape=(iend,jend,1))
  zdldt = np.ndarray(order="C", shape=(iend,jend,1))
  zdidt = np.ndarray(order="C", shape=(iend,jend,1))
  zqcrit = np.ndarray(order="C", shape=(iend,jend,1))
  zcovpclr = np.ndarray(order="C", shape=(iend,jend,1))
  zcovptot = np.ndarray(order="C", shape=(iend,jend,1))
  zdqsdtemp = np.ndarray(order="C", shape=(iend,jend,1))
  zcorqs = np.ndarray(order="C", shape=(iend,jend,1))
  zdtgdp = np.ndarray(order="C", shape=(iend,jend,1))
  zqold = np.ndarray(order="C", shape=(iend,jend,1))
  zpp = np.ndarray(order="C", shape=(iend,jend,1))
  zdq = np.ndarray(order="C", shape=(iend,jend,1))
  zqlim = np.ndarray(order="C", shape=(iend,jend,1))
  zscalm = np.ndarray(order="C", shape=(klev,))
  zqsat = np.ndarray(order="C", shape=(iend,jend,1))
  zfoeew = np.ndarray(order="C", shape=(iend,jend,1))
  zfwat = np.ndarray(order="C", shape=(iend,jend,1))
  
  
  ztrpaus = np.ndarray(order="C", shape=(iend,jend,1))
  
  
  llflag = np.ndarray(order="C", shape=(iend,jend,1))
  #REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
  
  #     ------------------------------------------------------------------
  ztp1 = np.ndarray(order="F", shape=(iend,jend,klev))
  zqp1 = np.ndarray(order="F", shape=(iend,jend,klev))
  zl = np.ndarray(order="F", shape=(iend,jend,klev))
  zi = np.ndarray(order="F", shape=(iend,jend,klev))
  zlude = np.ndarray(order="F", shape=(iend,jend,klev))
  zqc = np.ndarray(order="F", shape=(iend,jend,klev))
  zqlwc = np.ndarray(order="F", shape=(iend,jend,klev))
  zqiwc = np.ndarray(order="F", shape=(iend,jend,klev))
  zdp = np.ndarray(order="F", shape=(iend,jend,klev))
  zlsdcp = np.ndarray(order="F", shape=(iend,jend,klev))
  zlfdcp = np.ndarray(order="F", shape=(iend,jend,klev))
  zlvdcp = np.ndarray(order="F", shape=(iend,jend,klev))
  zrfreeze = np.ndarray(order="F", shape=(iend,jend,klev))
  zcondl = np.ndarray(order="F", shape=(iend,jend,klev))
  zcondi = np.ndarray(order="F", shape=(iend,jend,klev))
  zevapr = np.ndarray(order="F", shape=(iend,jend,klev))
  zevaps = np.ndarray(order="F", shape=(iend,jend,klev))
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
    for i, j in product(range(iend), range(jend)):
      ztp1[i,j,jk] = ptm1[i,j,jk] + ptsphy*pgtent[i,j,jk]
      zqp1[i,j,jk] = \
        pqm1[i,j,jk] + ptsphy*pgtenq[i,j,jk] + psupsat[i,j,jk]
      zl[i,j,jk] = pl[i,j,jk] + ptsphy*pgtenl[i,j,jk]
      zi[i,j,jk] = pi[i,j,jk] + ptsphy*pgteni[i,j,jk]
  
  for jk in range(klev):
    
    # Parameter for cloud formation

    zscalm[jk] = zscal*max((yrecld.ceta[jk] - 0.2), zeps1)**0.2
    
    for i, j in product(range(iend), range(jend)):
      
      # thermodynamic constants
      zdp[i,j,jk] = paphp1[i,j,jk+1] - paphp1[i,j,jk]
      zzz = 1.0 / (yrmcst.rcpd + yrmcst.rcpd*yrethf.rvtmp2*zqp1[i,j,jk])
      zlfdcp[i,j,jk] = yrmcst.rlmlt*zzz
      zlsdcp[i,j,jk] = yrmcst.rlstt*zzz
      zlvdcp[i,j,jk] = yrmcst.rlvtt*zzz
      llflag[i,j,0] = True
  
  #     ------------------------------------------------------------------
  
  #*         2.2    INITIALIZATION OF CLOUD AND PRECIPITATION ARRAYS
  #                 ------------------------------------------------
  
  #       Clear cloud and freezing arrays
  
  for jk in range(klev):
    for i, j in product(range(iend), range(jend)):
      pclc[i,j,jk] = 0.0
      zqc[i,j,jk] = 0.0
      zqlwc[i,j,jk] = 0.0
      zqiwc[i,j,jk] = 0.0
      zrfreeze[i,j,jk] = 0.0
      zcondl[i,j,jk] = 0.0
      zcondi[i,j,jk] = 0.0
      zevapr[i,j,jk] = 0.0
      zevaps[i,j,jk] = 0.0
      pcovptot[i,j,jk] = 0.0
  
  #       Set to zero precipitation fluxes at the top
  
  for i, j in product(range(iend), range(jend)):
    zrfl[i,j,0] = 0.0
    zsfl[i,j,0] = 0.0
    pfplsl[i,j,0] = 0.0
    pfplsn[i,j,0] = 0.0
    zcovptot[i,j,0] = 0.0
    zcovpclr[i,j,0] = 0.0
  
  # Eta value at tropopause
  for i, j in product(range(iend), range(jend)):
    ztrpaus[i,j,0] = 0.1
  for jk in range(klev - 1):
    for i, j in product(range(iend), range(jend)):
      llo1 = \
        yrecld.ceta[jk] > 0.1 and yrecld.ceta[jk] < 0.4 and ztp1[i,j,jk] > ztp1[i,j,jk+1]
      if llo1:
        ztrpaus[i,j,0] = yrecld.ceta[jk]
  
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
    
    for i, j in product(range(iend), range(jend)):
      
      #-----------------------------------
      # calculate dqs/dT correction factor
      #-----------------------------------
      
      if yrephli.lphylin or ldrain1d:
        zoealfaw = 0.545*(np.tanh(0.17*(ztp1[i,j,jk] - yrephli.rlptrc)) + 1.0)
        if ztp1[i,j,jk] < yrmcst.rtt:
          zfwat[i,j,0] = zoealfaw
          z3es = yrethf.r3ies
          z4es = yrethf.r4ies
        else:
          zfwat[i,j,0] = 1.0
          z3es = yrethf.r3les
          z4es = yrethf.r4les
        zfoeew[i,j,0] = \
          yrethf.r2es*np.exp((z3es*(ztp1[i,j,jk] - yrmcst.rtt)) / (ztp1[i,j,jk] - z4es))
        zesdp = zfoeew[i,j,0] / papp1[i,j,jk]
        if zesdp > zqmax:
          zesdp = zqmax
      else:
        zfwat[i,j,0] = foealfa(ztp1[i,j,jk], yrethf)
        zfoeew[i,j,0] = foeewm(ztp1[i,j,jk], yrethf, yrmcst)
        zesdp = zfoeew[i,j,0] / papp1[i,j,jk]
      zfacw = yrethf.r5les / ((ztp1[i,j,jk] - yrethf.r4les)**2)
      zfaci = yrethf.r5ies / ((ztp1[i,j,jk] - yrethf.r4ies)**2)
      zfac = zfwat[i,j,0]*zfacw + (1.0 - zfwat[i,j,0])*zfaci
      zcor = 1.0 / (1.0 - yrmcst.retv*zesdp)
      zdqsdtemp[i,j,0] = zfac*zcor*pqs[i,j,jk]
      zcorqs[i,j,0] = 1.0 + zcons3*zdqsdtemp[i,j,0]
      
      # use clipped state
      
      zqlim[i,j,0] = zqp1[i,j,jk]
      if zqp1[i,j,jk] > pqs[i,j,jk]:
        zqlim[i,j,0] = pqs[i,j,jk]
      
      # set up critical value of humidity
      
      zeta3 = ztrpaus[i,j,0]
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
      if ztp1[i,j,jk] < yrethf.rtice:
        zsupsat = 1.8 - 3.E-03*ztp1[i,j,jk]
      else:
        zsupsat = 1.0
      zqsat[i,j,0] = pqs[i,j,jk]*zsupsat
      zqcrit[i,j,0] = zcrh2*zqsat[i,j,0]
    
    # Simple UNIFORM distribution of total water from Letreut & Li (90)
    
    for i, j in product(range(iend), range(jend)):
      zqt = zqp1[i,j,jk] + zl[i,j,jk] + zi[i,j,jk]
      if zqt <= zqcrit[i,j,0]:
        pclc[i,j,jk] = 0.0
        zqc[i,j,jk] = 0.0
      elif zqt >= zqsat[i,j,0]:
        pclc[i,j,jk] = 1.0
        zqc[i,j,jk] = (1.0 - zscalm[jk])*(zqsat[i,j,0] - zqcrit[i,j,0])
      else:
        zqpd = zqsat[i,j,0] - zqt
        zqcd = zqsat[i,j,0] - zqcrit[i,j,0]
        pclc[i,j,jk] = \
          1.0 - np.sqrt(zqpd / (zqcd - zscalm[jk]*(zqt - zqcrit[i,j,0])))
        zqc[i,j,jk] = \
          (zscalm[jk]*zqpd + (1.0 - zscalm[jk])*zqcd)*pclc[i,j,jk]**2
    
    # Add convective component
    
    for i, j in product(range(iend), range(jend)):
      zgdp[i,j,0] = yrmcst.rg / (paphp1[i,j,jk+1] - paphp1[i,j,jk])
      zlude[i,j,jk] = plude[i,j,jk]*ptsphy*zgdp[i,j,0]
      if jk < klev:
        llo1 = zlude[i,j,jk] >= yrecldp.rlmin and plu[i,j,jk+1] >= zeps2
      else:
        llo1 = False
      if llo1:
        pclc[i,j,jk] = pclc[i,j,jk] + (1.0 - pclc[i,j,jk])*(1.0 - \
          np.exp(-zlude[i,j,jk] / plu[i,j,jk+1]))
        zqc[i,j,jk] = zqc[i,j,jk] + zlude[i,j,jk]
    
    # Add compensating subsidence component
    
    for i, j in product(range(iend), range(jend)):
      zfac1 = 1.0 / ((yrmcst.rd*ztp1[i,j,jk]))
      zrho = papp1[i,j,jk]*zfac1
      zfac2 = 1.0 / (papp1[i,j,jk] - yrmcst.retv*zfoeew[i,j,0])
      zrodqsdp = -zrho*pqs[i,j,jk]*zfac2
      zldcp = \
        zfwat[i,j,0]*zlvdcp[i,j,jk] + (1.0 - zfwat[i,j,0])*zlsdcp[i,j,jk]
      zfac3 = 1.0 / (1.0 + zldcp*zdqsdtemp[i,j,0])
      dtdzmo = yrmcst.rg*(1.0 / yrmcst.rcpd - zldcp*zrodqsdp)*zfac3
      zdqsdz = zdqsdtemp[i,j,0]*dtdzmo - yrmcst.rg*zrodqsdp
      zfac4 = 1.0 / zrho
      zdqc = min(zdqsdz*(pmfu[i,j,jk] + pmfd[i,j,jk])*ptsphy*zfac4, zqc[i,j,jk])
      zqc[i,j,jk] = zqc[i,j,jk] - zdqc
    
    # New cloud liquid/ice contents and condensation rates (liquid/ice)
    
    for i, j in product(range(iend), range(jend)):
      zqlwc[i,j,jk] = zqc[i,j,jk]*zfwat[i,j,0]
      zqiwc[i,j,jk] = zqc[i,j,jk]*(1.0 - zfwat[i,j,0])
      zcondl[i,j,jk] = (zqlwc[i,j,jk] - zl[i,j,jk])*zqtmst
      zcondi[i,j,jk] = (zqiwc[i,j,jk] - zi[i,j,jk])*zqtmst
    
    
    # Calculate precipitation overlap.
    # Simple form based on Maximum Overlap.
    
    for i, j in product(range(iend), range(jend)):
      if pclc[i,j,jk] > zcovptot[i,j,0]:
        # total rain frac
        zcovptot[i,j,0] = pclc[i,j,jk]
      zcovpclr[i,j,0] = zcovptot[i,j,0] - pclc[i,j,jk]        # clear sky frac
      zcovpclr[i,j,0] = max(zcovpclr[i,j,0], 0.0)
    
    #*         3.3    CALCULATE PRECIPITATION
    
    # Melting of incoming snow
    
    for i, j in product(range(iend), range(jend)):
      if zsfl[i,j,0] != 0.0:
        zcons = (zcons2*zdp[i,j,jk]) / zlfdcp[i,j,jk]
        zsnmlt = min(zsfl[i,j,0], zcons*max(0.0, (ztp1[i,j,jk] - zmeltp2)))
        zrfln[i,j,0] = zrfl[i,j,0] + zsnmlt
        zsfln[i,j,0] = zsfl[i,j,0] - zsnmlt
        ztp1[i,j,jk] = ztp1[i,j,jk] - zsnmlt / zcons
      else:
        zrfln[i,j,0] = zrfl[i,j,0]
        zsfln[i,j,0] = zsfl[i,j,0]
    
    for i, j in product(range(iend), range(jend)):
      
      #   Diagnostic calculation of rain production from cloud liquid water
      
      if pclc[i,j,jk] > zeps2:
        # if yrphnc.levapls2 or ldrain1d:
        if False or ldrain1d:
          zlcrit = 1.9*yrecldp.rclcrit
        else:
          zlcrit = yrecldp.rclcrit*2.
        zcldl = zqlwc[i,j,jk] / pclc[i,j,jk]          # in-cloud liquid
        zd = zckcodtl*(1.0 - np.exp(-(zcldl / zlcrit)**2))
        zlnew = pclc[i,j,jk]*zcldl*np.exp(-zd)
        zprr = zqlwc[i,j,jk] - zlnew
        zqlwc[i,j,jk] = zqlwc[i,j,jk] - zprr
      else:
        zprr = 0.0
      
      #   Diagnostic calculation of snow production from cloud ice
      
      if pclc[i,j,jk] > zeps2:
        # if yrphnc.levapls2 or ldrain1d:
        if False or ldrain1d:
          zlcrit = 1.E-04
        else:
          zlcrit = yrecldp.rclcrit*2.
        zcldi = zqiwc[i,j,jk] / pclc[i,j,jk]          # in-cloud ice
        zd = zckcodti*np.exp(0.025*(ztp1[i,j,jk] - yrmcst.rtt))*(1.0 - np.exp(-(zcldi / \
          zlcrit)**2))
        zinew = pclc[i,j,jk]*zcldi*np.exp(-zd)
        zprs = zqiwc[i,j,jk] - zinew
        zqiwc[i,j,jk] = zqiwc[i,j,jk] - zprs
      else:
        zprs = 0.0
      
      #   New precipitation (rain + snow)
      
      zdr = zcons2*zdp[i,j,jk]*(zprr + zprs)
      
      #   Rain fraction (different from cloud liquid water fraction!)
      
      if ztp1[i,j,jk] < yrmcst.rtt:
        zrfreeze[i,j,jk] = zcons2*zdp[i,j,jk]*zprr
        zfwatr = 0.0
      else:
        zfwatr = 1.0
      
      zrn = zfwatr*zdr
      zsn = (1.0 - zfwatr)*zdr
      zrfln[i,j,0] = zrfln[i,j,0] + zrn
      zsfln[i,j,0] = zsfln[i,j,0] + zsn
      
      #   Precip evaporation
      
      zprtot = zrfln[i,j,0] + zsfln[i,j,0]
      llo2 = zprtot > zeps2 and zcovpclr[i,j,0] > zeps2 and ldrain1d
      if llo2:
        
        zpreclr = (zprtot*zcovpclr[i,j,0]) / zcovptot[i,j,0]
        
        #     This is the humidity in the moistest zcovpclr region
        
        zqe = pqs[i,j,jk] - ((pqs[i,j,jk] - zqlim[i,j,0])*zcovpclr[i,j,0]) \
          / (1.0 - pclc[i,j,jk])**2
        zbeta = yrmcst.rg*yrecldp.rpecons*(((np.sqrt(papp1[i,j,jk] / paphp1[klev, jl]) / \
          5.09E-3)*zpreclr) / zcovpclr[i,j,0])**0.5777
        
        #     implicit solution:
        zb = (ptsphy*zbeta*(pqs[i,j,jk] - zqe)) / (1.0 + zbeta*ptsphy*zcorqs[i,j,0])
        
        #     exact solution:
        #     ZB=(PQS(JL,JK)-ZQE)*(_ONE_-EXP(-ZBETA*ZCORQS(JL)*PTSPHY))/ZCORQS(JL)
        
        zdtgdp[i,j,0] = (ptsphy*yrmcst.rg) / (paphp1[i,j,jk+1] - paphp1[i,j,jk])
        
        zdpr = (zcovpclr[i,j,0]*zb) / zdtgdp[i,j,0]
        zdpr = min(zdpr, zpreclr)
        zpreclr = zpreclr - zdpr          # take away from clr sky flux
        if zpreclr <= 0.0:
          zcovptot[i,j,0] = pclc[i,j,jk]
        #reset
        pcovptot[i,j,jk] = zcovptot[i,j,0]
        
        # warm proportion
        zevapr[i,j,jk] = (zdpr*zrfln[i,j,0]) / zprtot
        zrfln[i,j,0] = zrfln[i,j,0] - zevapr[i,j,jk]
        
        # ice proportion
        zevaps[i,j,jk] = (zdpr*zsfln[i,j,0]) / zprtot
        zsfln[i,j,0] = zsfln[i,j,0] - zevaps[i,j,jk]
      
    
    # Update of T and Q tendencies due to:
    #  - condensation/evaporation of cloud liquid water/ice
    #  - detrainment of convective cloud condensate
    #  - evaporation of precipitation
    #  - freezing of rain (impact on T only).
    
    for i, j in product(range(iend), range(jend)):
      zdqdt[i,j,0] = -(zcondl[i,j,jk] + zcondi[i,j,jk]) + (plude[i,j,jk] \
        + zevapr[i,j,jk] + zevaps[i,j,jk])*zgdp[i,j,0]
      
      zdtdt[i,j,0] = zlvdcp[i,j,jk]*zcondl[i,j,jk] + zlsdcp[i,j,jk] \
        *zcondi[i,j,jk] - (zlvdcp[i,j,jk]*zevapr[i,j,jk] + \
        zlsdcp[i,j,jk]*zevaps[i,j,jk] + plude[i,j,jk]*(zfwat[i,j,0] \
        *zlvdcp[i,j,jk] + (1.0 - zfwat[i,j,0])*zlsdcp[i,j,jk]) - \
        (zlsdcp[i,j,jk] - zlvdcp[i,j,jk])*zrfreeze[i,j,jk])*zgdp[i,j,0]
      
      # first guess T and Q
      ztp1[i,j,jk] = ztp1[i,j,jk] + ptsphy*zdtdt[i,j,0]
      zqp1[i,j,jk] = zqp1[i,j,jk] + ptsphy*zdqdt[i,j,0]
      
      zpp[i,j,0] = papp1[i,j,jk]
      zqold[i,j,0] = zqp1[i,j,jk]
    
    # clipping of final qv
    
    # -----------------------------------
    # IK=JK
    # ICALL=0
    # CALL CUADJTQS ( IEND, JEND, KLON, KLEV, IK,&
    #   & ZPP  , ZTP1  , ZQP1 , LLFLAG, ICALL  )
    # -----------------------------------
    # Manually inlined CUADJTQS
    # -----------------------------------
    zqmax = 0.5
    for i, j in product(range(iend), range(jend)):
      if ztp1[i,j,jk] > yrmcst.rtt:
        z3es = yrethf.r3les
        z4es = yrethf.r4les
        z5alcp = yrethf.r5alvcp
        zaldcp = yrethf.ralvdcp
      else:
        z3es = yrethf.r3ies
        z4es = yrethf.r4ies
        z5alcp = yrethf.r5alscp
        zaldcp = yrethf.ralsdcp
      
      zqp = 1.0 / zpp[i,j,0]
      ztarg = ztp1[i,j,jk]
      zfoeew[i,j,0] = yrethf.r2es*np.exp((z3es*(ztarg - yrmcst.rtt)) / (ztarg - z4es))
      zqsat[i,j,0] = zqp*zfoeew[i,j,0]
      if zqsat[i,j,0] > zqmax:
        zqsat[i,j,0] = zqmax
      zcor = 1.0 / (1.0 - yrmcst.retv*zqsat[i,j,0])
      zqsat[i,j,0] = zqsat[i,j,0]*zcor
      z2s = z5alcp / (ztarg - z4es)**2
      zcond1 = (zqp1[i,j,jk] - zqsat[i,j,0]) / (1.0 + zqsat[i,j,0]*zcor*z2s)
      ztp1[i,j,jk] = ztp1[i,j,jk] + zaldcp*zcond1
      zqp1[i,j,jk] = zqp1[i,j,jk] - zcond1
      ztarg = ztp1[i,j,jk]
      zfoeew[i,j,0] = yrethf.r2es*np.exp((z3es*(ztarg - yrmcst.rtt)) / (ztarg - z4es))
      zqsat[i,j,0] = zqp*zfoeew[i,j,0]
      if zqsat[i,j,0] > zqmax:
        zqsat[i,j,0] = zqmax
      zcor = 1.0 / (1.0 - yrmcst.retv*zqsat[i,j,0])
      zqsat[i,j,0] = zqsat[i,j,0]*zcor
      z2s = z5alcp / (ztarg - z4es)**2
      zcond1 = (zqp1[i,j,jk] - zqsat[i,j,0]) / (1.0 + zqsat[i,j,0]*zcor*z2s)
      ztp1[i,j,jk] = ztp1[i,j,jk] + zaldcp*zcond1
      zqp1[i,j,jk] = zqp1[i,j,jk] - zcond1
    # -----------------------------------
    
    for i, j in product(range(iend), range(jend)):
      zdq[i,j,0] = max(0.0, zqold[i,j,0] - zqp1[i,j,jk])
      zdr2 = zcons2*zdp[i,j,jk]*zdq[i,j,0]
      # Update rain fraction and freezing.
      # Note: impact of new temperature ZTP1 on ZFWAT is neglected here.
      if ztp1[i,j,jk] < yrmcst.rtt:
        zrfreeze2 = zfwat[i,j,0]*zdr2
        zfwatr = 0.0
      else:
        zrfreeze2 = 0.0
        zfwatr = 1.0
      zrn = zfwatr*zdr2
      zsn = (1.0 - zfwatr)*zdr2
      # Note: The extra condensation due to the adjustment goes directly to precipitation
      zcondl[i,j,jk] = zcondl[i,j,jk] + zfwatr*zdq[i,j,0]*zqtmst
      zcondi[i,j,jk] = zcondi[i,j,jk] + (1.0 - zfwatr)*zdq[i,j,0]*zqtmst
      zrfln[i,j,0] = zrfln[i,j,0] + zrn
      zsfln[i,j,0] = zsfln[i,j,0] + zsn
      zrfreeze[i,j,jk] = zrfreeze[i,j,jk] + zrfreeze2
    
    for i, j in product(range(iend), range(jend)):
      zdqdt[i,j,0] = -(zcondl[i,j,jk] + zcondi[i,j,jk]) + (plude[i,j,jk] \
        + zevapr[i,j,jk] + zevaps[i,j,jk])*zgdp[i,j,0]
      
      zdtdt[i,j,0] = zlvdcp[i,j,jk]*zcondl[i,j,jk] + zlsdcp[i,j,jk] * \
        zcondi[i,j,jk] - (zlvdcp[i,j,jk]*zevapr[i,j,jk] + \
        zlsdcp[i,j,jk]*zevaps[i,j,jk] + plude[i,j,jk]*(zfwat[i,j,0] * \
        zlvdcp[i,j,jk] + (1.0 - zfwat[i,j,0])*zlsdcp[i,j,jk]) - \
        (zlsdcp[i,j,jk] - zlvdcp[i,j,jk])*zrfreeze[i,j,jk])*zgdp[i,j,0]
      
      zdldt[i,j,0] = (zqlwc[i,j,jk] - zl[i,j,jk])*zqtmst
      
      zdidt[i,j,0] = (zqiwc[i,j,jk] - zi[i,j,jk])*zqtmst
      
      ptenq[i,j,jk] = zdqdt[i,j,0]
      ptent[i,j,jk] = zdtdt[i,j,0]
      ptenl[i,j,jk] = zdldt[i,j,0]
      pteni[i,j,jk] = zdidt[i,j,0]
      
      pfplsl[i,j,jk+1] = zrfln[i,j,0]
      pfplsn[i,j,jk+1] = zsfln[i,j,0]
    
    # record rain flux for next level
    
    for i, j in product(range(iend), range(jend)):
      zrfl[i,j,0] = zrfln[i,j,0]
      zsfl[i,j,0] = zsfln[i,j,0]
    
  #jk
  
  #*     ENTHALPY FLUXES DUE TO PRECIPITATION
  #      ------------------------------------
  
  for jk in range(klev + 1):
    for i, j in product(range(iend), range(jend)):
      pfhpsl[i,j,jk] = -pfplsl[i,j,jk]*yrmcst.rlvtt
      pfhpsn[i,j,jk] = -pfplsn[i,j,jk]*yrmcst.rlstt
  
  #     ------------------------------------------------------------------
  
  #IF (LHOOK) CALL DR_HOOK('CLOUDSC2',1,ZHOOK_HANDLE)
  return 
