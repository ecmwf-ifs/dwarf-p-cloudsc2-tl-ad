import numpy as np
from pathlib import Path
from collections import OrderedDict
from itertools import product

import gt4py
import gt4py.gtscript as gtscript
import gt4py.storage as gt_storage
from gt4py.gtscript import Field, IJK

from cloudsc2_inputs import load_input_parameters


__all__ = ['wrap_input_arrays', 'satur_py_gt4py', 'cloudsc2_py_gt4py', 'cloudsc2_gt4py']


backend = 'numpy'
dtype = np.float64
origin = (0, 0, 0)

# Get the constant parameters from file
rootpath = Path(__file__).resolve().parents[2]
input_path = rootpath/'config-files/input.h5'
yrecldp, yrmcst, yrethf, yrephli, yrecld = load_input_parameters(path=input_path)
klev = 137


#=======================================================================
#     TUNABLE CONSTANTS (to be moved to include files later)
#=======================================================================

# zscal is a scale factor that linearly reduces the variance between
# qv=qv-crit and qv-qsat
# 0 = No scaling
# 1 = full scaling (i.e. variance=0 when qv=qsat)
zscal = 0.9

#----------------------------------------------------------------------
#*    1.           DEFINE CONSTANTS
zqmax = 0.5
zeps1 = 1.E-12
zeps2 = 1.E-10


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
        pclc, pfplsl, pfplsn, pfhpsl, pfhpsn, pcovptot
):

    #=======================================================================

    zrfln = np.ndarray(order="C", shape=(iend,jend,klev))
    zsfln = np.ndarray(order="C", shape=(iend,jend,klev))
    zqold = np.ndarray(order="C", shape=(iend,jend,klev))
    zscalm = np.ndarray(order="C", shape=(iend,jend,klev,))
    ztrpaus = np.ndarray(order="C", shape=(iend,jend,1))
    #     ------------------------------------------------------------------
    paphp1_top = np.ndarray(order="C", shape=(iend,jend,klev))
    plu_p1 = np.ndarray(order="C", shape=(iend,jend,klev))
    #     ------------------------------------------------------------------
    ztp1 = np.ndarray(order="F", shape=(iend,jend,klev))
    zqp1 = np.ndarray(order="F", shape=(iend,jend,klev))

    # -------
    zeta3 = np.ndarray(order="F", shape=(iend,jend))

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
  
    #     --------------------------------------------------------------------
    #*         2.1    COMPUTE CRITICAL RELATIVE HUMIDITY AND RELATIVE HUMIDITY
    #                 --------------------------------------------------------

    # first guess values for T, q, ql and qi
    for jk in range(klev):
        for i, j in product(range(iend), range(jend)):
            ztp1[i,j,jk] = ptm1[i,j,jk] + ptsphy*pgtent[i,j,jk]
            zqp1[i,j,jk] = pqm1[i,j,jk] + ptsphy*pgtenq[i,j,jk] + psupsat[i,j,jk]
            # zl[i,j,jk] = pl[i,j,jk] + ptsphy*pgtenl[i,j,jk]
            # zi[i,j,jk] = pi[i,j,jk] + ptsphy*pgteni[i,j,jk]

    #     ------------------------------------------------------------------

    #*         2.2    INITIALIZATION OF CLOUD AND PRECIPITATION ARRAYS
    #                 ------------------------------------------------

    #       Set to zero precipitation fluxes at the top
    for jk in range(klev):
        for i, j in product(range(iend), range(jend)):
            zrfln[i,j,jk] = 0.0
            zsfln[i,j,jk] = 0.0
            pfplsl[i,j,jk] = 0.0
            pfplsn[i,j,jk] = 0.0

    # Eta value at tropopause
    for i, j in product(range(iend), range(jend)):
        ztrpaus[i,j,0] = 0.1
    for jk in range(klev - 1):
        for i, j in product(range(iend), range(jend)):
            if yrecld.ceta[jk] > 0.1 and yrecld.ceta[jk] < 0.4 and ztp1[i,j,jk] > ztp1[i,j,jk+1]:
                ztrpaus[i,j,0] = yrecld.ceta[jk]

    # set up critical value of humidity
    for jk in range(klev):
        for i, j in product(range(iend), range(jend)):
            zeta3[i,j] = ztrpaus[i,j,0]

    for jk in range(klev-1):
        for i, j in product(range(iend), range(jend)):
            plu_p1[i,j,jk] = plu[i,j,jk+1]

    #     ------------------------------------------------------------------

    #*        3. COMPUTE LAYER CLOUD AMOUNTS
    #            ---------------------------

    # Large loop over KLEV
    # Calculates
    #   1. diagnostic CC and QL
    #   2. Convective CC and QL
    #   3. Rainfall

    cloudsc2_kernel_gt4py(
        paphp1=paphp1, papp1=papp1, pqm1=pqm1, pqs=pqs, ptm1=ptm1,
        pl=pl, pi=pi, plude=plude, plu=plu, pmfu=pmfu, pmfd=pmfd,
        ptent=ptent, pgtent=pgtent, ptenq=ptenq, pgtenq=pgtenq,
        ptenl=ptenl, pgtenl=pgtenl, pteni=pteni, pgteni=pgteni,
        psupsat=psupsat, pclc=pclc, pfhpsl=pfhpsl, pfhpsn=pfhpsn,
        pfplsl=pfplsl, pfplsn=pfplsn, pcovptot=pcovptot,
        # ---------------------------
        zrfln=zrfln, zsfln=zsfln, ztp1=ztp1, zqp1=zqp1, zqold=zqold,
        # ---------------------------
        paphp1_top=paphp1_top, plu_p1=plu_p1, ceta=yrecld.ceta, zeta3=zeta3,
        # ---------------------------
        ptsphy=ptsphy, ldrain1d=ldrain1d,
        zckcodtl=zckcodtl, zckcodti=zckcodti, zcons2=zcons2,
        zmeltp2=zmeltp2, zqtmst=zqtmst, zcons3=zcons3,
    )


externals = {
    'yrmcst': yrmcst,
    'yrethf': yrethf,
    'yrephli': yrephli,
    'zscal': zscal,
}


@gtscript.function
def tanh(x):
    return (exp(x) - exp(-x)) / (exp(x) + exp(-x))

# Note that the use of compile-time constant structs kill the clean
# use of gtscript.function here!
# ------------
# @gtscript.function
# def foealfa_gt4py(ptare):
#   return min(1.0,((max(rtice,min(yrethf.rtwat,ptare))-yrethf.rtice)*yrethf.rtwat_rtice_r)**2)

# @gtscript.function
# def foeewm(ptare):
#   return yrethf.r2es*(foealfa(ptare, yrethf)*exp(yrethf.r3les*(ptare-yrmcst.rtt)/(ptare-yrethf.r4les)) + \
#                       (1.0-foealfa(ptare, yrethf))*exp(yrethf.r3ies*(ptare-yrmcst.rtt)/(ptare-yrethf.r4ies)))

@gtscript.stencil(backend=backend, externals=externals, rebuild=True)
def cloudsc2_kernel_gt4py(
        paphp1: Field[IJK, dtype],
        papp1: Field[IJK, dtype],
        pqm1: Field[IJK, dtype],
        pqs: Field[IJK, dtype],
        ptm1: Field[IJK,dtype],
        pl: Field[IJK,dtype],
        pi: Field[IJK,dtype],
        plude: Field[IJK,dtype],
        plu: Field[IJK, dtype],
        pmfu: Field[IJK, dtype],
        pmfd: Field[IJK, dtype],
        ptent: Field[IJK,dtype],
        pgtent: Field[IJK,dtype],
        ptenq: Field[IJK,dtype],
        pgtenq: Field[IJK,dtype],
        ptenl: Field[IJK,dtype],
        pgtenl: Field[IJK,dtype],
        pteni: Field[IJK,dtype],
        pgteni: Field[IJK,dtype],
        psupsat: Field[IJK, dtype],
        pclc: Field[IJK, dtype],
        pfhpsl: Field[IJK,dtype],
        pfhpsn: Field[IJK,dtype],
        pfplsl: Field[IJK,dtype],
        pfplsn: Field[IJK,dtype],
        pcovptot: Field[IJK, dtype],
        # ---------------------------
        zrfln: gtscript.Field[gtscript.IJK, dtype],
        zsfln: gtscript.Field[gtscript.IJK, dtype],
        ztp1: gtscript.Field[gtscript.IJK, dtype],
        zqp1: gtscript.Field[gtscript.IJK, dtype],
        # Note that removing this seems to upset things, not sure why...
        zqold: gtscript.Field[gtscript.IJK, dtype],
        # ---------------------------
        paphp1_top: gtscript.Field[gtscript.IJK, dtype],
        plu_p1: gtscript.Field[gtscript.IJK, dtype],
        ceta: gtscript.Field[gtscript.K, dtype],
        zeta3: gtscript.Field[gtscript.IJ, dtype],
        # ---------------------------
        ptsphy: np.float64,
        ldrain1d: bool,
        zckcodtl: np.float64,
        zckcodti: np.float64,
        zcons2: np.float64,
        zmeltp2: np.float64,
        zqtmst: np.float64,
        zcons3: np.float64,
):
    with computation(FORWARD), interval(...):

        # set up critical value of humidity
        zrh1 = 1.0
        zrh2 = 0.35 + 0.14*((zeta3[0,0] - 0.25) / 0.15)**2 + (0.04*min(zeta3[0,0] - 0.25, 0.0)) / 0.15
        zrh3 = 1.0
        zdeta2 = 0.3
        zdeta1 = 0.09 + (0.16*(0.4 - zeta3[0,0])) / 0.3

        if ceta[0] < zeta3[0,0]:
            zcrh2 = zrh3
        elif ceta[0] >= zeta3[0,0] and ceta[0] < (zeta3[0,0] + zdeta2):
            zcrh2 = zrh3 + (zrh2 - zrh3)*((ceta[0] - zeta3[0,0]) / zdeta2)
        elif ceta[0] >= (zeta3[0,0] + zdeta2) and ceta[0] < (1.0 - zdeta1):
            zcrh2 = zrh2
        elif ceta[0] >= (1.0 - zdeta1):
            zcrh2 = zrh1 + (zrh2 - zrh1)*((1.0 - ceta[0]) / zdeta1)**0.5

        # Parameter for cloud formation
        zscalm = zscal*max((ceta[0] - 0.2), zeps1)**0.2

        # Clear cloud and freezing arrays
        pclc[0,0,0] = 0.0
        pcovptot[0,0,0] = 0.0

        # first guess values for T, q, ql and qi
        # ztp1[0,0,0] = ptm1[0,0,0] + ptsphy*pgtent[0,0,0]
        # zqp1[0,0,0] = pqm1[0,0,0] + ptsphy*pgtenq[0,0,0] + psupsat[0,0,0]
        zl = pl[0,0,0] + ptsphy*pgtenl[0,0,0]
        zi = pi[0,0,0] + ptsphy*pgteni[0,0,0]

        # thermodynamic constants
        zdp = paphp1[0,0,+1] - paphp1[0,0,0]
        zzz = 1.0 / (yrmcst.rcpd + yrmcst.rcpd*yrethf.rvtmp2*zqp1[0,0,0])
        zlfdcp = yrmcst.rlmlt*zzz
        zlsdcp = yrmcst.rlstt*zzz
        zlvdcp = yrmcst.rlvtt*zzz

        #-----------------------------------
        # calculate dqs/dT correction factor
        #-----------------------------------

        zrfln[0,0,0] = zrfln[0,0,-1]
        zsfln[0,0,0] = zsfln[0,0,-1]
      
        if yrephli.lphylin or ldrain1d:
            rlptrc = yrephli.rlptrc
            zoealfaw = 0.545*(tanh(0.17*(ztp1[0,0,0] - rlptrc)) + 1.0)
            if ztp1[0,0,0] < yrmcst.rtt:
                zfwat = zoealfaw
                z3es = yrethf.r3ies
                z4es = yrethf.r4ies
            else:
                zfwat = 1.0
                z3es = yrethf.r3les
                z4es = yrethf.r4les
            zfoeew = yrethf.r2es*exp((z3es*(ztp1[0,0,0] - yrmcst.rtt)) / (ztp1[0,0,0] - z4es))
            zesdp = zfoeew / papp1[0,0,0]
            if zesdp > zqmax:
                zesdp = zqmax
        else:
            # Have to unroll these!
            # zfwat = foealfa(ztp1[0,0,0])
            zfwat = min(1.0,((max(yrethf.rtice,min(yrethf.rtwat,ztp1[0,0,0]))-yrethf.rtice)*yrethf.rtwat_rtice_r)**2)
            # zfoeew[0,0,0] = foeewm(ztp1[0,0,0], yrethf, yrmcst)
            # pfoealfa = foealfa(ztp1[0,0,0], yrethf)
            pfoealfa = min(1.0,((max(yrethf.rtice,min(yrethf.rtwat,ztp1[0,0,0]))-yrethf.rtice)*yrethf.rtwat_rtice_r)**2)
            zfoeew = yrethf.r2es*(pfoealfa*exp(yrethf.r3les*(ztp1[0,0,0]-yrmcst.rtt)/(ztp1[0,0,0]-yrethf.r4les)) \
                                       + (1.0-pfoealfa)*exp(yrethf.r3ies*(ztp1[0,0,0]-yrmcst.rtt)/(ztp1[0,0,0]-yrethf.r4ies)))

            zesdp = zfoeew / papp1[0,0,0]
        zfacw = yrethf.r5les / ((ztp1[0,0,0] - yrethf.r4les)**2)
        zfaci = yrethf.r5ies / ((ztp1[0,0,0] - yrethf.r4ies)**2)
        zfac = zfwat*zfacw + (1.0 - zfwat)*zfaci
        zcor = 1.0 / (1.0 - yrmcst.retv*zesdp)
        zdqsdtemp = zfac*zcor*pqs[0,0,0]
        zcorqs = 1.0 + zcons3*zdqsdtemp

        # use clipped state
      
        zqlim = zqp1[0,0,0]
        if zqp1[0,0,0] > pqs[0,0,0]:
            zqlim = pqs[0,0,0]
      
        # Allow ice supersaturation at cold temperatures
        if ztp1[0,0,0] < yrethf.rtice:
            zsupsat = 1.8 - 3.E-03*ztp1[0,0,0]
        else:
            zsupsat = 1.0
        zqsat = pqs[0,0,0]*zsupsat
        zqcrit = zcrh2*zqsat


        # simple UNIFORM distribution of total water from Letreut & Li (90)
    
        zqt = zqp1[0,0,0] + zl + zi
        if zqt <= zqcrit:
            pclc[0,0,0] = 0.0
            zqc = 0.0
        elif zqt >= zqsat:
            pclc[0,0,0] = 1.0
            zqc = (1.0 - zscalm)*(zqsat - zqcrit)
        else:
            zqpd = zqsat - zqt
            zqcd = zqsat - zqcrit
            pclc[0,0,0] = 1.0 - sqrt(zqpd / (zqcd - zscalm*(zqt - zqcrit)))
            zqc = (zscalm*zqpd + (1.0 - zscalm)*zqcd)*pclc[0,0,0]**2

        # Add convective component
    
        zgdp = yrmcst.rg / (paphp1[0,0,+1] - paphp1[0,0,0])
        zlude = plude[0,0,0]*ptsphy*zgdp
        if zlude >= yrecldp.rlmin and plu_p1[0,0,0] >= zeps2:
            pclc[0,0,0] = pclc[0,0,0] + (1.0 - pclc[0,0,0])*(1.0 - exp(-zlude / plu_p1[0,0,0]))
            zqc = zqc + zlude
    
        # Add compensating subsidence component
    
        zfac1 = 1.0 / ((yrmcst.rd*ztp1[0,0,0]))
        zrho = papp1[0,0,0]*zfac1
        zfac2 = 1.0 / (papp1[0,0,0] - yrmcst.retv*zfoeew)
        zrodqsdp = -zrho*pqs[0,0,0]*zfac2
        zldcp = zfwat*zlvdcp + (1.0 - zfwat)*zlsdcp
        zfac3 = 1.0 / (1.0 + zldcp*zdqsdtemp)
        dtdzmo = yrmcst.rg*(1.0 / yrmcst.rcpd - zldcp*zrodqsdp)*zfac3
        zdqsdz = zdqsdtemp*dtdzmo - yrmcst.rg*zrodqsdp
        zfac4 = 1.0 / zrho
        zdqc = min(zdqsdz*(pmfu[0,0,0] + pmfd[0,0,0])*ptsphy*zfac4, zqc)
        zqc = zqc - zdqc
    
        # New cloud liquid/ice contents and condensation rates (liquid/ice)
        zqlwc = zqc*zfwat
        zqiwc = zqc*(1.0 - zfwat)
        zcondl = (zqlwc - zl)*zqtmst
        zcondi = (zqiwc - zi)*zqtmst
    
        # Calculate precipitation overlap.
        # Simple form based on Maximum Overlap.
    
        if pclc[0,0,0] > pcovptot[0,0,0]:
            # total rain frac
            pcovptot[0,0,0] = pclc[0,0,0]
        zcovpclr = pcovptot[0,0,0] - pclc[0,0,0]        # clear sky frac
        zcovpclr = max(zcovpclr, 0.0)

    
    with computation(FORWARD), interval(...):

        #*         3.3    CALCULATE PRECIPITATION

        # Melting of incoming snow

        # if zsfln[0,0,0] != 0.0:
        zcons = (zcons2*zdp) / zlfdcp
        zsnmlt = min(zsfln[0,0,0], zcons*max(0.0, (ztp1[0,0,0] - zmeltp2)))
        zrfln[0,0,0] = zrfln[0,0,0] + zsnmlt
        zsfln[0,0,0] = zsfln[0,0,0] - zsnmlt
        ztp1[0,0,0] = ztp1[0,0,0] - zsnmlt / zcons


        #   Diagnostic calculation of rain production from cloud liquid water
      
        if pclc[0,0,0] > zeps2:
            # if yrphnc.levapls2 or ldrain1d:
            if False or ldrain1d:
                zlcrit = 1.9*yrecldp.rclcrit
            else:
                zlcrit = yrecldp.rclcrit*2.
            zcldl = zqlwc / pclc[0,0,0]          # in-cloud liquid
            zd = zckcodtl*(1.0 - exp(-(zcldl / zlcrit)**2))
            zlnew = pclc[0,0,0]*zcldl*exp(-zd)
            zprr = zqlwc - zlnew
            zqlwc = zqlwc - zprr
        else:
            zprr = 0.0
      
        #   Diagnostic calculation of snow production from cloud ice

        if pclc[0,0,0] > zeps2:
            # if yrphnc.levapls2 or ldrain1d:
            if False or ldrain1d:
                zlcrit = 1.E-04
            else:
                zlcrit = yrecldp.rclcrit*2.
            zcldi = zqiwc / pclc[0,0,0]          # in-cloud ice
            zd = zckcodti*exp(0.025*(ztp1[0,0,0] - yrmcst.rtt))*(1.0 - exp(-(zcldi / zlcrit)**2))
            zinew = pclc[0,0,0]*zcldi*exp(-zd)
            zprs = zqiwc - zinew
            zqiwc = zqiwc - zprs
        else:
            zprs = 0.0
      
        #   New precipitation (rain + snow)

        zdr = zcons2*zdp*(zprr + zprs)
      
        #   Rain fraction (different from cloud liquid water fraction!)
      
        if ztp1[0,0,0] < yrmcst.rtt:
            zrfreeze = zcons2*zdp*zprr
            zfwatr = 0.0
        else:
            zfwatr = 1.0
      
        zrn = zfwatr*zdr
        zsn = (1.0 - zfwatr)*zdr
        zrfln[0,0,0] = zrfln[0,0,0] + zrn
        zsfln[0,0,0] = zsfln[0,0,0] + zsn

    with computation(FORWARD), interval(...):
        #   Precip evaporation
        zprtot = zrfln[0,0,0] + zsfln[0,0,0]
        llo2 = zprtot > zeps2 and zcovpclr > zeps2 and ldrain1d
        if llo2:
            zpreclr = (zprtot*zcovpclr) / pcovptot[0,0,0]
        
            #     This is the humidity in the moistest zcovpclr region
        
            zqe = pqs[0,0,0] - ((pqs[0,0,0] - zqlim)*zcovpclr) \
                / (1.0 - pclc[0,0,0])**2
            zbeta = yrmcst.rg*yrecldp.rpecons*(((sqrt(papp1[0,0,0] / paphp1_top[0,0,0]) / \
                                                 5.09E-3)*zpreclr) / zcovpclr)**0.5777
        
            #     implicit solution:
            zb = (ptsphy*zbeta*(pqs[0,0,0] - zqe)) / (1.0 + zbeta*ptsphy*zcorqs)
        
            #     exact solution:
            #     ZB=(PQS(JL,JK)-ZQE)*(_ONE_-EXP(-ZBETA*ZCORQS(JL)*PTSPHY))/ZCORQS(JL)
        
            zdtgdp = (ptsphy*yrmcst.rg) / (paphp1[0,0,+1] - paphp1[0,0,0])
        
            zdpr = (zcovpclr*zb) / zdtgdp
            zdpr = min(zdpr, zpreclr)
            zpreclr = zpreclr - zdpr          # take away from clr sky flux
            if zpreclr <= 0.0:
                pcovptot[0,0,0] = pclc[0,0,0]
            #reset
            pcovptot[0,0,0] = pcovptot[0,0,0]
        
            # warm proportion
            zevapr = (zdpr*zrfln[0,0,0]) / zprtot
            zrfln[0,0,0] = zrfln[0,0,0] - zevapr
        
            # ice proportion
            zevaps = (zdpr*zsfln[0,0,0]) / zprtot
            zsfln[0,0,0] = zsfln[0,0,0] - zevaps

    
    # Update of T and Q tendencies due to:
    #  - condensation/evaporation of cloud liquid water/ice
    #  - detrainment of convective cloud condensate
    #  - evaporation of precipitation
    #  - freezing of rain (impact on T only).
    with computation(FORWARD), interval(...):
        zdqdt = -(zcondl + zcondi) + (plude[0,0,0] + \
            zevapr + zevaps)*zgdp
      
        zdtdt = zlvdcp*zcondl + zlsdcp * \
            zcondi - (zlvdcp*zevapr + \
            zlsdcp*zevaps + plude[0,0,0]*(zfwat * \
            zlvdcp + (1.0 - zfwat)*zlsdcp) - \
            (zlsdcp - zlvdcp)*zrfreeze)*zgdp
      
        # first guess T and Q
        ztp1[0,0,0] = ztp1[0,0,0] + ptsphy*zdtdt
        zqp1[0,0,0] = zqp1[0,0,0] + ptsphy*zdqdt

        zpp = papp1[0,0,0]
        zqold = zqp1[0,0,0]

    with computation(FORWARD), interval(...):
        if ztp1[0,0,0] > yrmcst.rtt:
            z3es = yrethf.r3les
            z4es = yrethf.r4les
            z5alcp = yrethf.r5alvcp
            zaldcp = yrethf.ralvdcp
        else:
            z3es = yrethf.r3ies
            z4es = yrethf.r4ies
            z5alcp = yrethf.r5alscp
            zaldcp = yrethf.ralsdcp
      
        zqp = 1.0 / zpp
        ztarg = ztp1[0,0,0]
        zfoeew = yrethf.r2es*exp((z3es*(ztarg - yrmcst.rtt)) / (ztarg - z4es))
        zqsat = zqp*zfoeew
        if zqsat > zqmax:
            zqsat = zqmax
        zcor = 1.0 / (1.0 - yrmcst.retv*zqsat)
        zqsat = zqsat*zcor
        z2s = z5alcp / (ztarg - z4es)**2
        zcond1 = (zqp1[0,0,0] - zqsat) / (1.0 + zqsat*zcor*z2s)
        ztp1[0,0,0] = ztp1[0,0,0] + zaldcp*zcond1
        zqp1[0,0,0] = zqp1[0,0,0] - zcond1
        ztarg = ztp1[0,0,0]
        zfoeew = yrethf.r2es*exp((z3es*(ztarg - yrmcst.rtt)) / (ztarg - z4es))
        zqsat = zqp*zfoeew
        if zqsat > zqmax:
            zqsat = zqmax
        zcor = 1.0 / (1.0 - yrmcst.retv*zqsat)
        zqsat = zqsat*zcor
        z2s = z5alcp / (ztarg - z4es)**2
        zcond1 = (zqp1[0,0,0] - zqsat) / (1.0 + zqsat*zcor*z2s)
        ztp1[0,0,0] = ztp1[0,0,0] + zaldcp*zcond1
        zqp1[0,0,0] = zqp1[0,0,0] - zcond1

    with computation(FORWARD), interval(...):
        zdq = max(0.0, zqold - zqp1[0,0,0])
        zdr2 = zcons2*zdp*zdq
        # Update rain fraction and freezing.
        # Note: impact of new temperature ZTP1 on ZFWAT is neglected here.
        if ztp1[0,0,0] < yrmcst.rtt:
            zrfreeze2 = zfwat*zdr2
            zfwatr = 0.0
        else:
            zrfreeze2 = 0.0
            zfwatr = 1.0
        zrn = zfwatr*zdr2
        zsn = (1.0 - zfwatr)*zdr2
        # Note: The extra condensation due to the adjustment goes directly to precipitation
        zcondl = zcondl + zfwatr*zdq*zqtmst
        zcondi = zcondi + (1.0 - zfwatr)*zdq*zqtmst
        zrfln[0,0,0] = zrfln[0,0,0] + zrn
        zsfln[0,0,0] = zsfln[0,0,0] + zsn
        zrfreeze = zrfreeze + zrfreeze2

    #*     ENTHALPY FLUXES DUE TO PRECIPITATION
    #      ------------------------------------

    with computation(PARALLEL), interval(...):
        zdqdt = -(zcondl + zcondi) + (plude[0,0,0] + \
            zevapr + zevaps)*zgdp

        zdtdt = zlvdcp*zcondl + zlsdcp * \
            zcondi - (zlvdcp*zevapr + \
            zlsdcp*zevaps + plude[0,0,0]*(zfwat * \
            zlvdcp + (1.0 - zfwat)*zlsdcp) - \
            (zlsdcp - zlvdcp)*zrfreeze)*zgdp
        
        ptenq = zdqdt
        ptent = zdtdt
        ptenl = (zqlwc - zl)*zqtmst
        pteni = (zqiwc - zi)*zqtmst
    
    with computation(PARALLEL), interval(1, klev+1):
        pfplsl = zrfln[0,0,-1]
        pfplsn = zsfln[0,0,-1]
    with computation(PARALLEL), interval(0, klev+1):
        pfhpsl = -pfplsl[0,0,0]*yrmcst.rlvtt
        pfhpsn = -pfplsn[0,0,0]*yrmcst.rlstt
