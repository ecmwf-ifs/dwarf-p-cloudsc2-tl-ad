# -*- coding: utf-8 -*-
from gt4py.cartesian import gtscript

from cloudsc2py.framework.stencil import function_collection
from cloudsc2py.physics.nonlinear.stencils.cuadjtqs import cuadjtqs_nl
from cloudsc2py.utils.f2py import ported_function


@ported_function(from_file="cloudsc2_ad/cloudsc2ad.F90", from_line=339, to_line=866)
@function_collection("cloudsc2_ad_0")
@gtscript.function
def cloudsc2_ad_0(
    in_eta,
    in_ap,
    in_aph,
    in_q,
    in_qsat,
    in_ql,
    in_qi,
    in_lu,
    in_lude,
    in_mfd,
    in_mfu,
    in_supsat,
    in_tnd_cml_q,
    in_tnd_cml_ql,
    in_tnd_cml_qi,
    tmp_aph_s,
    tmp_trpaus,
    dt,
    covptot,
    rfl,
    sfl,
    t,
):
    from __externals__ import (
        LDRAIN1D,
        LEVAPLS2,
        R2ES,
        R3IES,
        R3LES,
        R4IES,
        R4LES,
        R5IES,
        R5LES,
        RCLCRIT,
        RCPD,
        RD,
        RETV,
        RG,
        RKCONV,
        RLMIN,
        RLMLT,
        RLPTRC,
        RLSTT,
        RLVTT,
        RPECONS,
        RTICE,
        RTT,
        RVTMP2,
        ZEPS1,
        ZEPS2,
        ZQMAX,
        ZSCAL,
    )

    # set up constants required
    ckcodtl = 2 * RKCONV * dt
    ckcodti = 5 * RKCONV * dt
    cons2 = 1 / (RG * dt)
    cons3 = RLVTT / RCPD
    meltp2 = RTT + 2

    # first guess values for q, ql and qi
    q = in_q[0, 0, 0] + dt * in_tnd_cml_q[0, 0, 0] + in_supsat[0, 0, 0]
    ql = in_ql[0, 0, 0] + dt * in_tnd_cml_ql[0, 0, 0]
    qi = in_qi[0, 0, 0] + dt * in_tnd_cml_qi[0, 0, 0]
    # store trajectory arrays for adjoint
    q2 = q + 0

    # parameter for cloud formation
    scalm = ZSCAL * max(in_eta[0] - 0.2, ZEPS1) ** 0.2

    # thermodynamic constants
    dp = in_aph[0, 0, 1] - in_aph[0, 0, 0]
    zz = RCPD + RCPD * RVTMP2 * q
    lfdcp = RLMLT / zz
    lsdcp = RLSTT / zz
    lvdcp = RLVTT / zz

    # calculate dqs/dT correction factor
    # if __INLINED(LPHYLIN or LDRAIN1D):
    if t < RTT:
        fwat = 0.545 * (tanh(0.17 * (t - RLPTRC)) + 1)
        z3es = R3IES
        z4es = R4IES
    else:
        fwat = 1.0
        z3es = R3LES
        z4es = R4LES
    foeew = R2ES * exp(z3es * (t - RTT) / (t - z4es))
    esdp1 = foeew / in_ap[0, 0, 0]
    esdp = min(esdp1, ZQMAX)
    # else:
    #     fwat = foealfa(t)
    #     foeew = foeewm(t)
    #     esdp = foeew / in_ap[0, 0, 0]
    #     esdp1 = foeew / in_ap[0, 0, 0]
    facw = R5LES / ((t - R4LES) ** 2)
    faci = R5IES / ((t - R4IES) ** 2)
    fac = fwat * facw + (1 - fwat) * faci
    cor = 1 / (1 - RETV * esdp)
    dqsdtemp = fac * in_qsat[0, 0, 0] / (1 - RETV * esdp)
    corqs = 1 + cons3 * dqsdtemp

    # use clipped state
    qlim = min(q, in_qsat[0, 0, 0])

    # set up critical value of humidity
    rh1 = 1.0
    rh2 = (
        0.35 + 0.14 * ((tmp_trpaus[0, 0] - 0.25) / 0.15) ** 2 + 0.04 * min(tmp_trpaus[0, 0] - 0.25, 0.0) / 0.15
    )
    rh3 = 1.0
    if in_eta[0] < tmp_trpaus[0, 0]:
        crh2 = rh3
    else:
        deta2 = 0.3
        bound1 = tmp_trpaus[0, 0] + deta2
        if in_eta[0] < bound1:
            crh2 = rh3 + (rh2 - rh3) * (in_eta[0] - tmp_trpaus[0, 0]) / deta2
        else:
            deta1 = 0.09 + 0.16 * (0.4 - tmp_trpaus[0, 0]) / 0.3
            bound2 = 1 - deta1
            if in_eta[0] < bound2:
                crh2 = rh2
            else:
                crh2 = rh1 + (rh2 - rh1) * sqrt((1 - in_eta[0]) / deta1)

    # allow ice supersaturation at cold temperatures
    if t < RTICE:
        supsat = 1.8 - 0.003 * t
        qsat = in_qsat[0, 0, 0] * supsat
    else:
        supsat = 1.0
        qsat = in_qsat[0, 0, 0] + 0
    qcrit = crh2 * qsat

    # simple uniform distribution of total water from Leutreut & Li (1990)
    qt = q + ql + qi
    if qt < qcrit:
        clc = 0.0
        qc1 = 0.0
        qcd = 0.0
        qpd = 0.0
        tmp3 = 0.0
    elif qt >= qsat:
        clc = 1.0
        qc1 = (1 - scalm) * (qsat - qcrit)
        qcd = 0.0
        qpd = 0.0
        tmp3 = 0.0
    else:
        qcd = qsat - qcrit
        qpd = qsat - qt
        tmp3 = sqrt(qpd / (qcd - scalm * (qt - qcrit)))
        clc = 1 - tmp3
        qc1 = (scalm * qpd + (1 - scalm) * qcd) * (clc**2)

    # add convective component
    gdp = RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
    lude = dt * in_lude[0, 0, 0] * gdp
    lo1 = lude[0, 0, 0] >= RLMIN and in_lu[0, 0, 1] >= ZEPS2
    if lo1:
        out_clc = clc[0, 0, 0] + (1 - clc[0, 0, 0]) * (1 - exp(-lude[0, 0, 0] / in_lu[0, 0, 1]))
        qc2 = qc1 + lude
    else:
        out_clc = clc + 0
        qc2 = qc1 + 0

    # add compensating subsidence component
    fac1 = 1 / (RD * t)
    rho = in_ap[0, 0, 0] * fac1
    fac2 = 1 / (in_ap[0, 0, 0] - RETV * foeew)
    rodqsdp = -rho * in_qsat[0, 0, 0] * fac2
    ldcp = fwat * lvdcp + (1 - fwat) * lsdcp
    fac3 = 1 / (1 + ldcp * dqsdtemp)
    dtdzmo = RG * (1 / RCPD - ldcp * rodqsdp) * fac3
    dqsdz = dqsdtemp * dtdzmo - RG * rodqsdp
    lo3 = dt * dqsdz * (in_mfu[0, 0, 0] + in_mfd[0, 0, 0]) / rho < qc2
    dqc = min(dt * dqsdz * (in_mfu[0, 0, 0] + in_mfd[0, 0, 0]) / rho, qc2)
    qc3 = qc2 - dqc

    # new cloud liquid/ice contents and condensation rates (liquid/ice)
    qlwc1 = qc3 * fwat
    qiwc1 = qc3 * (1 - fwat)
    condl1 = (qlwc1 - ql) / dt
    condi1 = (qiwc1 - qi) / dt

    # calculate precipitation overlap
    # simple form based on Maximum Overlap
    covptot1 = max(covptot, out_clc)
    covptotn = covptot1 + 0
    covpclr1 = covptotn - out_clc
    covpclr = max(covpclr1, 0.0)

    # melting of incoming snow
    if sfl != 0:
        cons = cons2 * dp / lfdcp
        z2s = cons * max(t - meltp2, 0.0)
        snmlt = min(sfl, z2s)
        rfln = rfl + snmlt
        sfln = sfl - snmlt
        t -= snmlt / cons
    else:
        cons = 0.0
        rfln = rfl + 0
        sfln = sfl + 0
        snmlt = 0.0
        z2s = 0.0

    if out_clc > ZEPS2:
        # diagnostic calculation of rain production from cloud liquid water
        if __INLINED(LEVAPLS2 or LDRAIN1D):
            lcrit = 1.9 * RCLCRIT
        else:
            lcrit = 2.0 * RCLCRIT
        cldl = qlwc1 / out_clc
        ltmp1 = exp(-((cldl / lcrit) ** 2))
        ltmp2 = exp(-ckcodtl * (1 - ltmp1))
        prr = qlwc1 - out_clc * cldl * ltmp2
        qlwc = qlwc1 - prr

        # diagnostic calculation of snow production from cloud ice
        if __INLINED(LEVAPLS2 or LDRAIN1D):
            icrit = 0.0001
        else:
            icrit = 2 * RCLCRIT
        cldi = qiwc1 / out_clc
        itmp11 = exp(-((cldi / icrit) ** 2))
        itmp12 = exp(0.025 * (t - RTT))
        itmp2 = exp(-ckcodti * itmp12 * (1 - itmp11))
        prs = qiwc1 - out_clc * cldi * itmp2
        qiwc = qiwc1 - prs
    else:
        cldi = 0.0
        cldl = 0.0
        itmp11 = 0.0
        itmp12 = 0.0
        itmp2 = 0.0
        ltmp1 = 0.0
        ltmp2 = 0.0
        qiwc = qiwc1 + 0
        qlwc = qlwc1 + 0
        prr = 0.0
        prs = 0.0

    # new precipitation (rain + snow)
    dr1 = cons2 * dp * (prr + prs)

    # rain fraction (different from cloud liquid water fraction!)
    if t < RTT:
        rfreeze1 = cons2 * dp * prr
        fwatr1 = 0.0
    else:
        rfreeze1 = 0.0
        fwatr1 = 1.0
    rfln += fwatr1 * dr1
    sfln += (1 - fwatr1) * dr1
    # store trajectory for adjoint
    rfln2 = rfln + 0
    sfln2 = sfln + 0

    # precipitation evaporation
    prtot = rfln + sfln
    if prtot > ZEPS2 and covpclr > ZEPS2 and (LEVAPLS2 or LDRAIN1D):
        preclr1 = prtot * covpclr / covptot1

        # this is the humidity in the moisest zcovpclr region
        qe = in_qsat[0, 0, 0] - (in_qsat[0, 0, 0] - qlim) * covpclr / (1 - out_clc) ** 2
        beta = RG * RPECONS * (preclr1 * sqrt(in_ap[0, 0, 0] / tmp_aph_s[0, 0]) / (0.00509 * covpclr)) ** 0.5777

        # implicit solution
        b = dt * beta * (in_qsat[0, 0, 0] - qe) / (1 + dt * beta * corqs)

        dtgdp = dt * RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
        dpr1 = covpclr * b / dtgdp
        dpr = min(dpr1, preclr1)

        # take away from clear sky flux
        preclr = preclr1 - dpr
        if preclr <= 0:
            covptotn = out_clc
        out_covptot = covptotn

        # warm proportion
        evapr = dpr * rfln2 / prtot
        rfln -= evapr

        # ice proportion
        evaps = dpr * sfln2 / prtot
        sfln -= evaps
    else:
        b = 0.0
        beta = 0.0
        dpr = 0.0
        dpr1 = 0.0
        dtgdp = 0.0
        evapr = 0.0
        evaps = 0.0
        out_covptot = 0.0
        preclr = 0.0
        preclr1 = 0.0
        qe = 0.0

    # incrementation of t and q, and fluxes swap
    dqdt = -(condl1 + condi1) + (in_lude[0, 0, 0] + evapr + evaps) * gdp
    dtdt = (
        lvdcp * condl1
        + lsdcp * condi1
        - (
            lvdcp * evapr
            + lsdcp * evaps
            + in_lude[0, 0, 0] * (fwat * lvdcp + (1 - fwat) * lsdcp)
            - (lsdcp - lvdcp) * rfreeze1
        )
        * gdp
    )

    # first guess T and Q
    t3 = t + dt * dtdt
    q += dt * dqdt
    told = t3 + 0
    qold = q + 0
    qold1 = q + 0

    # clipping of final qv
    t3, q = cuadjtqs_nl(in_ap, t3, q)

    # update rain fraction and freezing
    dq = max(qold1 - q, 0.0)
    dr2 = cons2 * dp * dq
    if t3 < RTT:
        rfreeze2 = fwat * dr2
        fwatr2 = 0.0
    else:
        rfreeze2 = 0.0
        fwatr2 = 1.0
    condl2 = condl1 + fwatr2 * dq / dt
    condi2 = condi1 + (1 - fwatr2) * dq / dt
    rfln += fwatr2 * dr2
    sfln += (1 - fwatr2) * dr2
    rfreeze3 = rfreeze1 + rfreeze2

    # calculate output tendencies
    out_tnd_q = -(condl2 + condi2) + (in_lude[0, 0, 0] + evapr + evaps) * gdp
    out_tnd_t = (
        lvdcp * condl2
        + lsdcp * condi2
        - (
            lvdcp * evapr
            + lsdcp * evaps
            + in_lude[0, 0, 0] * (fwat * lvdcp + (1 - fwat) * lsdcp)
            - (lsdcp - lvdcp) * rfreeze3
        )
        * gdp
    )
    out_tnd_ql = (qlwc - ql) / dt
    out_tnd_qi = (qiwc - qi) / dt

    return (
        out_tnd_t,
        out_tnd_q,
        out_tnd_ql,
        out_tnd_qi,
        out_clc,
        out_covptot,
        b,
        beta,
        clc,
        cldi,
        cldl,
        condi1,
        condi2,
        condl1,
        condl2,
        cons,
        cor,
        corqs,
        covpclr,
        covpclr1,
        covptot1,
        covptotn,
        crh2,
        dp,
        dpr,
        dpr1,
        dq,
        dqc,
        dqsdtemp,
        dqsdz,
        dr1,
        dr2,
        dtdzmo,
        dtgdp,
        esdp1,
        evapr,
        evaps,
        fac,
        fac1,
        fac2,
        fac3,
        faci,
        facw,
        foeew,
        fwat,
        fwatr1,
        fwatr2,
        gdp,
        itmp11,
        itmp12,
        itmp2,
        ldcp,
        lfdcp,
        lo3,
        lsdcp,
        ltmp1,
        ltmp2,
        lude,
        lvdcp,
        preclr,
        preclr1,
        prr,
        prs,
        prtot,
        q,
        q2,
        qc3,
        qcd,
        qcrit,
        qe,
        qiwc1,
        qlim,
        qlwc1,
        qold,
        qold1,
        qpd,
        qsat,
        qt,
        rho,
        rfln,
        rfln2,
        rfreeze1,
        rfreeze3,
        rodqsdp,
        scalm,
        sfln,
        sfln2,
        snmlt,
        supsat,
        t,
        t3,
        tmp3,
        told,
        z2s,
    )
