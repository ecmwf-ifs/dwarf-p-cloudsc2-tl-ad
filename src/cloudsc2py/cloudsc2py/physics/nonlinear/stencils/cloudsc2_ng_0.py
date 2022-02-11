# -*- coding: utf-8 -*-
from gt4py import gtscript

from cloudsc2py.framework.stencil import function_collection


@function_collection(
    "cloudsc2_nl_0",
    external_names=[
        "ICALL",
        "LDRAIN1D",
        "LEVAPLS2",
        "LPHYLIN",
        "R2ES",
        "R3IES",
        "R3LES",
        "R4IES",
        "R4LES",
        "R5ALSCP",
        "R5ALVCP",
        "R5IES",
        "R5LES",
        "RALSDCP",
        "RALVDCP",
        "RCLCRIT",
        "RCPD",
        "RD",
        "RETV",
        "RG",
        "RKCONV",
        "RLMIN",
        "RLMLT",
        "RLPTRC",
        "RLSTT",
        "RLVTT",
        "RPECONS",
        "RTICE",
        "RTICECU",
        "RTT",
        "RTWAT",
        "RTWAT_RTICE_R",
        "RTWAT_RTICECU_R",
        "RVTMP2",
        "ZEPS1",
        "ZEPS2",
        "ZQMAX",
        "ZSCAL",
        "cuadjtqs_nl",
    ],
)
@gtscript.function
def cloudsc2_nl_0(
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
    tmp_trpaus,
    covptot,
    rfl,
    sfl,
    t,
    dt,
):
    from __externals__ import ext

    # first guess values for q, ql and qi
    q = in_q + dt * in_tnd_cml_q + in_supsat
    ql = in_ql + dt * in_tnd_cml_ql
    qi = in_qi + dt * in_tnd_cml_qi

    # set up constants required
    ckcodtl = 2 * ext.RKCONV * dt
    ckcodti = 5 * ext.RKCONV * dt
    cons2 = 1 / (ext.RG * dt)
    cons3 = ext.RLVTT / ext.RCPD
    meltp2 = ext.RTT + 2

    # parameter for cloud formation
    scalm = ext.ZSCAL * max(in_eta - 0.2, ext.ZEPS1) ** 0.2

    # thermodynamic constants
    dp = in_aph[0, 0, 1] - in_aph[0, 0, 0]
    zz = ext.RCPD + ext.RCPD * ext.RVTMP2 * q
    lfdcp = ext.RLMLT / zz
    lsdcp = ext.RLSTT / zz
    lvdcp = ext.RLVTT / zz

    # calculate dqs/dT correction factor
    if __INLINED(ext.LPHYLIN or ext.LDRAIN1D):
        if t < ext.RTT:
            fwat = 0.545 * (tanh(0.17 * (t - ext.RLPTRC)) + 1)
            z3es = ext.R3IES
            z4es = ext.R4IES
        else:
            fwat = 1.0
            z3es = ext.R3LES
            z4es = ext.R4LES
        foeew = ext.R2ES * exp(z3es * (t - ext.RTT) / (t - z4es))
        esdp = min(foeew / in_ap, ext.ZQMAX)
    else:
        fwat = ext.foealfa(t)
        foeew = ext.foeewm(t)
        esdp = foeew / in_ap
    facw = ext.R5LES / ((t - ext.R4LES) ** 2)
    faci = ext.R5IES / ((t - ext.R4IES) ** 2)
    fac = fwat * facw + (1 - fwat) * faci
    dqsdtemp = fac * in_qsat / (1 - ext.RETV * esdp)
    corqs = 1 + cons3 * dqsdtemp

    # use clipped state
    qlim = min(q, in_qsat)

    # set up critical value of humidity
    rh1 = 1.0
    rh2 = (
        0.35
        + 0.14 * ((tmp_trpaus - 0.25) / 0.15) ** 2
        + 0.04 * min(tmp_trpaus - 0.25, 0.0) / 0.15
    )
    rh3 = 1.0
    if in_eta < tmp_trpaus:
        crh2 = rh3
    else:
        deta2 = 0.3
        bound1 = tmp_trpaus + deta2
        if in_eta < bound1:
            crh2 = rh3 + (rh2 - rh3) * (in_eta - tmp_trpaus) / deta2
        else:
            deta1 = 0.09 + 0.16 * (0.4 - tmp_trpaus) / 0.3
            bound2 = 1 - deta1
            if in_eta < bound2:
                crh2 = rh2
            else:
                crh2 = rh1 + (rh2 - rh1) * sqrt((1 - in_eta) / deta1)

    # allow ice supersaturation at cold temperatures
    if t < ext.RTICE:
        qsat = in_qsat * (1.8 - 0.003 * t)
    else:
        qsat = in_qsat
    qcrit = crh2 * qsat

    # simple uniform distribution of total water from Leutreut & Li (1990)
    qt = q + ql + qi
    if qt < qcrit:
        out_clc = 0.0
        qc = 0.0
    elif qt >= qsat:
        out_clc = 1.0
        qc = (1 - scalm) * (qsat - qcrit)
    else:
        qpd = qsat - qt
        qcd = qsat - qcrit
        out_clc = 1 - sqrt(qpd / (qcd - scalm * (qt - qcrit)))
        qc = (scalm * qpd + (1 - scalm) * qcd) * (out_clc ** 2)

    # add convective component
    gdp = ext.RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
    lude = dt * in_lude * gdp
    lo1 = lude[0, 0, 0] >= ext.RLMIN and in_lu[0, 0, 1] >= ext.ZEPS2
    if lo1:
        out_clc += (1 - out_clc[0, 0, 0]) * (
            1 - exp(-lude[0, 0, 0] / in_lu[0, 0, 1])
        )
        qc += lude

    # add compensating subsidence component
    rho = in_ap / (ext.RD * t)
    rodqsdp = -rho * in_qsat / (in_ap - ext.RETV * foeew)
    ldcp = fwat * lvdcp + (1 - fwat) * lsdcp
    dtdzmo = ext.RG * (1 / ext.RCPD - ldcp * rodqsdp) / (1 + ldcp * dqsdtemp)
    dqsdz = dqsdtemp * dtdzmo - ext.RG * rodqsdp
    dqc = min(dt * dqsdz * (in_mfu + in_mfd) / rho, qc)
    qc -= dqc

    # new cloud liquid/ice contents and condensation rates (liquid/ice)
    qlwc = qc * fwat
    qiwc = qc * (1 - fwat)
    condl = (qlwc - ql) / dt
    condi = (qiwc - qi) / dt

    # calculate precipitation overlap
    # simple form based on Maximum Overlap
    covptotn = max(covptot, out_clc)
    covpclr = max(covptotn - out_clc, 0.0)

    # melting of incoming snow
    if sfl != 0:
        cons = cons2 * dp / lfdcp
        snmlt = min(sfl, cons * max(t - meltp2, 0.0))
        rfln = rfl + snmlt
        sfln = sfl - snmlt
        t -= snmlt / cons
    else:
        rfln = rfl
        sfln = sfl

    # diagnostic calculation of rain production from cloud liquid water
    if out_clc > ext.ZEPS2:
        if __INLINED(ext.LEVAPLS2 or ext.LDRAIN1D):
            lcrit = 1.9 * ext.RCLCRIT
        else:
            lcrit = 2.0 * ext.RCLCRIT
        cldl = qlwc / out_clc
        dl = ckcodtl * (1 - exp(-((cldl / lcrit) ** 2)))
        prr = qlwc - out_clc * cldl * exp(-dl)
        qlwc -= prr
    else:
        prr = 0.0

    # diagnostic calculation of snow production from cloud ice
    if out_clc > ext.ZEPS2:
        if __INLINED(ext.LEVAPLS2 or ext.LDRAIN1D):
            icrit = 0.0001
        else:
            icrit = 2 * ext.RCLCRIT
        cldi = qiwc / out_clc
        di = (
            ckcodti
            * exp(0.025 * (t - ext.RTT))
            * (1 - exp(-((cldi / icrit) ** 2)))
        )
        prs = qiwc - out_clc * cldi * exp(-di)
        qiwc -= prs
    else:
        prs = 0.0

    # new precipitation (rain + snow)
    dr = cons2 * dp * (prr + prs)

    # rain fraction (different from cloud liquid water fraction!)
    if t < ext.RTT:
        rfreeze = cons2 * dp * prr
        fwatr = 0.0
    else:
        rfreeze = 0.0
        fwatr = 1.0
    rfln += fwatr * dr
    sfln += (1 - fwatr) * dr

    # precipitation evaporation
    prtot = rfln + sfln
    if (
        prtot > ext.ZEPS2
        and covpclr > ext.ZEPS2
        and (ext.LEVAPLS2 or ext.LDRAIN1D)
    ):
        preclr = prtot * covpclr / covptotn

        # this is the humidity in the moisest zcovpclr region
        qe = in_qsat - (in_qsat - qlim) * covpclr / ((1 - out_clc) ** 2)
        beta = (
            ext.RG
            * ext.RPECONS
            * (
                sqrt(in_ap[0, 0, 0] / in_aph[0, 0, 1])
                / 0.00509
                * preclr
                / covpclr
            )
            ** 0.5777
        )

        # implicit solution
        b = dt * beta * (in_qsat - qe) / (1 + dt * beta * corqs)

        dtgdp = dt * ext.RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
        dpr = min(covpclr * b / dtgdp, preclr)
        preclr -= dpr
        if preclr <= 0:
            covptotn = out_clc
        out_covptot = covptotn

        # warm proportion
        evapr = dpr * rfln / prtot
        rfln -= evapr

        # ice proportion
        evaps = dpr * sfln / prtot
        sfln -= evaps
    else:
        out_covptot = 0.0
        evapr = 0.0
        evaps = 0.0

    # update of T and Q tendencies due to:
    # - condensation/evaporation of cloud water/ice
    # - detrainment of convective cloud condensate
    # - evaporation of precipitation
    # - freezing of rain (impact on T only)
    dqdt = -(condl + condi) + (in_lude + evapr + evaps) * gdp
    dtdt = (
        lvdcp * condl
        + lsdcp * condi
        - (
            lvdcp * evapr
            + lsdcp * evaps
            + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
            - (lsdcp - lvdcp) * rfreeze
        )
        * gdp
    )

    # first guess T and Q
    t += dt * dtdt
    q += dt * dqdt
    qold = q

    # clipping of final qv
    t, q = ext.cuadjtqs_nl(in_ap, t, q)

    # update rain fraction and freezing
    dq = max(qold - q, 0.0)
    dr2 = cons2 * dp * dq
    if t < ext.RTT:
        rfreeze2 = fwat * dr2
        fwatr = 0.0
    else:
        rfreeze2 = 0.0
        fwatr = 1.0
    rn = fwatr * dr2
    sn = (1 - fwatr) * dr2
    condl += fwatr * dq / dt
    condi += (1 - fwatr) * dq / dt
    rfln += rn
    sfln += sn
    rfreeze += rfreeze2

    # calculate output tendencies
    out_tnd_q = -(condl + condi) + (in_lude + evapr + evaps) * gdp
    out_tnd_t = (
        lvdcp * condl
        + lsdcp * condi
        - (
            lvdcp * evapr
            + lsdcp * evaps
            + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
            - (lsdcp - lvdcp) * rfreeze
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
        covptotn,
        rfln,
        sfln,
    )
