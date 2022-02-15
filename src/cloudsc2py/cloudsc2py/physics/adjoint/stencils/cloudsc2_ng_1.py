# -*- coding: utf-8 -*-
from gt4py import gtscript

from cloudsc2py.framework.stencil import function_collection
from cloudsc2py.utils.f2py import ported_function


@ported_function(
    from_file="cloudsc2_ad/cloudsc2ad.F90", from_line=339, to_line=866
)
@function_collection(
    "cloudsc2_ad_1",
    external_names=[
        "LDRAIN1D",
        "LEVAPLS2",
        "LREGCL",
        "R3IES",
        "R3LES",
        "R4IES",
        "R4LES",
        "R5IES",
        "RCLCRIT",
        "RETV",
        "RG",
        "RLMIN",
        "RLPTRC",
        "RLSTT",
        "RLVTT",
        "RPECONS",
        "RTICE",
        "RTT",
        "ZEPS2",
        "ZQMAX",
        "cuadjtqs_ad",
    ],
)
@gtscript.function
def cloudsc2_ad_1(
    in_ap,
    in_aph,
    in_qsat,
    in_lu,
    in_lude,
    in_mfd,
    in_mfu,
    in_tnd_t_i,
    in_tnd_q_i,
    in_tnd_ql_i,
    in_tnd_qi_i,
    in_clc_i,
    in_covptot_i,
    in_fhpsl_i,
    in_fhpsn_i,
    in_fplsl_i,
    in_fplsn_i,
    out_clc,
    dt,
    aph_s,
    aph_i_s,
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
    covptot,
    covptot_i,
    covptot1,
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
    fplsl_i,
    fplsn_i,
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
    on_last_level,
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
    rfl_i,
    rfln2,
    rfreeze1,
    rfreeze3,
    rho,
    rodqsdp,
    scalm,
    sfl,
    sfl_i,
    sfln2,
    snmlt,
    supsat,
    t,
    t2,
    t3,
    tmp3,
    told,
    z2s,
):
    from __externals__ import ext

    # set up constants required
    ckcodtl = 2 * ext.RKCONV * dt
    ckcodtla = ckcodtl / 100
    ckcodti = 5 * ext.RKCONV * dt
    ckcodtia = ckcodti / 100
    cons2 = 1 / (ext.RG * dt)
    cons3 = ext.RLVTT / ext.RCPD
    meltp2 = ext.RTT + 2

    # incrementation of T and q, and fluxes swap
    rfln_i = rfl_i + fplsl_i[0, 0, 1]
    sfln_i = sfl_i + fplsn_i[0, 0, 1]

    # qice tendency
    out_qi_i = -in_tnd_qi_i / dt
    qiwc_i = in_tnd_qi_i / dt

    # qliq tendency
    out_ql_i = -in_tnd_ql_i / dt
    qlwc_i = in_tnd_ql_i / dt

    # T and q tendency
    gdp_i = -in_tnd_t_i * (
        lvdcp * evapr
        + lsdcp * evaps
        + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
        - (lsdcp - lvdcp) * rfreeze3
    ) + in_tnd_q_i * (in_lude + evapr + evaps)
    condl_i = in_tnd_t_i * lvdcp - in_tnd_q_i
    condi_i = in_tnd_t_i * lsdcp - in_tnd_q_i
    evapr_i = -in_tnd_t_i * lvdcp * gdp + in_tnd_q_i * gdp
    evaps_i = -in_tnd_t_i * lsdcp * gdp + in_tnd_q_i * gdp
    lvdcp_i = in_tnd_t_i * (condl2 - gdp * (evapr + in_lude * fwat + rfreeze3))
    lsdcp_i = in_tnd_t_i * (
        condi2 - gdp * (evaps + in_lude * (1 - fwat) - rfreeze3)
    )
    out_lude_i = gdp * (
        -in_tnd_t_i * (fwat * lvdcp + (1 - fwat) * lsdcp) + in_tnd_q_i
    )
    fwat_i = -in_tnd_t_i * in_lude * gdp * (lvdcp - lsdcp)
    rfreeze_i = in_tnd_t_i * (lsdcp - lvdcp) * gdp

    # clipping of final qv
    rn_i = rfln_i + 0
    sn_i = sfln_i + 0

    # note: the extra condensation due to the adjustment goes directly
    # to precipitation
    dq_i = (fwatr2 * condl_i + (1 - fwatr2) * condi_i) / dt
    # fwatr_i = (condl_i - condi_i) * dq / dt + dr2 * (rn_i - sn_i)
    dr2_i = fwatr2 * rn_i + (1 - fwatr2) * sn_i

    # update rain fraction adn freezing
    # note: impact of new temperature out_t_i on fwat_i is neglected here
    if t3 < ext.RTT:
        fwat_i += dr2 * rfreeze_i
        dr2_i += fwat * rfreeze_i
    fwatr_i = 0.0

    dq_i += cons2 * dp * dr2_i
    dp_i = cons2 * dq * dr2_i

    if qold1 >= q:
        # regularization
        if __INLINED(ext.LREGCL):
            dq_i *= 0.7
        qold_i = dq_i
        out_q_i = -dq_i
    else:
        qold_i = 0.0
        out_q_i = 0.0

    out_ap_i = 0.0
    out_t_i = 0.0
    out_ap_i, told, out_t_i, qold, out_q_i = ext.cuadjtqs_ad(
        in_ap, out_ap_i, told, out_t_i, qold, out_q_i
    )

    # first guess T and q
    out_q_i += qold_i
    dtdt_i = dt * out_t_i
    dqdt_i = dt * out_q_i

    # incrementation of T and q
    gdp_i += -dtdt_i * (
        lvdcp * evapr
        + lsdcp * evaps
        + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
        - (lsdcp - lvdcp) * rfreeze1
    ) + dqdt_i * (in_lude + evapr + evaps)
    condl_i += dtdt_i * lvdcp - dqdt_i
    condi_i += dtdt_i * lsdcp - dqdt_i
    evapr_i += -dtdt_i * lvdcp * gdp + dqdt_i * gdp
    evaps_i += -dtdt_i * lsdcp * gdp + dqdt_i * gdp
    lvdcp_i += dtdt_i * (condl1 - gdp * (evapr + fwat * in_lude + rfreeze1))
    lsdcp_i += dtdt_i * (
        condi1 - gdp * (evaps + (1 - fwat) * in_lude - rfreeze1)
    )
    out_lude_i += (
        -dtdt_i * gdp * (fwat * lvdcp + (1 - fwat) * lsdcp) + dqdt_i * gdp
    )
    fwat_i -= dtdt_i * in_lude * gdp * (lvdcp - lsdcp)
    rfreeze_i += dtdt_i * (lsdcp - lvdcp) * gdp

    if (
        prtot > ext.ZEPS2
        and covpclr > ext.ZEPS2
        and (ext.LEVAPLS2 or ext.LDRAIN1D)
    ):
        # ice proportion
        evaps_i -= sfln_i
        sfln_i += dpr * evaps_i / prtot
        dpr_i = sfln2 * evaps_i / prtot
        prtot_i = dpr * sfln2 * evaps_i / prtot ** 2

        # warm proportion
        evapr_i -= rfln_i
        rfln_i += dpr * evapr_i / prtot
        dpr_i += rfln2 * evapr_i / prtot
        prtot_i -= dpr * rfln2 * evapr_i / prtot ** 2

        # take away clear sky flux
        covptotn_i = covptot_i + in_covptot_i
        if preclr <= 0:
            clc_i = in_clc_i + covptotn_i
            covptotn_i = 0.0
        else:
            clc_i = in_clc_i + 0

        if dpr1 > preclr1:
            preclr_i = dpr_i
            dpr_i = 0.0
        else:
            preclr_i = 0.0

        b_i = covpclr * dpr_i / dtgdp
        covpclr_i = b * dpr_i / dtgdp
        dtgdp_i = -covpclr * b * dpr_i / dtgdp ** 2
        out_aph_i = dt * ext.RG * dtgdp_i / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
        daph_i = dt * ext.RG * dtgdp_i / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
        if on_last_level:
            aphn_i_s = aph_i_s - daph_i
        else:
            aphn_i_s = aph_i_s + 0

        # implicit solution
        tmp1 = 1 + dt * beta * corqs
        beta_i = (
            dt * (in_qsat - qe) * b_i / tmp1
            - (dt ** 2) * beta * (in_qsat - qe) * corqs * b_i / tmp1 ** 2
        )
        out_qsat_i = dt * beta * b_i / tmp1
        qe_i = -dt * beta * b_i / tmp1
        corqs_i = -(dt ** 2) * beta * (in_qsat - qe) * beta * b_i / tmp1 ** 2

        # this is the humidity in the moistest covpclr region
        xx = (
            0.5777
            * (ext.RG * ext.RPECONS / 0.00509)
            * (0.00509 * covpclr / (preclr1 * sqrt(in_ap / aph_s))) ** 0.4223
        )
        preclr_i += xx * sqrt(in_ap / aph_s) * beta_i / covpclr
        out_ap_i += (
            0.5 * xx * preclr1 * beta_i / (covpclr * sqrt(in_ap * aph_s))
        )
        aphn_i_s -= (
            0.5
            * xx
            * preclr1
            * sqrt(in_ap / aph_s)
            * beta_i
            / (covpclr * aph_s)
        )
        covpclr_i -= (
            (xx * preclr1 * sqrt(in_ap / aph_s) * beta_i / covpclr ** 2)
            + (in_qsat - qlim) * qe_i / (1 - out_clc) ** 2
        ) - prtot * preclr_i / covptot1
        out_qsat_i += qe_i - covpclr * qe_i / (1 - out_clc) ** 2
        qlim_i = covpclr * qe_i / (1 - out_clc) ** 2
        clc_i -= 2 * (in_qsat - qlim) * covpclr * qe_i / (1 - out_clc) ** 3
        prtot_i += covpclr * preclr_i / covptot1
        covptotn_i -= prtot * covpclr * preclr_i / covptot1 ** 2
    else:
        clc_i = in_clc_i + 0
        corqs_i = 0.0
        covpclr_i = 0.0
        covptotn_i = covptot_i + 0
        daph_i = 0.0
        out_aph_i = 0.0
        aphn_i_s = aph_i_s + 0
        out_qsat_i = 0.0
        prtot_i = 0.0
        qlim_i = 0.0

    # new precipitation
    rfln_i += prtot_i
    sfln_i += prtot_i
    fwatr_i += dr1 * (rfln_i - sfln_i)
    dr_i = fwatr1 * rfln_i + (1 - fwatr1) * sfln_i

    # update rain fraction and freezing
    # note: impact of new temperature out_t_i on fwat_i in neglected here
    if t < ext.RTT:
        dp_i += rfreeze_i * cons2 * prr
        prr_i = rfreeze_i * cons2 * dp
    else:
        prr_i = 0.0
    prr_i += cons2 * dp * dr_i
    prs_i = cons2 * dp * dr_i
    dp_i += cons2 * (prr + prs) * dr_i

    if out_clc > ext.ZEPS2:
        # diagnostic calculation of rain production from cloud ice
        if __INLINED(ext.LEVAPLS2 or ext.LDRAIN1D):
            icrit = 0.0001
            lcrit = 1.9 * ext.RCLCRIT
        else:
            icrit = 2 * ext.RCLCRIT
            lcrit = 2 * ext.RCLCRIT
        prs_i -= qiwc_i
        qiwc_i += prs_i
        qinew_i = -prs_i
        clc_i += qinew_i * cldi * itmp2
        cldi_i = qinew_i * out_clc * itmp2
        di_i = -qinew_i * out_clc * cldi * itmp2

        # regularization
        if __INLINED(ext.LREGCL):
            itmp4 = ckcodtia
        else:
            itmp4 = ckcodti
        out_t_i += 0.025 * itmp4 * itmp12 * (1 - itmp11) * di_i
        cldi_i += 2 * ckcodti * itmp12 * itmp11 * cldi * di_i / icrit ** 2

        qiwc_i += cldi_i / out_clc
        clc_i -= qiwc1 * cldi_i / out_clc ** 2

        # diagnostic calculation of rain production from cloud liquid water
        prr_i -= qlwc_i
        qlwc_i += prr_i
        qlnew_i = -prr_i
        clc_i += qlnew_i * cldl * ltmp2
        cldl_i = qlnew_i * out_clc * ltmp2
        dl_i = -qlnew_i * out_clc * cldl * ltmp2

        # regularization
        if __INLINED(ext.LREGCL):
            ltmp4 = ckcodtla
        else:
            ltmp4 = ckcodtl
        cldl_i += 2 * ltmp4 * ltmp1 * cldl * dl_i / lcrit ** 2

        qlwc_i += cldl_i / out_clc
        clc_i -= qlwc1 * cldl_i / out_clc ** 2

    if sfl != 0.0:
        snmlt_i = -out_t_i / cons + rfln_i - sfln_i
        cons_i = out_t_i * snmlt / cons ** 2

        if sfl <= z2s:
            sfln_i += snmlt_i
            z2s_i = 0.0
        else:
            z2s_i = snmlt_i + 0

        if t2 > meltp2:
            out_t_i += cons * z2s_i
            cons_i += (t2 - meltp2) * z2s_i

        dp_i += cons_i * cons2 / lfdcp
        lfdcp_i = -cons2 * dp * cons_i / lfdcp ** 2
    else:
        lfdcp_i = 0.0

    if covpclr1 < 0:
        covpclr_i = 0.0
    covptotn_i += covpclr_i
    clc_i -= covpclr_i

    if out_clc > covptot:
        clc_i += covptotn_i

    # new cloud liquid/ice contents and condensation rates (liquid/ice)
    qiwc_i += condi_i / dt
    out_qi_i -= condi_i / dt
    qlwc_i += condl_i / dt
    out_ql_i -= condl_i / dt
    qc_i = fwat * qlwc_i + (1 - fwat) * qiwc_i
    fwat_i += qc3 * (qlwc_i - qiwc_i)

    # add compensating subsidence component
    dqc_i = -qc_i
    if lo3:
        if __INLINED(ext.LREGCL):
            # regularization
            dqc_i *= 0.1
        dqsdz_i = dt * dqc_i * (in_mfd + in_mfu) / rho
        out_mfd_i = dt * dqc_i * dqsdz / rho
        out_mfu_i = dt * dqc_i * dqsdz / rho
        rho_i = -dqc_i * dqc / rho
    else:
        qc_i += dqc_i
        dqsdz_i = 0.0
        out_mfd_i = 0.0
        out_mfu_i = 0.0
        rho_i = 0.0

    dtdzmo_i = dqsdz_i * dqsdtemp
    dqsdtemp_i = dqsdz_i * dtdzmo - dtdzmo * dtdzmo_i * ldcp * fac3
    rodqsdp_i = -ext.RG * (dqsdz_i + dtdzmo_i * ldcp * fac3)
    ldcp_i = -dtdzmo_i * (ext.RG * rodqsdp + dtdzmo * dqsdtemp) * fac3
    fwat_i += ldcp_i * (lvdcp - lsdcp)
    lvdcp_i += fwat * ldcp_i
    lsdcp_i += (1 - fwat) * ldcp_i
    rho_i -= rodqsdp_i * in_qsat * fac2
    out_qsat_i -= rodqsdp_i * rho * in_qsat * fac2 ** 2
    out_ap_i += rodqsdp_i * rho * in_qsat * fac2 ** 2 + rho_i * fac1
    foeew_i = -ext.RETV * rodqsdp_i * rho * in_qsat * fac2 ** 2
    out_t_i -= rho_i * in_ap * fac1 / t2

    if lude >= ext.RLMIN and in_lu[0, 0, 1] >= ext.ZEPS2:
        lude_i = (
            qc_i[0, 0, 0]
            + (1 - clc[0, 0, 0])
            * exp(-lude[0, 0, 0] / in_lu[0, 0, 1])
            * clc_i[0, 0, 0]
            / in_lu[0, 0, 1]
        )
        dlu_i = (
            (1 - clc[0, 0, 0])
            * lude[0, 0, 0]
            * exp(-lude[0, 0, 0] / in_lu[0, 0, 1])
            * clc_i[0, 0, 0]
            / in_lu[0, 0, 1] ** 2
        )
        clc_i *= 1 - (1 - exp(-lude[0, 0, 0] / in_lu[0, 0, 1]))
    else:
        dlu_i = 0.0
        lude_i = 0.0

    out_lude_i += dt * gdp * lude_i
    gdp_i += dt * in_lude * lude_i
    out_aph_i += (
        ext.RG * gdp_i[0, 0, 0] / (in_aph[0, 0, 1] - in_aph[0, 0, 0]) ** 2
    )
    daph_i += (
        ext.RG * gdp_i[0, 0, 0] / (in_aph[0, 0, 1] - in_aph[0, 0, 0]) ** 2
    )

    if qt < qcrit:
        qsat_i = 0.0
        qcrit_i = 0.0
        qt_i = 0.0
    elif qt >= qsat:
        qsat_i = (1 - scalm) * qc_i
        qcrit_i = -(1 - scalm) * qc_i
        qt_i = 0.0
    else:
        qpd_i = scalm * qc_i * clc ** 2
        qcd_i = (1 - scalm) * qc_i * clc ** 2
        clc_i += 2 * (scalm * qpd + (1 - scalm) * qcd) * clc * qc_i

        if __INLINED(ext.LREGCL):
            # regularization of cloud fraction
            rat = qpd / qcd
            yyy = min(
                0.3,
                3.5 * sqrt(rat * (1 - scalm * (1 - rat)) ** 3) / (1 - scalm),
            )
            clc_i *= yyy

        qpd_i -= 0.5 / tmp3 * clc_i / (qcd - scalm * (qt - qcrit))
        qcd_i += 0.5 / tmp3 * qpd * clc_i / (qcd - scalm * (qt - qcrit)) ** 2
        qt_i = (
            -0.5
            / tmp3
            * (qpd * scalm * clc_i)
            / (qcd - scalm * (qt - qcrit)) ** 2
        ) - qpd_i
        qcrit_i = (
            0.5
            / tmp3
            * (qpd * scalm * clc_i)
            / (qcd - scalm * (qt - qcrit)) ** 2
        ) - qcd_i
        qsat_i = qcd_i + qpd_i

    out_q_i += qt_i
    out_ql_i += qt_i
    out_qi_i += qt_i

    # set up critical value of humidity
    qsat_i += qcrit_i * crh2
    out_qsat_i += qsat_i * supsat
    supsat_i = qsat_i * in_qsat

    # allow ice supersaturation at cold temperatures
    if t2 < ext.RTICE:
        out_t_i -= 0.003 * supsat_i

    # use clipped state
    if q2 > in_qsat:
        out_qsat_i += qlim_i
    else:
        out_q_i += qlim_i

    # calculate dqs/dT correction factor
    dqsdtemp_i += cons3 * corqs_i
    out_qsat_i += fac * cor * dqsdtemp_i
    cor_i = fac * in_qsat * dqsdtemp_i
    fac_i = cor * in_qsat * dqsdtemp_i
    esdp_i = ext.RETV * cor_i * cor ** 2
    facw_i = fwat * fac_i
    faci_i = (1 - fwat) * fac_i
    fwat_i += (facw - faci) * fac_i
    out_t_i -= 2 * (
        ext.R5IES * faci_i / (t2 - ext.R4IES) ** 3
        + ext.R5LES * facw_i / (t2 - ext.R4LES) ** 3
    )

    if esdp1 > ext.ZQMAX:
        esdp_i = 0.0
    foeew_i += esdp_i / in_ap
    out_ap_i -= esdp_i * foeew / in_ap ** 2

    if t2 < ext.RTT:
        z3es = ext.R3IES
        z4es = ext.R4IES
    else:
        z3es = ext.R3LES
        z4es = ext.R4LES
    out_t_i += z3es * (ext.RTT - z4es) * foeew_i * foeew / (t2 - z4es) ** 2

    if t2 < ext.RTT:
        out_t_i += 0.545 * 0.17 * fwat_i / cosh(0.17 * (t2 - ext.RLPTRC)) ** 2

    return (
        out_ap_i,
        out_aph_i,
        out_t_i,
        out_q_i,
        out_qsat_i,
        out_ql_i,
        out_qi_i,
        out_lude_i,
        out_mfd_i,
        out_mfu_i,
        aph_s,
        aphn_i_s,
        covptotn_i,
        daph_i,
        dlu_i,
        dp_i,
        lfdcp_i,
        lsdcp_i,
        lvdcp_i,
        rfln_i,
        sfln_i,
    )
