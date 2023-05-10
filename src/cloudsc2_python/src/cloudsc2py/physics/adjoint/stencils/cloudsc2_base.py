# -*- coding: utf-8 -*-
from gt4py.cartesian import gtscript

from cloudsc2py.framework.stencil import stencil_collection
from cloudsc2py.physics.adjoint.stencils.cuadjtqs import cuadjtqs_ad
from cloudsc2py.physics.nonlinear.stencils.cuadjtqs import cuadjtqs_nl
from src.cloudsc2_python.utils import ported_function


@ported_function(from_file="cloudsc2_ad/cloudsc2ad.F90", from_line=346, to_line=1740)
@stencil_collection("cloudsc2_ad")
def cloudsc2_ad_def(
    in_eta: gtscript.Field[gtscript.K, "float"],
    in_ap: gtscript.Field["float"],
    in_aph: gtscript.Field["float"],
    in_aph_s: gtscript.Field[gtscript.IJ, "float"],
    in_t: gtscript.Field["float"],
    in_q: gtscript.Field["float"],
    in_qsat: gtscript.Field["float"],
    in_ql: gtscript.Field["float"],
    in_qi: gtscript.Field["float"],
    in_lu: gtscript.Field["float"],
    in_lude: gtscript.Field["float"],
    in_mfd: gtscript.Field["float"],
    in_mfu: gtscript.Field["float"],
    in_supsat: gtscript.Field["float"],
    in_tnd_cml_t: gtscript.Field["float"],
    in_tnd_t_i: gtscript.Field["float"],
    in_tnd_cml_q: gtscript.Field["float"],
    in_tnd_q_i: gtscript.Field["float"],
    in_tnd_cml_ql: gtscript.Field["float"],
    in_tnd_ql_i: gtscript.Field["float"],
    in_tnd_cml_qi: gtscript.Field["float"],
    in_tnd_qi_i: gtscript.Field["float"],
    in_clc_i: gtscript.Field["float"],
    in_fhpsl_i: gtscript.Field["float"],
    in_fhpsn_i: gtscript.Field["float"],
    in_fplsl_i: gtscript.Field["float"],
    in_fplsn_i: gtscript.Field["float"],
    in_covptot_i: gtscript.Field["float"],
    tmp_aph_s: gtscript.Field[gtscript.IJ, "float"],
    tmp_rfl: gtscript.Field[gtscript.IJ, "float"],
    tmp_sfl: gtscript.Field[gtscript.IJ, "float"],
    tmp_trpaus: gtscript.Field[gtscript.IJ, "float"],
    out_ap_i: gtscript.Field["float"],
    out_aph_i: gtscript.Field["float"],
    out_t_i: gtscript.Field["float"],
    out_q_i: gtscript.Field["float"],
    out_qsat_i: gtscript.Field["float"],
    out_ql_i: gtscript.Field["float"],
    out_qi_i: gtscript.Field["float"],
    out_lu_i: gtscript.Field["float"],
    out_lude_i: gtscript.Field["float"],
    out_mfd_i: gtscript.Field["float"],
    out_mfu_i: gtscript.Field["float"],
    out_supsat_i: gtscript.Field["float"],
    out_tnd_t: gtscript.Field["float"],
    out_tnd_cml_t_i: gtscript.Field["float"],
    out_tnd_q: gtscript.Field["float"],
    out_tnd_cml_q_i: gtscript.Field["float"],
    out_tnd_ql: gtscript.Field["float"],
    out_tnd_cml_ql_i: gtscript.Field["float"],
    out_tnd_qi: gtscript.Field["float"],
    out_tnd_cml_qi_i: gtscript.Field["float"],
    out_clc: gtscript.Field["float"],
    out_fhpsl: gtscript.Field["float"],
    out_fhpsn: gtscript.Field["float"],
    out_fplsl: gtscript.Field["float"],
    out_fplsn: gtscript.Field["float"],
    out_covptot: gtscript.Field["float"],
    *,
    dt: "float",
):
    from __externals__ import (
        LDRAIN1D,
        LEVAPLS2,
        LREGCL,
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

    # set to zero precipitation fluxes at the top
    with computation(FORWARD), interval(0, 1):
        covptot = 0
        out_fplsl[0, 0, 0] = 0
        out_fplsn[0, 0, 0] = 0

    with computation(PARALLEL), interval(0, -1):
        # first guess values for T
        t = in_t + dt * in_tnd_cml_t
        # store trajectory arrays for adjoint
        t2 = t + 0
        t3 = t + 0

    # eta value at tropopause
    with computation(FORWARD), interval(0, 1):
        tmp_trpaus = 0.1
    with computation(FORWARD), interval(0, -2):
        if in_eta[0] > 0.1 and in_eta[0] < 0.4 and t[0, 0, 0] > t[0, 0, 1]:
            tmp_trpaus = in_eta[0]

    with computation(FORWARD), interval(0, -1):
        # first guess values for q, ql and qi
        q = in_q + dt * in_tnd_cml_q + in_supsat
        ql = in_ql + dt * in_tnd_cml_ql
        qi = in_qi + dt * in_tnd_cml_qi
        # store trajectory arrays for adjoint
        q2 = q + 0

        # set up constants required
        ckcodtl = 2 * RKCONV * dt
        ckcodti = 5 * RKCONV * dt
        cons2 = 1 / (RG * dt)
        cons3 = RLVTT / RCPD
        meltp2 = RTT + 2

        # parameter for cloud formation
        scalm = ZSCAL * max(in_eta - 0.2, ZEPS1) ** 0.2

        # thermodynamic constants
        dp = in_aph[0, 0, 1] - in_aph[0, 0, 0]
        zz = RCPD + RCPD * RVTMP2 * q
        lfdcp = RLMLT / zz
        lsdcp = RLSTT / zz
        lvdcp = RLVTT / zz

        # clear cloud and freezing arrays
        out_clc[0, 0, 0] = 0.0
        out_covptot[0, 0, 0] = 0.0

        # calculate dqs/dT correction factor
        if t < RTT:
            fwat = 0.545 * (tanh(0.17 * (t - RLPTRC)) + 1)
            z3es = R3IES
            z4es = R4IES
        else:
            fwat = 1.0
            z3es = R3LES
            z4es = R4LES
        foeew = R2ES * exp(z3es * (t - RTT) / (t - z4es))
        esdp1 = foeew / in_ap
        esdp = min(esdp1, ZQMAX)
        facw = R5LES / ((t - R4LES) ** 2)
        faci = R5IES / ((t - R4IES) ** 2)
        fac = fwat * facw + (1 - fwat) * faci
        dqsdtemp = fac * in_qsat / (1 - RETV * esdp)
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
        if t < RTICE:
            supsat = 1.8 - 0.003 * t
            qsat = in_qsat * supsat
        else:
            supsat = 1.0
            qsat = in_qsat + 0
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
            qc1 = (scalm * qpd + (1 - scalm) * qcd) * (clc ** 2)

        # add convective component
        gdp = RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
        lude = dt * in_lude * gdp
        lo1 = lude[0, 0, 0] >= RLMIN and in_lu[0, 0, 1] >= ZEPS2
        if lo1:
            out_clc[0, 0, 0] = clc[0, 0, 0] + (1 - clc[0, 0, 0]) * (
                1 - exp(-lude[0, 0, 0] / in_lu[0, 0, 1])
            )
            qc2 = qc1 + lude
        else:
            out_clc[0, 0, 0] = clc + 0
            qc2 = qc1 + 0

        # add compensating subsidence component
        fac1 = 1 / (RD * t)
        rho = in_ap * fac1
        fac2 = 1 / (in_ap - RETV * foeew)
        rodqsdp = -rho * in_qsat * fac2
        ldcp = fwat * lvdcp + (1 - fwat) * lsdcp
        fac3 = 1 / (1 + ldcp * dqsdtemp)
        dtdzmo = RG * (1 / RCPD - ldcp * rodqsdp) * fac3
        dqsdz = dqsdtemp * dtdzmo - RG * rodqsdp
        fac4 = 1 / rho
        lo3 = dt * dqsdz * (in_mfu + in_mfd) * fac4 < qc2
        dqc = min(dt * dqsdz * (in_mfu + in_mfd) * fac4, qc2)
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
        if tmp_sfl != 0:
            cons = cons2 * dp / lfdcp
            snmlt = min(tmp_sfl, cons * max(t - meltp2, 0.0))
            rfln = tmp_rfl + snmlt
            sfln = tmp_sfl - snmlt
            t -= snmlt / cons
        else:
            rfln = tmp_rfl
            sfln = tmp_sfl

        # diagnostic calculation of rain production from cloud liquid water
        if out_clc > ZEPS2:
            if __INLINED(LEVAPLS2 or LDRAIN1D):
                lcrit = 1.9 * RCLCRIT
            else:
                lcrit = 2.0 * RCLCRIT
            cldl = qlwc1 / out_clc
            ltmp1 = exp(-((cldl / lcrit) ** 2))
            dl = ckcodtl * (1 - ltmp1)
            ltmp2 = exp(-dl)
            qlnew = out_clc * cldl * ltmp2
            prr = qlwc1 - qlnew
            qlwc = qlwc1 - prr
        else:
            prr = 0.0
            qlwc = qlwc1 + 0

        # diagnostic calculation of snow production from cloud ice
        if out_clc > ZEPS2:
            if __INLINED(LEVAPLS2 or LDRAIN1D):
                icrit = 0.0001
            else:
                icrit = 2 * RCLCRIT
            cldi = qiwc1 / out_clc
            itmp11 = exp(-((cldi / icrit) ** 2))
            itmp12 = exp(0.025 * (t - RTT))
            di = ckcodti * itmp12 * (1 - itmp11)
            itmp2 = exp(-di)
            qinew = out_clc * cldi * itmp2
            prs = qiwc1 - qinew
            qiwc = qiwc1 - prs
        else:
            prs = 0.0
            qiwc = qiwc1 + 0

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

            # this is the humidity in the moistest zcovpclr region
            qe = in_qsat - (in_qsat - qlim) * covpclr / ((1 - out_clc) ** 2)
            beta = RG * RPECONS * (sqrt(in_ap / tmp_aph_s) / 0.00509 * preclr1 / covpclr) ** 0.5777

            # implicit solution
            b = dt * beta * (in_qsat - qe) / (1 + dt * beta * corqs)

            dtgdp = dt * RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
            dpr1 = covpclr * b / dtgdp
            dpr = min(dpr1, preclr1)

            # take away clear sky flux
            preclr = preclr1 - dpr
            if preclr <= 0:
                covptotn = out_clc
            out_covptot[0, 0, 0] = covptotn

            # warm proportion
            evapr = dpr * rfln2 / prtot
            rfln -= evapr

            # ice proportion
            evaps = dpr * sfln2 / prtot
            sfln -= evaps
        else:
            evapr = 0.0
            evaps = 0.0

        # update of T and Q tendencies due to:
        # - condensation/evaporation of cloud water/ice
        # - detrainment of convective cloud condensate
        # - evaporation of precipitation
        # - freezing of rain (impact on T only)
        dqdt = -(condl1 + condi1) + (in_lude + evapr + evaps) * gdp
        dtdt = (
            lvdcp * condl1
            + lsdcp * condi1
            - (
                lvdcp * evapr
                + lsdcp * evaps
                + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
                - (lsdcp - lvdcp) * rfreeze1
            )
            * gdp
        )

        # first guess T and Q
        t3 = t + dt * dtdt
        q = q2 + dt * dqdt
        told = t3
        qold = q
        qold1 = q

        # clipping of final qv
        t, q = cuadjtqs_nl(in_ap, t, q)

        # update rain fraction and freezing
        dq = max(qold1 - q, 0.0)
        dr2 = cons2 * dp * dq
        if t3 < RTT:
            rfreeze2 = fwat * dr2
            fwatr2 = 0.0
        else:
            rfreeze2 = 0.0
            fwatr2 = 1.0
        rn = fwatr2 * dr2
        sn = (1 - fwatr2) * dr2
        condl2 = condl1 + fwatr2 * dq / dt
        condi2 = condi1 + (1 - fwatr2) * dq / dt
        rfln += rn
        sfln += sn
        rfreeze3 = rfreeze1 + rfreeze2

        # calculate output tendencies
        out_tnd_q[0, 0, 0] = -(condl2 + condi2) + (in_lude + evapr + evaps) * gdp
        out_tnd_t[0, 0, 0] = (
            lvdcp * condl2
            + lsdcp * condi2
            - (
                lvdcp * evapr
                + lsdcp * evaps
                + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
                - (lsdcp - lvdcp) * rfreeze3
            )
            * gdp
        )
        out_tnd_ql[0, 0, 0] = (qlwc - ql) / dt
        out_tnd_qi[0, 0, 0] = (qiwc - qi) / dt

        # these fluxes will later be shifted one level downward
        fplsl = rfln
        fplsn = sfln

        # record rain flux for next level
        tmp_rfl[0, 0] = rfln
        tmp_sfl[0, 0] = sfln

        covptoto = covptot
        covptot = covptotn

    # enthalpy fluxes due to precipitation
    with computation(FORWARD):
        with interval(0, 1):
            covptot = covptoto
            out_fhpsl[0, 0, 0] = 0.0
            out_fhpsn[0, 0, 0] = 0.0
        with interval(1, None):
            covptot = covptotn[0, 0, -1]
            out_fplsl[0, 0, 0] = fplsl[0, 0, -1]
            out_fplsn[0, 0, 0] = fplsn[0, 0, -1]
            out_fhpsl[0, 0, 0] = -out_fplsl * RLVTT
            out_fhpsn[0, 0, 0] = -out_fplsn * RLSTT

    with computation(BACKWARD):
        with interval(-1, None):
            # enthalpy fluxes due to precipitation
            fplsn_i = in_fplsn_i - in_fhpsn_i * RLSTT
            fplsl_i = in_fplsl_i - in_fhpsl_i * RLVTT
        with interval(-2, -1):
            # set up constants required
            ckcodtl = 2 * RKCONV * dt
            ckcodtla = ckcodtl / 100
            ckcodti = 5 * RKCONV * dt
            ckcodtia = ckcodti / 100
            cons2 = 1 / (RG * dt)
            cons3 = RLVTT / RCPD
            meltp2 = RTT + 2

            # enthalpy fluxes due to precipitation
            fplsn_i = in_fplsn_i - in_fhpsn_i * RLSTT
            fplsl_i = in_fplsl_i - in_fhpsl_i * RLVTT

            # incrementation of T and q, and fluxes swap
            rfl_i = fplsl_i[0, 0, 1]
            sfl_i = fplsn_i[0, 0, 1]

            # qice tendency
            out_qi_i = -in_tnd_qi_i / dt
            qiwc_i = in_tnd_qi_i / dt

            # qliq tendency
            out_ql_i = -in_tnd_ql_i / dt
            qlwc_i = in_tnd_ql_i / dt

            # # T tendency
            # gdp_i = -in_tnd_t_i * (
            #     lvdcp * evapr
            #     + lsdcp * evaps
            #     + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
            #     - (lsdcp - lvdcp) * rfreeze3
            # )
            # condl_i = in_tnd_t_i * lvdcp
            # condi_i = in_tnd_t_i * lsdcp
            # evapr_i = -in_tnd_t_i * lvdcp * gdp
            # evaps_i = -in_tnd_t_i * lsdcp * gdp
            # lvdcp_i = in_tnd_t_i * (condl2 - evapr * gdp)
            # lsdcp_i = in_tnd_t_i * (condi2 - evaps * gdp)
            # out_lude_i = (
            #     -in_tnd_t_i * gdp * (fwat * lvdcp + (1 - fwat) * lsdcp)
            # )
            # lvdcp_i -= in_tnd_t_i * in_lude * gdp * fwat
            # lsdcp_i -= in_tnd_t_i * in_lude * gdp * (1 - fwat)
            # fwat_i = -in_tnd_t_i * in_lude * gdp * (lvdcp - lsdcp)
            # lvdcp_i -= in_tnd_t_i * rfreeze3 * gdp
            # lsdcp_i += in_tnd_t_i * rfreeze3 * gdp
            # rfreeze_i = in_tnd_t_i * (lsdcp - lvdcp) * gdp
            #
            # # q tendency
            # gdp_i += in_tnd_q_i * (in_lude + evapr + evaps)
            # out_lude_i += in_tnd_q_i * gdp
            # evapr_i += in_tnd_q_i * gdp
            # evaps_i += in_tnd_q_i * gdp
            # condl_i -= in_tnd_q_i
            # condi_i -= in_tnd_q_i

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
            lsdcp_i = in_tnd_t_i * (condi2 - gdp * (evaps + in_lude * (1 - fwat) - rfreeze3))
            out_lude_i = gdp * (-in_tnd_t_i * (fwat * lvdcp + (1 - fwat) * lsdcp) + in_tnd_q_i)
            fwat_i = -in_tnd_t_i * in_lude * gdp * (lvdcp - lsdcp)
            rfreeze_i = in_tnd_t_i * (lsdcp - lvdcp) * gdp

            # clipping of final qv
            rn_i = rfl_i + 0
            sn_i = sfl_i + 0

            # note: the extra condensation due to the adjustment goes directly
            # to precipitation
            dq_i = (fwatr2 * condl_i + (1 - fwatr2) * condi_i) / dt
            # fwatr_i = (condl_i - condi_i) * dq / dt + dr2 * (rn_i - sn_i)
            dr2_i = fwatr2 * rn_i + (1 - fwatr2) * sn_i

            # update rain fraction adn freezing
            # note: impact of new temperature out_t_i on fwat_i is neglected here
            if t3 < RTT:
                fwat_i += dr2 * rfreeze_i
                dr2_i += fwat * rfreeze_i
            fwatr_i = 0.0

            dq_i += cons2 * dp * dr2_i
            dp_i = cons2 * dq * dr2_i

            if qold1 >= q:
                # regularization
                if __INLINED(LREGCL):
                    dq_i *= 0.7
                qold_i = dq_i
                out_q_i = -dq_i
            else:
                qold_i = 0.0
                out_q_i = 0.0

            out_ap_i = 0.0
            out_t_i = 0.0
            out_ap_i, told, out_t_i, qold, out_q_i = cuadjtqs_ad(
                in_ap, out_ap_i, told, out_t_i, qold, out_q_i
            )

            # first guess T and q
            out_q_i += qold_i
            dtdt_i = dt * out_t_i
            dqdt_i = dt * out_q_i

            # incrementation of T and q
            # T tendency
            gdp_i -= dtdt_i * (
                lvdcp * evapr
                + lsdcp * evaps
                + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
                + (lsdcp - lvdcp) * rfreeze1
            )
            condl_i += dtdt_i * lvdcp
            condi_i += dtdt_i * lsdcp
            evapr_i -= dtdt_i * lvdcp * gdp
            evaps_i -= dtdt_i * lsdcp * gdp
            lvdcp_i += dtdt_i * (condl1 - evapr * gdp)
            lsdcp_i += dtdt_i * (condi1 - evaps * gdp)
            out_lude_i -= dtdt_i * gdp * (fwat * lvdcp + (1 - fwat) * lsdcp)
            lvdcp_i -= dtdt_i * in_lude * gdp * fwat
            lsdcp_i -= dtdt_i * in_lude * gdp * (1 - fwat)
            fwat_i -= dtdt_i * in_lude * gdp * (lvdcp - lsdcp)
            lvdcp_i -= dtdt_i * rfreeze1 * gdp
            lsdcp_i += dtdt_i * rfreeze1 * gdp
            rfreeze_i += dtdt_i * (lsdcp - lvdcp) * gdp

            # q tendency
            gdp_i += dqdt_i * (in_lude + evapr + evaps)
            out_lude_i += dqdt_i * gdp
            evapr_i += dqdt_i * gdp
            evaps_i += dqdt_i * gdp
            condl_i -= dqdt_i
            condi_i -= dqdt_i

            if prtot > ZEPS2 and covpclr > ZEPS2 and (LEVAPLS2 or LDRAIN1D):
                # ice proportion
                evaps_i -= sfl_i
                sfl_i += dpr * evaps_i / prtot
                dpr_i = sfln2 * evaps_i / prtot
                prtot_i = dpr * sfln2 * evaps_i / prtot ** 2

                # warm proportion
                evapr_i -= rfl_i
                rfl_i += dpr * evapr_i / prtot
                dpr_i += rfln2 * evapr_i / prtot
                prtot_i -= dpr * rfln2 * evapr_i / prtot ** 2

                # take away clear sky flux
                covptot_i = in_covptot_i + 0
                if preclr <= 0:
                    clc_i = in_clc_i + covptot_i
                    covptot_i = 0.0
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
                out_aph_i = dt * RG * dtgdp_i / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
                daph_i = dt * RG * dtgdp_i / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
                aph_i_s = -daph_i

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
                    * (RG * RPECONS / 0.00509)
                    * (0.00509 * covpclr / (preclr1 * sqrt(in_ap / in_aph_s))) ** 0.4223
                )
                preclr_i += xx * sqrt(in_ap / in_aph_s) * beta_i / covpclr
                out_ap_i += 0.5 * xx * preclr1 * beta_i / (covpclr * sqrt(in_ap * in_aph_s))
                aph_i_s -= (
                    0.5 * xx * preclr1 * sqrt(in_ap / in_aph_s) * beta_i / (covpclr * in_aph_s)
                )
                covpclr_i -= (
                    (xx * preclr1 * sqrt(in_ap / in_aph_s) * beta_i / covpclr ** 2)
                    + (in_qsat - qlim) * qe_i / (1 - out_clc) ** 2
                ) - prtot * preclr_i / covptot1
                out_qsat_i += qe_i - covpclr * qe_i / (1 - out_clc) ** 2
                qlim_i = covpclr * qe_i / (1 - out_clc) ** 2
                clc_i -= 2 * (in_qsat - qlim) * covpclr * qe_i / (1 - out_clc) ** 3
                prtot_i += covpclr * preclr_i / covptot1
                covptot_i -= prtot * covpclr * preclr_i / covptot1 ** 2
            else:
                clc_i = in_clc_i + 0
                corqs_i = 0.0
                covpclr_i = 0.0
                covptot_i = 0
                daph_i = 0.0
                out_aph_i = 0.0
                aph_i_s = 0
                out_qsat_i = 0.0
                prtot_i = 0.0
                qlim_i = 0.0

            # new precipitation
            rfl_i += prtot_i
            sfl_i += prtot_i
            fwatr_i += dr1 * (rfl_i - sfl_i)
            dr_i = fwatr1 * rfl_i + (1 - fwatr1) * sfl_i

            # update rain fraction and freezing
            # note: impact of new temperature out_t_i on fwat_i in neglected here
            if t < RTT:
                dp_i += rfreeze_i * cons2 * prr
                prr_i = rfreeze_i * cons2 * dp
            else:
                prr_i = 0.0
            prr_i += cons2 * dp * dr_i
            prs_i = cons2 * dp * dr_i
            dp_i += cons2 * (prr + prs) * dr_i

            if out_clc > ZEPS2:
                # diagnostic calculation of rain production from cloud ice
                if __INLINED(LEVAPLS2 or LDRAIN1D):
                    icrit = 0.0001
                    lcrit = 1.9 * RCLCRIT
                else:
                    icrit = 2 * RCLCRIT
                    lcrit = 2 * RCLCRIT
                prs_i -= qiwc_i
                qiwc_i += prs_i
                qinew_i = -prs_i
                clc_i += qinew_i * cldi * itmp2
                cldi_i = qinew_i * out_clc * itmp2
                di_i = -qinew_i * out_clc * cldi * itmp2

                # regularization
                if __INLINED(LREGCL):
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
                if __INLINED(LREGCL):
                    ltmp4 = ckcodtla
                else:
                    ltmp4 = ckcodtl
                cldl_i += 2 * ltmp4 * ltmp1 * cldl * dl_i / lcrit ** 2

                qlwc_i += cldl_i / out_clc
                clc_i -= qlwc1 * cldl_i / out_clc ** 2

            if sfl[0, 0, -1] != 0.0:
                snmlt_i = -out_t_i / cons + rfl_i - sfl_i

                if sfl[0, 0, -1] <= z2s:
                    sfl_i += snmlt_i
                    z2s_i = 0.0
                else:
                    z2s_i = snmlt_i + 0

                if t2 > meltp2:
                    cons_i = out_t_i * snmlt / cons ** 2 + (t2 - meltp2) * z2s_i
                    out_t_i += cons * z2s_i
                else:
                    cons_i = out_t_i * snmlt / cons ** 2

                dp_i += cons_i * cons2 / lfdcp
                lfdcp_i = -cons2 * dp * cons_i / lfdcp ** 2
            else:
                lfdcp_i = 0.0

            if covpclr1 < 0:
                covpclr_i = 0.0
            covptot_i += covpclr_i
            clc_i -= covpclr_i

            if out_clc > covptot:
                clc_i += covptot_i

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
                if __INLINED(LREGCL):
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
            rodqsdp_i = -RG * (dqsdz_i + dtdzmo_i * ldcp * fac3)
            ldcp_i = -dtdzmo_i * (RG * rodqsdp + dtdzmo * dqsdtemp) * fac3
            fwat_i += ldcp_i * (lvdcp - lsdcp)
            lvdcp_i += fwat * ldcp_i
            lsdcp_i += (1 - fwat) * ldcp_i
            rho_i -= rodqsdp_i * in_qsat * fac2
            out_qsat_i -= rodqsdp_i * rho * fac2
            out_ap_i += rodqsdp_i * rho * in_qsat * fac2 ** 2
            foeew_i = -RETV * rodqsdp_i * rho * in_qsat * fac2 ** 2
            out_ap_i += rho_i * fac1
            out_t_i -= rho_i * in_ap * fac1 / t2

            if lude >= RLMIN and in_lu[0, 0, 1] >= ZEPS2:
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
            out_aph_i += RG * gdp_i[0, 0, 0] / (in_aph[0, 0, 1] - in_aph[0, 0, 0]) ** 2
            daph_i += RG * gdp_i[0, 0, 0] / (in_aph[0, 0, 1] - in_aph[0, 0, 0]) ** 2

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

                if __INLINED(LREGCL):
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
                    -0.5 / tmp3 * (qpd * scalm * clc_i) / (qcd - scalm * (qt - qcrit)) ** 2
                ) - qpd_i
                qcrit_i = (
                    0.5 / tmp3 * (qpd * scalm * clc_i) / (qcd - scalm * (qt - qcrit)) ** 2
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
            if t2 < RTICE:
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
            esdp_i = RETV * cor_i * cor ** 2
            facw_i = fwat * fac_i
            faci_i = (1 - fwat) * fac_i
            fwat_i += (facw - faci) * fac_i
            out_t_i -= 2 * (R5IES * faci_i / (t2 - R4IES) ** 3 + R5LES * facw_i / (t2 - R4LES) ** 3)

            if esdp1 > ZQMAX:
                esdp_i = 0.0
            foeew_i += esdp_i / in_ap
            out_ap_i -= esdp_i * foeew / in_ap ** 2

            if t2 < RTT:
                z3es = R3IES
                z4es = R4IES
            else:
                z3es = R3LES
                z4es = R4LES
            out_t_i += z3es * (RTT - z4es) * foeew_i * foeew / (t2 - z4es) ** 2

            if t2 < RTT:
                out_t_i += 0.545 * 0.17 * fwat_i / cosh(0.17 * (t2 - RLPTRC)) ** 2
        with interval(1, -2):
            # set up constants required
            ckcodtl = 2 * RKCONV * dt
            ckcodtla = ckcodtl / 100
            ckcodti = 5 * RKCONV * dt
            ckcodtia = ckcodti / 100
            cons2 = 1 / (RG * dt)
            cons3 = RLVTT / RCPD
            meltp2 = RTT + 2

            # enthalpy fluxes due to precipitation
            fplsn_i = in_fplsn_i - in_fhpsn_i * RLSTT
            fplsl_i = in_fplsl_i - in_fhpsl_i * RLVTT

            # incrementation of t and q, and fluxes swap
            rfl_i = rfl_i[0, 0, 1] + fplsl_i[0, 0, 1]
            sfl_i = sfl_i[0, 0, 1] + fplsn_i[0, 0, 1]

            # qice tendency
            out_qi_i = -in_tnd_qi_i / dt
            qiwc_i = in_tnd_qi_i / dt

            # qliq tendency
            out_ql_i = -in_tnd_ql_i / dt
            qlwc_i = in_tnd_ql_i / dt

            # T tendency
            gdp_i = -in_tnd_t_i * (
                lvdcp * evapr
                + lsdcp * evaps
                + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
                - (lsdcp - lvdcp) * rfreeze3
            )
            condl_i = in_tnd_t_i * lvdcp
            condi_i = in_tnd_t_i * lsdcp
            evapr_i = -in_tnd_t_i * lvdcp * gdp
            evaps_i = -in_tnd_t_i * lsdcp * gdp
            lvdcp_i = in_tnd_t_i * (condl2 - evapr * gdp)
            lsdcp_i = in_tnd_t_i * (condi2 - evaps * gdp)
            out_lude_i = -in_tnd_t_i * gdp * (fwat * lvdcp + (1 - fwat) * lsdcp)
            lvdcp_i -= in_tnd_t_i * in_lude * gdp * fwat
            lsdcp_i -= in_tnd_t_i * in_lude * gdp * (1 - fwat)
            fwat_i = -in_tnd_t_i * in_lude * gdp * (lvdcp - lsdcp)
            lvdcp_i -= in_tnd_t_i * rfreeze3 * gdp
            lsdcp_i += in_tnd_t_i * rfreeze3 * gdp
            rfreeze_i = in_tnd_t_i * (lsdcp - lvdcp) * gdp

            # q tendency
            gdp_i += in_tnd_q_i * (in_lude + evapr + evaps)
            out_lude_i += in_tnd_q_i * gdp
            evapr_i += in_tnd_q_i * gdp
            evaps_i += in_tnd_q_i * gdp
            condl_i -= in_tnd_q_i
            condi_i -= in_tnd_q_i

            # clipping of final qv
            rn_i = rfl_i + 0
            sn_i = sfl_i + 0

            # note: the extra condensation due to the adjustment goes directly
            # to precipitation
            dq_i = (fwatr2 * condl_i + (1 - fwatr2) * condi_i) / dt
            # fwatr_i = (condl_i - condi_i) * dq / dt + dr2 * (rn_i - sn_i)
            dr2_i = fwatr2 * rn_i + (1 - fwatr2) * sn_i

            # update rain fraction adn freezing
            # note: impact of new temperature out_t_i on fwat_i is neglected here
            if t3 < RTT:
                fwat_i += dr2 * rfreeze_i
                dr2_i += fwat * rfreeze_i
            fwatr_i = 0.0

            dq_i += cons2 * dp * dr2_i
            dp_i = cons2 * dq * dr2_i

            if qold1 >= q:
                # regularization
                if __INLINED(LREGCL):
                    dq_i *= 0.7
                qold_i = dq_i
                out_q_i = -dq_i
            else:
                qold_i = 0.0
                out_q_i = 0.0

            out_ap_i = 0.0
            out_t_i = 0.0
            out_ap_i, told, out_t_i, qold, out_q_i = cuadjtqs_ad(
                in_ap, out_ap_i, told, out_t_i, qold, out_q_i
            )

            # first guess t and q
            out_q_i += qold_i
            dtdt_i = dt * out_t_i
            dqdt_i = dt * out_q_i

            # incrementation of T and q
            # T tendency
            gdp_i -= dtdt_i * (
                lvdcp * evapr
                + lsdcp * evaps
                + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
                + (lsdcp - lvdcp) * rfreeze1
            )
            condl_i += dtdt_i * lvdcp
            condi_i += dtdt_i * lsdcp
            evapr_i -= dtdt_i * lvdcp * gdp
            evaps_i -= dtdt_i * lsdcp * gdp
            lvdcp_i += dtdt_i * (condl1 - evapr * gdp)
            lsdcp_i += dtdt_i * (condi1 - evaps * gdp)
            out_lude_i -= dtdt_i * gdp * (fwat * lvdcp + (1 - fwat) * lsdcp)
            lvdcp_i -= dtdt_i * in_lude * gdp * fwat
            lsdcp_i -= dtdt_i * in_lude * gdp * (1 - fwat)
            fwat_i -= dtdt_i * in_lude * gdp * (lvdcp - lsdcp)
            lvdcp_i -= dtdt_i * rfreeze1 * gdp
            lsdcp_i += dtdt_i * rfreeze1 * gdp
            rfreeze_i += dtdt_i * (lsdcp - lvdcp) * gdp

            # q tendency
            gdp_i += dqdt_i * (in_lude + evapr + evaps)
            out_lude_i += dqdt_i * gdp
            evapr_i += dqdt_i * gdp
            evaps_i += dqdt_i * gdp
            condl_i -= dqdt_i
            condi_i -= dqdt_i

            if prtot > ZEPS2 and covpclr > ZEPS2 and (LEVAPLS2 or LDRAIN1D):
                # ice proportion
                evaps_i -= sfl_i
                sfl_i += dpr * evaps_i / prtot
                dpr_i = sfln2 * evaps_i / prtot
                prtot_i = dpr * sfln2 * evaps_i / prtot ** 2

                # warm proportion
                evapr_i -= rfl_i
                rfl_i += dpr * evapr_i / prtot
                dpr_i += rfln2 * evapr_i / prtot
                prtot_i -= dpr * rfln2 * evapr_i / prtot ** 2

                # take away clear sky flux
                covptot_i = covptot_i[0, 0, 1] + in_covptot_i
                if preclr <= 0:
                    clc_i = in_clc_i + covptot_i
                    covptot_i = 0.0
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
                out_aph_i = dt * RG * dtgdp_i / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
                daph_i = dt * RG * dtgdp_i / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
                aph_i_s = aph_i_s[0, 0, 1]

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
                    * (RG * RPECONS / 0.00509)
                    * (0.00509 * covpclr / (preclr1 * sqrt(in_ap / in_aph_s))) ** 0.4223
                )
                preclr_i += xx * sqrt(in_ap / in_aph_s) * beta_i / covpclr
                out_ap_i += 0.5 * xx * preclr1 * beta_i / (covpclr * sqrt(in_ap * in_aph_s))
                aph_i_s -= (
                    0.5 * xx * preclr1 * sqrt(in_ap / in_aph_s) * beta_i / (covpclr * in_aph_s)
                )
                covpclr_i -= (
                    (xx * preclr1 * sqrt(in_ap / in_aph_s) * beta_i / covpclr ** 2)
                    + (in_qsat - qlim) * qe_i / (1 - out_clc) ** 2
                ) - prtot * preclr_i / covptot1
                out_qsat_i += qe_i - covpclr * qe_i / (1 - out_clc) ** 2
                qlim_i = covpclr * qe_i / (1 - out_clc) ** 2
                clc_i -= 2 * (in_qsat - qlim) * covpclr * qe_i / (1 - out_clc) ** 3
                prtot_i += covpclr * preclr_i / covptot1
                covptot_i -= prtot * covpclr * preclr_i / covptot1 ** 2
            else:
                clc_i = in_clc_i + 0
                corqs_i = 0.0
                covpclr_i = 0.0
                covptot_i = covptot_i[0, 0, 1]
                daph_i = 0.0
                out_aph_i = 0.0
                aph_i_s = aph_i_s[0, 0, 1]
                out_qsat_i = 0.0
                prtot_i = 0.0
                qlim_i = 0.0

            # new precipitation
            rfl_i += prtot_i
            sfl_i += prtot_i
            fwatr_i += dr1 * (rfl_i - sfl_i)
            dr_i = fwatr1 * rfl_i + (1 - fwatr1) * sfl_i

            # update rain fraction and freezing
            # note: impact of new temperature out_t_i on fwat_i in neglected here
            if t < RTT:
                dp_i += rfreeze_i * cons2 * prr
                prr_i = rfreeze_i * cons2 * dp
                rfreeze_i = 0.0
            else:
                prr_i = 0.0
            prr_i += cons2 * dp * dr_i
            prs_i = cons2 * dp * dr_i
            dp_i += cons2 * (prr + prs) * dr_i

            if out_clc > ZEPS2:
                # diagnostic calculation of rain production from cloud ice
                if __INLINED(LEVAPLS2 or LDRAIN1D):
                    icrit = 0.0001
                    lcrit = 1.9 * RCLCRIT
                else:
                    icrit = 2 * RCLCRIT
                    lcrit = 2 * RCLCRIT
                prs_i -= qiwc_i
                qiwc_i += prs_i
                qinew_i = -prs_i
                clc_i += qinew_i * cldi * itmp2
                cldi_i = qinew_i * out_clc * itmp2
                di_i = -qinew_i * out_clc * cldi * itmp2

                # regularization
                if __INLINED(LREGCL):
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
                if __INLINED(LREGCL):
                    ltmp4 = ckcodtla
                else:
                    ltmp4 = ckcodtl
                cldl_i += 2 * ltmp4 * ltmp1 * cldl * dl_i / lcrit ** 2

                qlwc_i += cldl_i / out_clc
                clc_i -= qlwc1 * cldl_i / out_clc ** 2

            if sfl[0, 0, -1] != 0.0:
                snmlt_i = -out_t_i / cons + rfl_i - sfl_i

                if sfl[0, 0, -1] <= z2s:
                    sfl_i += snmlt_i
                    z2s_i = 0.0
                else:
                    z2s_i = snmlt_i + 0

                if t2 > meltp2:
                    cons_i = out_t_i * snmlt / cons ** 2 + (t2 - meltp2) * z2s_i
                    out_t_i += cons * z2s_i
                else:
                    cons_i = out_t_i * snmlt / cons ** 2

                dp_i += cons_i * cons2 / lfdcp
                lfdcp_i = -cons2 * dp * cons_i / lfdcp ** 2
            else:
                lfdcp_i = 0.0

            if covpclr1 < 0:
                covpclr_i = 0.0
            covptot_i += covpclr_i
            clc_i -= covpclr_i

            if out_clc > covptot:
                clc_i += covptot_i

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
                if __INLINED(LREGCL):
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
            rodqsdp_i = -RG * (dqsdz_i + dtdzmo_i * ldcp * fac3)
            ldcp_i = -dtdzmo_i * (RG * rodqsdp + dtdzmo * dqsdtemp) * fac3
            fwat_i += ldcp_i * (lvdcp - lsdcp)
            lvdcp_i += fwat * ldcp_i
            lsdcp_i += (1 - fwat) * ldcp_i
            rho_i -= rodqsdp_i * in_qsat * fac2
            out_qsat_i -= rodqsdp_i * rho * fac2
            out_ap_i += rodqsdp_i * rho * in_qsat * fac2 ** 2
            foeew_i = -RETV * rodqsdp_i * rho * in_qsat * fac2 ** 2
            out_ap_i += rho_i * fac1
            out_t_i -= rho_i * in_ap * fac1 / t2

            if lude >= RLMIN and in_lu[0, 0, 1] >= ZEPS2:
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
            out_aph_i += RG * gdp_i[0, 0, 0] / (in_aph[0, 0, 1] - in_aph[0, 0, 0]) ** 2
            daph_i += RG * gdp_i[0, 0, 0] / (in_aph[0, 0, 1] - in_aph[0, 0, 0]) ** 2

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

                if __INLINED(LREGCL):
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
                    -0.5 / tmp3 * (qpd * scalm * clc_i) / (qcd - scalm * (qt - qcrit)) ** 2
                ) - qpd_i
                qcrit_i = (
                    0.5 / tmp3 * (qpd * scalm * clc_i) / (qcd - scalm * (qt - qcrit)) ** 2
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
            if t2 < RTICE:
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
            esdp_i = RETV * cor_i * cor ** 2
            facw_i = fwat * fac_i
            faci_i = (1 - fwat) * fac_i
            fwat_i += (facw - faci) * fac_i
            out_t_i -= 2 * (R5IES * faci_i / (t2 - R4IES) ** 3 + R5LES * facw_i / (t2 - R4LES) ** 3)

            if esdp1 > ZQMAX:
                esdp_i = 0.0
            foeew_i += esdp_i / in_ap
            out_ap_i -= esdp_i * foeew / in_ap ** 2

            if t2 < RTT:
                z3es = R3IES
                z4es = R4IES
            else:
                z3es = R3LES
                z4es = R4LES
            out_t_i += z3es * (RTT - z4es) * foeew_i * foeew / (t2 - z4es) ** 2

            if t2 < RTT:
                out_t_i += 0.545 * 0.17 * fwat_i / cosh(0.17 * (t2 - RLPTRC)) ** 2
        with interval(0, 1):
            # set up constants required
            ckcodtl = 2 * RKCONV * dt
            ckcodtla = ckcodtl / 100
            ckcodti = 5 * RKCONV * dt
            ckcodtia = ckcodti / 100
            cons2 = 1 / (RG * dt)
            cons3 = RLVTT / RCPD
            meltp2 = RTT + 2

            # enthalpy fluxes due to precipitation
            fplsn_i = in_fplsn_i - in_fhpsn_i * RLSTT
            fplsl_i = in_fplsl_i - in_fhpsl_i * RLVTT

            # incrementation of t and q, and fluxes swap
            rfl_i = rfl_i[0, 0, 1] + fplsl_i[0, 0, 1]
            sfl_i = sfl_i[0, 0, 1] + fplsn_i[0, 0, 1]

            # qice tendency
            out_qi_i = -in_tnd_qi_i / dt
            qiwc_i = in_tnd_qi_i / dt

            # qliq tendency
            out_ql_i = -in_tnd_ql_i / dt
            qlwc_i = in_tnd_ql_i / dt

            # T tendency
            gdp_i = -in_tnd_t_i * (
                lvdcp * evapr
                + lsdcp * evaps
                + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
                - (lsdcp - lvdcp) * rfreeze3
            )
            condl_i = in_tnd_t_i * lvdcp
            condi_i = in_tnd_t_i * lsdcp
            evapr_i = -in_tnd_t_i * lvdcp * gdp
            evaps_i = -in_tnd_t_i * lsdcp * gdp
            lvdcp_i = in_tnd_t_i * (condl2 - evapr * gdp)
            lsdcp_i = in_tnd_t_i * (condi2 - evaps * gdp)
            out_lude_i = -in_tnd_t_i * gdp * (fwat * lvdcp + (1 - fwat) * lsdcp)
            lvdcp_i -= in_tnd_t_i * in_lude * gdp * fwat
            lsdcp_i -= in_tnd_t_i * in_lude * gdp * (1 - fwat)
            fwat_i = -in_tnd_t_i * in_lude * gdp * (lvdcp - lsdcp)
            lvdcp_i -= in_tnd_t_i * rfreeze3 * gdp
            lsdcp_i += in_tnd_t_i * rfreeze3 * gdp
            rfreeze_i = in_tnd_t_i * (lsdcp - lvdcp) * gdp

            # q tendency
            gdp_i += in_tnd_q_i * (in_lude + evapr + evaps)
            out_lude_i += in_tnd_q_i * gdp
            evapr_i += in_tnd_q_i * gdp
            evaps_i += in_tnd_q_i * gdp
            condl_i -= in_tnd_q_i
            condi_i -= in_tnd_q_i

            # clipping of final qv
            rn_i = rfl_i + 0
            sn_i = sfl_i + 0

            # note: the extra condensation due to the adjustment goes directly
            # to precipitation
            dq_i = (fwatr2 * condl_i + (1 - fwatr2) * condi_i) / dt
            # fwatr_i = (condl_i - condi_i) * dq / dt + dr2 * (rn_i - sn_i)
            dr2_i = fwatr2 * rn_i + (1 - fwatr2) * sn_i

            # update rain fraction adn freezing
            # note: impact of new temperature out_t_i on fwat_i is neglected here
            if t3 < RTT:
                fwat_i += dr2 * rfreeze_i
                dr2_i += fwat * rfreeze_i
            fwatr_i = 0.0

            dq_i += cons2 * dp * dr2_i
            dp_i = cons2 * dq * dr2_i

            if qold1 >= q:
                # regularization
                if __INLINED(LREGCL):
                    dq_i *= 0.7
                qold_i = dq_i
                out_q_i = -dq_i
            else:
                qold_i = 0.0
                out_q_i = 0.0

            out_ap_i = 0.0
            out_t_i = 0.0
            out_ap_i, told, out_t_i, qold, out_q_i = cuadjtqs_ad(
                in_ap, out_ap_i, told, out_t_i, qold, out_q_i
            )

            # first guess t and q
            out_q_i += qold_i
            dtdt_i = dt * out_t_i
            dqdt_i = dt * out_q_i

            # incrementation of T and q
            # T tendency
            gdp_i -= dtdt_i * (
                lvdcp * evapr
                + lsdcp * evaps
                + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
                + (lsdcp - lvdcp) * rfreeze1
            )
            condl_i += dtdt_i * lvdcp
            condi_i += dtdt_i * lsdcp
            evapr_i -= dtdt_i * lvdcp * gdp
            evaps_i -= dtdt_i * lsdcp * gdp
            lvdcp_i += dtdt_i * (condl1 - evapr * gdp)
            lsdcp_i += dtdt_i * (condi1 - evaps * gdp)
            out_lude_i -= dtdt_i * gdp * (fwat * lvdcp + (1 - fwat) * lsdcp)
            lvdcp_i -= dtdt_i * in_lude * gdp * fwat
            lsdcp_i -= dtdt_i * in_lude * gdp * (1 - fwat)
            fwat_i -= dtdt_i * in_lude * gdp * (lvdcp - lsdcp)
            lvdcp_i -= dtdt_i * rfreeze1 * gdp
            lsdcp_i += dtdt_i * rfreeze1 * gdp
            rfreeze_i += dtdt_i * (lsdcp - lvdcp) * gdp

            # q tendency
            gdp_i += dqdt_i * (in_lude + evapr + evaps)
            out_lude_i += dqdt_i * gdp
            evapr_i += dqdt_i * gdp
            evaps_i += dqdt_i * gdp
            condl_i -= dqdt_i
            condi_i -= dqdt_i

            if prtot > ZEPS2 and covpclr > ZEPS2 and (LEVAPLS2 or LDRAIN1D):
                # ice proportion
                evaps_i -= sfl_i
                sfl_i += dpr * evaps_i / prtot
                dpr_i = sfln2 * evaps_i / prtot
                prtot_i = dpr * sfln2 * evaps_i / prtot ** 2

                # warm proportion
                evapr_i -= rfl_i
                rfl_i += dpr * evapr_i / prtot
                dpr_i += rfln2 * evapr_i / prtot
                prtot_i -= dpr * rfln2 * evapr_i / prtot ** 2

                # take away clear sky flux
                covptot_i = covptot_i[0, 0, 1] + in_covptot_i
                if preclr <= 0:
                    clc_i = in_clc_i + covptot_i
                    covptot_i = 0.0
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
                out_aph_i = dt * RG * dtgdp_i / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
                daph_i = dt * RG * dtgdp_i / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
                aph_i_s = aph_i_s[0, 0, 1]

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
                    * (RG * RPECONS / 0.00509)
                    * (0.00509 * covpclr / (preclr1 * sqrt(in_ap / in_aph_s))) ** 0.4223
                )
                preclr_i += xx * sqrt(in_ap / in_aph_s) * beta_i / covpclr
                out_ap_i += 0.5 * xx * preclr1 * beta_i / (covpclr * sqrt(in_ap * in_aph_s))
                aph_i_s -= (
                    0.5 * xx * preclr1 * sqrt(in_ap / in_aph_s) * beta_i / (covpclr * in_aph_s)
                )
                covpclr_i -= (
                    (xx * preclr1 * sqrt(in_ap / in_aph_s) * beta_i / covpclr ** 2)
                    + (in_qsat - qlim) * qe_i / (1 - out_clc) ** 2
                ) - prtot * preclr_i / covptot1
                out_qsat_i += qe_i - covpclr * qe_i / (1 - out_clc) ** 2
                qlim_i = covpclr * qe_i / (1 - out_clc) ** 2
                clc_i -= 2 * (in_qsat - qlim) * covpclr * qe_i / (1 - out_clc) ** 3
                prtot_i += covpclr * preclr_i / covptot1
                covptot_i -= prtot * covpclr * preclr_i / covptot1 ** 2
            else:
                clc_i = in_clc_i + 0
                corqs_i = 0.0
                covpclr_i = 0.0
                covptot_i = covptot_i[0, 0, 1]
                daph_i = 0.0
                out_aph_i = 0.0
                aph_i_s = aph_i_s[0, 0, 1]
                out_qsat_i = 0.0
                prtot_i = 0.0
                qlim_i = 0.0

            # new precipitation
            rfl_i += prtot_i
            sfl_i += prtot_i
            fwatr_i += dr1 * (rfl_i - sfl_i)
            dr_i = fwatr1 * rfl_i + (1 - fwatr1) * sfl_i

            # update rain fraction and freezing
            # note: impact of new temperature out_t_i on fwat_i in neglected here
            if t < RTT:
                dp_i += rfreeze_i * cons2 * prr
                prr_i = rfreeze_i * cons2 * dp
                rfreeze_i = 0.0
            else:
                prr_i = 0.0
            prr_i += cons2 * dp * dr_i
            prs_i = cons2 * dp * dr_i
            dp_i += cons2 * (prr + prs) * dr_i

            if out_clc > ZEPS2:
                # diagnostic calculation of rain production from cloud ice
                if __INLINED(LEVAPLS2 or LDRAIN1D):
                    icrit = 0.0001
                    lcrit = 1.9 * RCLCRIT
                else:
                    icrit = 2 * RCLCRIT
                    lcrit = 2 * RCLCRIT
                prs_i -= qiwc_i
                qiwc_i += prs_i
                qinew_i = -prs_i
                clc_i += qinew_i * cldi * itmp2
                cldi_i = qinew_i * out_clc * itmp2
                di_i = -qinew_i * out_clc * cldi * itmp2

                # regularization
                if __INLINED(LREGCL):
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
                if __INLINED(LREGCL):
                    ltmp4 = ckcodtla
                else:
                    ltmp4 = ckcodtl
                cldl_i += 2 * ltmp4 * ltmp1 * cldl * dl_i / lcrit ** 2

                qlwc_i += cldl_i / out_clc
                clc_i -= qlwc1 * cldl_i / out_clc ** 2

            lfdcp_i = 0.0

            if covpclr1 < 0:
                covpclr_i = 0.0
            covptot_i += covpclr_i
            clc_i -= covpclr_i

            if out_clc > covptot:
                clc_i += covptot_i

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
                if __INLINED(LREGCL):
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
            rodqsdp_i = -RG * (dqsdz_i + dtdzmo_i * ldcp * fac3)
            ldcp_i = -dtdzmo_i * (RG * rodqsdp + dtdzmo * dqsdtemp) * fac3
            fwat_i += ldcp_i * (lvdcp - lsdcp)
            lvdcp_i += fwat * ldcp_i
            lsdcp_i += (1 - fwat) * ldcp_i
            rho_i -= rodqsdp_i * in_qsat * fac2
            out_qsat_i -= rodqsdp_i * rho * fac2
            out_ap_i += rodqsdp_i * rho * in_qsat * fac2 ** 2
            foeew_i = -RETV * rodqsdp_i * rho * in_qsat * fac2 ** 2
            out_ap_i += rho_i * fac1
            out_t_i -= rho_i * in_ap * fac1 / t2

            if lude >= RLMIN and in_lu[0, 0, 1] >= ZEPS2:
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
            out_aph_i += RG * gdp_i[0, 0, 0] / (in_aph[0, 0, 1] - in_aph[0, 0, 0]) ** 2
            daph_i += RG * gdp_i[0, 0, 0] / (in_aph[0, 0, 1] - in_aph[0, 0, 0]) ** 2

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

                if __INLINED(LREGCL):
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
                    -0.5 / tmp3 * (qpd * scalm * clc_i) / (qcd - scalm * (qt - qcrit)) ** 2
                ) - qpd_i
                qcrit_i = (
                    0.5 / tmp3 * (qpd * scalm * clc_i) / (qcd - scalm * (qt - qcrit)) ** 2
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
            if t2 < RTICE:
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
            esdp_i = RETV * cor_i * cor ** 2
            facw_i = fwat * fac_i
            faci_i = (1 - fwat) * fac_i
            fwat_i += (facw - faci) * fac_i
            out_t_i -= 2 * (R5IES * faci_i / (t2 - R4IES) ** 3 + R5LES * facw_i / (t2 - R4LES) ** 3)

            if esdp1 > ZQMAX:
                esdp_i = 0.0
            foeew_i += esdp_i / in_ap
            out_ap_i -= esdp_i * foeew / in_ap ** 2

            if t2 < RTT:
                z3es = R3IES
                z4es = R4IES
            else:
                z3es = R3LES
                z4es = R4LES
            out_t_i += z3es * (RTT - z4es) * foeew_i * foeew / (t2 - z4es) ** 2

            if t2 < RTT:
                out_t_i += 0.545 * 0.17 * fwat_i / cosh(0.17 * (t2 - RLPTRC)) ** 2

    with computation(BACKWARD):
        with interval(-1, None):
            out_aph_i += dp_i[0, 0, -1] - daph_i[0, 0, -1]
        with interval(1, -1):
            out_aph_i += -dp_i[0, 0, 0] + dp_i[0, 0, -1] - daph_i[0, 0, -1]
            out_lu_i = -dlu_i

            zz = RLVTT * lvdcp_i + RLSTT * lsdcp_i + RLMLT * lfdcp_i
            out_q_i -= zz * RCPD * RVTMP2 / (RCPD + RCPD * RVTMP2 * q) ** 2
            out_supsat_i = dt * out_q_i
            out_tnd_cml_t_i = dt * out_t_i
            out_tnd_cml_q_i = dt * out_q_i
            out_tnd_cml_ql_i = dt * out_ql_i
            out_tnd_cml_qi_i = dt * out_qi_i
        with interval(0, 1):
            out_aph_i += -dp_i[0, 0, 0]
            out_lu_i = 0.0

            zz = RLVTT * lvdcp_i + RLSTT * lsdcp_i + RLMLT * lfdcp_i
            out_q_i -= zz * RCPD * RVTMP2 / (RCPD + RCPD * RVTMP2 * q) ** 2
            out_supsat_i = dt * out_q_i
            out_tnd_cml_t_i = dt * out_t_i
            out_tnd_cml_q_i = dt * out_q_i
            out_tnd_cml_ql_i = dt * out_ql_i
            out_tnd_cml_qi_i = dt * out_qi_i
