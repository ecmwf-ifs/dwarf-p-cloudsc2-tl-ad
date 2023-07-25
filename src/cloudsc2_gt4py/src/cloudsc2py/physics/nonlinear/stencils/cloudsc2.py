# -*- coding: utf-8 -*-
from gt4py.cartesian import gtscript

from cloudsc2py.physics.common.stencils.fcttre import foealfa, foeewm
from cloudsc2py.physics.nonlinear.stencils.cuadjtqs import cuadjtqs_nl
from ifs_physics_common.framework.stencil import stencil_collection
from ifs_physics_common.utils.f2py import ported_function


@ported_function(from_file="cloudsc2_nl/cloudsc2.F90", from_line=235, to_line=735)
@stencil_collection("cloudsc2_nl")
def cloudsc2_nl_def(
    in_ap: gtscript.Field["float"],
    in_aph: gtscript.Field["float"],
    in_eta: gtscript.Field[gtscript.K, "float"],
    in_lu: gtscript.Field["float"],
    in_lude: gtscript.Field["float"],
    in_mfd: gtscript.Field["float"],
    in_mfu: gtscript.Field["float"],
    in_q: gtscript.Field["float"],
    in_qi: gtscript.Field["float"],
    in_ql: gtscript.Field["float"],
    in_qsat: gtscript.Field["float"],
    in_supsat: gtscript.Field["float"],
    in_t: gtscript.Field["float"],
    in_tnd_cml_q: gtscript.Field["float"],
    in_tnd_cml_qi: gtscript.Field["float"],
    in_tnd_cml_ql: gtscript.Field["float"],
    in_tnd_cml_t: gtscript.Field["float"],
    out_clc: gtscript.Field["float"],
    out_covptot: gtscript.Field["float"],
    out_fhpsl: gtscript.Field["float"],
    out_fhpsn: gtscript.Field["float"],
    out_fplsl: gtscript.Field["float"],
    out_fplsn: gtscript.Field["float"],
    out_tnd_q: gtscript.Field["float"],
    out_tnd_qi: gtscript.Field["float"],
    out_tnd_ql: gtscript.Field["float"],
    out_tnd_t: gtscript.Field["float"],
    tmp_aph_s: gtscript.Field[gtscript.IJ, "float"],
    tmp_covptot: gtscript.Field[gtscript.IJ, "float"],
    tmp_rfl: gtscript.Field[gtscript.IJ, "float"],
    tmp_sfl: gtscript.Field[gtscript.IJ, "float"],
    tmp_trpaus: gtscript.Field[gtscript.IJ, "float"],
    *,
    dt: "float",
):
    from __externals__ import (
        LDRAIN1D,
        LEVAPLS2,
        LPHYLIN,
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
        tmp_rfl[0, 0] = 0.0
        tmp_sfl[0, 0] = 0.0
        tmp_covptot[0, 0] = 0.0

    with computation(PARALLEL), interval(0, -1):
        # first guess values for T
        t = in_t + dt * in_tnd_cml_t

    # eta value at tropopause
    with computation(FORWARD), interval(0, 1):
        tmp_trpaus[0, 0] = 0.1
    with computation(FORWARD), interval(0, -2):
        if in_eta > 0.1 and in_eta < 0.4 and t[0, 0, 0] > t[0, 0, 1]:
            tmp_trpaus[0, 0] = in_eta

    with computation(FORWARD), interval(0, -1):
        # first guess values for q, ql and qi
        q = in_q + dt * in_tnd_cml_q + in_supsat
        ql = in_ql + dt * in_tnd_cml_ql
        qi = in_qi + dt * in_tnd_cml_qi

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
        if LPHYLIN or LDRAIN1D:
            if t < RTT:
                fwat = 0.545 * (tanh(0.17 * (t - RLPTRC)) + 1)
                z3es = R3IES
                z4es = R4IES
            else:
                fwat = 1.0
                z3es = R3LES
                z4es = R4LES
            foeew = R2ES * exp(z3es * (t - RTT) / (t - z4es))
            esdp = min(foeew / in_ap, ZQMAX)
        else:
            fwat = foealfa(t)
            foeew = foeewm(t)
            esdp = foeew / in_ap
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
            qsat = in_qsat * (1.8 - 0.003 * t)
        else:
            qsat = in_qsat
        qcrit = crh2 * qsat

        # simple uniform distribution of total water from Leutreut & Li (1990)
        qt = q + ql + qi
        if qt < qcrit:
            out_clc[0, 0, 0] = 0.0
            qc = 0.0
        elif qt >= qsat:
            out_clc[0, 0, 0] = 1.0
            qc = (1 - scalm) * (qsat - qcrit)
        else:
            qpd = qsat - qt
            qcd = qsat - qcrit
            out_clc[0, 0, 0] = 1 - sqrt(qpd / (qcd - scalm * (qt - qcrit)))
            qc = (scalm * qpd + (1 - scalm) * qcd) * (out_clc**2)

        # add convective component
        gdp = RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
        lude = dt * in_lude * gdp
        lo1 = lude >= RLMIN and in_lu[0, 0, 1] >= ZEPS2
        if lo1:
            out_clc[0, 0, 0] += (1 - out_clc) * (1 - exp(-lude / in_lu[0, 0, 1]))
            qc += lude

        # add compensating subsidence component
        rho = in_ap / (RD * t)
        rodqsdp = -rho * in_qsat / (in_ap - RETV * foeew)
        ldcp = fwat * lvdcp + (1 - fwat) * lsdcp
        dtdzmo = RG * (1 / RCPD - ldcp * rodqsdp) / (1 + ldcp * dqsdtemp)
        dqsdz = dqsdtemp * dtdzmo - RG * rodqsdp
        dqc = min(dt * dqsdz * (in_mfu + in_mfd) / rho, qc)
        qc -= dqc

        # new cloud liquid/ice contents and condensation rates (liquid/ice)
        qlwc = qc * fwat
        qiwc = qc * (1 - fwat)
        condl = (qlwc - ql) / dt
        condi = (qiwc - qi) / dt

        # calculate precipitation overlap
        # simple form based on Maximum Overlap
        tmp_covptot[0, 0] = max(tmp_covptot, out_clc)
        covpclr = max(tmp_covptot - out_clc, 0.0)

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
        if out_clc[0, 0, 0] > ZEPS2:
            if LEVAPLS2 or LDRAIN1D:
                lcrit = 1.9 * RCLCRIT
            else:
                lcrit = 2.0 * RCLCRIT
            cldl = qlwc / out_clc
            dl = ckcodtl * (1 - exp(-((cldl / lcrit) ** 2)))
            prr = qlwc - out_clc * cldl * exp(-dl)
            qlwc -= prr
        else:
            prr = 0.0

        # diagnostic calculation of snow production from cloud ice
        if out_clc > ZEPS2:
            if LEVAPLS2 or LDRAIN1D:
                icrit = 0.0001
            else:
                icrit = 2 * RCLCRIT
            cldi = qiwc / out_clc
            di = ckcodti * exp(0.025 * (t - RTT)) * (1 - exp(-((cldi / icrit) ** 2)))
            prs = qiwc - out_clc * cldi * exp(-di)
            qiwc -= prs
        else:
            prs = 0.0

        # new precipitation (rain + snow)
        dr = cons2 * dp * (prr + prs)

        # rain fraction (different from cloud liquid water fraction!)
        if t < RTT:
            rfreeze = cons2 * dp * prr
            fwatr = 0.0
        else:
            rfreeze = 0.0
            fwatr = 1.0
        rfln += fwatr * dr
        sfln += (1 - fwatr) * dr

        # precipitation evaporation
        prtot = rfln + sfln
        if prtot > ZEPS2 and covpclr > ZEPS2 and (LEVAPLS2 or LDRAIN1D):
            preclr = prtot * covpclr / tmp_covptot

            # this is the humidity in the moisest zcovpclr region
            qe = in_qsat - (in_qsat - qlim) * covpclr / ((1 - out_clc) ** 2)
            beta = RG * RPECONS * (sqrt(in_ap / tmp_aph_s) / 0.00509 * preclr / covpclr) ** 0.5777

            # implicit solution
            b = dt * beta * (in_qsat - qe) / (1 + dt * beta * corqs)

            dtgdp = dt * RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
            dpr = min(covpclr * b / dtgdp, preclr)
            preclr -= dpr
            if preclr <= 0:
                tmp_covptot[0, 0] = out_clc
            out_covptot[0, 0, 0] = tmp_covptot

            # warm proportion
            evapr = dpr * rfln / prtot
            rfln -= evapr

            # ice proportion
            evaps = dpr * sfln / prtot
            sfln -= evaps
        else:
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
        t, q = cuadjtqs_nl(in_ap, t, q)

        # update rain fraction and freezing
        dq = max(qold - q, 0.0)
        dr2 = cons2 * dp * dq
        if t < RTT:
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
        out_tnd_q[0, 0, 0] = -(condl + condi) + (in_lude + evapr + evaps) * gdp
        out_tnd_t[0, 0, 0] = (
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
        out_tnd_ql[0, 0, 0] = (qlwc - ql) / dt
        out_tnd_qi[0, 0, 0] = (qiwc - qi) / dt

        # these fluxes will later be shifted one level downward
        fplsl = rfln
        fplsn = sfln

        # record rain flux for next level
        tmp_rfl[0, 0] = rfln
        tmp_sfl[0, 0] = sfln

    # enthalpy fluxes due to precipitation
    with computation(FORWARD):
        with interval(0, 1):
            out_fhpsl[0, 0, 0] = 0.0
            out_fhpsn[0, 0, 0] = 0.0
        with interval(1, None):
            out_fplsl[0, 0, 0] = fplsl[0, 0, -1]
            out_fplsn[0, 0, 0] = fplsn[0, 0, -1]
            out_fhpsl[0, 0, 0] = -out_fplsl * RLVTT
            out_fhpsn[0, 0, 0] = -out_fplsn * RLSTT
