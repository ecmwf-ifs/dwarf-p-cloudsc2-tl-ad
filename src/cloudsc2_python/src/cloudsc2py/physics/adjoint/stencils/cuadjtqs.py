# -*- coding: utf-8 -*-
from gt4py.cartesian import gtscript

from ifs_physics_common.framework.stencil import function_collection
from ifs_physics_common.utils.f2py import ported_function


@ported_function(from_file="cloudsc2_ad/cuadjtqsad.F90", from_line=10, to_line=871)
@function_collection("cuadjtqs_ad")
@gtscript.function
def cuadjtqs_ad(ap, ap_i, t, t_i, q, q_i):
    from __externals__ import (
        ICALL,
        R2ES,
        R3IES,
        R3LES,
        R4IES,
        R4LES,
        R5ALSCP,
        R5ALVCP,
        RALSDCP,
        RALVDCP,
        RETV,
        RTT,
        ZQMAX,
    )

    if t > RTT:
        z3es = R3LES
        z4es = R4LES
        z5alcp = R5ALVCP
        zaldcp = RALVDCP
    else:
        z3es = R3IES
        z4es = R4IES
        z5alcp = R5ALSCP
        zaldcp = RALSDCP

    if ICALL == 0:
        targ = t + 0
        foeew = R2ES * exp(z3es * (targ - RTT) / (targ - z4es))
        foeew_b = foeew + 0
        qsat = foeew / ap
        ltest2 = qsat > ZQMAX
        if ltest2:
            qsat = ZQMAX
        cor = 1 / (1 - RETV * qsat)
        qsat_d = qsat + 0
        qsat *= cor
        targ_b = targ + 0
        z2s = z5alcp / (targ - z4es) ** 2
        qsat_b = qsat + 0
        cor_b = cor + 0
        z2s_b = z2s + 0
        q_b = q + 0
        cond1 = (q - qsat) / (1 + qsat * cor * z2s)
        t += zaldcp * cond1
        q -= cond1

        targ = t + 0
        foeew = R2ES * exp(z3es * (targ - RTT) / (targ - z4es))
        foeew_a = foeew + 0
        qsat = foeew / ap
        ltest1 = qsat > ZQMAX
        if ltest1:
            qsat = ZQMAX
        cor = 1 / (1 - RETV * qsat)
        qsat_c = qsat + 0
        qsat *= cor
        targ_a = targ + 0
        z2s = z5alcp / (targ - z4es) ** 2
        qsat_a = qsat + 0
        cor_a = cor + 0
        z2s_a = z2s + 0
        q_a = q + 0
        cond1 = (q - qsat) / (1 + qsat * cor * z2s)
        t += zaldcp * cond1
        q -= cond1

        cond1_i = -q_i + zaldcp * t_i
        qsat = qsat_a + 0
        cor = cor_a + 0
        z2s = z2s_a + 0
        q_i += cond1_i / (1 + qsat * cor * z2s)
        qsat_i = (
            -cond1_i / (1 + qsat * cor * z2s)
            - cond1_i * (q_a - qsat) * cor * z2s / (1 + qsat * cor * z2s) ** 2
        )
        cor_i = -cond1_i * (q_a - qsat) * qsat * z2s / (1 + qsat * cor * z2s) ** 2
        z2s_i = -cond1_i * (q_a - qsat) * qsat * cor / (1 + qsat * cor * z2s) ** 2
        targ = targ_a + 0
        targ_i = -2 * z2s_i * z5alcp / (targ - z4es) ** 3
        qsat = qsat_c + 0
        cor_i += qsat_i * qsat
        qsat_i *= cor
        qsat_i += cor_i * RETV / (1 - RETV * qsat) ** 2
        if ltest1:
            qsat_i = 0.0
        foeew_i = qsat_i / ap
        foeew = foeew_a + 0
        qp_i = qsat_i * foeew
        targ_i += (
            foeew_i
            * R2ES
            * z3es
            * (RTT - z4es)
            * exp(z3es * (targ - RTT) / (targ - z4es))
            / (targ - z4es) ** 2
        )
        t_i += targ_i

        cond1_i = -q_i + zaldcp * t_i
        qsat = qsat_b + 0
        cor = cor_b + 0
        z2s = z2s_b + 0
        q_i += cond1_i / (1 + qsat * cor * z2s)
        qsat_i = (
            -cond1_i / (1 + qsat * cor * z2s)
            - cond1_i * (q_b - qsat) * cor * z2s / (1 + qsat * cor * z2s) ** 2
        )
        cor_i = -cond1_i * (q_b - qsat) * qsat * z2s / (1 + qsat * cor * z2s) ** 2
        z2s_i = -cond1_i * (q_b - qsat) * qsat * cor / (1 + qsat * cor * z2s) ** 2
        targ = targ_b + 0
        targ_i = -2 * z2s_i * z5alcp / (targ - z4es) ** 3
        qsat = qsat_d + 0
        cor_i += qsat_i * qsat
        qsat_i *= cor
        qsat_i += cor_i * RETV / (1 - RETV * qsat) ** 2
        if ltest2:
            qsat_i = 0.0
        foeew_i = qsat_i / ap
        foeew = foeew_b + 0
        qp_i += qsat_i * foeew
        targ_i += (
            foeew_i
            * R2ES
            * z3es
            * (RTT - z4es)
            * exp(z3es * (targ - RTT) / (targ - z4es))
            / (targ - z4es) ** 2
        )
        t_i += targ_i
        ap_i -= qp_i / ap**2

        return ap_i, t, t_i, q, q_i
