# -*- coding: utf-8 -*-
from gt4py.cartesian import gtscript

from cloudsc2py.framework.stencil import function_collection
from cloudsc2py.utils.f2py import ported_function


@ported_function(from_file="cloudsc2_ad/cuadjtqsad.F90", from_line=318, to_line=340)
@function_collection("cuadjtqs_ad_00")
@gtscript.function
def cuadjtqs_ad_00(ap, t, q, z3es, z4es, z5alcp, zaldcp):
    from __externals__ import R2ES, RETV, RTT, ZQMAX

    foeew = R2ES * exp(z3es * (t - RTT) / (t - z4es))
    qsat = foeew / ap
    test = qsat > ZQMAX
    qsat = min(qsat, ZQMAX)
    cor = 1 / (1 - RETV * qsat)
    qsat_ = qsat + 0
    qsat *= cor
    z2s = z5alcp / (t - z4es) ** 2
    cond = (q - qsat) / (1 + qsat * cor * z2s)
    t += zaldcp * cond
    q -= cond

    return t, q, cor, foeew, qsat, qsat_, test, z2s


@ported_function(from_file="cloudsc2_as/cuadjtqsad.F90", from_line=546, to_line=590)
@function_collection("cuadjtqs_ad_01")
@gtscript.function
def cuadjtqs_ad_01(
    ap,
    ap_i,
    t_i,
    q_i,
    z3es,
    z4es,
    z5alcp,
    zaldcp,
    cor_a,
    cor_b,
    foeew_a,
    foeew_b,
    q_a,
    q_b,
    qsat_a,
    qsat_b,
    qsat_c,
    qsat_d,
    targ_a,
    targ_b,
    test1,
    test2,
    z2s_a,
    z2s_b,
):
    from __externals__ import R2ES, RETV, RTT

    cond = zaldcp * t_i - q_i
    q_i += cond / (1 + qsat_a * cor_a * z2s_a)
    qsat_i = (
        -cond / (1 + qsat_a * cor_a * z2s_a)
        - cond * (q_a - qsat_a) * cor_a * z2s_a / (1 + qsat_a * cor_a * z2s_a) ** 2
    )
    cor_i = -cond * (q_a - qsat_a) * qsat_a * z2s_a / (1 + qsat_a * cor_a * z2s_a) ** 2
    z2s_i = -cond * (q_a - qsat_a) * qsat_a * cor_a / (1 + qsat_a * cor_a * z2s_a) ** 2
    targ_i = -2 * z2s_i * z5alcp / (targ_a - z4es) ** 3
    cor_i += qsat_i * qsat_c
    qsat_i *= cor_a
    qsat_i += cor_i * RETV / (1 - RETV * qsat_c) ** 2
    if test1:
        qsat_i = 0.0
    foeew_i = qsat_i / ap
    qp = qsat_i * foeew_a
    targ_i += (
        foeew_i
        * R2ES
        * z3es
        * (RTT - z4es)
        * exp(z3es * (targ_a - RTT) / (targ_a - z4es))
        / (targ_a - z4es) ** 2
    )
    t_i += targ_i

    cond = zaldcp * t_i - q_i
    q_i += cond / (1 + qsat_b * cor_b * z2s_b)
    qsat_i = (
        -cond / (1 + qsat_b * cor_b * z2s_b)
        - cond * (q_b - qsat_b) * cor_b * z2s_b / (1 + qsat_b * cor_b * z2s_b) ** 2
    )
    cor_i = -cond * (q_b - qsat_b) * qsat_b * z2s_b / (1 + qsat_b * cor_b * z2s_b) ** 2
    z2s_i = -cond * (q_b - qsat_b) * qsat_b * cor_b / (1 + qsat_b * cor_b * z2s_b) ** 2
    targ_i = -2 * z2s_i * z5alcp / (targ_b - z4es) ** 3
    cor_i += qsat_i * qsat_d
    qsat_i *= cor_b
    qsat_i += cor_i * RETV / (1 - RETV * qsat_d) ** 2
    if test2:
        qsat_i = 0.0
    foeew_i = qsat_i / ap
    qp += qsat_i * foeew_b
    targ_i += (
        foeew_i
        * R2ES
        * z3es
        * (RTT - z4es)
        * exp(z3es * (targ_b - RTT) / (targ_b - z4es))
        / (targ_b - z4es) ** 2
    )
    t_i += targ_i
    ap_i -= qp / ap**2

    return ap_i, t_i, q_i


@ported_function(from_file="cloudsc2_ad/cuadjtqsad.F90", from_line=10, to_line=871)
@function_collection("cuadjtqs_ad")
@gtscript.function
def cuadjtqs_ad(ap, ap_i, t, t_i, q, q_i):
    from __externals__ import (
        ICALL,
        R3IES,
        R3LES,
        R4IES,
        R4LES,
        R5ALSCP,
        R5ALVCP,
        RALSDCP,
        RALVDCP,
        RTT,
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

    if __INLINED(ICALL == 0):
        targ_b = t + 0
        q_b = q + 0
        (
            targ_a,
            q_a,
            cor_b,
            foeew_b,
            qsat_b,
            qsat_d,
            test2,
            z2s_b,
        ) = cuadjtqs_ad_00(ap, t, q, z3es, z4es, z5alcp, zaldcp)
        (
            t,
            q,
            cor_a,
            foeew_a,
            qsat_a,
            qsat_c,
            test1,
            z2s_a,
        ) = cuadjtqs_ad_00(ap, targ_a, q_a, z3es, z4es, z5alcp, zaldcp)
        ap_i, t_i, q_i = cuadjtqs_ad_01(
            ap,
            ap_i,
            t_i,
            q_i,
            z3es,
            z4es,
            z5alcp,
            zaldcp,
            cor_a,
            cor_b,
            foeew_a,
            foeew_b,
            q_a,
            q_b,
            qsat_a,
            qsat_b,
            qsat_c,
            qsat_d,
            targ_a,
            targ_b,
            test1,
            test2,
            z2s_a,
            z2s_b,
        )
        return ap_i, t, t_i, q, q_i
