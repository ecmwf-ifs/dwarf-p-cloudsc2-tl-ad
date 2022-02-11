# -*- coding: utf-8 -*-
from gt4py import gtscript

from cloudsc2py.framework.stencil import function_collection
from cloudsc2py.utils.f2py import ported_function


@ported_function(
    from_file="cloudsc2_tl/cuadjtqstl.F90", from_line=337, to_line=369
)
@function_collection("cuadjtqs_tl_0", external_names=["R2ES", "RETV", "ZQMAX"])
@gtscript.function
def cuadjtqs_tl_0(ap, ap_i, t, t_i, q, q_i, z3es, z4es, z5alcp, zaldcp):
    from __externals__ import ext

    foeew = ext.R2ES * exp(z3es * (t - ext.RTT) / (t - z4es))
    foeew_i = foeew * z3es * t_i * (ext.RTT - z4es) / (t - z4es) ** 2
    qsat = foeew / ap
    qsat_i = foeew_i / ap - foeew * ap_i / ap ** 2
    if qsat > ext.ZQMAX:
        qsat = ext.ZQMAX
        qsat_i = 0.0

    cor = 1 / (1 - ext.RETV * qsat)
    cor_i = ext.RETV * qsat_i / (1 - ext.RETV * qsat) ** 2
    qsat_i = qsat_i * cor + qsat * cor_i
    qsat *= cor
    z2s = z5alcp / (t - z4es) ** 2
    z2s_i = -2 * z5alcp * t_i / (t - z4es) ** 3
    cond = (q - qsat) / (1 + qsat * cor * z2s)
    cond_i = (q_i - qsat_i) / (1 + qsat * cor * z2s) - (q - qsat) * (
        qsat_i * cor * z2s + qsat * cor_i * z2s + qsat * cor * z2s_i
    ) / (1 + qsat * cor * z2s) ** 2
    t += zaldcp * cond
    t_i += zaldcp * cond_i
    q -= cond
    q_i -= cond_i

    return t, t_i, q, q_i


@ported_function(
    from_file="cloudsc2_tl/cuadjtqstl.F90", from_line=10, to_line=482
)
@function_collection(
    "cuadjtqs_tl",
    external_names=[
        "ICALL",
        "R3IES",
        "R3LES",
        "R4IES",
        "R4LES",
        "R5ALSCP",
        "R5ALVCP",
        "RALSDCP",
        "RALVDCP",
        "cuadjtqs_tl_0",
    ],
)
@gtscript.function
def cuadjtqs_tl(ap, ap_i, t, t_i, q, q_i):
    from __externals__ import ext

    if t > ext.RTT:
        z3es = ext.R3LES
        z4es = ext.R4LES
        z5alcp = ext.R5ALVCP
        zaldcp = ext.RALVDCP
    else:
        z3es = ext.R3IES
        z4es = ext.R4IES
        z5alcp = ext.R5ALSCP
        zaldcp = ext.RALSDCP

    if __INLINED(ext.ICALL == 0):
        t, t_i, q, q_i = ext.cuadjtqs_tl_0(
            ap, ap_i, t, t_i, q, q_i, z3es, z4es, z5alcp, zaldcp
        )
        t, t_i, q, q_i = ext.cuadjtqs_tl_0(
            ap, ap_i, t, t_i, q, q_i, z3es, z4es, z5alcp, zaldcp
        )
        return t, t_i, q, q_i
