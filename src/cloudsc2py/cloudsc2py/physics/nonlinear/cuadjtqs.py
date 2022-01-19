# -*- coding: utf-8 -*-
from gt4py import gtscript

from cloudsc2py.framework.stencil import function_collection
from cloudsc2py.utils.f2py import ported_function


@staticmethod
@ported_function(
    from_file="cloudsc2_nl/cloudsc2.F90", from_line=644, to_line=656
)
@function_collection("cuadjtqs_nl_0", ["R2ES", "RETV", "ZQMAX"])
@gtscript.function
def cuadjtqs_nl_0(ap, t, q, z3es, z4es, z5alcp, zaldcp):
    from __externals__ import ext

    foeew = ext.R2ES * exp(z3es * (t - ext.RTT) / (t - z4es))
    qsat = min(foeew / ap, ext.ZQMAX)
    cor = 1 / (1 - ext.RETV * qsat)
    qsat *= cor
    z2s = z5alcp / (t - z4es) ** 2
    cond = (q - qsat) / (1 + qsat * cor * z2s)
    t += zaldcp * cond
    q -= cond
    return t, q


@staticmethod
@ported_function(
    from_file="cloudsc2_nl/cloudsc2.F90", from_line=630, to_line=669
)
@function_collection(
    "cuadjtqs_nl",
    [
        "ICALL",
        "R3IES",
        "R3LES",
        "R4IES",
        "R4LES",
        "R5ALSCP",
        "R5ALVCP",
        "RALSDCP",
        "RALVDCP",
        "cuadjtqs_nl_0",
    ],
)
@gtscript.function
def cuadjtqs_nl(ap, t, q):
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
        t, q = ext.cuadjtqs_nl_0(ap, t, q, z3es, z4es, z5alcp, zaldcp)
        t, q = ext.cuadjtqs_nl_0(ap, t, q, z3es, z4es, z5alcp, zaldcp)
        return t, q
