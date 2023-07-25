# -*- coding: utf-8 -*-
from gt4py.cartesian import gtscript

from ifs_physics_common.framework.stencil import function_collection
from ifs_physics_common.utils.f2py import ported_function


@ported_function(from_file="cloudsc2_nl/cloudsc2.F90", from_line=644, to_line=656)
@function_collection("cuadjtqs_nl_0")
@gtscript.function
def cuadjtqs_nl_0(ap, t, q, z3es, z4es, z5alcp, zaldcp):
    from __externals__ import R2ES, RETV, RTT, ZQMAX

    foeew = R2ES * exp(z3es * (t - RTT) / (t - z4es))
    qsat = min(foeew / ap, ZQMAX)
    cor = 1 / (1 - RETV * qsat)
    qsat *= cor
    z2s = z5alcp / (t - z4es) ** 2
    cond = (q - qsat) / (1 + qsat * cor * z2s)
    t += zaldcp * cond
    q -= cond
    return t, q


@ported_function(from_file="cloudsc2_nl/cloudsc2.F90", from_line=630, to_line=669)
@function_collection("cuadjtqs_nl")
@gtscript.function
def cuadjtqs_nl(ap, t, q):
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

    if ICALL == 0:
        t, q = cuadjtqs_nl_0(ap, t, q, z3es, z4es, z5alcp, zaldcp)
        t, q = cuadjtqs_nl_0(ap, t, q, z3es, z4es, z5alcp, zaldcp)
        return t, q
