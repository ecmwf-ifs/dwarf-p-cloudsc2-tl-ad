# -*- coding: utf-8 -*-
from gt4py.cartesian import gtscript

from ifs_physics_common.framework.stencil import function_collection
from ifs_physics_common.utils.f2py import ported_function


@ported_function(from_file="common/include/fcttre.func.h", from_line=73, to_line=75)
@function_collection("foealfa")
@gtscript.function
def foealfa(t):
    from __externals__ import RTICE, RTWAT, RTWAT_RTICE_R

    return min(1.0, ((max(RTICE, min(RTWAT, t)) - RTICE) * RTWAT_RTICE_R) ** 2)


@ported_function(from_file="common/include/fcttre.func.h", from_line=121, to_line=123)
@function_collection("foealfcu")
@gtscript.function
def foealfcu(t):
    from __externals__ import RTICECU, RTWAT, RTWAT_RTICECU_R

    return min(1.0, ((max(RTICECU, min(RTWAT, t)) - RTICECU) * RTWAT_RTICECU_R) ** 2)


@ported_function(from_file="common/include/fcttre.func.h", from_line=80, to_line=83)
@function_collection("foeewm")
@gtscript.function
def foeewm(t):
    from __externals__ import R2ES, R3IES, R3LES, R4IES, R4LES, RTT

    return R2ES * (
        foealfa(t) * exp(R3LES * (t - RTT) / (t - R4LES))
        + (1 - foealfa(t)) * (exp(R3IES * (t - RTT) / (t - R4IES)))
    )


@ported_function(from_file="common/include/fcttre.func.h", from_line=129, to_line=131)
@function_collection("foeewmcu")
@gtscript.function
def foeewmcu(t):
    from __externals__ import R2ES, R3IES, R3LES, R4IES, R4LES, RTT

    return R2ES * (
        foealfcu(t) * exp(R3LES * (t - RTT) / (t - R4LES))
        + (1 - foealfcu(t)) * (exp(R3IES * (t - RTT) / (t - R4IES)))
    )
