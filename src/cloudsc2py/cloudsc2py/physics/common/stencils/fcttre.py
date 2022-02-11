# -*- coding: utf-8 -*-
from gt4py import gtscript

from cloudsc2py.framework.stencil import function_collection
from cloudsc2py.utils.f2py import ported_function


@staticmethod
@ported_function(
    from_file="common/include/fcttre.func.h", from_line=73, to_line=75
)
@function_collection("foealfa", ["RTICE", "RTWAT", "RTWAT_RTICE_R"])
@gtscript.function
def foealfa(t):
    from __externals__ import ext

    return min(
        1.0,
        ((max(ext.RTICE, min(ext.RTWAT, t)) - ext.RTICE) * ext.RTWAT_RTICE_R)
        ** 2,
    )


@staticmethod
@ported_function(
    from_file="common/include/fcttre.func.h", from_line=121, to_line=123
)
@function_collection("foealfcu", ["RTICECU", "RTWAT", "RTWAT_RTICECU_R"])
@gtscript.function
def foealfcu(t):
    from __externals__ import ext

    return min(
        1.0,
        (
            (max(ext.RTICECU, min(ext.RTWAT, t)) - ext.RTICECU)
            * ext.RTWAT_RTICECU_R
        )
        ** 2,
    )


@staticmethod
@ported_function(
    from_file="common/include/fcttre.func.h", from_line=80, to_line=83
)
@function_collection(
    "foeewm", ["R2ES", "R3IES", "R3LES", "R4IES", "R4LES", "RTT", "foealfa"]
)
@gtscript.function
def foeewm(t):
    from __externals__ import ext

    return ext.R2ES * (
        ext.foealfa(t) * exp(ext.R3LES * (t - ext.RTT) / (t - ext.R4LES))
        + (1 - ext.foealfa(t))
        * (exp(ext.R3IES * (t - ext.RTT) / (t - ext.R4IES)))
    )


@staticmethod
@ported_function(
    from_file="common/include/fcttre.func.h", from_line=129, to_line=131
)
@function_collection(
    "foeewmcu",
    ["R2ES", "R3IES", "R3LES", "R4IES", "R4LES", "RTT", "foealfcu"],
)
@gtscript.function
def foeewmcu(t):
    from __externals__ import ext

    return ext.R2ES * (
        ext.foealfcu(t) * exp(ext.R3LES * (t - ext.RTT) / (t - ext.R4LES))
        + (1 - ext.foealfcu(t))
        * (exp(ext.R3IES * (t - ext.RTT) / (t - ext.R4IES)))
    )
