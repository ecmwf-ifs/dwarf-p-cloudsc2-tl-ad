# -*- coding: utf-8 -*-
from gt4py import gtscript

from cloudsc2py.framework.stencil import stencil_collection
from cloudsc2py.utils.f2py import ported_function


@ported_function(from_file="clouds2_nl/satur.F90", from_line=106, to_line=140)
@stencil_collection(
    "saturation",
    external_names=[
        "KFLAG",
        "LPHYLIN",
        "QMAX",
        "R2ES",
        "R3IES",
        "R3LES",
        "R4IES",
        "R4LES",
        "RETV",
        "RTT",
        "foealfa",
        "foeewm",
        "foeewmcu",
    ],
)
def saturation_def(
    in_ap: gtscript.Field["ftype"],
    in_t: gtscript.Field["ftype"],
    out_qsat: gtscript.Field["ftype"],
):
    from __externals__ import ext

    with computation(PARALLEL), interval(...):
        if __INLINED(ext.LPHYLIN):
            alfa = ext.foealfa(in_t)
            foeewl = ext.R2ES * exp(
                ext.R3LES * (in_t - ext.RTT) / (in_t - ext.R4LES)
            )
            foeewi = ext.R2ES * exp(
                ext.R3IES * (in_t - ext.RTT) / (in_t - ext.R4IES)
            )
            foeew = alfa * foeewl + (1 - alfa) * foeewi
            qs = min(foeew / in_ap, ext.QMAX)
        else:
            if __INLINED(ext.KFLAG == 1):
                ew = ext.foeewmcu(in_t)
            else:
                ew = ext.foeewm(in_t)
            qs = min(ew / in_ap, ext.QMAX)
        out_qsat = qs / (1.0 - ext.RETV * qs)
