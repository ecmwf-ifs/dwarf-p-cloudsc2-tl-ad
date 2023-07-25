# -*- coding: utf-8 -*-
from gt4py.cartesian import gtscript

from cloudsc2py.physics.common.stencils.fcttre import foealfa, foeewm, foeewmcu
from ifs_physics_common.framework.stencil import stencil_collection
from ifs_physics_common.utils.f2py import ported_function


@ported_function(from_file="clouds2_nl/satur.F90", from_line=106, to_line=140)
@stencil_collection("saturation")
def saturation_def(
    in_ap: gtscript.Field["float"], in_t: gtscript.Field["float"], out_qsat: gtscript.Field["float"]
):
    from __externals__ import KFLAG, LPHYLIN, QMAX, R2ES, R3IES, R3LES, R4IES, R4LES, RETV, RTT

    with computation(PARALLEL), interval(...):
        if LPHYLIN:
            alfa = foealfa(in_t)
            foeewl = R2ES * exp(R3LES * (in_t - RTT) / (in_t - R4LES))
            foeewi = R2ES * exp(R3IES * (in_t - RTT) / (in_t - R4IES))
            foeew = alfa * foeewl + (1 - alfa) * foeewi
            qs = min(foeew / in_ap, QMAX)
        else:
            if KFLAG == 1:
                ew = foeewmcu(in_t)
            else:
                ew = foeewm(in_t)
            qs = min(ew / in_ap, QMAX)
        out_qsat[0, 0, 0] = qs / (1.0 - RETV * qs)
