# -*- coding: utf-8 -*-
from gt4py import gtscript

from cloudsc2py.framework.stencil import stencil_collection
from cloudsc2py.physics.nonlinear.stencils.cloudsc2_ng_0 import cloudsc2_nl_0
from src.cloudsc2_python.utils import ported_function


@ported_function(from_file="cloudsc2_nl/cloudsc2.F90", from_line=235, to_line=735)
@stencil_collection("cloudsc2_nl_ng")
def cloudsc2_nl_ng_def(
    in_eta: gtscript.Field[gtscript.K, "float"],
    in_ap: gtscript.Field["float"],
    in_aph: gtscript.Field["float"],
    in_t: gtscript.Field["float"],
    in_q: gtscript.Field["float"],
    in_qsat: gtscript.Field["float"],
    in_ql: gtscript.Field["float"],
    in_qi: gtscript.Field["float"],
    in_lu: gtscript.Field["float"],
    in_lude: gtscript.Field["float"],
    in_mfd: gtscript.Field["float"],
    in_mfu: gtscript.Field["float"],
    in_supsat: gtscript.Field["float"],
    in_tnd_cml_t: gtscript.Field["float"],
    in_tnd_cml_q: gtscript.Field["float"],
    in_tnd_cml_ql: gtscript.Field["float"],
    in_tnd_cml_qi: gtscript.Field["float"],
    tmp_aph_s: gtscript.Field[gtscript.IJ, "float"],
    tmp_trpaus: gtscript.Field[gtscript.IJ, "float"],
    out_tnd_t: gtscript.Field["float"],
    out_tnd_q: gtscript.Field["float"],
    out_tnd_ql: gtscript.Field["float"],
    out_tnd_qi: gtscript.Field["float"],
    out_clc: gtscript.Field["float"],
    out_fhpsl: gtscript.Field["float"],
    out_fhpsn: gtscript.Field["float"],
    out_fplsl: gtscript.Field["float"],
    out_fplsn: gtscript.Field["float"],
    out_covptot: gtscript.Field["float"],
    *,
    dt: "float",
):
    from __externals__ import RLSTT, RLVTT

    with computation(PARALLEL), interval(0, -1):
        # first guess values for T
        t = in_t + dt * in_tnd_cml_t

    # eta value at tropopause
    with computation(FORWARD), interval(0, 1):
        tmp_trpaus = 0.1
    with computation(FORWARD), interval(0, -2):
        if in_eta[0] > 0.1 and in_eta[0] < 0.4 and t[0, 0, 0] > t[0, 0, 1]:
            tmp_trpaus = in_eta[0]

    with computation(FORWARD):
        with interval(0, 1):
            (
                out_tnd_t,
                out_tnd_q,
                out_tnd_ql,
                out_tnd_qi,
                out_clc,
                out_covptot,
                covptot,
                rfl,
                sfl,
            ) = cloudsc2_nl_0(
                in_eta,
                in_ap,
                in_aph,
                in_q,
                in_qsat,
                in_ql,
                in_qi,
                in_lu,
                in_lude,
                in_mfd,
                in_mfu,
                in_supsat,
                in_tnd_cml_q,
                in_tnd_cml_ql,
                in_tnd_cml_qi,
                tmp_aph_s,
                tmp_trpaus,
                dt,
                covptot=0.0,
                rfl=0.0,
                sfl=0.0,
                t=t,
            )
            out_fhpsl = 0.0
            out_fhpsn = 0.0
        with interval(1, -1):
            (
                out_tnd_t,
                out_tnd_q,
                out_tnd_ql,
                out_tnd_qi,
                out_clc,
                out_covptot,
                covptot,
                rfl,
                sfl,
            ) = cloudsc2_nl_0(
                in_eta,
                in_ap,
                in_aph,
                in_q,
                in_qsat,
                in_ql,
                in_qi,
                in_lu,
                in_lude,
                in_mfd,
                in_mfu,
                in_supsat,
                in_tnd_cml_q,
                in_tnd_cml_ql,
                in_tnd_cml_qi,
                tmp_aph_s,
                tmp_trpaus,
                dt,
                covptot=covptot[0, 0, -1],
                rfl=rfl[0, 0, -1],
                sfl=sfl[0, 0, -1],
                t=t,
            )
            out_fplsl = rfl[0, 0, -1]
            out_fplsn = sfl[0, 0, -1]
            out_fhpsl = -out_fplsl * RLVTT
            out_fhpsn = -out_fplsn * RLSTT
        with interval(-1, None):
            out_fplsl = rfl[0, 0, -1]
            out_fplsn = sfl[0, 0, -1]
            out_fhpsl = -out_fplsl * RLVTT
            out_fhpsn = -out_fplsn * RLSTT
