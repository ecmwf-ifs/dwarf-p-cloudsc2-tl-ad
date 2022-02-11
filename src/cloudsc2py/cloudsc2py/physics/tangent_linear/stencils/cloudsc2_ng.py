# -*- coding: utf-8 -*-
from gt4py import gtscript

from cloudsc2py.framework.stencil import stencil_collection
from cloudsc2py.utils.f2py import ported_function


@ported_function(
    from_file="cloudsc2_tl/cloudsc2tl.F90", from_line=321, to_line=1113
)
@stencil_collection(
    "cloudsc2_tl_ng", external_names=["RLSTT", "RLVTT", "cloudsc2_tl_0"]
)
def cloudsc2_tl_ng_def(
    in_eta: gtscript.Field[gtscript.K, "ftype"],
    in_ap: gtscript.Field["ftype"],
    in_ap_i: gtscript.Field["ftype"],
    in_aph: gtscript.Field["ftype"],
    in_aph_i: gtscript.Field["ftype"],
    in_t: gtscript.Field["ftype"],
    in_t_i: gtscript.Field["ftype"],
    in_q: gtscript.Field["ftype"],
    in_q_i: gtscript.Field["ftype"],
    in_qsat: gtscript.Field["ftype"],
    in_qsat_i: gtscript.Field["ftype"],
    in_ql: gtscript.Field["ftype"],
    in_ql_i: gtscript.Field["ftype"],
    in_qi: gtscript.Field["ftype"],
    in_qi_i: gtscript.Field["ftype"],
    in_lu: gtscript.Field["ftype"],
    in_lu_i: gtscript.Field["ftype"],
    in_lude: gtscript.Field["ftype"],
    in_lude_i: gtscript.Field["ftype"],
    in_mfd: gtscript.Field["ftype"],
    in_mfd_i: gtscript.Field["ftype"],
    in_mfu: gtscript.Field["ftype"],
    in_mfu_i: gtscript.Field["ftype"],
    in_supsat: gtscript.Field["ftype"],
    in_supsat_i: gtscript.Field["ftype"],
    in_tnd_cml_t: gtscript.Field["ftype"],
    in_tnd_cml_t_i: gtscript.Field["ftype"],
    in_tnd_cml_q: gtscript.Field["ftype"],
    in_tnd_cml_q_i: gtscript.Field["ftype"],
    in_tnd_cml_ql: gtscript.Field["ftype"],
    in_tnd_cml_ql_i: gtscript.Field["ftype"],
    in_tnd_cml_qi: gtscript.Field["ftype"],
    in_tnd_cml_qi_i: gtscript.Field["ftype"],
    tmp_trpaus: gtscript.Field[gtscript.IJ, "ftype"],
    out_tnd_t: gtscript.Field["ftype"],
    out_tnd_t_i: gtscript.Field["ftype"],
    out_tnd_q: gtscript.Field["ftype"],
    out_tnd_q_i: gtscript.Field["ftype"],
    out_tnd_ql: gtscript.Field["ftype"],
    out_tnd_ql_i: gtscript.Field["ftype"],
    out_tnd_qi: gtscript.Field["ftype"],
    out_tnd_qi_i: gtscript.Field["ftype"],
    out_clc: gtscript.Field["ftype"],
    out_clc_i: gtscript.Field["ftype"],
    out_fhpsl: gtscript.Field["ftype"],
    out_fhpsl_i: gtscript.Field["ftype"],
    out_fhpsn: gtscript.Field["ftype"],
    out_fhpsn_i: gtscript.Field["ftype"],
    out_fplsl: gtscript.Field["ftype"],
    out_fplsl_i: gtscript.Field["ftype"],
    out_fplsn: gtscript.Field["ftype"],
    out_fplsn_i: gtscript.Field["ftype"],
    out_covptot: gtscript.Field["ftype"],
    out_covptot_i: gtscript.Field["ftype"],
    *,
    dt: "ftype",
):
    from __externals__ import ext

    with computation(PARALLEL), interval(0, -1):
        # first guess values for T
        t = in_t + dt * in_tnd_cml_t
        t_i = in_t_i + dt * in_tnd_cml_t_i

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
                out_tnd_t_i,
                out_tnd_q,
                out_tnd_q_i,
                out_tnd_ql,
                out_tnd_ql_i,
                out_tnd_qi,
                out_tnd_qi_i,
                out_clc,
                out_clc_i,
                out_covptot,
                out_covptot_i,
                covptot,
                covptot_i,
                rfl,
                rfl_i,
                sfl,
                sfl_i,
            ) = ext.cloudsc2_tl_0(
                in_eta,
                in_ap,
                in_ap_i,
                in_aph,
                in_aph_i,
                in_q,
                in_q_i,
                in_qsat,
                in_qsat_i,
                in_ql,
                in_ql_i,
                in_qi,
                in_qi_i,
                in_lu,
                in_lu_i,
                in_lude,
                in_lude_i,
                in_mfd,
                in_mfd_i,
                in_mfu,
                in_mfu_i,
                in_supsat,
                in_supsat_i,
                in_tnd_cml_q,
                in_tnd_cml_q_i,
                in_tnd_cml_ql,
                in_tnd_cml_ql_i,
                in_tnd_cml_qi,
                in_tnd_cml_qi_i,
                tmp_trpaus,
                dt,
                covptot=0.0,
                covptot_i=0.0,
                rfl=0.0,
                rfl_i=0.0,
                sfl=0.0,
                sfl_i=0.0,
                t=t,
                t_i=t_i,
            )
            out_fhpsl = 0.0
            out_fhpsl_i = 0.0
            out_fhpsn = 0.0
            out_fhpsn_i = 0.0
        with interval(1, -1):
            (
                out_tnd_t,
                out_tnd_t_i,
                out_tnd_q,
                out_tnd_q_i,
                out_tnd_ql,
                out_tnd_ql_i,
                out_tnd_qi,
                out_tnd_qi_i,
                out_clc,
                out_clc_i,
                out_covptot,
                out_covptot_i,
                covptot,
                covptot_i,
                rfl,
                rfl_i,
                sfl,
                sfl_i,
            ) = ext.cloudsc2_tl_0(
                in_eta,
                in_ap,
                in_ap_i,
                in_aph,
                in_aph_i,
                in_q,
                in_q_i,
                in_qsat,
                in_qsat_i,
                in_ql,
                in_ql_i,
                in_qi,
                in_qi_i,
                in_lu,
                in_lu_i,
                in_lude,
                in_lude_i,
                in_mfd,
                in_mfd_i,
                in_mfu,
                in_mfu_i,
                in_supsat,
                in_supsat_i,
                in_tnd_cml_q,
                in_tnd_cml_q_i,
                in_tnd_cml_ql,
                in_tnd_cml_ql_i,
                in_tnd_cml_qi,
                in_tnd_cml_qi_i,
                tmp_trpaus,
                dt,
                covptot=covptot[0, 0, -1],
                covptot_i=covptot_i[0, 0, -1],
                rfl=rfl[0, 0, -1],
                rfl_i=rfl_i[0, 0, -1],
                sfl=sfl[0, 0, -1],
                sfl_i=sfl_i[0, 0, -1],
                t=t,
                t_i=t_i,
            )
            out_fplsl = rfl[0, 0, -1]
            out_fplsl_i = rfl_i[0, 0, -1]
            out_fplsn = sfl[0, 0, -1]
            out_fplsn_i = sfl_i[0, 0, -1]
            out_fhpsl = -out_fplsl * ext.RLVTT
            out_fhpsl_i = -out_fplsl_i * ext.RLVTT
            out_fhpsn = -out_fplsn * ext.RLSTT
            out_fhpsn_i = -out_fplsn_i * ext.RLSTT
        with interval(-1, None):
            out_fplsl = rfl[0, 0, -1]
            out_fplsl_i = rfl_i[0, 0, -1]
            out_fplsn = sfl[0, 0, -1]
            out_fplsn_i = sfl_i[0, 0, -1]
            out_fhpsl = -out_fplsl * ext.RLVTT
            out_fhpsl_i = -out_fplsl_i * ext.RLVTT
            out_fhpsn = -out_fplsn * ext.RLSTT
            out_fhpsn_i = -out_fplsn_i * ext.RLSTT
