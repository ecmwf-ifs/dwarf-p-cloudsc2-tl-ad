# -*- coding: utf-8 -*-
from gt4py import gtscript

from cloudsc2py.framework.stencil import stencil_collection


@stencil_collection("state_increment")
def state_increment_def(
    in_aph: gtscript.Field["ftype"],
    in_ap: gtscript.Field["ftype"],
    in_q: gtscript.Field["ftype"],
    in_qsat: gtscript.Field["ftype"],
    in_t: gtscript.Field["ftype"],
    in_ql: gtscript.Field["ftype"],
    in_qi: gtscript.Field["ftype"],
    in_lude: gtscript.Field["ftype"],
    in_lu: gtscript.Field["ftype"],
    in_mfu: gtscript.Field["ftype"],
    in_mfd: gtscript.Field["ftype"],
    in_tnd_cml_t: gtscript.Field["ftype"],
    in_tnd_cml_q: gtscript.Field["ftype"],
    in_tnd_cml_ql: gtscript.Field["ftype"],
    in_tnd_cml_qi: gtscript.Field["ftype"],
    in_supsat: gtscript.Field["ftype"],
    out_aph_i: gtscript.Field["ftype"],
    out_ap_i: gtscript.Field["ftype"],
    out_q_i: gtscript.Field["ftype"],
    out_qsat_i: gtscript.Field["ftype"],
    out_t_i: gtscript.Field["ftype"],
    out_ql_i: gtscript.Field["ftype"],
    out_qi_i: gtscript.Field["ftype"],
    out_lude_i: gtscript.Field["ftype"],
    out_lu_i: gtscript.Field["ftype"],
    out_mfu_i: gtscript.Field["ftype"],
    out_mfd_i: gtscript.Field["ftype"],
    out_tnd_cml_t_i: gtscript.Field["ftype"],
    out_tnd_cml_q_i: gtscript.Field["ftype"],
    out_tnd_cml_ql_i: gtscript.Field["ftype"],
    out_tnd_cml_qi_i: gtscript.Field["ftype"],
    out_supsat_i: gtscript.Field["ftype"],
    *,
    f: "ftype",
):
    with computation(PARALLEL), interval(...):
        out_aph_i = f * in_aph
        out_ap_i = f * in_ap
        out_q_i = f * in_q
        out_qsat_i = f * in_qsat
        out_t_i = f * in_t
        out_ql_i = f * in_ql
        out_qi_i = f * in_qi
        out_lude_i = f * in_lude
        out_lu_i = f * in_lu
        out_mfu_i = f * in_mfu
        out_mfd_i = f * in_mfd
        out_tnd_cml_t_i = f * in_tnd_cml_t
        out_tnd_cml_q_i = f * in_tnd_cml_q
        out_tnd_cml_ql_i = f * in_tnd_cml_ql
        out_tnd_cml_qi_i = f * in_tnd_cml_qi
        out_supsat_i = f * in_supsat
