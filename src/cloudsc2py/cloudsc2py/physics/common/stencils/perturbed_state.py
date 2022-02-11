# -*- coding: utf-8 -*-
from gt4py import gtscript

from cloudsc2py.framework.stencil import stencil_collection


@stencil_collection("perturbed_state")
def pertubed_state_def(
    in_aph: gtscript.Field["ftype"],
    in_aph_i: gtscript.Field["ftype"],
    in_ap: gtscript.Field["ftype"],
    in_ap_i: gtscript.Field["ftype"],
    in_q: gtscript.Field["ftype"],
    in_q_i: gtscript.Field["ftype"],
    in_qsat: gtscript.Field["ftype"],
    in_qsat_i: gtscript.Field["ftype"],
    in_t: gtscript.Field["ftype"],
    in_t_i: gtscript.Field["ftype"],
    in_ql: gtscript.Field["ftype"],
    in_ql_i: gtscript.Field["ftype"],
    in_qi: gtscript.Field["ftype"],
    in_qi_i: gtscript.Field["ftype"],
    in_lude: gtscript.Field["ftype"],
    in_lude_i: gtscript.Field["ftype"],
    in_lu: gtscript.Field["ftype"],
    in_lu_i: gtscript.Field["ftype"],
    in_mfu: gtscript.Field["ftype"],
    in_mfu_i: gtscript.Field["ftype"],
    in_mfd: gtscript.Field["ftype"],
    in_mfd_i: gtscript.Field["ftype"],
    in_tnd_cml_t: gtscript.Field["ftype"],
    in_tnd_cml_t_i: gtscript.Field["ftype"],
    in_tnd_cml_q: gtscript.Field["ftype"],
    in_tnd_cml_q_i: gtscript.Field["ftype"],
    in_tnd_cml_ql: gtscript.Field["ftype"],
    in_tnd_cml_ql_i: gtscript.Field["ftype"],
    in_tnd_cml_qi: gtscript.Field["ftype"],
    in_tnd_cml_qi_i: gtscript.Field["ftype"],
    in_supsat: gtscript.Field["ftype"],
    in_supsat_i: gtscript.Field["ftype"],
    out_aph: gtscript.Field["ftype"],
    out_ap: gtscript.Field["ftype"],
    out_q: gtscript.Field["ftype"],
    out_qsat: gtscript.Field["ftype"],
    out_t: gtscript.Field["ftype"],
    out_ql: gtscript.Field["ftype"],
    out_qi: gtscript.Field["ftype"],
    out_lude: gtscript.Field["ftype"],
    out_lu: gtscript.Field["ftype"],
    out_mfu: gtscript.Field["ftype"],
    out_mfd: gtscript.Field["ftype"],
    out_tnd_cml_t: gtscript.Field["ftype"],
    out_tnd_cml_q: gtscript.Field["ftype"],
    out_tnd_cml_ql: gtscript.Field["ftype"],
    out_tnd_cml_qi: gtscript.Field["ftype"],
    out_supsat: gtscript.Field["ftype"],
    *,
    f: "ftype",
):
    with computation(PARALLEL), interval(...):
        out_aph = in_aph + f * in_aph_i
        out_ap = in_ap + f * in_ap_i
        out_q = in_q + f * in_q_i
        out_qsat = in_qsat + f * in_qsat_i
        out_t = in_t + f * in_t_i
        out_ql = in_ql + f * in_ql_i
        out_qi = in_qi + f * in_qi_i
        out_lude = in_lude + f * in_lude_i
        out_lu = in_lu + f * in_lu_i
        out_mfu = in_mfu + f * in_mfu_i
        out_mfd = in_mfd + f * in_mfd_i
        out_tnd_cml_t = in_tnd_cml_t + f * in_tnd_cml_t_i
        out_tnd_cml_q = in_tnd_cml_q + f * in_tnd_cml_q_i
        out_tnd_cml_ql = in_tnd_cml_ql + f * in_tnd_cml_ql_i
        out_tnd_cml_qi = in_tnd_cml_qi + f * in_tnd_cml_qi_i
        out_supsat = in_supsat + f * in_supsat_i
