# -*- coding: utf-8 -*-
from gt4py.cartesian import gtscript

from cloudsc2py.framework.stencil import stencil_collection


@stencil_collection("perturbed_state")
def pertubed_state_def(
    in_aph: gtscript.Field["float"],
    in_aph_i: gtscript.Field["float"],
    in_ap: gtscript.Field["float"],
    in_ap_i: gtscript.Field["float"],
    in_q: gtscript.Field["float"],
    in_q_i: gtscript.Field["float"],
    in_qsat: gtscript.Field["float"],
    in_qsat_i: gtscript.Field["float"],
    in_t: gtscript.Field["float"],
    in_t_i: gtscript.Field["float"],
    in_ql: gtscript.Field["float"],
    in_ql_i: gtscript.Field["float"],
    in_qi: gtscript.Field["float"],
    in_qi_i: gtscript.Field["float"],
    in_lude: gtscript.Field["float"],
    in_lude_i: gtscript.Field["float"],
    in_lu: gtscript.Field["float"],
    in_lu_i: gtscript.Field["float"],
    in_mfu: gtscript.Field["float"],
    in_mfu_i: gtscript.Field["float"],
    in_mfd: gtscript.Field["float"],
    in_mfd_i: gtscript.Field["float"],
    in_tnd_cml_t: gtscript.Field["float"],
    in_tnd_cml_t_i: gtscript.Field["float"],
    in_tnd_cml_q: gtscript.Field["float"],
    in_tnd_cml_q_i: gtscript.Field["float"],
    in_tnd_cml_ql: gtscript.Field["float"],
    in_tnd_cml_ql_i: gtscript.Field["float"],
    in_tnd_cml_qi: gtscript.Field["float"],
    in_tnd_cml_qi_i: gtscript.Field["float"],
    in_supsat: gtscript.Field["float"],
    in_supsat_i: gtscript.Field["float"],
    out_aph: gtscript.Field["float"],
    out_ap: gtscript.Field["float"],
    out_q: gtscript.Field["float"],
    out_qsat: gtscript.Field["float"],
    out_t: gtscript.Field["float"],
    out_ql: gtscript.Field["float"],
    out_qi: gtscript.Field["float"],
    out_lude: gtscript.Field["float"],
    out_lu: gtscript.Field["float"],
    out_mfu: gtscript.Field["float"],
    out_mfd: gtscript.Field["float"],
    out_tnd_cml_t: gtscript.Field["float"],
    out_tnd_cml_q: gtscript.Field["float"],
    out_tnd_cml_ql: gtscript.Field["float"],
    out_tnd_cml_qi: gtscript.Field["float"],
    out_supsat: gtscript.Field["float"],
    *,
    f: "float",
):
    with computation(PARALLEL), interval(...):
        out_aph[0, 0, 0] = in_aph + f * in_aph_i
        out_ap[0, 0, 0] = in_ap + f * in_ap_i
        out_q[0, 0, 0] = in_q + f * in_q_i
        out_qsat[0, 0, 0] = in_qsat + f * in_qsat_i
        out_t[0, 0, 0] = in_t + f * in_t_i
        out_ql[0, 0, 0] = in_ql + f * in_ql_i
        out_qi[0, 0, 0] = in_qi + f * in_qi_i
        out_lude[0, 0, 0] = in_lude + f * in_lude_i
        out_lu[0, 0, 0] = in_lu + f * in_lu_i
        out_mfu[0, 0, 0] = in_mfu + f * in_mfu_i
        out_mfd[0, 0, 0] = in_mfd + f * in_mfd_i
        out_tnd_cml_t[0, 0, 0] = in_tnd_cml_t + f * in_tnd_cml_t_i
        out_tnd_cml_q[0, 0, 0] = in_tnd_cml_q + f * in_tnd_cml_q_i
        out_tnd_cml_ql[0, 0, 0] = in_tnd_cml_ql + f * in_tnd_cml_ql_i
        out_tnd_cml_qi[0, 0, 0] = in_tnd_cml_qi + f * in_tnd_cml_qi_i
        out_supsat[0, 0, 0] = in_supsat + f * in_supsat_i
