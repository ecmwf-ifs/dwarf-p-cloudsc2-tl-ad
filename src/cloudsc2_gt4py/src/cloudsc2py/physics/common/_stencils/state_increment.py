# -*- coding: utf-8 -*-
from gt4py.cartesian import gtscript

from ifs_physics_common.framework.stencil import stencil_collection


@stencil_collection("state_increment")
def state_increment_def(
    in_aph: gtscript.Field["float"],
    in_ap: gtscript.Field["float"],
    in_q: gtscript.Field["float"],
    in_qsat: gtscript.Field["float"],
    in_t: gtscript.Field["float"],
    in_ql: gtscript.Field["float"],
    in_qi: gtscript.Field["float"],
    in_lude: gtscript.Field["float"],
    in_lu: gtscript.Field["float"],
    in_mfu: gtscript.Field["float"],
    in_mfd: gtscript.Field["float"],
    in_tnd_cml_t: gtscript.Field["float"],
    in_tnd_cml_q: gtscript.Field["float"],
    in_tnd_cml_ql: gtscript.Field["float"],
    in_tnd_cml_qi: gtscript.Field["float"],
    in_supsat: gtscript.Field["float"],
    out_aph_i: gtscript.Field["float"],
    out_ap_i: gtscript.Field["float"],
    out_q_i: gtscript.Field["float"],
    out_qsat_i: gtscript.Field["float"],
    out_t_i: gtscript.Field["float"],
    out_ql_i: gtscript.Field["float"],
    out_qi_i: gtscript.Field["float"],
    out_lude_i: gtscript.Field["float"],
    out_lu_i: gtscript.Field["float"],
    out_mfu_i: gtscript.Field["float"],
    out_mfd_i: gtscript.Field["float"],
    out_tnd_cml_t_i: gtscript.Field["float"],
    out_tnd_cml_q_i: gtscript.Field["float"],
    out_tnd_cml_ql_i: gtscript.Field["float"],
    out_tnd_cml_qi_i: gtscript.Field["float"],
    out_supsat_i: gtscript.Field["float"],
    *,
    f: "float",
):
    from __externals__ import IGNORE_SUPSAT

    with computation(PARALLEL), interval(...):
        out_aph_i[0, 0, 0] = f * in_aph
        out_ap_i[0, 0, 0] = f * in_ap
        out_q_i[0, 0, 0] = f * in_q
        out_qsat_i[0, 0, 0] = f * in_qsat
        out_t_i[0, 0, 0] = f * in_t
        out_ql_i[0, 0, 0] = f * in_ql
        out_qi_i[0, 0, 0] = f * in_qi
        out_lude_i[0, 0, 0] = f * in_lude
        out_lu_i[0, 0, 0] = f * in_lu
        out_mfu_i[0, 0, 0] = f * in_mfu
        out_mfd_i[0, 0, 0] = f * in_mfd
        out_tnd_cml_t_i[0, 0, 0] = f * in_tnd_cml_t
        out_tnd_cml_q_i[0, 0, 0] = f * in_tnd_cml_q
        out_tnd_cml_ql_i[0, 0, 0] = f * in_tnd_cml_ql
        out_tnd_cml_qi_i[0, 0, 0] = f * in_tnd_cml_qi
        if not IGNORE_SUPSAT:
            out_supsat_i[0, 0, 0] = f * in_supsat
        else:
            out_supsat_i[0, 0, 0] = 0.0
