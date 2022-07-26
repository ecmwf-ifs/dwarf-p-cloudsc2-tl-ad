# -*- coding: utf-8 -*-
from gt4py import gtscript

from cloudsc2py.framework.stencil import stencil_collection
from cloudsc2py.utils.f2py import ported_function


@ported_function(from_file="cloudsc2_nl/dwarf_cloudsc.F90", from_line=100, to_line=102)
@stencil_collection(name="diagnose_eta")
def diagnose_eta_def(
    in_ap: gtscript.Field["float"],
    out_eta: gtscript.Field[gtscript.K, "float"],
    *,
    ap_top: "float",
):
    with computation(FORWARD), interval(...):
        out_eta = in_ap / ap_top
