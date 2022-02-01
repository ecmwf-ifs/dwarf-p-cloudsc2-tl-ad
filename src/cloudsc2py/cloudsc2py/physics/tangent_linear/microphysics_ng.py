# -*- coding: utf-8 -*-
from functools import partial
from typing import Dict, Optional, TYPE_CHECKING

from gt4py import gtscript

from cloudsc2py.framework.components import ImplicitTendencyComponent
from cloudsc2py.framework.stencil import (
    function_collection,
    stencil_collection,
)
from cloudsc2py.utils.f2py import ported_method
from cloudsc2py.utils.storage import get_array

if TYPE_CHECKING:
    from datetime import timedelta

    from sympl._core.typingx import PropertyDict

    from cloudsc2py.framework.grid import Grid
    from cloudsc2py.framework.options import BackendOptions, StorageOptions
    from cloudsc2py.utils.typingx import ArrayDict, ParameterDict


class Cloudsc2TL(ImplicitTendencyComponent):
    def __init__(
        self,
        grid: "Grid",
        ldphylin: bool,
        ldrain1d: bool,
        yoethf_parameters: Optional["ParameterDict"] = None,
        yomcst_parameters: Optional["ParameterDict"] = None,
        yrecld_parameters: Optional["ParameterDict"] = None,
        yrecldp_parameters: Optional["ParameterDict"] = None,
        yrephli_parameters: Optional["ParameterDict"] = None,
        yrncl_parameters: Optional["ParameterDict"] = None,
        yrphnc_parameters: Optional["ParameterDict"] = None,
        *,
        enable_checks: bool = True,
        backend: str = "numpy",
        backend_options: Optional["BackendOptions"] = None,
        storage_options: Optional["StorageOptions"] = None,
    ) -> None:
        super().__init__(
            grid,
            enable_checks=enable_checks,
            backend=backend,
            backend_options=backend_options,
            storage_options=storage_options,
        )

        externals = {
            "ICALL": 0,
            "LDPHYLIN": ldphylin,
            "LDRAIN1D": ldrain1d,
            "ZEPS1": 1e-12,
            "ZEPS2": 1e-10,
            "ZQMAX": 0.5,
            "ZSCAL": 0.9,
        }
        externals.update(yoethf_parameters or {})
        externals.update(yomcst_parameters or {})
        externals.update(yrecld_parameters or {})
        externals.update(yrecldp_parameters or {})
        externals.update(yrephli_parameters or {})
        externals.update(yrncl_parameters or {})
        externals.update(yrphnc_parameters or {})
        self.bo.externals.update(externals)

        self.cloudsc = self.compile_stencil("cloudsc2_tl")

        # allocate temporary 2d arrays
        allocate_f = partial(
            get_array,
            self.grid,
            (self.grid.dims_x, self.grid.dims_y),
            backend=backend,
            dtype=self.so.dtypes.float,
            storage_options=self.so,
        )
        self.temporary_fields = {"tmp_trpaus": allocate_f()}

    @property
    @ported_method(
        from_file="cloudsc2_tl/cloudsc2tl.F90", from_line=53, to_line=78
    )
    def input_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        dims_zh = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_zh)
        return {
            "f_eta": {"dims": dims[2:], "units": ""},
            "f_aph": {"dims": dims_zh, "units": "Pa"},
            "f_aph_i": {"dims": dims_zh, "units": "Pa"},
            "f_ap": {"dims": dims, "units": "Pa"},
            "f_ap_i": {"dims": dims, "units": "Pa"},
            "f_q": {"dims": dims, "units": "g g^-1"},
            "f_q_i": {"dims": dims, "units": "g g^-1"},
            "f_qsat": {"dims": dims, "units": "g g^-1"},
            "f_qsat_i": {"dims": dims, "units": "g g^-1"},
            "f_t": {"dims": dims, "units": "K"},
            "f_t_i": {"dims": dims, "units": "K"},
            "f_ql": {"dims": dims, "units": "g g^-1"},
            "f_ql_i": {"dims": dims, "units": "g g^-1"},
            "f_qi": {"dims": dims, "units": "g g^-1"},
            "f_qi_i": {"dims": dims, "units": "g g^-1"},
            "f_lude": {"dims": dims, "units": "kg m^-3 s^-1"},
            "f_lude_i": {"dims": dims, "units": "kg m^-3 s^-1"},
            "f_lu": {"dims": dims, "units": "g g^-1"},
            "f_lu_i": {"dims": dims, "units": "g g^-1"},
            "f_mfu": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_mfu_i": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_mfd": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_mfd_i": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_cml_tnd_t": {"dims": dims, "units": "K s^-1"},
            "f_cml_tnd_t_i": {"dims": dims, "units": "K s^-1"},
            "f_cml_tnd_q": {"dims": dims, "units": "K s^-1"},
            "f_cml_tnd_q_i": {"dims": dims, "units": "K s^-1"},
            "f_cml_tnd_ql": {"dims": dims, "units": "K s^-1"},
            "f_cml_tnd_ql_i": {"dims": dims, "units": "K s^-1"},
            "f_cml_tnd_qi": {"dims": dims, "units": "K s^-1"},
            "f_cml_tnd_qi_i": {"dims": dims, "units": "K s^-1"},
            "f_supsat": {"dims": dims, "units": "g g^-1"},
            "f_supsat_i": {"dims": dims, "units": "g g^-1"},
        }

    @property
    @ported_method(
        from_file="cloudsc2_tl/cloudsc2tl.F90", from_line=80, to_line=90
    )
    def tendency_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        return {
            "f_t": {"dims": dims, "units": "K s^-1"},
            "f_t_i": {"dims": dims, "units": "K s^-1"},
            "f_q": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_q_i": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_ql": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_ql_i": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_qi": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_qi_i": {"dims": dims, "units": "g g^-1 s^-1"},
        }

    @property
    @ported_method(
        from_file="cloudsc2_tl/cloudsc2tl.F90", from_line=92, to_line=106
    )
    def diagnostic_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        dims_zh = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_zh)
        return {
            "f_clc": {"dims": dims, "units": ""},
            "f_clc_i": {"dims": dims, "units": ""},
            "f_fhpsl": {"dims": dims_zh, "units": "J m^-2 s^-1"},
            "f_fhpsl_i": {"dims": dims_zh, "units": "J m^-2 s^-1"},
            "f_fhpsn": {"dims": dims_zh, "units": "J m^-2 s^-1"},
            "f_fhpsn_i": {"dims": dims_zh, "units": "J m^-2 s^-1"},
            "f_fplsl": {"dims": dims_zh, "units": "Kg m^-2 s^-1"},
            "f_fplsl_i": {"dims": dims_zh, "units": "Kg m^-2 s^-1"},
            "f_fplsn": {"dims": dims_zh, "units": "Kg m^-2 s^-1"},
            "f_fplsn_i": {"dims": dims_zh, "units": "Kg m^-2 s^-1"},
            "f_covptot": {"dims": dims, "units": ""},
            "f_covptot_i": {"dims": dims, "units": ""},
        }

    def array_call(
        self,
        state: "ArrayDict",
        timestep: "timedelta",
        out_tendencies: "ArrayDict",
        out_diagnostics: "ArrayDict",
        overwrite_tendencies: Dict[str, bool],
    ) -> None:
        self.cloudsc(
            in_eta=state["f_eta"],
            in_ap=state["f_ap"],
            in_ap_i=state["f_ap_i"],
            in_aph=state["f_aph"],
            in_aph_i=state["f_aph_i"],
            in_t=state["f_t"],
            in_t_i=state["f_t_i"],
            in_q=state["f_q"],
            in_q_i=state["f_q_i"],
            in_qsat=state["f_qsat"],
            in_qsat_i=state["f_qsat_i"],
            in_ql=state["f_ql"],
            in_ql_i=state["f_ql_i"],
            in_qi=state["f_qi"],
            in_qi_i=state["f_qi_i"],
            in_lu=state["f_lu"],
            in_lu_i=state["f_lu_i"],
            in_lude=state["f_lude"],
            in_lude_i=state["f_lude_i"],
            in_mfd=state["f_mfd"],
            in_mfd_i=state["f_mfd_i"],
            in_mfu=state["f_mfu"],
            in_mfu_i=state["f_mfu_i"],
            in_supsat=state["f_supsat"],
            in_supsat_i=state["f_supsat_i"],
            in_tnd_t=state["f_cml_tnd_t"],
            in_tnd_t_i=state["f_cml_tnd_t_i"],
            in_tnd_q=state["f_cml_tnd_q"],
            in_tnd_q_i=state["f_cml_tnd_q_i"],
            in_tnd_ql=state["f_cml_tnd_ql"],
            in_tnd_ql_i=state["f_cml_tnd_ql_i"],
            in_tnd_qi=state["f_cml_tnd_qi"],
            in_tnd_qi_i=state["f_cml_tnd_qi_i"],
            **self.temporary_fields,
            out_tnd_t=out_tendencies["f_t"],
            out_tnd_t_i=out_tendencies["f_t_i"],
            out_tnd_q=out_tendencies["f_q"],
            out_tnd_q_i=out_tendencies["f_q_i"],
            out_tnd_ql=out_tendencies["f_ql"],
            out_tnd_ql_i=out_tendencies["f_ql_i"],
            out_tnd_qi=out_tendencies["f_qi"],
            out_tnd_qi_i=out_tendencies["f_qi_i"],
            out_clc=out_diagnostics["f_clc"],
            out_clc_i=out_diagnostics["f_clc_i"],
            out_fhpsl=out_diagnostics["f_fhpsl"],
            out_fhpsl_i=out_diagnostics["f_fhpsl_i"],
            out_fhpsn=out_diagnostics["f_fhpsn"],
            out_fhpsn_i=out_diagnostics["f_fhpsn_i"],
            out_fplsl=out_diagnostics["f_fplsl"],
            out_fplsl_i=out_diagnostics["f_fplsl_i"],
            out_fplsn=out_diagnostics["f_fplsn"],
            out_fplsn_i=out_diagnostics["f_fplsn_i"],
            out_covptot=out_diagnostics["f_covptot"],
            out_covptot_i=out_diagnostics["f_covptot_i"],
            dt=timestep.total_seconds(),
            origin=(0, 0, 0),
            domain=(self.grid.nx, self.grid.ny, self.grid.nz + 1),
            validate_args=self.bo.validate_args,
            exec_info=self.bo.exec_info,
        )

    @staticmethod
    @ported_method(
        from_file="cloudsc2_tl/cloudsc2tl.F90", from_line=321, to_line=1113
    )
    @stencil_collection(
        "cloudsc2_tl", external_names=["RLSTT", "RLVTT", "cloudsc2_tl_0"]
    )
    def cloudsc2_tl_def(
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
        in_tnd_t: gtscript.Field["ftype"],
        in_tnd_t_i: gtscript.Field["ftype"],
        in_tnd_q: gtscript.Field["ftype"],
        in_tnd_q_i: gtscript.Field["ftype"],
        in_tnd_ql: gtscript.Field["ftype"],
        in_tnd_ql_i: gtscript.Field["ftype"],
        in_tnd_qi: gtscript.Field["ftype"],
        in_tnd_qi_i: gtscript.Field["ftype"],
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
            t = in_t + dt * in_tnd_t
            t_i = in_t_i + dt * in_tnd_t_i

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
                    in_tnd_q,
                    in_tnd_q_i,
                    in_tnd_ql,
                    in_tnd_ql_i,
                    in_tnd_qi,
                    in_tnd_qi_i,
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
                    in_tnd_q,
                    in_tnd_q_i,
                    in_tnd_ql,
                    in_tnd_ql_i,
                    in_tnd_qi,
                    in_tnd_qi_i,
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

    @staticmethod
    @ported_method(
        from_file="cloudsc2_tl/cloudsc2tl.F90", from_line=321, to_line=1113
    )
    @function_collection(
        "cloudsc2_tl_0",
        external_names=[
            "LDPHYLIN",
            "LDRAIN1D",
            "LEVAPLS2",
            "LREGCL",
            "R2ES",
            "R3IES",
            "R3LES",
            "R4IES",
            "R4LES",
            "R5IES",
            "R5LES",
            "RCLCRIT",
            "RCPD",
            "RD",
            "RETV",
            "RG",
            "RKCONV",
            "RLMIN",
            "RLMLT",
            "RLPTRC",
            "RLSTT",
            "RLVTT",
            "RPECONS",
            "RTICE",
            "RTT",
            "RVTMP2",
            "ZEPS1",
            "ZEPS2",
            "ZQMAX",
            "ZSCAL",
            "cuadjtqs_tl",
        ],
    )
    @gtscript.function
    def cloudsc2_tl_0(
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
        in_tnd_q,
        in_tnd_q_i,
        in_tnd_ql,
        in_tnd_ql_i,
        in_tnd_qi,
        in_tnd_qi_i,
        tmp_trpaus,
        dt,
        covptot,
        covptot_i,
        rfl,
        rfl_i,
        sfl,
        sfl_i,
        t,
        t_i,
    ):
        from __externals__ import ext

        # first guess values for q, ql and qi
        q = in_q + dt * in_tnd_q + in_supsat
        q_i = in_q_i + dt * in_tnd_q_i + in_supsat_i
        ql = in_ql + dt * in_tnd_ql
        ql_i = in_ql_i + dt * in_tnd_ql_i
        qi = in_qi + dt * in_tnd_qi
        qi_i = in_qi_i + dt * in_tnd_qi_i

        # set up constants required
        ckcodtl = 2 * ext.RKCONV * dt
        ckcodti = 5 * ext.RKCONV * dt
        ckcodtla = ckcodtl / 100
        ckcodtia = ckcodti / 100
        cons2 = 1 / (ext.RG * dt)
        cons3 = ext.RLVTT / ext.RCPD
        meltp2 = ext.RTT + 2

        # parameter for cloud formation
        scalm = ext.ZSCAL * max(in_eta - 0.2, ext.ZEPS1) ** 0.2

        # thermodynamic constants
        dp = in_aph[0, 0, 1] - in_aph[0, 0, 0]
        dp_i = in_aph_i[0, 0, 1] - in_aph_i[0, 0, 0]
        zz = ext.RCPD + ext.RCPD * ext.RVTMP2 * q
        zz_i = (
            -ext.RCPD
            * ext.RVTMP2
            * q_i
            / (ext.RCPD + ext.RCPD * ext.RVTMP2 * q) ** 2
        )
        lfdcp = ext.RLMLT / zz
        lfdcp_i = ext.RLMLT * zz_i
        lsdcp = ext.RLSTT / zz
        lsdcp_i = ext.RLSTT * zz_i
        lvdcp = ext.RLVTT / zz
        lvdcp_i = ext.RLVTT * zz_i

        # calculate dqs/dT correction factor
        if t < ext.RTT:
            fwat = 0.545 * (tanh(0.17 * (t - ext.RLPTRC)) + 1)
            fwat_i = 0.545 * 0.17 * t_i / (cosh(0.17 * (t - ext.RLPTRC)) ** 2)
            z3es = ext.R3IES
            z4es = ext.R4IES
        else:
            fwat = 1.0
            fwat_i = 0.0
            z3es = ext.R3LES
            z4es = ext.R4LES
        foeew = ext.R2ES * exp(z3es * (t - ext.RTT) / (t - z4es))
        foeew_i = z3es * (ext.RTT - z4es) * t_i * foeew / (t - z4es) ** 2
        esdp = foeew / in_ap
        esdp_i = foeew_i / in_ap - foeew * in_ap_i / (in_ap ** 2)
        if esdp > ext.ZQMAX:
            esdp = ext.ZQMAX
            esdp_i = 0.0

        facw = ext.R5LES / (t - ext.R4LES) ** 2
        facw_i = -2 * ext.R5LES * t_i / (t - ext.R4LES) ** 3
        faci = ext.R5IES / (t - ext.R4IES) ** 2
        faci_i = -2 * ext.R5IES * t_i / (t - ext.R4IES) ** 3
        fac = fwat * facw + (1 - fwat) * faci
        fac_i = fwat_i * (facw - faci) + fwat * facw_i + (1 - fwat) * faci_i
        cor = 1 / (1 - ext.RETV * esdp)
        cor_i = ext.RETV * esdp_i / (1 - ext.RETV * esdp) ** 2
        dqsdtemp = fac * cor * in_qsat
        dqsdtemp_i = (
            fac_i * cor * in_qsat
            + fac * cor_i * in_qsat
            + fac * cor * in_qsat_i
        )
        corqs = 1 + cons3 * dqsdtemp
        corqs_i = cons3 * dqsdtemp_i

        # use clipped state
        if q > in_qsat:
            qlim = in_qsat
            qlim_i = in_qsat_i
        else:
            qlim = q
            qlim_i = q_i

        # set up critical value of humidity
        rh1 = 1.0
        rh2 = (
            0.35
            + 0.14 * ((tmp_trpaus - 0.25) / 0.15) ** 2
            + 0.04 * min(tmp_trpaus - 0.25, 0.0) / 0.15
        )
        rh3 = 1.0
        if in_eta < tmp_trpaus:
            crh2 = rh3
        else:
            deta2 = 0.3
            bound1 = tmp_trpaus + deta2
            if in_eta < bound1:
                crh2 = rh3 + (rh2 - rh3) * (in_eta - tmp_trpaus) / deta2
            else:
                deta1 = 0.09 + 0.16 * (0.4 - tmp_trpaus) / 0.3
                bound2 = 1 - deta1
                if in_eta < bound2:
                    crh2 = rh2
                else:
                    crh2 = rh1 + (rh2 - rh1) * ((1 - in_eta) / deta1) ** 0.5

        # allow ice supersaturation at cold temperatures
        if t < ext.RTICE:
            qsat = in_qsat * (1.8 - 0.003 * t)
            qsat_i = in_qsat_i * (1.8 - 0.003 * t) - 0.003 * in_qsat * t_i
        else:
            qsat = in_qsat
            qsat_i = in_qsat_i
        qcrit = crh2 * qsat
        qcrit_i = crh2 * qsat_i

        # simple uniform distribution of total water from Leutreut & Li (1990)
        qt = q + ql + qi
        qt_i = q_i + ql_i + qi_i
        if qt < qcrit:
            out_clc = 0.0
            out_clc_i = 0.0
            qc = 0.0
            qc_i = 0.0
        elif qt >= qsat:
            out_clc = 1.0
            out_clc_i = 0.0
            qc = (1 - scalm) * (qsat - qcrit)
            qc_i = (1 - scalm) * (qsat_i - qcrit_i)
        else:
            qpd = qsat - qt
            qpd_i = qsat_i - qt_i
            qcd = qsat - qcrit
            qcd_i = qsat_i - qcrit_i
            tmp1 = sqrt(qpd / (qcd - scalm * (qt - qcrit)))
            out_clc = 1 - tmp1
            out_clc_i = (
                -(0.5 / tmp1)
                * (
                    qpd_i * (qcd - scalm * (qt - qcrit))
                    - qpd * (qcd_i - scalm * (qt_i - qcrit_i))
                )
                / (qcd - scalm * (qt - qcrit)) ** 2
            )

            # regularization of cloud fraction perturbation
            if __INLINED(ext.LREGCL):
                rat = qpd / qcd
                yyy = min(
                    0.3,
                    3.5
                    * sqrt(rat * (1 - scalm * (1 - rat)) ** 3)
                    / (1 - scalm),
                )
                out_clc_i *= yyy

            qc = (scalm * qpd + (1 - scalm) * qcd) * out_clc ** 2
            qc_i = (scalm * qpd_i + (1 - scalm) * qcd_i) * out_clc ** 2 + 2 * (
                scalm * qpd + (1 - scalm) * qcd
            ) * out_clc * out_clc_i

        # add convective component
        gdp = ext.RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
        gdp_i = (
            -ext.RG
            * (in_aph_i[0, 0, 1] - in_aph_i[0, 0, 0])
            / (in_aph[0, 0, 1] - in_aph[0, 0, 0]) ** 2
        )
        lude = dt * in_lude * gdp
        lude_i = dt * (in_lude_i * gdp + in_lude * gdp_i)
        lo1 = lude[0, 0, 0] >= ext.RLMIN and in_lu[0, 0, 1] >= ext.ZEPS2
        if lo1:
            tmp2 = exp(-lude[0, 0, 0] / in_lu[0, 0, 1])
            out_clc_i += -out_clc_i[0, 0, 0] * (1 - tmp2) + (
                1 - out_clc[0, 0, 0]
            ) * tmp2 * (
                lude_i[0, 0, 0] / in_lu[0, 0, 1]
                - lude[0, 0, 0] * in_lu_i[0, 0, 1] / in_lu[0, 0, 1] ** 2
            )
            out_clc += (1 - out_clc[0, 0, 0]) * (1 - tmp2)
            qc += lude
            qc_i += lude_i

        # add compensating subsidence component
        rho = in_ap / (ext.RD * t)
        rho_i = (in_ap_i - in_ap * t_i / t) / (ext.RD * t)

        fac2 = in_ap - ext.RETV * foeew
        rodqsdp = -rho * in_qsat / fac2
        rodqsdp_i = (
            -rho_i * in_qsat
            - rho * in_qsat_i
            + rho * in_qsat * (in_ap_i - ext.RETV * foeew_i) / fac2
        ) / fac2

        ldcp = fwat * lvdcp + (1 - fwat) * lsdcp
        ldcp_i = (
            fwat_i * (lvdcp - lsdcp) + fwat * lvdcp_i + (1 - fwat) * lsdcp_i
        )

        fac3 = 1 + ldcp * dqsdtemp
        dtdzmo = ext.RG * (1 / ext.RCPD - ldcp * rodqsdp) / fac3
        dtdzmo_i = (
            -(
                ext.RG * (ldcp_i * rodqsdp + ldcp * rodqsdp_i)
                + dtdzmo * (ldcp_i * dqsdtemp + ldcp * dqsdtemp_i)
            )
            / fac3
        )

        dqsdz = dqsdtemp * dtdzmo - ext.RG * rodqsdp
        dqsdz_i = (
            dqsdtemp_i * dtdzmo + dqsdtemp * dtdzmo_i - ext.RG * rodqsdp_i
        )

        tmp3 = dt * dqsdz * (in_mfu + in_mfd) / rho
        if tmp3 < qc:
            dqc = tmp3
            dqc_i = (
                dt
                * (dqsdz_i * (in_mfu + in_mfd) + dqsdz * (in_mfu_i + in_mfd_i))
                - dqc * rho_i
            ) / rho
            if __INLINED(ext.LREGCL):
                dqc_i *= 0.1
        else:
            dqc = qc
            dqc_i = qc_i
        qc -= dqc
        qc_i -= dqc_i

        # new cloud liquid/ice contents and condensation rates (liquid/ice)
        qlwc = qc * fwat
        qlwc_i = qc_i * fwat + qc * fwat_i

        qiwc = qc * (1 - fwat)
        qiwc_i = qc_i * (1 - fwat) - qc * fwat_i

        condl = (qlwc - ql) / dt
        condl_i = (qlwc_i - ql_i) / dt

        condi = (qiwc - qi) / dt
        condi_i = (qiwc_i - qi_i) / dt

        # calculate precipitation overlap
        # simple form based on Maximum Overlap
        if out_clc > covptot:
            covptotn = out_clc
            covptotn_i = out_clc_i
        else:
            covptotn = covptot
            covptotn_i = covptot_i
        if covptotn > out_clc:
            covpclr = covptotn - out_clc
            covpclr_i = covptotn_i - out_clc_i
        else:
            covpclr = 0.0
            covpclr_i = 0.0

        # melting of incoming snow
        if sfl != 0:
            cons = cons2 * dp / lfdcp
            cons_i = cons2 * (dp_i * lfdcp - dp * lfdcp_i) / lfdcp ** 2
            if t > meltp2:
                z2s = cons * (t - meltp2)
                z2s_i = cons_i * (t - meltp2) + cons * t_i
            else:
                z2s = 0.0
                z2s_i = 0.0

            if sfl < z2s:
                snmlt = sfl
                snmlt_i = sfl_i
            else:
                snmlt = z2s
                snmlt_i = z2s_i

            rfln = rfl + snmlt
            rfln_i = rfl_i + snmlt_i
            sfln = sfl - snmlt
            sfln_i = sfl_i - snmlt_i
            t -= snmlt / cons
            t_i -= (snmlt_i * cons - snmlt * cons_i) / cons ** 2
        else:
            rfln = rfl
            rfln_i = rfl_i
            sfln = sfl
            sfln_i = sfl_i

        if out_clc > ext.ZEPS2:
            # diagnostic calculation of rain production from cloud liquid water
            if __INLINED(ext.LEVAPLS2 or ext.LDRAIN1D):
                lcrit = 1.9 * ext.RCLCRIT
            else:
                lcrit = 2.0 * ext.RCLCRIT

            # in-cloud liquid
            cldl = qlwc / out_clc
            cldl_i = qlwc_i / out_clc - qlwc * out_clc_i / out_clc ** 2

            ltmp4 = exp(-((cldl / lcrit) ** 2))
            dl = ckcodtl * (1 - ltmp4)
            if __INLINED(ext.LREGCL):
                dl_i = (2 * ckcodtla / lcrit ** 2) * ltmp4 * cldl * cldl_i
            else:
                dl_i = (2 * ckcodtl / lcrit ** 2) * ltmp4 * cldl * cldl_i

            ltmp5 = exp(-dl)
            prr = qlwc - out_clc * cldl * ltmp5
            prr_i = qlwc_i - ltmp5 * (
                out_clc_i * cldl + out_clc * cldl_i - out_clc * cldl * dl_i
            )
            qlwc -= prr
            qlwc_i -= prr_i

            # diagnostic calculation of snow production from cloud ice
            if __INLINED(ext.LEVAPLS2 or ext.LDRAIN1D):
                icrit = 0.0001
            else:
                icrit = 2.0 * ext.RCLCRIT

            cldi = qiwc / out_clc
            cldi_i = qiwc_i / out_clc - qiwc * out_clc_i / out_clc ** 2

            itmp41 = exp(-((cldi / icrit) ** 2))
            itmp42 = exp(0.025 * (t - ext.RTT))
            di = ckcodti * itmp42 * (1 - itmp41)
            if __INLINED(ext.LREGCL):
                di_i = (
                    ckcodtia
                    * itmp42
                    * (
                        itmp41 * (2 * cldi * cldi_i / icrit ** 2 - 0.025 * t_i)
                        + 0.025 * t_i
                    )
                )
            else:
                di_i = (
                    ckcodti
                    * itmp42
                    * (
                        0.025 * t_i * (1 - itmp41)
                        + itmp41 * (2 * cldi * cldi_i / icrit ** 2)
                    )
                )

            itmp5 = exp(-di)
            prs = qiwc - out_clc * cldi * itmp5
            prs_i = qiwc_i - itmp5 * (
                out_clc_i * cldi + out_clc * cldi_i - out_clc * cldi * di_i
            )
            qiwc -= prs
            qiwc_i -= prs_i
        else:
            prr = 0.0
            prr_i = 0.0
            prs = 0.0
            prs_i = 0.0

        # new precipitation
        dr = cons2 * dp * (prr + prs)
        dr_i = cons2 * (dp_i * (prr + prs) + dp * (prr_i + prs_i))

        # rain fraction (different from cloud liquid water fraction!)
        if t < ext.RTT:
            rfreeze = cons2 * dp * prr
            rfreeze_i = cons2 * (dp_i * prr + dp * prr_i)
            fwatr = 0.0
            fwatr_i = 0.0
        else:
            rfreeze = 0.0
            rfreeze_i = 0.0
            fwatr = 1.0
            fwatr_i = 0.0
        rfln += fwatr * dr
        rfln_i += fwatr_i * dr + fwatr * dr_i
        sfln += (1 - fwatr) * dr
        sfln_i += -fwatr_i * dr + (1 - fwatr) * dr_i

        # precipitation evaporation
        prtot = rfln + sfln
        prtot_i = rfln_i + sfln_i
        if (
            prtot > ext.ZEPS2
            and covpclr > ext.ZEPS2
            and (ext.LEVAPLS2 or ext.LDRAIN1D)
        ):
            # note: the code never enters this branch when input data
            # are retrieved from input.h5
            preclr = prtot * covpclr / covptotn
            preclr_i = (
                prtot_i * covpclr + prtot * covpclr_i
            ) / covptotn - prtot * covpclr * covptotn_i / covptotn ** 2

            # this is the humidity in the moisest zcovpclr region
            qe = in_qsat - (in_qsat - qlim) * covpclr / (1 - out_clc) ** 2
            qe_i = (
                in_qsat_i
                - (
                    in_qsat_i * covpclr
                    - qlim_i * covpclr
                    + (in_qsat - qlim) * covpclr_i
                )
                / (1 - out_clc) ** 2
                - 2
                * (in_qsat - qlim)
                * covpclr
                * out_clc_i
                / (1 - out_clc) ** 3
            )

            tmp6 = sqrt(in_ap[0, 0, 0] / in_aph[0, 0, 1])
            beta = (
                ext.RG
                * ext.RPECONS
                * (tmp6 * preclr / (0.00509 * covpclr)) ** 0.5777
            )
            beta_i = (
                0.5777
                * ext.RG
                * ext.RPECONS
                / (0.00509 ** 0.5777)
                * (covpclr / (tmp6 * preclr)) ** 0.4223
                * (
                    (
                        0.5
                        / tmp6
                        * (
                            in_ap_i[0, 0, 0] * in_aph[0, 0, 1]
                            - in_ap[0, 0, 0] * in_aph_i[0, 0, 1]
                        )
                        / in_aph_i[0, 0, 1] ** 2
                        * preclr
                        + tmp6 * preclr_i
                    )
                    / covpclr
                    - tmp6 * preclr * covpclr_i / covpclr ** 2
                )
            )

            # implicit solution
            b = dt * beta * (in_qsat - qe) / (1 + dt * beta * corqs)
            b_i = (
                dt * (beta_i * (in_qsat - qe) + beta * (in_qsat_i - qe_i))
                - dt * b * (beta_i * corqs + beta * corqs_i)
            ) / (1 + dt * beta * corqs)

            dtgdp = dt * ext.RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
            dtgdp_i = (
                dt
                * ext.RG
                * (in_aph_i[0, 0, 0] - in_aph_i[0, 0, 1])
                / (in_aph[0, 0, 1] - in_aph[0, 0, 0]) ** 2
            )
            dpr = covpclr * b / dtgdp
            dpr_i = ((covpclr_i * b + covpclr * b_i) - dpr * dtgdp_i) / dtgdp
            if dpr > preclr:
                dpr = preclr
                dpr_i = preclr_i
            preclr -= dpr
            preclr_i -= dpr_i
            if preclr <= 0:
                covptotn = out_clc
                covptotn_i = out_clc_i
            out_covptot = covptotn
            out_covptot_i = covptotn_i

            # warm proportion
            evapr = dpr * rfln / prtot
            evapr_i = (dpr_i * rfln + dpr * rfln_i - evapr * prtot_i) / prtot
            rfln -= evapr
            rfln_i -= evapr_i

            # ice proportion
            evaps = dpr * sfln / prtot
            evaps_i = (dpr_i * sfln + dpr * sfln_i - evaps * prtot_i) / prtot
            sfln -= evaps
            sfln_i -= evaps_i
        else:
            evapr = 0.0
            evapr_i = 0.0
            evaps = 0.0
            evaps_i = 0.0
            out_covptot = 0.0
            out_covptot_i = 0.0

        # incrementation of T and Q fluxes' swap
        dqdt = -(condl + condi) + (in_lude + evapr + evaps) * gdp
        dqdt_i = (
            -(condl_i + condi_i)
            + (in_lude_i + evapr_i + evaps_i) * gdp
            + (in_lude + evapr + evaps) * gdp_i
        )

        tmp7 = (
            lvdcp * evapr
            + lsdcp * evaps
            + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
            - (lsdcp - lvdcp) * rfreeze
        )
        dtdt = lvdcp * condl + lsdcp * condi - tmp7 * gdp
        dtdt_i = (
            lvdcp_i * condl
            + lvdcp * condl_i
            + lsdcp_i * condi
            + lsdcp * condi_i
            - (
                lvdcp_i * evapr
                + lvdcp * evapr_i
                + lsdcp_i * evaps
                + lsdcp * evaps_i
                + in_lude_i * (fwat * lvdcp + (1 - fwat) * lsdcp)
                + in_lude
                * (
                    fwat_i * (lvdcp - lsdcp)
                    + fwat * lvdcp_i
                    + (1 - fwat) * lsdcp_i
                )
                - (lsdcp_i - lvdcp_i) * rfreeze
                - (lsdcp - lvdcp) * rfreeze_i
            )
            * gdp
            - tmp7 * gdp_i
        )

        # first guess T and Q
        t += dt * dtdt
        t_i += dt * dtdt_i
        q += dt * dqdt
        q_i += dt * dqdt_i
        qold = q
        qold_i = q_i

        # clipping of final qv
        t, t_i, q, q_i = ext.cuadjtqs_tl(in_ap, in_ap_i, t, t_i, q, q_i)

        if qold > q:
            dq = qold - q
            dq_i = qold_i - q_i
            if __INLINED(ext.LREGCL):
                dq_i *= 0.7
        else:
            dq = 0.0
            dq_i = 0.0
        dr2 = cons2 * dp * dq
        dr2_i = cons2 * (dp_i * dq + dp * dq_i)

        # update rain fraction and freezing
        # note: impact of new temperature t_i on fwat_i is neglected here
        if t < ext.RTT:
            rfreeze2 = fwat * dr2
            rfreeze2_i = fwat_i * dr2 + fwat * dr2_i
            fwatr = 0.0
            fwatr_i = 0.0
        else:
            rfreeze2 = 0.0
            rfreeze2_i = 0.0
            fwatr = 1.0
            fwatr_i = 0.0

        rn = fwatr * dr2
        rn_i = fwatr_i * dr2 + fwatr * dr2_i
        sn = (1 - fwatr) * dr2
        sn_i = -fwatr_i * dr2 + (1 - fwatr) * dr2_i

        # note: the extra condensation due to the adjustment goes directly
        # to precipitation
        condl += fwatr * dq / dt
        condl_i += (fwatr_i * dq + fwatr * dq_i) / dt
        condi += (1 - fwatr) * dq / dt
        condi_i += (-fwatr_i * dq + (1 - fwatr) * dq_i) / dt
        rfln += rn
        rfln_i += rn_i
        sfln += sn
        sfln_i += sn_i
        rfreeze += rfreeze2
        rfreeze_i += rfreeze2_i

        # calculate output tendencies
        out_tnd_q = -(condl + condi) + (in_lude + evapr + evaps) * gdp
        out_tnd_q_i = (
            -(condl_i + condi_i)
            + (in_lude_i + evapr_i + evaps_i) * gdp
            + (in_lude + evapr + evaps) * gdp_i
        )
        tmp8 = (
            lvdcp * evapr
            + lsdcp * evaps
            + in_lude * (fwat * lvdcp + (1 - fwat) * lsdcp)
            - (lsdcp - lvdcp) * rfreeze
        )
        out_tnd_t = lvdcp * condl + lsdcp * condi - tmp8 * gdp
        out_tnd_t_i = (
            lvdcp_i * condl
            + lvdcp * condl_i
            + lsdcp_i * condi
            + lsdcp * condi_i
            - (
                lvdcp_i * evapr
                + lvdcp * evapr_i
                + lsdcp_i * evaps
                + lsdcp * evaps_i
                + in_lude_i * (fwat * lvdcp + (1 - fwat) * lsdcp)
                + in_lude
                * (
                    fwat_i * (lvdcp - lsdcp)
                    + fwat * lvdcp_i
                    + (1 - fwat) * lsdcp_i
                )
                - (lsdcp_i - lvdcp_i) * rfreeze
                - (lsdcp - lvdcp) * rfreeze_i
            )
            * gdp
            - tmp8 * gdp_i
        )
        out_tnd_ql = (qlwc - ql) / dt
        out_tnd_ql_i = (qlwc_i - ql_i) / dt
        out_tnd_qi = (qiwc - qi) / dt
        out_tnd_qi_i = (qiwc_i - qi_i) / dt

        return (
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
            covptotn,
            covptotn_i,
            rfln,
            rfln_i,
            sfln,
            sfln_i,
        )
