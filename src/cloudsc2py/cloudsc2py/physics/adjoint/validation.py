# -*- coding: utf-8 -*-
import numpy as np
import sys
from typing import Optional, TYPE_CHECKING, Tuple

import gt4py

from cloudsc2py.physics.adjoint.microphysics_ng import Cloudsc2AD
from cloudsc2py.physics.common.increment import StateIncrement
from cloudsc2py.physics.common.saturation import Saturation
from cloudsc2py.physics.tangent_linear.microphysics_ng import Cloudsc2TL
from cloudsc2py.utils.f2py import ported_method

if TYPE_CHECKING:
    from datetime import timedelta

    from sympl._core.typingx import DataArrayDict

    from cloudsc2py.framework.grid import Grid
    from cloudsc2py.framework.options import BackendOptions, StorageOptions
    from cloudsc2py.utils.typingx import ParameterDict


class SymmetryTest:
    def __init__(
        self,
        grid: "Grid",
        factor: float,
        kflag: int,
        lphylin: bool,
        ldrain1d: bool,
        yoethf_parameters: Optional["ParameterDict"] = None,
        yomcst_parameters: Optional["ParameterDict"] = None,
        yrecld_parameters: Optional["ParameterDict"] = None,
        yrecldp_parameters: Optional["ParameterDict"] = None,
        yrephli_parameters: Optional["ParameterDict"] = None,
        yrncl_parameters: Optional["ParameterDict"] = None,
        yrphnc_parameters: Optional["ParameterDict"] = None,
        enable_checks: bool = True,
        backend: str = "numpy",
        backend_options: Optional["BackendOptions"] = None,
        storage_options: Optional["StorageOptions"] = None,
    ) -> None:
        self.f = factor

        # saturation
        self.saturation = Saturation(
            grid,
            kflag,
            lphylin,
            yoethf_parameters,
            yomcst_parameters,
            enable_checks=enable_checks,
            backend=backend,
            backend_options=backend_options,
            storage_options=storage_options,
        )

        # microphysics
        self.cloudsc2_tl = Cloudsc2TL(
            grid,
            lphylin,
            ldrain1d,
            yoethf_parameters,
            yomcst_parameters,
            yrecld_parameters,
            yrecldp_parameters,
            yrephli_parameters,
            yrncl_parameters,
            yrphnc_parameters,
            enable_checks=enable_checks,
            backend=backend,
            backend_options=backend_options,
            storage_options=storage_options,
        )
        self.cloudsc2_ad = Cloudsc2AD(
            grid,
            lphylin,
            ldrain1d,
            yoethf_parameters,
            yomcst_parameters,
            yrecld_parameters,
            yrecldp_parameters,
            yrephli_parameters,
            yrncl_parameters,
            yrphnc_parameters,
            enable_checks=enable_checks,
            backend=backend,
            backend_options=backend_options,
            storage_options=storage_options,
        )

        # perturbation
        self.state_increment = StateIncrement(
            grid,
            factor,
            enable_checks=enable_checks,
            backend=backend,
            backend_options=backend_options,
            storage_options=storage_options,
        )

    def __call__(self, state: "DataArrayDict", timestep: "timedelta") -> None:
        self.validate(*self.run(state, timestep))
        # self.run(state, timestep)

    @ported_method(
        from_file="cloudsc2_ad/cloudsc_driver_ad_mod.F90",
        from_line=22,
        to_line=297,
    )
    def run(
        self, state: "DataArrayDict", timestep: "timedelta"
    ) -> Tuple[
        "DataArrayDict",
        "DataArrayDict",
        "DataArrayDict",
        "DataArrayDict",
        "DataArrayDict",
    ]:
        diags_sat = self.saturation(state)
        state.update(diags_sat)

        state_i = self.state_increment(state)
        state.update(state_i)
        tends_tl, diags_tl = self.cloudsc2_tl(state, timestep)

        self.add_tendencies_to_state(state, tends_tl)
        state.update(diags_tl)
        tends_ad, diags_ad = self.cloudsc2_ad(state, timestep)

        # field = tends_ad["f_cml_ql_i"].data
        # for k in range(137, -1, -1):
        #     print(k+1, " ", field[0, 0, k])

        return state_i, tends_tl, diags_tl, tends_ad, diags_ad

    @ported_method(
        from_file="cloudsc2_ad/cloudsc_driver_ad_mod.F90",
        from_line=260,
        to_line=294,
    )
    def validate(
        self,
        state_i: "DataArrayDict",
        tends_tl: "DataArrayDict",
        diags_tl: "DataArrayDict",
        tends_ad: "DataArrayDict",
        diags_ad: "DataArrayDict",
    ) -> None:
        norm1 = self.get_norm1(tends_tl, diags_tl)
        norm2 = self.get_norm2(state_i, tends_ad, diags_ad)

        # norm3 = np.where(
        #     norm2 == 0,
        #     abs((norm1 - norm2) / sys.float_info.epsilon),
        #     abs((norm1 - norm2) / (sys.float_info.epsilon * norm2)),
        # )
        # if norm3.max() < 1e4:
        #     print("The symmetry test passed. HOORAY!")
        # else:
        #     print("The symmetry test failed.")
        # print(f"The error is {norm3.max()} times the machine epsilon.")

        out = np.where(np.isclose(norm1, norm2) == False)[0]
        if out.size == 0:
            print("The symmetry test passed. HOORAY!")
        else:
            print(
                f"The symmetry test failed on the following columns: "
                f"{', '.join([str(idx+1) for idx in out])}."
            )

    @ported_method(
        from_file="cloudsc2_ad/cloudsc_driver_ad_mod.F90",
        from_line=183,
        to_line=195,
    )
    def get_norm1(
        self, tends_tl: "DataArrayDict", diags_tl: "DataArrayDict"
    ) -> np.ndarray:
        out = None

        tend_names = ("f_t_i", "f_q_i", "f_ql_i", "f_qi_i")
        for name in tend_names:
            field = self.get_field(name, tends_tl)
            out = np.zeros(field.shape[0]) if out is None else out
            out += np.sum(field * field, axis=1)

        diag_names = (
            "f_clc_i",
            "f_fhpsl_i",
            "f_fhpsn_i",
            "f_fplsl_i",
            "f_fplsn_i",
            "f_covptot_i",
        )
        for name in diag_names:
            field = self.get_field(name, diags_tl)
            out += np.sum(field ** 2, axis=1)

        return out

    @ported_method(
        from_file="cloudsc2_ad/cloudsc_driver_ad_mod.F90",
        from_line=239,
        to_line=256,
    )
    def get_norm2(
        self,
        state_i: "DataArrayDict",
        tends_ad: "DataArrayDict",
        diags_ad: "DataArrayDict",
    ) -> np.ndarray:
        out = None

        tend_names = ("f_cml_t_i", "f_cml_q_i", "f_cml_ql_i", "f_cml_qi_i")
        for name in tend_names:
            field_a = self.get_field("f_tnd_" + name[2:], state_i)
            field_b = self.get_field(name, tends_ad)
            out = np.zeros(field_a.shape[0]) if out is None else out
            out += np.sum(field_a * field_b, axis=1)

        diag_names = (
            "f_ap_i",
            "f_aph_i",
            "f_t_i",
            "f_q_i",
            "f_qsat_i",
            "f_ql_i",
            "f_qi_i",
            "f_lu_i",
            "f_lude_i",
            "f_mfd_i",
            "f_mfu_i",
            "f_supsat_i",
        )
        for name in diag_names:
            field_a = self.get_field(name, state_i)
            field_b = self.get_field(name, diags_ad)
            out += np.sum(field_a * field_b, axis=1)

        return out

    @staticmethod
    def get_field(name: str, dct: "DataArrayDict") -> np.ndarray:
        gt4py.storage.restore_numpy()
        field = dct[name].data
        field.synchronize()
        field = np.asarray(field)[:, 0, :]
        gt4py.storage.prepare_numpy()
        return field

    @staticmethod
    def add_tendencies_to_state(
        state: "DataArrayDict", tends_tl: "DataArrayDict"
    ) -> None:
        state["f_tnd_t"] = tends_tl["f_t"]
        state["f_tnd_t_i"] = tends_tl["f_t_i"]
        state["f_tnd_q"] = tends_tl["f_q"]
        state["f_tnd_q_i"] = tends_tl["f_q_i"]
        state["f_tnd_ql"] = tends_tl["f_ql"]
        state["f_tnd_ql_i"] = tends_tl["f_ql_i"]
        state["f_tnd_qi"] = tends_tl["f_qi"]
        state["f_tnd_qi_i"] = tends_tl["f_qi_i"]
