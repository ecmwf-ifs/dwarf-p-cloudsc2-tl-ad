# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
import sys
from typing import TYPE_CHECKING

import gt4py

from cloudsc2py.physics.common.increment import PerturbedState, StateIncrement
from cloudsc2py.physics.common.saturation import Saturation
from cloudsc2py.physics.nonlinear.microphysics import Cloudsc2NL
from cloudsc2py.physics.tangent_linear.microphysics import Cloudsc2TL
from cloudsc2py.utils.f2py import ported_method

if TYPE_CHECKING:
    from datetime import timedelta
    from typing import Optional

    from sympl._core.typingx import DataArrayDict

    from cloudsc2py.framework.config import GT4PyConfig
    from cloudsc2py.framework.grid import ComputationalGrid
    from cloudsc2py.utils.typingx import ParameterDict


class TaylorTest:
    def __init__(
        self,
        computational_grid: ComputationalGrid,
        factor1: float,
        factor2s: tuple[float, ...],
        kflag: int,
        lphylin: bool,
        ldrain1d: bool,
        yoethf_parameters: Optional[ParameterDict] = None,
        yomcst_parameters: Optional[ParameterDict] = None,
        yrecld_parameters: Optional[ParameterDict] = None,
        yrecldp_parameters: Optional[ParameterDict] = None,
        yrephli_parameters: Optional[ParameterDict] = None,
        yrncl_parameters: Optional[ParameterDict] = None,
        yrphnc_parameters: Optional[ParameterDict] = None,
        *,
        enable_checks: bool = True,
        gt4py_config: GT4PyConfig,
    ) -> None:
        self.f1 = factor1
        self.f2s = factor2s

        # saturation
        self.saturation = Saturation(
            computational_grid,
            kflag,
            lphylin,
            yoethf_parameters,
            yomcst_parameters,
            enable_checks=enable_checks,
            gt4py_config=gt4py_config,
        )

        # microphysics
        self.cloudsc2_nl = Cloudsc2NL(
            computational_grid,
            lphylin,
            ldrain1d,
            yoethf_parameters,
            yomcst_parameters,
            yrecld_parameters,
            yrecldp_parameters,
            yrephli_parameters,
            yrphnc_parameters,
            enable_checks=enable_checks,
            gt4py_config=gt4py_config,
        )
        self.cloudsc2_tl = Cloudsc2TL(
            computational_grid,
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
            gt4py_config=gt4py_config,
        )

        # perturbation
        self.state_increment = StateIncrement(
            computational_grid, factor1, enable_checks=enable_checks, gt4py_config=gt4py_config
        )
        self.perturbed_states = [
            PerturbedState(
                computational_grid, factor2, enable_checks=enable_checks, gt4py_config=gt4py_config
            )
            for factor2 in factor2s
        ]

        # ausiliary dicts
        self.diags_nl = None
        self.diags_nl_p = None
        self.diags_tl = None
        self.state_p = None
        self.tends_nl = None
        self.tends_nl_p = None
        self.tends_tl = None

    def __call__(self, state: DataArrayDict, timestep: timedelta) -> None:
        self.validate(self.run(state, timestep))

    @ported_method(
        from_file="cloudsc2_tl/cloudsc_driver_tl_mod.F90",
        from_line=126,
        to_line=231,
    )
    def run(self, state: DataArrayDict, timestep: timedelta) -> np.ndarray:
        diags_sat = self.saturation(state)
        state.update(diags_sat)

        self.tends_nl, self.diags_nl = self.cloudsc2_nl(state, timestep)

        state_i = self.state_increment(state)
        state.update(state_i)
        self.tends_tl, self.diags_tl = self.cloudsc2_tl(state, timestep)

        norms = np.zeros(len(self.f2s))
        for i, perturbed_state in enumerate(self.perturbed_states):
            state_p = perturbed_state(state, out=self.state_p)
            state_p["time"] = state["time"]
            state_p["f_eta"] = state["f_eta"]
            self.tends_nl_p, self.diags_nl_p = self.cloudsc2_nl(
                state_p, timestep, out_tendencies=self.tends_nl_p, out_diagnostics=self.diags_nl_p
            )
            norms[i] = self.get_norm(i)

        return norms

    @ported_method(from_file="cloudsc2_tl/cloudsc_driver_tl_mod.F90", from_line=275, to_line=313)
    def validate(self, norms: np.ndarray) -> None:
        print(">>> Taylor test: Start")
        start = -1
        for i in range(norms.size):
            print(
                f"  factor1 = {self.f1:.3e}, factor2 = {self.f2s[i]:.3e}, "
                f"norm = {norms[i]:.10f}"
            )
            norms[i] = np.abs(1 - norms[i])
            if start == -1 and norms[i] < 0.5:
                start = i

        if start == -1 or start > 3:
            log = "The test failed with error 13."
        else:
            test = -10
            negat = 1
            for i in range(start, norms.size - 1):
                tmp_negat = (norms[i + 1] / norms[i]) < 1
                if negat > tmp_negat:
                    test += 10
                negat = tmp_negat
            if test == -10:
                test = 11
            if np.min(norms[start:]) > 1e-5:
                test += 7
            if np.min(norms[start:]) > 1e-6:
                test += 5
            if test > 5:
                log = f"The test failed with error {test}."
            else:
                log = f"The test passed with penalty {test}. HOORAY!"

        print("<<< Taylor test: End")
        print(log)

    @ported_method(from_file="cloudsc2_tl/cloudsc_driver_tl_mod.F90", from_line=233, to_line=245)
    def get_norm(self, i: int) -> float:
        assert self.diags_nl is not None, "Did you execute the run() method?"

        total_count = 0
        total_norm = 0.0

        tend_names = ("f_t", "f_q", "f_ql", "f_qi")
        for name in tend_names:
            field_nl, field_nl_p, field_tl = self.get_fields(name, "tends")
            norm = self.get_field_norm(i, field_nl, field_nl_p, field_tl)
            total_count += norm > 0
            total_norm += norm

        diag_names = ("f_clc", "f_fhpsl", "f_fhpsn", "f_fplsl", "f_fplsn", "f_covptot")
        for name in diag_names:
            field_nl, field_nl_p, field_tl = self.get_fields(name, "diags")
            norm = self.get_field_norm(i, field_nl, field_nl_p, field_tl)
            total_count += norm > 0
            total_norm += norm

        return total_norm / total_count if total_count > 0 else 0

    def get_fields(self, name: str, dct_name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        gt4py.storage.restore_numpy()

        dct_nl = getattr(self, dct_name + "_nl")
        field_nl = dct_nl[name].data
        field_nl.synchronize()
        field_nl = np.asarray(field_nl)[:, 0, :]

        dct_nl_p = getattr(self, dct_name + "_nl_p")
        field_nl_p = dct_nl_p[name].data
        field_nl_p.synchronize()
        field_nl_p = np.asarray(field_nl_p)[:, 0, :]

        dct_tl = getattr(self, dct_name + "_tl")
        field_tl = dct_tl[name + "_i"].data
        field_tl.synchronize()
        field_tl = np.asarray(field_tl)[:, 0, :]

        gt4py.storage.prepare_numpy()

        return field_nl, field_nl_p, field_tl

    @ported_method(from_file="cloudsc2_tl/cloudsc_driver_tl_mod.F90", from_line=21, to_line=31)
    def get_field_norm(
        self, i: int, field_nl: np.ndarray, field_nl_p: np.ndarray, field_tl: np.ndarray
    ) -> float:
        den = np.abs(self.f2s[i] * np.sum(field_tl))
        if den > sys.float_info.epsilon:
            norm = np.abs(np.sum(field_nl_p - field_nl)) / den
        else:
            norm = 0
        return norm
