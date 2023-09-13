# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING

from cloudsc2py.utils.iox import HDF5Reader
from ifs_physics_common.utils.numpyx import to_numpy

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from typing import Optional

    from sympl._core.typingx import DataArrayDict

    from ifs_physics_common.framework.config import DataTypes
    from ifs_physics_common.utils.typingx import ArrayLike


class Validator:
    atol: float
    reader: HDF5Reader
    rtol: float

    def __init__(
        self,
        reference_filename: str,
        data_types: DataTypes,
        atol: float = 1e-16,
        rtol: float = 1e-12,
    ) -> None:
        self.reader = HDF5Reader(reference_filename, data_types)
        self.atol = atol
        self.rtol = rtol

    def __call__(self, tendencies: DataArrayDict, diagnostics: DataArrayDict) -> list[str]:
        failing_fields = []

        # tendencies
        trg_keys = {
            "f_q": "TENDENCY_LOC_Q",
            "f_qi": "TENDENCY_LOC_CLD",
            "f_ql": "TENDENCY_LOC_CLD",
            "f_t": "TENDENCY_LOC_T",
        }
        data_indices = {"f_q": None, "f_qi": 1, "f_ql": 0, "f_t": None}
        for src_key, trg_key in trg_keys.items():
            if not self.validate_field(tendencies, src_key, trg_key, data_indices[src_key]):
                failing_fields.append(src_key)

        # diagnostics
        trg_keys = {
            "f_covptot": "PCOVPTOT",
            "f_fhpsl": "PFHPSL",
            "f_fhpsn": "PFHPSN",
            "f_fplsl": "PFPLSL",
            "f_fplsn": "PFPLSN",
        }
        for src_key, trg_key in trg_keys.items():
            if not self.validate_field(diagnostics, src_key, trg_key):
                failing_fields.append(src_key)

        return failing_fields

    def validate_field(
        self,
        src_dict: DataArrayDict,
        src_key: str,
        trg_key: str,
        trg_data_index: Optional[int] = None,
    ) -> bool:
        src = src_dict[src_key].data
        trg = self.reader.get_field(trg_key)
        if trg.ndim == 3 and trg_data_index is not None:
            trg = trg[..., trg_data_index]
        return self.compare_field(src, trg)

    def compare_field(self, src: ArrayLike, trg: NDArray) -> bool:
        src = to_numpy(src)[:, 0, :]
        mi = min(src.shape[0], trg.shape[0])
        mk = min(src.shape[1], trg.shape[1])
        out = np.allclose(
            src[:mi, :mk], trg[:mi, :mk], atol=self.atol, rtol=self.rtol, equal_nan=True
        )
        return out
