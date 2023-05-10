# -*- coding: utf-8 -*-
import numpy as np
from typing import Dict, TypeVar, Union

from sympl import DataArray as SymplDataArray

try:
    import cupy as cp
except (ImportError, ModuleNotFoundError):
    cp = np


Array = Union[np.ndarray, cp.ndarray]
ArrayDict = Dict[str, Array]
DataArray = SymplDataArray
DataArrayDict = Dict[str, DataArray]
ParameterDict = Dict[str, Union[bool, float, int]]
Range = TypeVar("Range")
