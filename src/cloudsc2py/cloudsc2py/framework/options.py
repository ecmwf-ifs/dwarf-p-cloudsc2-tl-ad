# -*- coding: utf-8 -*-
from collections import namedtuple
from dataclasses import dataclass, field
import numpy as np
from typing import Any, Mapping, Optional, Sequence, TYPE_CHECKING, Type, Union


DataTypes = namedtuple("DataTypes", ("bool", "float", "integer"))


@dataclass
class BackendOptions:
    backend_opts: Mapping[str, Any] = field(default_factory=dict)
    build_info: Optional[Mapping[str, Any]] = None
    device_sync: bool = True
    dtypes: Mapping[str, Type] = field(default_factory=dict)
    exec_info: Optional[Mapping[str, Any]] = None
    external_functions: Mapping[str, Any] = field(default_factory=dict)
    external_parameters: Mapping[str, Any] = field(default_factory=dict)
    rebuild: bool = False
    validate_args: bool = False
    verbose: bool = True


@dataclass
class StorageOptions:
    default_origin: Sequence[int] = (0, 0, 0)
    dtypes: DataTypes = DataTypes(bool=bool, float=float, integer=int)
    halo: Optional[Sequence[int]] = None
    managed: Union[bool, str] = "gt4py"


def fill_dtypes(bo: BackendOptions, so: StorageOptions) -> None:
    bo.dtypes = {
        "btype": so.dtypes.bool,
        "ftype": so.dtypes.float,
        "itype": so.dtypes.integer,
    }
