# -*- coding: utf-8 -*-
from __future__ import annotations
from typing import TYPE_CHECKING

from gt4py import gtscript

from cloudsc2py.framework.options import fill_dtypes

if TYPE_CHECKING:
    from typing import Any, Optional

    from gt4py import StencilObject

    from cloudsc2py.framework.options import BackendOptions, StorageOptions


FUNCTION_COLLECTION = {}
STENCIL_COLLECTION = {}


def function_collection(name: str):
    if name in FUNCTION_COLLECTION:
        raise RuntimeError(f"Another function called `{name}` found.")

    def core(definition):
        FUNCTION_COLLECTION[name] = {"definition": definition}
        return definition

    return core


def stencil_collection(name: str):
    if name in STENCIL_COLLECTION:
        raise RuntimeError(f"Another stencil called `{name}` found.")

    def core(definition):
        STENCIL_COLLECTION[name] = {"definition": definition}
        return definition

    return core


def compile_stencil(
    name: str,
    backend: str,
    backend_options: BackendOptions,
    externals: dict[str, Any] = None,
    storage_options: Optional[StorageOptions] = None,
) -> StencilObject:
    stencil_info = STENCIL_COLLECTION.get(name, None)
    if stencil_info is None:
        raise RuntimeError(f"Unknown stencil `{name}`.")
    definition = stencil_info["definition"]

    if storage_options:
        fill_dtypes(backend_options, storage_options)

    externals = externals or {}

    kwargs = backend_options.backend_opts.copy()
    if backend not in ("debug", "numpy", "gtc:numpy"):
        kwargs["verbose"] = backend_options.verbose

    return gtscript.stencil(
        backend,
        definition,
        name=name,
        build_info=backend_options.build_info,
        dtypes=backend_options.dtypes,
        externals=externals,
        rebuild=backend_options.rebuild,
        **kwargs,
    )
