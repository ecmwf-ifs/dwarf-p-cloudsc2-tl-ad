# -*- coding: utf-8 -*-
from typing import Optional, Sequence, TYPE_CHECKING

from gt4py import gtscript

from cloudsc2py.framework.options import fill_dtypes

if TYPE_CHECKING:
    from gt4py import StencilObject

    from cloudsc2py.framework.components import GridComponent
    from cloudsc2py.framework.options import BackendOptions, StorageOptions


FUNCTION_COLLECTION = {}
PARAMETER_COLLECTION = {}
STENCIL_COLLECTION = {}


def function_collection(name: str):
    if name in FUNCTION_COLLECTION:
        raise RuntimeError(f"Another function called {name} found.")

    def core(definition):
        FUNCTION_COLLECTION[name] = definition
        return definition

    return core


def stencil_collection(name: str):
    if name in STENCIL_COLLECTION:
        raise RuntimeError(f"Another stencil called {name} found.")

    def core(definition):
        STENCIL_COLLECTION[name] = definition
        return definition

    return core


def compile_stencil(
    name: str,
    backend: str,
    backend_options: "BackendOptions",
    storage_options: Optional["StorageOptions"] = None,
    used_externals: Optional[Sequence[str]] = None,
) -> "StencilObject":
    definition = STENCIL_COLLECTION.get(name, None)
    if definition is None:
        raise RuntimeError(f"Unknown stencil {name}.")

    if storage_options:
        fill_dtypes(backend_options, storage_options)

    if used_externals:
        collected_externals = PARAMETER_COLLECTION.copy()
        collected_externals.update(backend_options.external_parameters)
        collected_externals.update(FUNCTION_COLLECTION)
        collected_externals.update(backend_options.external_functions)
        externals = {key: collected_externals[key] for key in used_externals}
    else:
        externals = {}

    kwargs = backend_options.backend_opts.copy()
    if backend not in ("debug", "numpy", "gtc:numpy"):
        kwargs["verbose"] = backend_options.verbose

    return gtscript.stencil(
        backend,
        definition,
        build_info=backend_options.build_info,
        dtypes=backend_options.dtypes,
        externals=externals,
        rebuild=backend_options.rebuild,
        **kwargs,
    )
