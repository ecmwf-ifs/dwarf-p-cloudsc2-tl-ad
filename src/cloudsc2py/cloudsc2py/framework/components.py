# -*- coding: utf-8 -*-
import abc
from typing import Optional, Sequence, TYPE_CHECKING

from sympl._core.core_components import (
    DiagnosticComponent as SymplDiagnosticComponent,
    ImplicitTendencyComponent as SymplImplicitTendencyComponent,
)

from cloudsc2py.framework.options import (
    BackendOptions,
    StorageOptions,
    fill_dtypes,
)
from cloudsc2py.framework.stencil import compile_stencil
from cloudsc2py.utils.storage import (
    get_array,
    get_data_shape_from_name,
    get_dtype_from_name,
)

if TYPE_CHECKING:
    from sympl._core.typingx import PropertyDict

    from gt4py import StencilObject

    from cloudsc2py.framework.grid import Grid
    from cloudsc2py.utils.typingx import Array


class GridComponent:
    def __init__(
        self,
        grid: "Grid",
        *,
        backend: str = "numpy",
        backend_options: Optional[BackendOptions] = None,
        storage_options: Optional[StorageOptions] = None,
    ) -> None:
        self.grid = grid
        self.backend = backend
        self.bo = backend_options or BackendOptions()
        self.so = storage_options or StorageOptions()
        fill_dtypes(self.bo, self.so)

    def compile_stencil(self, name: str) -> "StencilObject":
        return compile_stencil(
            name, self.backend, self.bo, used_externals=self.used_externals
        )

    def get_field_shape(
        self, name: str, properties: "PropertyDict"
    ) -> Sequence[int]:
        grid_shape = tuple(
            self.grid.dims_to_shape[dim]
            for dim in properties[name]["dims"]
            if dim in self.grid.dims_to_shape
        )
        data_shape = get_data_shape_from_name(name)
        return *grid_shape, *data_shape

    def allocate(self, name: str, properties: "PropertyDict") -> "Array":
        shape = self.get_field_shape(name, properties)
        dtype = get_dtype_from_name(name, self.so)
        return get_array(
            shape, backend=self.backend, dtype=dtype, storage_options=self.so
        )

    @abc.abstractmethod
    def used_externals(self) -> Sequence[str]:
        pass


class DiagnosticComponent(GridComponent, SymplDiagnosticComponent):
    def __init__(
        self,
        grid: "Grid",
        *,
        enable_checks: bool = True,
        backend: str = "numpy",
        backend_options: Optional[BackendOptions] = None,
        storage_options: Optional[StorageOptions] = None,
    ) -> None:
        super().__init__(
            grid,
            backend=backend,
            backend_options=backend_options,
            storage_options=storage_options,
        )
        super(GridComponent, self).__init__(enable_checks=enable_checks)

    def allocate_diagnostic(self, name: str) -> "Array":
        return self.allocate(name, self.diagnostic_properties)


class ImplicitTendencyComponent(GridComponent, SymplImplicitTendencyComponent):
    def __init__(
        self,
        grid: "Grid",
        *,
        enable_checks: bool = True,
        backend: str = "numpy",
        backend_options: Optional[BackendOptions] = None,
        storage_options: Optional[StorageOptions] = None,
    ) -> None:
        super().__init__(
            grid,
            backend=backend,
            backend_options=backend_options,
            storage_options=storage_options,
        )
        super(GridComponent, self).__init__(enable_checks=enable_checks)

    def allocate_tendency(self, name: str) -> "Array":
        return self.allocate(name, self.tendency_properties)

    def allocate_diagnostic(self, name: str) -> "Array":
        return self.allocate(name, self.diagnostic_properties)
