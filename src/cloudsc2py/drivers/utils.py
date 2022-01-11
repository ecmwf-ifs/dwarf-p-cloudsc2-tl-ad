# -*- coding: utf-8 -*-
from typing import Dict, Sequence, TYPE_CHECKING, Type, Union

if TYPE_CHECKING:
    from cloudsc2py.utils.timing import Timer


def log_performance(
    nruns: int,
    timer: Type["Timer"],
    exec_info: Dict[str, Union[bool, float]],
    stencil_names: Sequence[str],
) -> None:
    total_time = timer.get_time("run", units="ms")
    cpp_time = 0.0
    call_time = 0.0
    for name in stencil_names:
        for key in exec_info:
            if key.startswith(name):
                cpp_time += exec_info[key]["total_run_cpp_time"] * 1e3
                call_time += exec_info[key]["total_call_time"] * 1e3
                break

    print(
        f"Average run time ({nruns} runs): {total_time / nruns:.3f} ms\n"
        f"  - GT4Py (stencil calculations): {cpp_time / nruns:.3f} ms\n"
        f"  - GT4Py (bindings overhead): {(call_time - cpp_time) / nruns:.3f} ms\n"
        f"  - Framework: {(total_time - call_time) / nruns:.3f} ms"
    )
