# `cloudsc2py`: GT4Py-based implementation of the CLOUDSC2 dwarf

This folder contains a Python implementation of the CLOUDSC2 microphysics scheme in its non-linear,
tangent-linear and adjoint formulations. The implementation is based on
[GT4Py](https://github.com/GridTools/gt4py.git). The code is bundled as an installable
package called `cloudsc2py`, whose source code is placed under `src/`.


## Installation

We strongly recommend installing the package in an isolated virtual environment:

```shell
# create virtual environment under `venv/`
$ python -m venv venv

# activate the virtual environment
$ . venv/bin/activate

# upgrade base packages
(venv) $ pip install --upgrade pip setuptools wheel
```

The package `cloudsc2py` can be installed using the Python package manager [pip](https://pip.pypa.io/en/stable/):

```shell
# install cloudsc2py along with the minimal set of requirements
(venv) $ pip install .

# create a fully reproducible development environment with an editable installation of cloudsc4py
(venv) $ pip install -r requirements-dev.txt
(venv) $ pip install -e .

# enable the GPU backends of GT4Py on an NVIDIA GPU using CUDA 11.x
(venv) $ pip install .[gpu-cuda11x]

# enable the GPU backends of GT4Py on an NVIDIA GPU using CUDA 12.x
(venv) $ pip install .[gpu-cuda12x]

# enable the GPU backends of GT4Py on an AMD GPU using ROCm/HIP
(venv) $ export CUPY_INSTALL_USE_HIP=1
(venv) $ export ROCM_HOME=<path to ROCm installation>
(venv) $ export HCC_AMDGPU_TARGET=<string denoting the Instruction Set Architecture (ISA) supported by the target GPU>
(venv) $ pip install .[gpu]
```

## Usage

The easiest way to run the dwarf is through the driver scripts contained in `drivers/`:

* `drivers/run_nonlinear.py` executes the non-linear formulation of CLOUDSC2;
* `drivers/run_taylor_test.py` performs the Taylor test for the tangent-linear formulation of CLOUDSC2;
* `drivers/run_symmetry_test.py` performs the symmetry test for the adjoint formulation of CLOUDSC2.

For the sake of convenience, we provide the driver `drivers/run_fortran.py` to invoke one of the
FORTRAN variants of the dwarf from within Python.

Run the two scripts with the `--help` option to get the full list of command-line options.
