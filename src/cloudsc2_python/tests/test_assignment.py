import gt4py

from src.cloudsc2_python.utils import Timer


Field3d = gt4py.gtscript.Field[gt4py.gtscript.IJK, float]


def copy_def(a: Field3d, b: Field3d):
    with computation(PARALLEL), interval(...):
        b[0, 0, 0] = a[0, 0, 0]


def main():
    backend = "gt:gpu"
    shape = (201, 31, 61)

    copy = gt4py.gtscript.stencil(backend, copy_def)

    a = gt4py.storage.ones(backend, (0, 0, 0), shape, float)
    a.synchronize()
    b = gt4py.storage.zeros(backend, (0, 0, 0), shape, float)
    b.synchronize()

    Timer.start("direct assignment")
    for _ in range(100):
        b[...] = a[...]
    Timer.stop()

    Timer.start("copy stencil")
    for _ in range(100):
        copy(a, b, origin=(0, 0, 0), domain=shape, validate_args=False)
    Timer.stop()

    print(f"backend: {backend}")
    Timer.print("direct assignment")
    Timer.print("copy stencil")


if __name__ == "__main__":
    main()
