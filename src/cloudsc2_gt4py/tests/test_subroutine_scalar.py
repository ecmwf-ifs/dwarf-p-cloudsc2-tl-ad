from collections import namedtuple
import gt4py


Field3d = gt4py.gtscript.Field[gt4py.gtscript.IJK, float]


@gt4py.gtscript.function
def core(a, b):
    return a + b


def foo_def(a: Field3d, b: Field3d, c: Field3d, d: Field3d):
    from __externals__ import core

    with computation(PARALLEL), interval(0, -1):
        c = core(a, b[0, 0, 1])
        d = core(a, 0.0)


def main():
    foo = gt4py.gtscript.stencil("gtc:gt:cpu_ifirst", foo_def, externals={"core": core})

    a = gt4py.storage.ones("gtc:gt:cpu_ifirst", (0, 0, 0), (30, 30, 30), float)
    b = gt4py.storage.zeros("gtc:gt:cpu_ifirst", (0, 0, 0), (30, 30, 30), float)
    c = gt4py.storage.zeros("gtc:gt:cpu_ifirst", (0, 0, 0), (30, 30, 30), float)
    d = gt4py.storage.zeros("gtc:gt:cpu_ifirst", (0, 0, 0), (30, 30, 30), float)
    foo(a, b, c, d, origin=(0, 0, 0), domain=(30, 30, 30))

    print("done")


if __name__ == "__main__":
    main()
