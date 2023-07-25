from collections import namedtuple
import gt4py


Field3d = gt4py.gtscript.Field[gt4py.gtscript.IJK, float]


TFACTOR = namedtuple("TFACTOR", "value")


@gt4py.gtscript.function
def core(a, b):
    from __externals__ import FACTOR

    return a + FACTOR.value * b


def apfb_def(a: Field3d, b: Field3d, c: Field3d):
    from __externals__ import core

    with computation(PARALLEL), interval(...):
        c = core(a, b)


def main():
    apfb = gt4py.gtscript.stencil(
        "numpy", apfb_def, externals={"FACTOR": TFACTOR(value=1.0), "core": core}
    )

    # a = gt4py.storage.ones("gtmc", (0, 0, 0), (30, 30, 30), float)
    # b = gt4py.storage.zeros("gtmc", (0, 0, 0), (30, 30, 30), float)
    # c = gt4py.storage.zeros("gtmc", (0, 0, 0), (30, 30, 30), float)
    # foo(a, b, c[:, 0, :], origin=(0, 0, 0), domain=(30, 30, 30))
    #
    # print("done")


if __name__ == "__main__":
    main()
