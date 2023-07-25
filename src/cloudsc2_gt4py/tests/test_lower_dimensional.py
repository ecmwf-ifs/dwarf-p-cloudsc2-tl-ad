import gt4py


Field2d = gt4py.gtscript.Field[gt4py.gtscript.IK, float]
Field3d = gt4py.gtscript.Field[gt4py.gtscript.IJK, float]


def foo_def(a: Field3d, b: Field3d, c: Field2d):
    # with computation(FORWARD), interval(...):
    #     c = a if a > 0.5 else -a
    with computation(PARALLEL), interval(...):
        b = a + c


def main():
    foo = gt4py.gtscript.stencil("gtmc", foo_def)

    a = gt4py.storage.ones("gtmc", (0, 0, 0), (30, 30, 30), float)
    b = gt4py.storage.zeros("gtmc", (0, 0, 0), (30, 30, 30), float)
    c = gt4py.storage.zeros("gtmc", (0, 0, 0), (30, 30, 30), float)
    foo(a, b, c[:, 0, :], origin=(0, 0, 0), domain=(30, 30, 30))

    print("done")


if __name__ == "__main__":
    main()
