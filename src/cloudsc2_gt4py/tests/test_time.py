# -*- coding: utf-8 -*-
import time

from src.cloudsc2_python.utils import Timer


def main():
    with Timer.timing("a"):
        time.sleep(1)

    Timer.start("b")
    time.sleep(1)
    Timer.stop()

    Timer.print("a")
    Timer.print("b")

if __name__ == "__main__":
    main()


