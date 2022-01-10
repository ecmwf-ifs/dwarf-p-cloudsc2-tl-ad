# -*- coding: utf-8 -*-
from sympl._core.time import Timer as SymplTimer


class TimerContextManager:
    def __init__(self, label: str) -> None:
        self.label = label

    def __enter__(self) -> "Timer":
        Timer.start(self.label)
        return Timer

    def __exit__(self, exc_type, exc_value, exc_tb) -> None:
        Timer.stop()


class Timer(SymplTimer):
    timing = TimerContextManager
