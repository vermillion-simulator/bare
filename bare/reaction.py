from typing import Callable
import numpy as np
import numba


class reaction():
    def __init__(self, update: list[int], propensity: Callable[[list[int]], float]):
        self.update = np.array(update)
        self.prop_func = propensity 

    @numba.jit
    def propensity(self, state: list[int]) -> float:
        return self.prop_func(state)