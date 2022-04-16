# Next Reaction Method
import numpy as np
from math import log
from typing import Callable

class reaction():

    def __init__(self, update: list[int], propensity: Callable[[list[int]], float]):
        self.update = update
        self.prop_func = propensity 

    def propensity(self, state: list[int]) -> float:
        return self.prop_func(state)



def simulate(start_state: list[int], end_time: float, reactions: list[reaction]) -> None:
    
    # Initialize. Set the initial number of molecules for each species and set t = 0.
    t = 0.0
    # initial_state = start_state.copy()
    state = start_state.copy()
    propensities = np.array([0.] * len(reactions))
    probabilities = np.array([0.] * len(reactions))


    while t < end_time:
        # Calculate the propensity function, a_k, for each reaction.
        for idx, r in enumerate(reactions):
            propensities[idx] = r.prop_func

        # Set A = sum_{k=1}^{M} a_k
        A = np.sum(propensities)

        # Calculate the probabilities of the reactions
        for idx, val in enumerate(propensities):
            probabilities[idx] = val / A

        # Generate two independent uniform (0,1) random numbers r_1 and r_2.
        r1 = np.random.rand()
        r2 = np.random.rand()

        # Set delta = 1 / A * ln( 1 / r1)
        delta = (1 / A) * log(1 / r1)

        # Find mu such that sum_{k=1}^{mu-1} a_k < r2 * A <= sum_{k=1}^{mu} a_k
        j = 0;
        p_sum = 0.;
        while p_sum < r2:
            p_sum += probabilities[j];
            j += 1
        j -= 1;

        # Set t = t + delta and update the number of each molecular species according to reaction mu.
        state += reactions[j].update
        t += delta
        print(f"Time: {t} \t State: {state}")
        # Return to step 2 or quit.
    return
