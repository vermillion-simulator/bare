# Direct Method
import numpy as np
from math import log
from reaction import reaction
import numba

@numba.njit
def fast_sum(propensities):
    return np.sum(propensities)

@numba.njit
def draw(probabilities, r2):
    j = 0
    p_sum = 0.
    while p_sum < r2:
        p_sum += probabilities[j]
        j += 1
    j -= 1
    return j

def simulate(start_state: list[int], end_time: float, reactions: list[reaction]) -> list[tuple[float, int]]:
    
    # Initialize. Set the initial number of molecules for each species and set t = 0.
    r_MAX = 101
    t = 0.0
    #initial_state = start_state.copy()
    results: list[tuple[float,int]] = []
    state = start_state.copy()
    propensities = np.array([0.] * len(reactions))
    probabilities = np.array([0.] * len(reactions))
    r1_lst = np.random.rand(r_MAX)
    r2_lst = np.random.rand(r_MAX)
    r_idx = 0

    while t < end_time:
        # Calculate the propensity function, a_k, for each reaction.
        for idx, r in enumerate(reactions):
            propensities[idx] = r.prop_func(state)

        # Set A = sum_{k=1}^{M} a_k
        A = fast_sum(propensities)
        if A == 0:
            break

        # Calculate the probabilities of the reactions
        probabilities = propensities / A

        # Generate two independent uniform (0,1) random numbers r_1 and r_2.
        if r_idx >= r_MAX:
            r1_lst = np.random.rand(r_MAX)
            r2_lst = np.random.rand(r_MAX)
            r_idx = 0
        r1 = r1_lst[r_idx]
        r2 = r2_lst[r_idx]
        r_idx += 1

        # Set delta = 1 / A * ln( 1 / r1)
        delta = (1 / A) * log(1 / r1)

        # Find mu such that sum_{k=1}^{mu-1} a_k < r2 * A <= sum_{k=1}^{mu} a_k
        j = draw(probabilities, r2)

        # Set t = t + delta and update the number of each molecular species according to reaction mu.
        state += reactions[j].update
        t += delta
        results.append((t, j))
        # Return to step 2 or quit.
    #df = reconstruct(results, reactions, initial_state)
    return results
