# Direct Method
import numpy as np
from math import log
from reaction import reaction
import numba
import pathos
from tqdm import tqdm


@numba.njit
def fast_sum(propensities):
    return np.sum(propensities)


@numba.njit
def draw(probabilities, r2):
    j = 0
    p_sum = 0.0
    while p_sum < r2:
        p_sum += probabilities[j]
        j += 1
    j -= 1
    return j


def simulate_wrapper(args):
    return simulate_run(*args)


def simulate(
    start_state: list[int],
    end_time: float,
    reactions: list[reaction],
    n_repeats: int,
    n_threads: int,
):

    with pathos.multiprocessing.Pool(n_threads) as pool:
        arguments = [(start_state, end_time, reactions, i) for i in range(n_repeats)]
        results = pool.map(simulate_wrapper, arguments)

    return results


def simulate_run(
    start_state: list[int], end_time: float, reactions: list[reaction], seed: int
) -> list[tuple[float, int]]:

    # Initialize. Set the initial number of molecules for each species and set t = 0.
    np.random.seed(seed)
    t = 0.0
    # initial_state = start_state.copy()
    results: list[tuple[float, int]] = []
    state = start_state.copy()
    propensities = np.array([0.0] * len(reactions))
    probabilities = np.array([0.0] * len(reactions))

    pbar = tqdm(total=100)
    chunks = np.linspace(0, end_time, 21)
    chunk_idx = 0

    while t < end_time:

        if t > chunks[chunk_idx]:
            pbar.update(5)
            chunk_idx += 1

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
        r1 = np.random.rand()
        r2 = np.random.rand()

        # Set delta = 1 / A * ln( 1 / r1)
        delta = (1 / A) * log(1 / r1)

        # Find mu such that sum_{k=1}^{mu-1} a_k < r2 * A <= sum_{k=1}^{mu} a_k
        j = draw(probabilities, r2)

        # Set t = t + delta and update the number of each molecular species according to reaction mu.
        state += reactions[j].update
        t += delta
        results.append((t, j))
        # Return to step 2 or quit.
    pbar.close()
    return results
