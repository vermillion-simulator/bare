import numpy as np
from typing import Tuple
from reaction import reaction
from math import log
import numba
from tqdm import tqdm


def simulate_run(
    start_state: list[int], end_time: float, reactions: list[reaction], seed: int
) -> list[Tuple[float, int]]:

    M = len(reactions)

    # Initialize. Set the initial number of molecules of each species and set t=0.
    np.random.seed(seed)
    t = 0.0
    results: list[tuple[float, int]] = []
    state = start_state.copy()

    propensities = np.array([0.0] * M)
    propensities_bar = np.array([0.0] * M)
    tau_k = np.array([0.0] * M)

    pbar = tqdm(total=100)
    chunks = np.linspace(0, end_time, 21)
    chunk_idx = 0

    # Calculate the propensity function, ak, for each reaction.
    for idx, r in enumerate(reactions):
        propensities[idx] = r.prop_func(state)

    # Generate M independent, uniform(0,1) random numbers r_k.
    r_k = np.random.rand(M)

    # For each k, set tau_k= 1/a_k * ln(1/r_k)
    for idx in range(M):
        tau_k[idx] = (1 / propensities[idx]) * log(1 / r_k[idx])

    while t < end_time:

        if t > chunks[chunk_idx]:
            pbar.update(5)
            chunk_idx += 1
        # Set t=min_k(tau_k) and let tau_mu be the time where the minimum is realized.
        T, j = min((val, idx) for (idx, val) in enumerate(tau_k))
        t = T

        # Update the number of each molecular species according to reaction mu
        state += reactions[j].update
        results.append((t, j))

        # Recalculate the propensity functions for each reaction and denote by a_kbar
        for idx, r in enumerate(reactions):
            propensities_bar[idx] = r.prop_func(state)

        # For each k != mu, set tau_k=(a_k /a_kbar)(tau_k−t)+t.
        for idx in range(M):
            if idx != j and tau_k[idx] != np.inf:
                tau_k[idx] = (propensities[idx] / propensities_bar[idx]) * (
                    tau_k[idx] - t
                ) + t
            else:
                # For reaction mu, let r be uniform(0,1) and set tau_k[mu] = 1 / a_kbar[mu] * ln(1 / r) + t
                r = np.random.rand()
                tau_k[idx] = (1 / propensities_bar[idx]) * log(1 / r) + t

        # For each k, set a_k = a_kbar
        propensities = np.copy(propensities_bar)

    return results
