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

    pbar = tqdm(total=100)
    chunks = np.linspace(0, end_time, 21)
    chunk_idx = 0

    # Initialize. Set the initial number of molecules of each species. Set t= 0.
    np.random.seed(seed)
    t = 0.0
    results: list[tuple[float, int]] = []
    state = start_state.copy()

    # For each k, set P_k = 0 and T_k = 0.
    P_k = np.array([0.0] * M)
    T_k = np.array([0.0] * M)
    deltaT = np.array([0.0] * M)
    propensities = np.array([0.0] * M)

    # Calculate the propensity function, a_k, for each reaction.
    for idx, r in enumerate(reactions):
        propensities[idx] = r.prop_func(state)

    # Generate M independent, uniform(0,1) random numbers r_k.
    r_k = np.random.rand(M)

    # For each k, set P_k=ln(1/r_k)
    for idx, val in enumerate(r_k):
        P_k[idx] = log(1 / val)

    while t < end_time:

        if t > chunks[chunk_idx]:
            pbar.update(5)
            chunk_idx += 1

        # For each k, set delta t_k =(P_k−T_k/a_k).
        for idx in range(M):
            deltaT[idx] = (
                (P_k[idx] - T_k[idx]) / propensities[idx]
                if propensities[idx] > 0
                else np.inf
            )

        # Set delta=min_k(delta t_k) and let delta t_mu be the time where the minimum is realized.
        delta, j = min((val, idx) for (idx, val) in enumerate(deltaT))

        # Set t=t+delta and update the number of each molecular species according to reaction mu.
        t += delta
        state += reactions[j].update
        results.append((t, j))

        # For each k, set T_k=T_k+a_k * delta
        for idx, val in enumerate(T_k):
            T_k[idx] = val + propensities[idx] * delta

        # For reaction mu, let r be uniform(0,1) and set P=P+ln(1/r).
        P_k[j] = P_k[j] + log(1 / np.random.rand())

        # Recalculate the propensity functions, a_k.
        for idx, r in enumerate(reactions):
            propensities[idx] = r.prop_func(state)

    return results
