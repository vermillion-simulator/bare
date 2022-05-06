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

    pbar = tqdm(total=100)
    chunks = np.linspace(0, end_time, 21)
    chunk_idx = 0

    # Initialize. Set the initial number of molecules of each species. Set t= 0.
    np.random.seed(seed)
    t = 0.0
    results: list[tuple[float, int]] = []
    state = start_state.copy()

    # For each k <= M, set Pk= 0 and Tk= 0.
    P_k = np.array([0.0] * M)
    T_k = np.array([0.0] * M)
    deltaT = np.array([0.0] * M)
    propensities = np.array([0.0] * M)

    # For each delayed reaction channel set s_k = [inf]
    S_k = {0: [np.inf], 2: [np.inf]}

    # Generate M independent, uniform (0,1) random numbers, r_k.
    r_k = np.random.rand(M)

    # For each k set P_k=ln(1/r_k)
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

        # Set delta=min_k(delta t_k, S_k(1)-t)
        delta = j = np.inf
        for idx, val in enumerate(deltaT):
            if idx in S_k:
                temp = min(deltaT[idx], S_k[idx][0])
            else:
                temp = deltaT[idx]
            if temp < delta:
                delta = temp
                j = idx

        # Set t=t+delta
        t += delta
