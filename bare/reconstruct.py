import pandas as pd
import numpy as np
from reaction import reaction


def identify(reactions, mapping, variables) -> list[int]:
    interest = []
    for idx, r in enumerate(reactions):
        changes = [idx for idx, val in enumerate(r.update) if val != 0]
        if any(mapping[x] in changes for x in variables):
            interest.append(idx)
    return interest


def reconstruct(
    results: list[tuple[float, int]],
    reactions: list[reaction],
    initial_state: list[int],
    mapping: dict[str, int],
    variables: list[str],
) -> pd.DataFrame:
    state = np.array(initial_state)
    interest = identify(reactions, mapping, variables)

    _ = {k: state[v] for k, v in mapping.items()}
    _.update({"time": 0.0})
    total = [_.copy()]

    for time, j in results:
        if j in interest:
            state += reactions[j].update
            _ = {k: state[v] for k, v in mapping.items()}
            _.update({"time": time})
            total.append(_.copy())

    # df = pd.DataFrame(total, columns = ["Time"] + [str(x) for x in list(range(len(reactions)))])
    df = pd.DataFrame(total)
    return df
