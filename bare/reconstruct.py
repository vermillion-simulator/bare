import pandas as pd
import numpy as np
from reaction import reaction

def reconstruct(results: list[tuple[float, int]], reactions: list[reaction], initial_state: list[int], mapping: dict[str, int]) -> pd.DataFrame:
    state = np.array(initial_state)
    _ = {k: state[v] for k,v in mapping.items()}
    _.update({"time": 0.0})
    total = [_.copy()]

    for time, j in results:
        state += reactions[j].update
        _ = {k: state[v] for k,v in mapping.items()} 
        _.update({"time": time})
        total.append(_.copy())

    # df = pd.DataFrame(total, columns = ["Time"] + [str(x) for x in list(range(len(reactions)))])
    df = pd.DataFrame(total)
    return df
