import pandas as pd
import numpy as np
from reaction import reaction

def reconstruct(results: list[tuple[float, int]], reactions: list[reaction], initial_state: list[int]) -> pd.DataFrame:
    state = np.array(initial_state)
    total = [{"time": 0.0, "state": state.copy()}]

    for time, j in results:
        state += reactions[j].update
        total.append({"time": time, "state": state.copy()})

    # df = pd.DataFrame(total, columns = ["Time"] + [str(x) for x in list(range(len(reactions)))])
    df = pd.DataFrame(total)
    return df
