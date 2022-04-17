from lib2to3.pytree import convert
from typing import Callable
from typing import Tuple
import numpy as np
import numba

class reaction():
    def __init__(self, update: list[int], propensity: Callable[[list[int]], float]):
        self.update = np.array(update)
        self.prop_func = propensity 

    @numba.jit
    def propensity(self, state: list[int]) -> float:
        return self.prop_func(state)

def split_reaction(reaction: str) -> Tuple[list[str], list[str], float]:
    rstr, rate = reaction.split(",")
    rate = float(rate)
    reac, prod = rstr.split("->")
    reac: list[str] = list(filter(lambda x: x!= '', [x.strip() for x in reac.split("+")]))
    prod: list[str] = list(filter(lambda x: x != '',[x.strip() for x in prod.split("+")]))
    return (reac, prod, rate)

def create_reactants_mapping(reactions: list[Tuple[list[str], list[str], float]]):
    reactants_mapping: dict[str, int] = {}
    r_idx = 0
    for r in reactions:
        reac, prod, _ = r
        for reactant in reac:
            if reactant not in reactants_mapping:
               reactants_mapping[reactant] = r_idx
               r_idx += 1
        for reactant in prod:
            if reactant not in reactants_mapping:
               reactants_mapping[reactant] = r_idx
               r_idx += 1 
    return reactants_mapping 

def convert_reaction(r: Tuple[list[str], list[str], float], mapping: dict[str, int]):
    reac, prod, rate = r
    state = [0] * len(mapping)

    if "None" not in reac:
        a = [mapping[i] for i in reac if i != ""]
        for i in a:
            state[i] -= 1
    
    if "None" not in prod:
        a = [mapping[i] for i in prod if i != ""]
        for i in a:
            state[i] += 1   

    # Zero order reaction
    if "None" in reac or "None" in prod:
        return reaction(state, lambda x: rate)

    # First order reaction
    if len(reac) == 1:
        i = mapping[reac[0]]
        return reaction(state, lambda x: x[i] * rate)

    # Second order reaction
    if len(reac) == 2:
        if len(set(reac)) == 2:
            i: int = mapping[reac[0]]
            j: int  = mapping[reac[1]]
            return reaction(state, lambda x: x[i] * x[j] * rate)
        else:
            i: int = mapping[reac[0]]
            return reaction(state, lambda x: x[i] * (x[i] - 1)/2 * rate)

def create_model(reactions: list[str], quantity: dict[str, int]):
    reactions = list(map(split_reaction, reactions))  # type: ignore
    reactant_map = create_reactants_mapping(reactions) # type: ignore
    for idx, val in enumerate(reactions):
        reactions[idx] = convert_reaction(val, reactant_map)
    state = np.array([0] * len(reactant_map))
    for item in quantity.keys():
        state[reactant_map[item]] = quantity[item]
    return reactions, state, reactant_map
