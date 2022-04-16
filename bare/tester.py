from direct import *
from reconstruct import *
import numpy as np

def main():
    f: Callable[[list[int]], float] = lambda x: x[0] * 1.0
    a = reaction(update = [-1], propensity = f)
    results = simulate([100], 10.0, [a])
    print(f"A update is {a.update} with type {type(a.update)}")
    df = reconstruct(results, [a], [100])
    print(df)

if __name__ == "__main__":
    main()