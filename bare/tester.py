from direct import *
from reconstruct import *

def main():
    f: Callable[[list[int]], float] = lambda x: x[0] * 1.0
    a = reaction(update = [-1], propensity = f)
    for x in range(10000):
        simulate([100], 10.0, [a])
    # results = simulate([100], 10.0, [a])
    # print(results)
    return

if __name__ == "__main__":
    main()