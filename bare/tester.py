from direct import *
from reconstruct import *

def main():
    f: Callable[[list[int]], float] = lambda x: x[0] * 1.0
    a = reaction(update = [-1], propensity = f)
    simulate([100], 10.0, [a], 1000, 8)
    # results = simulate([100], 10.0, [a], 2, 2)
    # print(results[0])
    # print('------------')
    # print(results[1])
    return

if __name__ == "__main__":
    main()