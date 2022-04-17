from venv import create
from reaction import *
from direct import *
from reconstruct import *
import matplotlib.pyplot as plt

def main():
    reactions = ["a -> b, 1.0", "b -> c, 0.5"]
    quantity = {"a": 100, "b": 0, "c": 0}
    reactions, initial_state, mapping = create_model(reactions, quantity)
    results = simulate_run(initial_state, 15.0, reactions, 0)
    results = reconstruct(results, reactions, initial_state, mapping)
    plt.figure()
    results.plot(x="time")
    plt.savefig("test.png")
    return

if __name__ == "__main__":
    main()