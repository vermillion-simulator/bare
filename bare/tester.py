from venv import create
from reaction import *
from direct import *
from reconstruct import *

def main():
    reactions = ["a -> b, 1.0"]
    quantity = {"a": 100, "b": 0}
    reactions, initial_state = create_model(reactions, quantity)
    results = simulate_run(initial_state, 10.0, reactions, 0)
    results = reconstruct(results, reactions, initial_state)
    print(results)
    return

if __name__ == "__main__":
    main()