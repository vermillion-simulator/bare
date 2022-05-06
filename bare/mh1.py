from reaction import *
import direct
from reconstruct import *
import matplotlib.pyplot as plt


def main():
    # reactions = [
    #     "mh10 -> mh11, 48.3",
    #     "mh11 -> mh12, 1.0",
    #     "mh12 -> mh1, 1.0",
    #     "mh1 ->, 0.3229",
    #     "mh1 -> mh1 + ph10, 49.9",
    #     "ph10 -> ph11, 3.0",
    #     "ph11 -> ph12, 3.0",
    #     "ph12 -> ph1, 3.0",
    #     "ph1 ->, 0.3495",
    # ]
    # quantity = {
    #     "mh10": 100,
    #     "mh11": 0,
    #     "mh12": 0,
    #     "mh1": 0,
    #     "ph11": 0,
    #     "ph12": 0,
    #     "ph1": 0,
    # }
    reactions = [
        "-> mh10, 400.3",
        "mh10 -> mh1, 200.0",
        "mh1 ->, 0.15",
        "mh1 -> mh1 + ph10, 1.5",
        "ph10 -> ph11, 0.5",
        "ph11 -> ph1, 0.5",
        "ph1 ->, 10.0",
    ]
    quantity = {"mh1": 50, "mh10": 0, "ph10": 0, "ph11": 0, "ph1": 0}
    reactions, initial_state, mapping = create_model(reactions, quantity)
    reactions[0].prop_func = lambda x: 1 / (1 + (x[mapping.get("ph1")]) ** 2)
    results = direct.simulate_run(initial_state, 30000.0, reactions, 5)
    results = reconstruct(results, reactions, initial_state, mapping, mapping.keys())
    plt.figure()
    results.plot(x="time", y=["mh1"])
    plt.savefig("mh1.png")
    return


if __name__ == "__main__":
    main()
