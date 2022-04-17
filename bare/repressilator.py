from reaction import *
from direct import *
from reconstruct import *
import matplotlib.pyplot as plt

def main():
    # kmu = 0.5, kmo = 0.0005, kp = 0.167, gamma_m = 0.005776, gamma_p = 0.001155, kr = 1, kul = 224, ku2 = 9
    reactions = [
        "-> m1, 0.5",
        "-> m2, 0.5",
        "-> m3, 0.5",
        "m1 -> p1, 0.167", # Translation of gene1 protein
        "m2 -> p2, 0.167", # Translation of gene 2 protein
        "m3 -> p3, 0.167", # Translation of gene 3 protein
        "m1 ->, 0.005776", # Degradation of gene 1 mRNA
        "m2 ->, 0.005776", # Degradation of gene 2 mRNA
        "m3 ->, 0.005776", # Degradation of gene 3 mRNA
        "p1 ->, 0.001155", # Degradation of unbound gene 1 protein
        "p2 ->, 0.001155", # Degradation of unbound gene 2 protein
        "p3 ->, 0.001155", # Degradation of unbound gene 3 protein
        "n1 ->, 0.001155",  # degradation of bound gene 1 protein
        "n2 ->, 0.001155", # degradation of bound gene 2 protein
        "n3 ->, 0.001155", # degradation of bound gene 3 protein
        "p3 -> n3, 1.0", # binding of protein to gene 1 operator
        "p1 -> n1, 1.0",
        "p2 -> n2, 1.0",
        "n3 -> p3, 224.0",
        "n1 -> p1, 224.0",
        "n2 -> p2, 224.0"
    ]
    quantity = {"m1": 10, "m2": 10, "m3": 10, "p1": 10, "p2": 10, "p3": 10, "n1": 0, "n2": 0, "n3": 0}
    reactions, initial_state, mapping = create_model(reactions, quantity)
    reactions[0].prop_func = lambda x: 0.5 if x[8] == 0 else 0.0005
    reactions[1].prop_func = lambda x: 0.5 if x[6] == 0 else 0.0005
    reactions[2].prop_func = lambda x: 0.5 if x[7] == 0 else 0.0005
    reactions[3].update = np.array([0,0,0,1,0,0,0,0,0])
    reactions[4].update = np.array([0,0,0,0,1,0,0,0,0])
    reactions[5].update = np.array([0,0,0,0,0,1,0,0,0])
    reactions[15].prop_func = lambda x: 1 * x[5] * (x[8] < 2)
    reactions[16].prop_func = lambda x: 1 * x[3] * (x[6] < 2)
    reactions[17].prop_func = lambda x: 1 * x[4] * (x[7] < 2)
    reactions[18].prop_func = lambda x: 224.0 * (x[8] == 1) + 2 * 9 * (x[8] == 2) 
    reactions[19].prop_func = lambda x: 224.0 * (x[6] == 1) + 2 * 9 * (x[6] == 2) 
    reactions[20].prop_func = lambda x: 224.0 * (x[7] == 1) + 2 * 9 * (x[7] == 2)
    results = simulate_run(initial_state, 60000.0, reactions, 5)
    results = reconstruct(results, reactions, initial_state, mapping, ["p1", "p2", "p3"])
    # results.to_csv("repressilator_results_1000.csv")
    plt.figure()
    results.plot(x="time", y=["p1", "p2", "p3"])
    plt.savefig("test.png")
    return

if __name__ == "__main__":
    main()