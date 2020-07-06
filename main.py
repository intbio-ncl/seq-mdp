
from MEnzDP\-SimpleTabuSearch.lib.uniprot_ec_dict import uniprot_ec_dict
from MEnzDP-SimpleTabuSearch.lib.max_min_diversity import compute_diverse_subset
from collections import defaultdict
from MEnzDP-SimpleTabuSearch.lib.gini_simpson import gini_simpson_dict, gini_simpson_value
from tabuSearch.simple_tabu_imp import compute_MDP_tabu

import matplotlib.pyplot as plt


#########################################

import argparse

parser = argparse.ArgumentParser(description="MDP Filler stuff") #TODO: Change desc
parser.add_argument("-a", "--annotation", help="Path to the annotation", type=str, required=True)
parser.add_argument("-hd", "--heading", help="Path to heading file", type=str, required=True)
parser.add_argument("-d", "--distance", help="Path to distance file", type=str, required=True)
parser.add_argument("-k", "--subset", help="Subset size", type=int, required=True)


args = parser.parse_args()
_ANNPATH = args.annotation
_HEADPATH = args.heading
_DISTPATH = args.distance
_K = args.subset


#########################################

def get_ec_subset(subset, ac_to_ec):

    ec_dict = defaultdict(int)

    for ac in subset:
        ec_str = ac_to_ec[ac]
        ec_lst = ec_str.split(';')

        for ec in ec_lst:
            ec = ec.strip()
            ec_dict[ec] += 1

    return ec_dict


def main():

    ac_to_ec, ec_to_ac = uniprot_ec_dict(_ANNPATH, 2)
    ec_num = len(set(ac_to_ec.values()))
    total_ec_count = len(ec_to_ac.keys())


    if _K != 0:
        ac_subset, binary = compute_diverse_subset(_DISTPATH, _HEADPATH, _K)
        ec_subset = get_ec_subset(ac_subset, ac_to_ec)

        print(ec_to_ac.keys())
        print(len(ec_to_ac.keys()))
        print(len(ec_subset.keys()))

        print(sorted(ec_subset))
        print(ac_subset)
        print(binary)


    else:
        node_num  = len(ac_to_ec)
        results = []
        k_arr = []
        gs_results = []

        for i in range(1, int(0.2*node_num)):
            ac_subset, binary = compute_diverse_subset(_DISTPATH, _HEADPATH, i)
            ec_subset = get_ec_subset(ac_subset, ac_to_ec)
            subset_ec_count = len(ec_subset.keys())

            print(f'K={i}\nTotal EC={total_ec_count}\nSubset EC={subset_ec_count}\n')
            results += [subset_ec_count/total_ec_count]
            k_arr += [i]

            if subset_ec_count == total_ec_count:
                print(f'Full Coverage at K={i}!')

            gs_dict = gini_simpson_dict(ac_subset, ac_to_ec)
            gs_val = gini_simpson_value(gs_dict)
            gs_results += [gs_val]

        print(results)
        plt.plot(k_arr, results, label='EC Coverage')
        plt.plot(k_arr, gs_results, label='Gini-Simpson Index')
        plt.legend(loc="upper left")
        plt.show()

    compute_MDP_tabu('../Datasets/SSNMatrices/trans1074_identities.npy', '../Datasets/SSNMatrices/trans1074_headings.json', 100)


if __name__ == '__main__':
    main()
