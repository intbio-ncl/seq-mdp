
from diversityStats.lib.uniprot_ec_dict import uniprot_ec_dict
from diversityStats.lib.max_min_diversity import compute_diverse_subset
from collections import defaultdict
from diversityStats.lib.gini_simpson import gini_simpson_dict, gini_simpson_value
from simpleTabuSearch.simple_tabu_imp import compute_MDP_tabu

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
        ac_subset, binary, solutions = compute_diverse_subset(_DISTPATH, _HEADPATH, _K)

        maxmin_subset = get_ec_subset(ac_subset, ac_to_ec)

        tabu_subset = compute_MDP_tabu('../Datasets/SSNMatrices/trans1074_identities.npy', '../Datasets/SSNMatrices/trans1074_headings.json', _K, 2)

        print(len(ec_to_ac.keys()))
        print(len(maxmin_subset.keys()))

        print(sorted(maxmin_subset))
        print(ac_subset)

        gs_dict = gini_simpson_dict(ac_subset, ac_to_ec)
        gs_val = gini_simpson_value(gs_dict)

        print(gs_val)

        best_tabu_subset = []
        best_tabu_index = 0
        best_score = 0

        for i in range(len(tabu_subset)):

            ecs = get_ec_subset(tabu_subset[i], ac_to_ec)
            ec_num = len(ecs.keys())

            if ec_num > best_score:
                best_score = ec_num
                best_tabu_index = i

        best_tabu = tabu_subset[best_tabu_index]
        print(get_ec_subset(best_tabu, ac_to_ec))
        print(best_score)
        print(best_tabu)
        gs_dict = gini_simpson_dict(best_tabu, ac_to_ec)
        gs_val = gini_simpson_value(gs_dict)
        print(gs_val)



    else:
        node_num  = len(ac_to_ec)
        results = []
        k_arr = []
        gs_results = []

        flag = True

        threshold = int(0.35 * node_num)
        sol, binary, solutions = compute_diverse_subset(_DISTPATH, _HEADPATH, threshold)

# for i in range(1, int(0.35 * node_num)):

        counter = len(solutions[0])
        for ac_subset in solutions:
            # ac_subset, binary = compute_diverse_subset(_DISTPATH, _HEADPATH, i)
            # print(ac_subset)
            ec_subset = get_ec_subset(ac_subset, ac_to_ec)
            subset_ec_count = len(ec_subset.keys())

            print(f'K={counter}\nTotal EC={total_ec_count}\nSubset EC={subset_ec_count}\n')
            results += [subset_ec_count/total_ec_count]
            k_arr += [counter]

            gs_dict = gini_simpson_dict(ac_subset, ac_to_ec)
            gs_val = gini_simpson_value(gs_dict)
            gs_results += [gs_val]

            if flag and subset_ec_count == total_ec_count:
                print(f'Full Coverage at K={counter}!')
                plt.axvline(counter, c="red", alpha=0.25, linestyle="--")
                flag = False
            counter += 1

        print(results)
        plt.plot(k_arr, results, label='EC Coverage')
        plt.plot(k_arr, gs_results, label='Gini-Simpson Index')
        plt.legend(loc='upper left')
        plt.title('Subset Diversity vs Subset Size')
        plt.xlabel('Subset Size')
        plt.ylabel('Subset Diversity')
        plt.show()


if __name__ == '__main__':
    main()
