#!/usr/bin/env python

from diversityStats.lib.uniprot_ec_dict import uniprot_ec_dict
from diversityStats.lib.max_min_diversity import compute_diverse_subset
from collections import defaultdict
from diversityStats.lib.gini_simpson import gini_simpson_dict, gini_simpson_value
from simpleTabuSearch.simple_tabu_imp import compute_MDP_tabu

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from copy import deepcopy
#########################################

import argparse

parser = argparse.ArgumentParser(description="MDP Filler stuff") #TODO: Change desc
parser.add_argument("-a", "--annotation", help="Path to the annotation", type=str, required=True)
parser.add_argument("-hd", "--heading", help="Path to heading file", type=str, required=True)
parser.add_argument("-d", "--distance", help="Path to distance file", type=str, required=True)
parser.add_argument("-k", "--subset", help="Subset size", type=int, required=True)
parser.add_argument("-g", "--graph", help="Path to Similarity Network (.gml) to output graph-based results", type=str, required=True)

args = parser.parse_args()
_ANNPATH = args.annotation
_HEADPATH = args.heading
_DISTPATH = args.distance
_K = args.subset
_GRAPHPATH = args.graph


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


def sep_indices(val):

    set_indices = []
    non_indices = []

    for i in range(0, len(val)):
        if val[i] == 1:
            set_indices += [i]
        else:
            non_indices += [i]

    return (set_indices, non_indices)


def score(mat, sol, use_clan=False):

    score = 0
    (set_indices, non_indices) = sep_indices(sol)

    clan_threshold = 0.60
    clan_penalty = 0

    sim_list = []

    for i in range(0, len(set_indices)):
        for j in range(i + 1, len(set_indices)):
            sim = ((1 - mat[set_indices[i], set_indices[j]])) - (mat[set_indices[i], set_indices[j]]/(1 - mat[set_indices[i], set_indices[j]]))
            score += sim

            sim_list.append(mat[set_indices[i], set_indices[j]])

            if mat[set_indices[i], set_indices[j]] >= clan_threshold:
                clan_penalty += 1

    # if use_clan:
    #     clan_penalty = 1 - (clan_penalty / len(set_indices))
    #     score *= clan_penalty

    return score, sim_list


def main():

    ac_to_ec, ec_to_ac = uniprot_ec_dict(_ANNPATH, 2)
    ec_num = len(set(ac_to_ec.values()))
    total_ec_count = len(ec_to_ac.keys())

    if _K != 0:
        ac_subset, binary, solutions, mat = compute_diverse_subset(_DISTPATH, _HEADPATH, _K)

        G = nx.read_gml(_GRAPHPATH)
        G.remove_edges_from(nx.selfloop_edges(G))

        maxmin_distscore, sim_list = score(mat, binary, True)

        fig, axs = plt.subplots(2)

        axs.flat[0].hist(sim_list, bins=np.linspace(0, 1, num=20), histtype='step')
        axs.flat[0].set_ylim(0, 5000)
        axs.flat[0].set_title("Max-Min MDP Solver", fontsize=10, fontweight='bold')

        mean = round(np.mean(sim_list), 3)
        std = round(np.std(sim_list), 3)

        axs.flat[0].text(0, 4000, f"{mean}±{std}")

        print(maxmin_distscore)
        maxmin_subset = get_ec_subset(ac_subset, ac_to_ec)
        maxmin_edgecount = len(G.subgraph(ac_subset).edges)
        maxmin_ECcount = len(maxmin_subset.keys()) / float(len(ec_to_ac.keys()))

        tabu_subset, score_list, best_progress, sim_list2 = compute_MDP_tabu(_DISTPATH, _HEADPATH, _K, 0, binary)
        axs.flat[1].hist(sim_list2, bins=np.linspace(0, 1, num=20), histtype='step')
        axs.flat[1].set_ylim(0, 5000)
        axs.flat[1].set_title("Tabu Search MDP Solver", fontsize=10, fontweight='bold')

        mean = round(np.mean(sim_list2), 3)
        std = round(np.std(sim_list2), 3)

        axs.flat[1].text(0, 4000, f"{mean}±{std}")

        fig.suptitle("Sequence Similarity Distributions for PF00171 ")

        print(len(ec_to_ac.keys()))
        print(len(maxmin_subset.keys()))

        greedy_coverage = len(maxmin_subset.keys()) / len(ec_to_ac.keys())

        print(f"Greedy Coverage: {greedy_coverage}")

        print(sorted(maxmin_subset))
        print(ac_subset)
        print(len(G.subgraph(ac_subset).edges))

        gs_dict = gini_simpson_dict(ac_subset, ac_to_ec)
        gs_val = gini_simpson_value(gs_dict)

        print(f"Greedy Gini: {gs_val}")

        best_tabu_subset = []
        best_tabu_index = 0
        best_score = 0

        ecnum_list = []
        edgenum_list = []

        for i in range(len(tabu_subset)):

            ecs = get_ec_subset(tabu_subset[i], ac_to_ec)
            ec_num = len(ecs.keys()) / float(total_ec_count)
            ecnum_list += [ec_num]
            print(ec_num, len(G.subgraph(tabu_subset[i]).edges), score_list[i])
            edgenum_list += [len(G.subgraph(tabu_subset[i]).edges)]

            if ec_num > best_score:
                best_score = ec_num
                best_tabu_index = i

        best_tabu = tabu_subset[best_tabu_index]
        print(get_ec_subset(best_tabu, ac_to_ec).keys())
        print(best_score)
        print(best_tabu)
        gs_dict = gini_simpson_dict(best_tabu, ac_to_ec)
        gs_val = gini_simpson_value(gs_dict)
        print(f"Tabu Gini: {gs_val}")

        fig1, ax1 = plt.subplots()

        ax1.plot(best_progress)
        ax1.set_xlabel("Epochs")
        ax1.set_ylabel("Distance Score")
        ax1.set_title("Distance Score Over Time for K=100 ")
        ax1.axhline(maxmin_distscore, c="red", alpha=0.25, linestyle="--")

        fig, axs = plt.subplots(2)
        fig.suptitle('EC Coverage vs Distance Score/Number of Edges for K=100')

        temp_ecnum_list = deepcopy(ecnum_list)
        print(edgenum_list, ecnum_list)
        edgenum_list, ecnum_list = (list(t) for t in zip(*sorted(zip(edgenum_list, ecnum_list))))
        print(edgenum_list, ecnum_list)
        print(score_list, temp_ecnum_list)
        score_list, temp_ecnum_list = (list(t) for t in zip(*sorted(zip(score_list, temp_ecnum_list))))
        print(score_list, temp_ecnum_list)

        axs[0].scatter(edgenum_list, ecnum_list, alpha=0.25)
        axs.flat[0].set(xlabel='Number of Edges')
        axs.flat[0].set(ylabel='EC Coverage')
        axs.flat[0].set_ylim(0, 1)
        axs[0].scatter(maxmin_edgecount, maxmin_ECcount)
        axs[0].plot(edgenum_list, ecnum_list, alpha=0.25)

        axs[1].scatter(score_list, temp_ecnum_list, alpha=0.25)
        axs.flat[1].set(xlabel='Distance Score')
        axs.flat[1].set(ylabel='EC Coverage')
        axs.flat[1].set_ylim(0, 1)
        axs[1].scatter(maxmin_distscore, maxmin_ECcount)
        axs[1].plot(score_list, temp_ecnum_list, alpha=0.25)

        plt.show()

        fw_greedy = open("./greedy_subset.txt", "w")
        fw_TS = open("./ts_subset.txt", "w")

        for entry in ac_subset:
            fw_greedy.write(entry + "\n")
        for entry in best_tabu:
            fw_TS.write(entry + "\n")

        fw_greedy.close()
        fw_TS.close()

    else:
        node_num = len(ac_to_ec)
        results = []
        k_arr = []
        gs_results = []

        flag = True

        threshold = int(0.1 * node_num)
        sol, binary, solutions, mat = compute_diverse_subset(_DISTPATH, _HEADPATH, threshold)

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
        plt.legend(loc='upper right')
        plt.title('Subset Diversity vs Subset Size for PF04055')
        plt.xlabel('Subset Size')
        plt.ylabel('Subset Diversity')
        plt.show()
        plt.save('./PF04055_EC_Diversity.pdf')


if __name__ == '__main__':
    main()
