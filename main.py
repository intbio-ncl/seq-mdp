#!/usr/bin/env python

from diversityStats.lib.uniprot_ec_dict import uniprot_ec_dict
from diversityStats.lib.max_min_diversity import compute_diverse_subset
from collections import defaultdict
from diversityStats.lib.gini_simpson import gini_simpson_dict, gini_simpson_value
from simpleTabuSearch.simple_tabu_imp import compute_MDP_tabu

import matplotlib.pyplot as plt
import numpy as np

#########################################

import argparse

parser = argparse.ArgumentParser(description="MDP Filler stuff") #TODO: Change desc
parser.add_argument("-a", "--annotation", help="Path to the annotation", type=str, required=False)
parser.add_argument("-hd", "--heading", help="Path to heading file", type=str, required=True)
parser.add_argument("-d", "--distance", help="Path to distance file", type=str, required=True)
parser.add_argument("-k", "--subset", help="Subset size", type=int, required=True)
parser.add_argument("-s", "--solver", choices=["greedy", "ts", "all"], help="Solver(s) to use (Choices are 'greedy', 'ts', or both with 'all').", required=True)


args = parser.parse_args()
_ANNPATH = args.annotation
_HEADPATH = args.heading
_DISTPATH = args.distance
_K = args.subset
_SOLVER = args.solver


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


def greedy_mdp(ac_to_ec=None, ec_to_ac=None):

    maxmin_subset = None
    greedy_coverage = None

    ac_subset, binary, solutions, mat = compute_diverse_subset(_DISTPATH, _HEADPATH, _K)

    maxmin_distscore, sim_list = score(mat, binary, True)

    if ac_to_ec is not None:
        maxmin_subset = get_ec_subset(ac_subset, ac_to_ec)
        greedy_coverage = len(maxmin_subset.keys()) / len(ec_to_ac.keys())

        print(f"Greedy Coverage: {greedy_coverage}")
        print(sorted(maxmin_subset))

        gs_dict = gini_simpson_dict(ac_subset, ac_to_ec)
        gs_val = gini_simpson_value(gs_dict)
        print(f"Greedy Gini: {gs_val}")

        fw_greedy = open("greedy_vals.txt", "w")
        fw_greedy.write(f"Greedy Coverage:{greedy_coverage}\nGreedy Gini-Simpson: {gs_val}")
        fw_greedy.close()

    print(ac_subset)

    fw_greedy = open("./greedy_subset.txt", "w")
    for entry in ac_subset:
        fw_greedy.write(entry + "\n")

    fw_greedy.close()

    plot_res(sim_list, 'greedy')


def ts_mdp(ac_to_ec=None, ec_to_ac=None, preset=None):

    tabu_subset, score_list, best_progress, sim_list = compute_MDP_tabu(_DISTPATH, _HEADPATH, _K, 0, preset)

    best_tabu_index = 0
    best_score = 0
    ecnum_list = []
    best_tabu = tabu_subset[0]

    for i in range(len(tabu_subset)):

        if ac_to_ec is not None:
            ecs = get_ec_subset(tabu_subset[i], ac_to_ec)
            ec_num = len(ecs.keys()) / float(len(ec_to_ac.keys()))
            ecnum_list += [ec_num]

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

            fw_ts = open("ts_vals.txt", "w")
            fw_ts.write(f"Tabu Search Coverage: {ec_num}\nTabu Search Gini-Simpson: {gs_val}")
            fw_ts.close()

    fw_TS = open("./ts_subset.txt", "w")
    for entry in best_tabu:
        fw_TS.write(entry + "\n")

    fw_TS.close()

    plot_res(sim_list, 'tabu')


def plot_res(sim_list, solver):

    fig, axs = plt.subplots()

    solver_title = ''
    file_title = ''

    if solver == 'tabu':
        solver_title = "Tabu Search"
        file_title = 'ts_distribution.pdf'

    elif solver == 'greedy':
        solver_title = 'Greedy'
        file_title = 'greedy_distribution.pdf'

    plt.hist(sim_list, bins=np.linspace(0, 1, num=20), histtype='step')
    plt.ylim(0, 5000)
    plt.title(f"Sequence Identity Distribution: K={_K}, Solver={solver_title}", fontsize=10, fontweight='bold')

    mean = round(np.mean(sim_list), 3)
    std = round(np.std(sim_list), 3)

    plt.text(0, 4000, f"{mean}±{std}")

    plt.savefig(file_title)


def main():

    ac_to_ec = None
    ec_to_ac = None
    total_ec_count = None

    if _ANNPATH is not None:
        ac_to_ec, ec_to_ac = uniprot_ec_dict(_ANNPATH, 2)
        total_ec_count = len(ec_to_ac.keys())

    if _K != 0:

        if _SOLVER in ('all', 'greedy'):
            if _ANNPATH is not None:
                greedy_mdp(ac_to_ec, ec_to_ac)
            else:
                greedy_mdp()

        if _SOLVER in ('all', 'ts'):
            if _ANNPATH is not None:
                ts_mdp(ac_to_ec, ec_to_ac)
            else:
                ts_mdp()

    else:
        node_num = len(ac_to_ec)
        results = []
        k_arr = []
        gs_results = []

        flag = True

        threshold = int(0.1 * node_num)
        sol, binary, solutions, mat = compute_diverse_subset(_DISTPATH, _HEADPATH, threshold)

        counter = len(solutions[0])
        for ac_subset in solutions:

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
