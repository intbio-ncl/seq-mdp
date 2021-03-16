
from diversityStats.lib.uniprot_ec_dict import uniprot_ec_dict
from diversityStats.lib.gini_simpson import gini_simpson_dict, gini_simpson_value

import numpy as np
import json
from collections import defaultdict
import matplotlib.pyplot as plt
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
import seaborn as sns
from copy import deepcopy


_MAT = "./PF04055_mat.npy"
_HEAD = "./PF04055_headings.json"
_ANN = "./PF04055_annotation.tab"


def sep_indices(val):

    set_indices = []
    non_indices = []

    for i in range(0, len(val)):
        if val[i] == 1:
            set_indices += [i]
        else:
            non_indices += [i]

    return (set_indices, non_indices)


def score_sol(mat, sol):

    score = 0
    (set_indices, non_indices) = sep_indices(sol)

    sim_list = []

    for i in range(0, len(set_indices)):
        for j in range(i + 1, len(set_indices)):
            sim = ((1 - mat[set_indices[i], set_indices[j]])) - (mat[set_indices[i], set_indices[j]]/(1 - mat[set_indices[i], set_indices[j]]))
            score += sim

            sim_list.append(mat[set_indices[i], set_indices[j]])

    return score, sim_list


def get_ec_subset(subset, ac_to_ec):

    ec_dict = defaultdict(int)

    for ac in subset:
        ec_str = ac_to_ec[ac]
        ec_lst = ec_str.split(';')

        for ec in ec_lst:
            ec = ec.strip()
            ec_dict[ec] += 1

    return ec_dict


def initialise_matrix(dist_file):

    dist = np.load(dist_file)

    for i in range(0, len(dist)):
        dist[i][i] = np.nan
        for j in range(i+1, len(dist)):
            if dist[i][j] == 1:
                dist[i][j] = 0.99
                dist[j][i] = 0.99

    return np.asmatrix(dist)


def initialise_headings(heading_file):

    with open(heading_file) as heading:
        headings = json.loads(heading.read())

    temp_dict = {}

    for i in range(0, len(headings)):
        temp_dict[i] = headings[i]

    return temp_dict


def random_solution(length, num_picked):

    arr = np.array([0] * (length - num_picked) + [1] * num_picked)
    np.random.shuffle(arr)
    return list(arr)


def main():

    mat = initialise_matrix(_MAT)
    head = initialise_headings(_HEAD)
    ac_to_ec, ec_to_ac = uniprot_ec_dict(_ANN, 2)

    # pf00171
    # custom_subset = ['Q6F9F7', 'A5GSH0', 'P0DPF0', 'B8I6T0', 'Q2FWX9', 'P43503', 'Q8GAI8', 'P46329', 'Q55167', 'P55653', 'Q9A777', 'H2IFE7', 'O05619', 'O32507', 'B2V9F3', 'Q5E2G9', 'Q985M6', 'P86808', 'B8DCT8', 'Q6F9G0', 'Q738L2', 'Q9I6C8', 'B7MVM5', 'Q2SKP1', 'P80668', 'Q3B2U3', 'Q55585', 'B3QPW0', 'P25553', 'P28810', 'Q2BN77', 'Q92YD2', 'Q8GAK7', 'P39616', 'Q5X4K4', 'Q2YV11', 'A5YBJ3', 'Q7VBM1', 'A6TUA0', 'Q59702', 'P42269', 'P23883', 'P23105', 'C0R0B8', 'B2VJX8', 'Q6G4Z0', 'A1WGI4', 'H8ZPX2', 'Q84DC3', 'P94358', 'Q2G9T9', 'Q7VI05', 'P0C1E0', 'Q3IC91', 'P76149', 'Q8CNI5', 'H1ZV37', 'P38947', 'B0S9A5', 'P12693', 'B2IZ89', 'Q47UQ0', 'E1V7V8', 'Q1JUP4', 'Q4JWT3', 'Q6LTX2', 'A1VYR7', 'A5CXP4', 'Q0K845', 'O86447', 'Q1GV29', 'A8EVN0', 'Q7N2G9', 'A5FYS4', 'A2C148', 'Q9I702', 'Q79EM7', 'Q9AHG1', 'P33008', 'Q165Y8', 'B2S2U7', 'O69497', 'Q9KWS5', 'B1XDF5', 'Q03ZF1', 'P0A391', 'Q5WH11', 'Q59931', 'B0U208', 'A2BQ71', 'B0T8I8', 'Q8CUQ4', 'Q72NQ9', 'Q4QJW6', 'A3DC22', 'B4U8A0', 'Q92UV7', 'Q04FB2', 'Q12SX8', 'Q4L919']
    #pf00155
    # custom_subset = ['Q8DTM1', 'Q058A6', 'C6C2Z3', 'Q6LT75', 'B1N009', 'P9WQ88', 'P39643', 'Q9I468', 'B8ERL9', 'P43089', 'P72173', 'C3KVX5', 'P95468', 'Q6MDE0', 'Q253K9', 'B6YRL2', 'P9WPZ4', 'B9KDN6', 'P43336', 'P04693', 'Q8DM42', 'Q08432', 'B0TY45', 'Q5SHW0', 'Q8XZC3', 'Q68XV9', 'Q9CJU0', 'Q92G23', 'Q9KM65', 'Q9Z856', 'Q8ABA8', 'Q7UZZ3', 'Q81V80', 'B1WY56', 'Q5F6R6', 'A6L8U2', 'P00509', 'B4RFX5', 'Q8KZM9', 'P0A959', 'B7VH15', 'A0L3L7', 'Q4FP52', 'Q3S8P9', 'P28735', 'Q02636', 'C0QFJ4', 'Q8R5U4', 'P96847', 'O25320', 'P63503', 'Q5ZW88', 'P77806', 'P74770', 'B8D707', 'P97084', 'Q72LL6', 'Q9RRM7', 'A4SPR6', 'Q6D3C0', 'Q84CG1', 'B4U9L1', 'Q7VL09', 'P09053', 'A5N7Q7', 'Q93QC6', 'Q7W9I4', 'P44425', 'Q492K2', 'B2A250', 'Q3ARM7', 'A7N6R9', 'Q55128', 'O07587', 'P36692', 'Q3BYN0', 'Q9L6I2', 'Q89AX7', 'B5YFU5', 'Q06965', 'A8MEH2', 'P37419', 'Q9CBM9', 'Q5FRR4', 'Q89AK6', 'Q72PG3', 'P16524', 'P21633', 'P63499', 'A5INE2', 'O87320', 'B8CX89', 'O52815', 'Q0BBD6', 'P0A4X5', 'A0PXP5', 'Q02135', 'A5CVR5', 'P36570', 'B8J3V0']

    #pf04055
    custom_subset = ['B3E599', 'A5G2D2', 'Q057Q1', 'B8FS78', 'Q9K864', 'Q7NCE3', 'Q8CJT5', 'Q185C5', 'Q97L63', 'Q7NIT2', 'Q2SWB9', 'A9FD89', 'P71011', 'A4XGB8', 'Q82K95', 'A1W1T3', 'O25376', 'Q057G5', 'Q9EYN8', 'B5QX73', 'B7GQG0', 'P39409', 'Q4FNN5', 'Q38HX2', 'Q057Q7', 'A0A1C7D1B7', 'B0CDZ6', 'P43751', 'P09825', 'P32131', 'Q0RD46', 'Q55373', 'P10390', 'Q8K9D9', 'B2RH08', 'P74132', 'Q8KBK9', 'A0Q2E1', 'A6H1N2', 'B8ENI9', 'A7ZE07', 'O83293', 'O67826', 'B0VQD7', 'Q8KC85', 'Q8YR77', 'B2GLQ7', 'B9M4F4', 'Q2RSY6', 'Q49573', 'Q89ZC3', 'P17434', 'A9EPV3', 'P20714', 'A8Z642', 'Q17XY7', 'Q8KCU0', 'Q9K0Q5', 'O34162', 'Q8KFK8', 'B1Y6D6', 'Q3J561', 'Q8EUX4', 'P9WJ78', 'B1IL14', 'Q81G67', 'Q02550', 'P51008', 'P73667', 'A0A384LP51', 'Q55914', 'Q53U14', 'P45097', 'Q8D1Y5', 'A0LV48', 'O33506', 'Q2JRI4', 'P55477', 'Q72DS4', 'A0RIB6', 'Q0TTH1', 'Q44634', 'A8Z609', 'B1ZVM5', 'P69848', 'C9XIS7', 'Q1IHK7', 'P75794', 'A9A0B5', 'Q6MED6', 'A6LSR6', 'P24427', 'Q8RHX4', 'Q30XT6', 'A0A069AMK2', 'O87941', 'Q9S498', 'Q8DII8', 'Q1GV98', 'A9CF16']

    score_lst = []
    cov_lst = []
    gsi_lst = []
    simscore_lst = []

    for i in range(1):
        sol = random_solution(len(head), 1)
        sol_subset = [list(ac_to_ec.keys())[i] for i in range(len(sol)) if sol[i] == 1]

        score, sim_list = score_sol(mat, sol)
        score_lst.append(score)
        simscore_lst.append(np.mean(sim_list))

        ecs = get_ec_subset(sol_subset, ac_to_ec)
        cov = len(ecs.keys()) / float(len(ec_to_ac.keys()))

        cov_lst.append(cov)

        gs_dict = gini_simpson_dict(sol_subset, ac_to_ec)
        gsi = gini_simpson_value(gs_dict)

        gsi_lst.append(gsi)

    print(f"Score:{np.mean(score_lst)}±{np.std(score_lst)}\n\
            Cov: {np.mean(cov_lst)}±{np.std(cov_lst)}\n\
            GSI: {np.mean(gsi_lst)}±{np.std(gsi_lst)}\n\
            SeqSim: {np.mean(simscore_lst)}±{np.std(simscore_lst)}")

    gs_dict = gini_simpson_dict(list(ac_to_ec.keys()), ac_to_ec)
    gsi = gini_simpson_value(gs_dict)
    full_ec_lst = set(sum([ec_to_ac[ec] for ec in ec_to_ac.keys() if '.-' not in ec],[]))

    print(f"Num Seq:{len(ac_to_ec.keys())}\nGSI: {gsi}\nNum Ann:{len(ec_to_ac.keys())}")
    print(gs_dict)

    gsi_dict1 = gs_dict
    gsi_dict2 = gini_simpson_dict(custom_subset, ac_to_ec)
    gsi2 = gini_simpson_value(gsi_dict2)

    print(gsi_dict2)
    new_dict = defaultdict(float)
    diff_lst = []
    prop_lst = []
    gsi_vals1 = []
    gsi_vals2 = []

    for sig in gsi_dict1.keys():
        if sig not in gsi_dict2.keys():
            gsi_dict2[sig] = 0
        gsi_vals1.append(gsi_dict1[sig])
        gsi_vals2.append(gsi_dict2[sig])

    for sig in gsi_dict1.keys():
        diff = gsi_dict2[sig] - gsi_dict1[sig]
        diff_lst.append(diff)

        new_dict[sig] = diff
        sig_prop = len(ec_to_ac[sig]) / len(ac_to_ec.keys())
        prop_lst.append(sig_prop)
        print(f"{sig}\t{diff}\t{sig_prop}")


    # plt.hist(diff_lst, bins=np.linspace(-0.25, 0.2, num=100))

    # plt.ylim(-0.1, 0.1)
    # plt.xlim(0, 1)
    # sns.regplot(prop_lst, diff_lst, truncate=False)


    temp_prop_lst = deepcopy(prop_lst)
    print(prop_lst)
    print(gsi_vals1, gsi_vals2)
    prop_lst, gsi_vals1 = (list(t) for t in zip(*sorted(zip(prop_lst, gsi_vals1))))
    # print(prop_lst, gsi_vals1)
    prop_lst, gsi_vals2 = (list(t) for t in zip(*sorted(zip(temp_prop_lst, gsi_vals2))))
    # print(prop_lst, gsi_vals2)

    print(prop_lst)
    print(gsi_vals1, gsi_vals2)

    #
    # prop_lst, gsi_vals1, gsi_vals2 = map(list, zip(*sorted(zip(prop_lst, gsi_vals1, gsi_vals2), reverse=True)))
    #

    fig, ax = plt.subplots()
    plt.title(r"Change in relative abundance for $\bf{SAM}$ dataset", fontsize=16)

    ax.scatter(prop_lst, gsi_vals1, label="Superset")
    ax.scatter(prop_lst, gsi_vals2, label="Subset")
    plt.xlabel("$p_i$", fontsize=14)
    plt.ylabel("$p^2_i$", fontsize=14)
    plt.text(0.025, 0.95, "$\overline{p^2_i}$, $i \in (U\cup Z)$=" + f"{round(np.mean(gsi_vals1), 4)}", ha='left', va='top', transform=ax.transAxes, fontsize=13)
    plt.text(0.025, 0.85, "$\overline{p^2_i}$, $i \in U$=" + f"{round(np.mean(gsi_vals2), 4)}", ha='left', va='top', transform=ax.transAxes, fontsize=13)

    plt.legend(fontsize=12)
    # print(prop_lst, gsi_vals1)
    plt.fill_between(prop_lst, gsi_vals1, gsi_vals2, step='mid', color='grey', alpha='0.3')


    print(len(gsi_vals1), len(gsi_vals2))
    plt.show()



if __name__ == '__main__':
    main()
