
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

_THR = 1e-15


def extract_scores(path, thr):

    fr = open(path, "r")

    e_vals = []
    bit_scores = []
    ids = defaultdict(float)

    for line in fr:
        line = line.strip()

        if line[0] == '#':
            continue

        temp_lst = line.split(' ')
        new_lst = []

        for i in temp_lst:
            if i != '':
                new_lst.append(i)
            elif new_lst[-1] == " ":
                new_lst[-1] = "\t"
            elif new_lst[-1] == "\t":
                continue

        if float(new_lst[4]) <= thr:
            e_vals.append(float(new_lst[4]))
            bit_scores.append(float(new_lst[5]))
            ids[new_lst[0].split("|")[1]] = new_lst[4]

    return ids, e_vals, bit_scores


def print_metrics(e_vals, bit_scores, label, flag=False):

    print(f"Average e-val: {np.mean(e_vals)}±{np.std(e_vals)}")
    print(f"Average bit score: {np.mean(bit_scores)}±{np.std(bit_scores)}")

    parsed_evals = []

    for x in e_vals:
        x = str(x)
        if 'e-' in x:
            parsed_evals.append(float(str(x).split('e-')[1]))
        else:
            parsed_evals.append(0.0)

    if isinstance(flag, bool):
        fig, ax = plt.subplots(2)
        ax[0].boxplot(parsed_evals, positions=[1], patch_artist=True)
        ax[1].boxplot(bit_scores, positions=[1], patch_artist=True)

    else:
        ax = flag
        ax[0].boxplot(parsed_evals, positions=[2], patch_artist=True)
        ax[1].boxplot(bit_scores, positions=[2], patch_artist=True)

    return ax


def check_overlap(ids_pfam, ids_mdp):

    overlap = [(x, ids_pfam[x]) for x in ids_pfam if x in ids_mdp]
    unique_pfam = [(x, ids_pfam[x]) for x in ids_pfam if x not in [y[0] for y in overlap]]
    unique_mdp = [(x, ids_mdp[x]) for x in ids_mdp if x not in [y[0] for y in overlap]]

    print(f"PFAM Count: {len(ids_pfam)}")
    print(f"MDP Count: {len(ids_mdp)}\n")

    print(f"Search Overlap Count: {len(overlap)}")
    print(f"Unique PFAM Count: {len(unique_pfam)}")
    print(f"Unique MDP Count: {len(unique_mdp)}")

    all_matches = [] + unique_mdp + unique_pfam + overlap

    save_ids(unique_pfam, "./hmmsearch_results/hmmsearch_pfam_sub_e15.tab")
    save_ids(unique_mdp, "./hmmsearch_results/hmmsearch_mdp_sub_e15.tab")
    save_ids(overlap, "./hmmsearch_results/hmmsearch_overlap_sub_e15.tab")


def save_ids(tups, out):

    fr = open(out, "w")

    for tup in tups:
        fr.write(f"{tup[0]}\t{tup[1]}\n")


def main():

    print(f"Threshold picked: {_THR}\n")

    ids_pfam, e_vals_pfam, bit_scores_pfam = extract_scores("./hmmsearch_results/hmmsearch_res_pfam.tab", _THR)
    ids_custom, e_vals_custom, bit_scores_custom = extract_scores("./hmmsearch_results/hmmsearch_res_mdp.tab", _THR)

    print("PFAM HMM:")
    ax = print_metrics(e_vals_pfam, bit_scores_pfam, 'PFAM')

    print("")

    print("CUSTOM HMM:")
    print_metrics(e_vals_custom, bit_scores_custom, 'Custom', ax)

    ax[0].set_title("e-val Distribution")
    ax[1].set_title("Bit score Distribution")
    ax[0].set_xticklabels(['PFAM', 'MDP'])
    ax[1].set_xticklabels(['PFAM', 'MDP'])

    # plt.show()

    check_overlap(ids_pfam, ids_custom)

    # save_ids(ids_pfam, "pfam_hmmsearch_results.tab")
    # save_ids(ids_custom, "mdp_hmmsearch_results.tab")

if __name__ == '__main__':
    main()
