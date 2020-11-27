
import json
from collections import defaultdict
import matplotlib.pyplot as plt

from download_dataset import SPARQL_query_families, SPARQL_query_ec_count


def gini_simpson_dict(key_subset, sig_dict):

    p_dict = {}

    for key in key_subset:

        # print key_subset
        sig = str(sig_dict[key])

        if sig in p_dict.keys():
            p_dict[sig] += 1.0
        else:
            p_dict[sig] = 1.0

    max = len(key_subset)

    for key, value in p_dict.items():
        p_dict[key] = (value/max)**2
        # print key, value
    return p_dict


def gini_simpson_value(gini_simps_dict):

    val = 1 - sum(gini_simps_dict.values())

    return val


def read_json_file(path):

    final_data = []

    with open(path) as f:
        data = json.load(f)

    for fam in data:
        final_data.append(fam['metadata']['accession'])

    return final_data


def process_ec_results(res):

    ec_set = set()
    ac_to_ec = defaultdict(str)
    full_ec_count = 0

    for seq in res:

        ac = seq['protein']['value'].split("/")[-1]
        ec = seq['enz']['value'].split("/")[-1]
        ec_set.add(ec)
        ac_to_ec[ac] = ec

        if '-' not in ec:
            full_ec_count += 1

    full_ec_perc = full_ec_count / len(res)

    gs_dict = gini_simpson_dict(ac_to_ec.keys(), ac_to_ec)
    gs_val = gini_simpson_value(gs_dict)

    return len(ec_set), full_ec_perc, gs_val


def main():

    fams = SPARQL_query_families()
    top100fams = [fam['fam']['value'] for fam in fams][0:200]
    seq_count_dict = {fam['fam']['value'].split('/')[-1]:fam['famcount']['value'] for fam in fams}

    ecs_dict = defaultdict(int)

    for fam in top100fams:
        fam = fam.split('/')[-1]
        ec_results = SPARQL_query_ec_count(fam)
        ec_set_len, full_ec_perc, gs_val = process_ec_results(ec_results)

        if gs_val < 0.6 or full_ec_perc < 0.9:
            continue

        if len(ec_results) < 900:
            break

        ecs_dict[fam] += ec_set_len
        print(fam, len(ec_results), ec_set_len, full_ec_perc, gs_val)

    top15fams = sorted(ecs_dict.items(), key=lambda item: item[1], reverse=True)[0:15]
    ecs_dict = {tup[0]:tup[1] for tup in top15fams}
    plt.bar(ecs_dict.keys(), ecs_dict.values())
    plt.title("Number of Unique EC Classes for the Top 15 PFAM Families")
    plt.ylabel("Unique EC Classes")
    plt.show()

    print(ecs_dict)


if __name__ == '__main__':
    main()
