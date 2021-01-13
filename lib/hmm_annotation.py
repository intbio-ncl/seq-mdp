
from collections import defaultdict
from tqdm import tqdm
import requests
import numpy as np

_PFAMPATH = "./hmmsearch_results/hmmsearch_pfam_sub.tab"
_MDPPATH = "./hmmsearch_results/hmmsearch_mdp_sub.tab"
_OVERLAPPATH = "./hmmsearch_results/hmmsearch_overlap_sub.tab"
_URL = "https://www.ebi.ac.uk/proteins/api/proteins/"


def read_file(path):

    res = defaultdict(tuple)
    fr = open(path, 'r')

    for line in fr:
        line = line.strip()
        temp_arr = line.split("\t")
        res[temp_arr[0]] = float(temp_arr[1])

    return res


def query_api(res, out):

    ann_dict = defaultdict(list)
    fw = open(out, 'w')

    for id in tqdm(res.keys()):
        info = requests.get(_URL + id).json()
        refs = info['dbReferences']

        flag = True
        for ref in refs:
            if ref['type'] == 'Pfam':
                ann_dict[ref['id']] += [res[id]]
                flag = False
                fw.write(f"{id}\t{ref['id']}\t{res[id]}\n")

        if flag:
            ann_dict['None'] += [res[id]]
                # TODO: track avg±std eval per signature (make it a list instead of a count)
                #TODO: for overlap, calculate the difference for each accession

    return ann_dict


def write_sig_res(ann, out):

    fr = open(out, 'w')
    for sig in ann:
        fr.write(f"{sig}\t{len(ann[sig])}\t{np.mean(ann[sig])}±{np.std(ann[sig])}\n")


def main():

    pfam_res = read_file(_PFAMPATH)
    mdp_res = read_file(_MDPPATH)
    overlap_res = read_file(_OVERLAPPATH)

    pfam_ann = query_api(pfam_res, "./hmmsearch_results/pfam_evals.tab")
    mdp_ann = query_api(mdp_res, "./hmmsearch_results/mdp_evals.tab")
    overlap_ann = query_api(overlap_res, "./hmmsearch_results/overlap_evals.tab")

    write_sig_res(pfam_ann, "./hmmsearch_results/pfam_sigs.tab")
    write_sig_res(mdp_ann, "./hmmsearch_results/mdp_sigs.tab")
    write_sig_res(overlap_ann, "./hmmsearch_results/overlap_sigs.tab")

    # [print(x, len(pfam_ann[x])) for x in pfam_ann]
    # [print(x, len(mdp_ann[x])) for x in mdp_ann]
    # [print(x, len(overlap_ann[x])) for x in overlap_ann]


if __name__ == '__main__':
    main()
