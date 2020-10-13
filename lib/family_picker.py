
import json
from collections import defaultdict

from download_dataset import SPARQL_query_families, SPARQL_query_ec_count


def read_json_file(path):

    final_data = []

    with open(path) as f:
        data = json.load(f)

    for fam in data:
        final_data.append(fam['metadata']['accession'])

    return final_data


def main():

    fams = SPARQL_query_families()
    top100fams = [fam['fam']['value'] for fam in fams][0:100]
    seq_count_dict = {fam['fam']['value'].split('/')[-1]:fam['famcount']['value'] for fam in fams}

    ecs_dict = defaultdict(int)

    for fam in top100fams:
        fam = fam.split('/')[-1]
        ecs = SPARQL_query_ec_count(fam)
        print(fam, len(ecs), seq_count_dict[fam])
        ecs_dict[fam] += len(ecs)

    print(ecs_dict)


if __name__ == '__main__':
    main()
