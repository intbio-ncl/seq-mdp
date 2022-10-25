
import networkx as nx
from tqdm import tqdm
#########################################

import argparse

parser = argparse.ArgumentParser(description="MDP Filler stuff") #TODO: Change desc
parser.add_argument("-i", "--input", help="Path to the identities file", type=str, required=True)
parser.add_argument("-t", "--threshold", help="Threshold for SSN", type=str, required=True)
parser.add_argument("-o", "--output", help="Path to output file", type=str, required=True)

args = parser.parse_args()
_INPUT = args.input
_THRESHOLD = args.threshold
_OUTPUT = args.output

#########################################


def main():

    fr = open(_INPUT, 'r')
    G = nx.Graph()

    for line in tqdm(list(fr)):

        line = line.strip()
        temp_arr = line.split(' ')

        source = temp_arr[0]
        sink = temp_arr[1]
        thr = temp_arr[2]

        if source not in list(G.nodes()):
            G.add_node(source)
        if sink not in list(G.nodes()):
            G.add_node(sink)

        if thr >= _THRESHOLD:
            G.add_edge(source, sink, identity=thr)

    nx.write_graphml(G, _OUTPUT)


if __name__ == '__main__':
    main()
