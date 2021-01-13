# Reads the tsv file produced by InterproScan and creates an Interpro network
# out of it

import networkx as nx


_IPRPATH = "./ipr/PF00171_seqs.fasta.tsv"
_OUTPATH = "./ipr/PF00171_iprnet.gml"


def create_net():
    fr = open(_IPRPATH, "r")
    G = nx.Graph()

    curr_entry = ""
    curr_sigs = []
    flag = False

    for line in fr:
        line = line.strip()
        line_lst = line.split("\t")
        col_num = len(line_lst)

        if curr_entry != line_lst[0]:  # Reached the signatures of a new entry
            if flag:
                G.add_node(curr_entry, type="Entry", graphics={
                                                    "fill": "#e30909",
                                                    "shape": "ellipse"})

                for sig in curr_sigs:
                    if sig not in G.nodes():
                        G.add_node(sig, type="IPR", graphics={
                                                    "fill": "#1a04de",
                                                    "shape": "ellipse"})

                    G.add_edge(curr_entry, sig)

            flag = True  # Flag to add nodes to graph after the first entry
            curr_entry = line_lst[0]
            curr_sigs[:] = []

        if col_num == 13:  # Rows with IPR signature have 13 columns
            sig = line_lst[11]  # Index 11 contains the IPR signature
            curr_sigs.append(sig)

        elif line_lst[3] == "PANTHER":
            sig = line_lst[4]  # Index 4 contains the Panther signature
            curr_sigs.append(sig)

    G.add_node(curr_entry, type="Entry", graphics={
                                        "fill": "#e30909",
                                        "shape": "ellipse"})

    for sig in curr_sigs:
        if sig not in G.nodes():
            G.add_node(sig, type="IPR", graphics={
                                        "fill": "#1a04de",
                                        "shape": "ellipse"})

        G.add_edge(curr_entry, sig)

    nx.write_gml(G, _OUTPATH)

    return G


def main():
    create_net()


if __name__ == '__main__':
    main()
