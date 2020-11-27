
from SPARQLWrapper import SPARQLWrapper, JSON
import matplotlib.pyplot as plt
import numpy as np
from random import seed


_ENDPOINT = "https://sparql.uniprot.org/sparql"
sig = 'PF00155'

def SPARQL_query_trembl(sig, swiss_len):

    trembl_len = 10000 - swiss_len

    sparql = SPARQLWrapper(_ENDPOINT)

    query = f"""
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

        SELECT DISTINCT ?protein ?taxon ?aa_seq
        WHERE
        {{
            ?protein a up:Protein .
            ?protein rdfs:seeAlso <http://purl.uniprot.org/pfam/{sig}> .
          	?protein up:reviewed false .
            ?protein up:organism ?taxon .
          	?taxon rdfs:subClassOf taxon:2 .
          	?protein up:sequence ?seq .
            ?seq rdf:value ?aa_seq .
            ?protein up:enzyme ?enz

        }}
            """

    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    seqs = results["results"]["bindings"]

    print(len(seqs))
    all_seqs = []
    keep_seqs = []

    for i in range(len(seqs)):
        aa_seq = seqs[i]['aa_seq']['value']
        all_seqs.append(aa_seq)

        if len(aa_seq) <= 500 and len(aa_seq) >=150:
            keep_seqs.append(seqs[i])

    print(len(keep_seqs))

    len_seqs = [len(seq) for seq in all_seqs]

    fig, ax = plt.subplots()
    ax.hist(len_seqs, bins=np.linspace(0, 1500, num=20))
    ax.set_title(f'TrEMBL Bacterial Sequence Length Distribution for {sig}')
    ax.set_xlabel('Length')
    ax.set_ylabel('Frequency')

    np.random.seed(6102019)

    sample = np.random.choice(keep_seqs, trembl_len, replace=False)

    return sample


def SPARQL_query_swissprot(sig):

    sparql = SPARQLWrapper(_ENDPOINT)

    query = f"""
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

        SELECT DISTINCT ?protein ?taxon ?aa_seq
        WHERE
        {{
            ?protein a up:Protein .
            ?protein rdfs:seeAlso <http://purl.uniprot.org/pfam/{sig}> .
          	?protein up:reviewed true .
            ?protein up:organism ?taxon .
          	?taxon rdfs:subClassOf taxon:2 .
          	?protein up:sequence ?seq .
            ?seq rdf:value ?aa_seq .
            ?protein up:enzyme ?enz

        }}
            """

    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    seqs = results["results"]["bindings"]

    print(len(seqs))
    all_seqs = []
    keep_seqs = []

    for i in range(len(seqs)):
        aa_seq = seqs[i]['aa_seq']['value']
        all_seqs.append(aa_seq)

        if len(aa_seq) <= 500 and len(aa_seq) >=150:
            keep_seqs.append(seqs[i])

    print(len(keep_seqs))

    len_seqs = [len(seq) for seq in all_seqs]

    fig, ax = plt.subplots()
    ax.hist(len_seqs, bins=np.linspace(0, 1500, num=20))
    ax.set_title(f'Swiss-Prot Bacterial Sequence Length Distribution for {sig}')
    ax.set_xlabel('Length')
    ax.set_ylabel('Frequency')

    return keep_seqs


def SPARQL_query_families():

    sparql = SPARQLWrapper(_ENDPOINT)

    query = """
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

        SELECT DISTINCT ?fam (COUNT(?fam) AS ?famcount)
        WHERE
            {{
            ?protein a up:Protein .
            ?protein up:reviewed true .
            ?protein rdfs:seeAlso ?fam .
            ?protein up:organism ?taxon .
            ?taxon rdfs:subClassOf taxon:2 .
            ?protein up:enzyme ?enz .
            FILTER (regex(?fam, "purl.uniprot.org/pfam/"))
            }}
        GROUP BY ?fam ORDER BY DESC(?famcount)
    """

    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    fams = results["results"]["bindings"]

    return fams


def SPARQL_query_ec_count(sig):

    sparql = SPARQLWrapper(_ENDPOINT)

    query = f"""
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        SELECT ?protein ?enz
        WHERE
        {{
            ?protein a up:Protein .
            ?protein rdfs:seeAlso <http://purl.uniprot.org/pfam/{sig}> .
          	?protein up:reviewed true .
            ?protein up:organism ?taxon .
          	?taxon rdfs:subClassOf taxon:2 .
          	?protein up:sequence ?seq .
            ?seq rdf:value ?aa_seq .
            ?protein up:enzyme ?enz.
        }}
    """

    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    ecs = results["results"]["bindings"]

    return ecs


def write_fasta(seqs):

    fw = open(f"{sig}_seqs_trembl.fasta", "w")

    for seq in seqs:
        id = seq['protein']['value'].split('/')[-1]
        seq = seq['aa_seq']['value']

        fw.write(f">{id}\n{seq}\n")


def download_dataset(sig):

    swiss_seqs = SPARQL_query_swissprot(sig)
    trembl_seqs = SPARQL_query_trembl(sig, len(swiss_seqs))

    print(len(swiss_seqs))
    print(len(trembl_seqs))

    combined_seqs = list(swiss_seqs) + list(trembl_seqs)
    write_fasta(combined_seqs)
    # plt.show()


def main():
    download_dataset(sig)


if __name__ == '__main__':
    main()
