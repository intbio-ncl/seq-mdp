#!/usr/bin/env nextflow

/*
Nextflow workflow that solves the Maximum Diversity Problem (MDP) for a set
of protein sequences using sequence identity as a distance metric.

Required parameters are either a FASTA file with sequences OR the matrix/heading
pair of files from a previous run, a solution subset size K, an output directory,
and the solver to be used (greedy, ts, or all).

Can optionally submit an annotation file, for which Coverage and Gini-Simpson
values will be computed.
*/



// Initialise Parameters
params.seqs = null
params.mat = null
params.head = null
params.ann = 'NO_FILE' // annotation file is optional
params.k = null
params.outdir = null
params.solver = null

flag = false

/*
In case the matrix and heading files have been created in a past run, then
they can be reused as inputs to avoid recomputing the identities. this bit makes
that happen.
*/

if (params.seqs == null && params.mat != null & params.head !=null){
    mat_1 = Channel.fromPath( params.mat )
    head_1 = Channel.fromPath( params.head )
    seqs = Channel.empty()
}
else {
    seqs = file ( params.seqs )
    mat_1 = Channel.empty()
    head_1 = Channel.empty()
    flag = true
}

ann = file( params.ann )
k = params.k
solver = params.solver

// Create output directory
outdir = file( params.outdir )
outdir.mkdirs()


process produce_identities{
    // Use needleall to produce all-vs-all global alignment sequence identities

    publishDir outdir, mode : "copy"

    input:
    file seqs

    output:
    file "./identities.txt" into identities

    when:
    flag

    """
    nextflow run ravenlocke/nf-needleall-ava --infile ${seqs}  --outdir needle_out --threshold 0.00000001 --cpu 4
    cp ./needle_out/identities.txt ./
    """
}

process make_matrix{
    // Restructure result of needleall to be used in MDP solving

    publishDir outdir, mode : "copy"

    input:
    file identities

    output:
    file "mat.npy" into mat_2
    file "headings.json" into head_2

    when:
    flag

    """
    python3 /seq-mdp/lib/make_matrix.py -i ${identities}
    """
}

// Merge the two matrix and heading channels so reusing previous runs can work
mat = mat_2.mix(mat_1).first()
head = head_2.mix(head_1).first()

process solve_mdp{
    // Solve the MDP

    publishDir outdir, mode : "copy"

    input:
    file mat
    file head
    file ann

    output:
    file "*.txt"
    file "*.pdf"

    script:
    def ann_filter = ann.name != 'NO_FILE' ? "-a $ann" : ''

    """
    python3 /seq-mdp/main.py $ann_filter -hd ${head} -d ${mat} -k ${k} -s ${solver}
    """
}
