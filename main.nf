#!/usr/bin/env nextflow

// Initialise Parameters
params.seqs = null
params.mat = null
params.head = null
params.ann = null
params.graph = null
params.k = null
params.outdir = null
params.solver = null

flag = false

if (params.seqs == null && params.mat != null & params.head !=null){

    mat_1 = Channel.fromPath( params.mat )
    head_1 = Channel.fromPath( params.head )
    seqs = Channel.empty()

}
else{
    seqs = file ( params.seqs )
    mat_1 = Channel.empty()
    head_1 = Channel.empty()
    flag = true

}

ann = file( params.ann )
graph = file ( params.graph )
k = params.k
solver = params.solver

// Create output directory and subdirectories
outdir = file( params.outdir )
outdir.mkdirs()
// subsets = file( outdir.name + '/subsets' )
// subsets.mkdirs()
// plots = file( outdir.name + '/plots' )
// plots.mkdirs()


process produce_identities{

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

    publishDir outdir, mode : "copy"
    container 'mdp-kit'

    input:
    file identities

    output:
    file "mat.npy" into mat_2
    file "headings.json" into head_2

    when:
    flag

    """
    python3 /code/lib/make_matrix.py -i ${identities}
    """
}

mat = mat_2.mix(mat_1).first()
head = head_2.mix(head_1).first()

process solve_mdp{

    publishDir outdir, mode : "copy"
    container 'mdp-kit'

    input:
    file mat
    file head
    file ann
    file graph

    output:
    file "*.txt"
    file "*.pdf"

    """
    python3 /code/main.py -a ${ann} -hd ${head} -d ${mat} -k ${k} -g ${graph} -s ${solver}
    """
}
