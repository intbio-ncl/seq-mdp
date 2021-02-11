#!/usr/bin/env nextflow

// Initialise Parameters
params.mat = null
params.head = null
params.ann = null
params.graph = null
params.k = null
params.outdir = null
params.solver = null


mat = file( params.mat )
head = file( params.head )
ann = file( params.ann )
graph = file ( params.graph )
k = params.k
solver = params.solver

// Create output directory and subdirectories
outdir = file( params.outdir )
outdir.mkdirs()
subsets = file(outdir.name + '/subsets')
subsets.mkdirs()

process solve_mdp{

    publishDir subsets, mode : "copy"

    input:
    file mat
    file head
    file ann
    file graph

    output:
    file "*.txt"

    """
    python3 /code/main.py -a ${ann} -hd ${head} -d ${mat} -k ${k} -g ${graph} -s ${solver}
    """
}
