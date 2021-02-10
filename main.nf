#!/usr/bin/env nextflow

params.mat = null
params.head = null
params.ann = null
params.graph = null
params.k = null

mat = file( params.mat )
head = file( params.head )
ann = file( params.ann )
graph = file ( params.graph )
k = params.k

process test{

    input:
    file mat
    file head
    file ann
    file graph

    """
    python3 /code/main.py -a ${ann} -hd ${head} -d ${mat} -k ${k} -g ${graph}
    """
}
