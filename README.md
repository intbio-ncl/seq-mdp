# seq-mdp
Nextflow workflow that solves the Maximum Diversity Problem (MDP) for protein sequence sampling.

## Requirements
The only requirements are [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/).

## Parameters

```
--seqs :    File containing protein sequences in FASTA format.

--k :       Subset size K for which the MDP should be solved.

--solver :  Which algorithm should be used to solve the MDP.
            Choices are 'greedy' for the Greedy Max-Min Solver,
            'ts' for our Tabu Search implementation, and 'all' for both.

--outdir :  Directory to store the results in.

--mat :     Numpy matrix containing the sequence identities from a previous run 'from scratch'.
            Should be paired with --head as input to rerun the workflow with the same sequences faster.

--head :    JSON file containing the headings from a previous run 'from scratch'.
            Should be paired with --mat as input to rerun the workflow with the same sequences faster.

--ann :     Table file containing class annotation about the sequences. Use to produce diversity assessment of solution.
            Format should just be two columns - Sequence ID in the first, classes in the second.
            If a sequence can have more than one class, separate them using ';'.
            Sequence IDs *MUST* match the sequence IDs in the FASTA headers/headings JSON file.
```

## How to run

This workflow can be run with two different sets of inputs:

### From Scratch

If you run this workflow from scratch, one of your inputs should be a FASTA file containing the protein sequences you want to sample from. If so, you want to use the `--seqs` input.

```
nextflow run intbio-ncl/seq-mdp --seqs test_seqs.fasta --k 20 --outdir Test --solver all
```
### Rerun

The sequence identity matrix and its headings can be reused in a new run. As the identity matrix can take a long time to produce, we recommend you use the 'rerun' set of inputs if you want to solve the MDP again on the same sequences, but with different other inputs.
In this case, the parameter `--seqs` should be replaced with two inputs: `--mat` for the Numpy matrix and `--head` for the JSON headings.

```
nextflow run intbio-ncl/seq-mdp --mat test_mat.npy --head test_headings.json --k 20 --outdir Test --solver all
```
You can find the .npy and .json files necessary for reruns in the output directory after your initial run.

### Diversity Assessment

If given an annotation file as described in the Parameters section, the parameter `--ann` can be used to produce metrics that give an idea of the overall diversity of the solution chosen by the workflow.

For example, if you used GO terms as the annotation, it would summarise the diversity in GO terms present in the solution.
Specifically, it produces:

* __Coverage__ : The proportion of the total classes present in the solution e.g. how many of GO classes contained in the annotation file are present in the subset chosen.
* __Gini-Simpson Index__ : The proportional abundance of classes present in the solution e.g. of all the GO classes present in the subset chosen, how equally abundant are they in their proportion.

## Citation

If you use this tool in your research, please cite our latest paper:

[Atallah, Christian, et al. "Automatic Diverse Subset Selection From Enzyme Families by Solving the Maximum Diversity Problem." 2022 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB). IEEE, 2022.](https://ieeexplore.ieee.org/abstract/document/9863021)
