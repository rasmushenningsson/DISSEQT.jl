# DISSEQT.jl
The DISSEQT.jl package is an implementation of the pipeline described in the paper [DISSEQT – DIStribution based modeling of SEQuence Space Time dynamics](https://www.biorxiv.org/content/10.1101/327338v1).



# Installation
Start Julia (1.1 or later) and enter the Package REPL by pressing ]. Then install DISSEQT.jl by entering:
```
add https://github.com/rasmushenningsson/SynapseClient.jl.git
add https://github.com/rasmushenningsson/DISSEQT.jl.git
```
Also see installation instructions for [SynapseClient.jl](https://github.com/rasmushenningsson/SynapseClient.jl) if you want to enable the [Synapse](https://www.synapse.org) features in DISSEQT.


# Examples
The complete analysis of deep sequencing data from the [DISSEQT paper](https://www.biorxiv.org/content/10.1101/327338v1) is available in Synapse [here](https://www.synapse.org/#!Synapse:syn11639899). 
Note how the Provenance system in Synapse makes it possible to trace the steps used to produce every result in Synapse, showing how the analysis was done (which script was called) and listing all input files.
All example scripts below upload their results to Synapse. The scripts themselves are also automatically uploaded to ensure that all analyses can be rerun elsewhere.
The steps in the DISSET pipeline are outlined below:

## Alignment
If you already have BAM files, you can start the DISSEQT pipeline from the next step. Note however that DISSEQT performs iterative alignment. That is, if the consensus sequence of an aligned mutant swarm is different from the reference sequence used during alignment, it is realigned using the new consensus sequence as the reference. The process is repeated until the reference does not change. Iterative alignment improves inference of codon frequencies close to consensus changes.

To run alignment locally, you need to have [bwa](https://github.com/lh3/bwa), [samtools](http://www.htslib.org) and [fastq-mcf](https://expressionanalysis.github.io/ea-utils/) installed and available in your path.

Example scripts and other relevant files for running alignment using DISSEQT can be found [here](https://www.synapse.org/#!Synapse:syn18694207). It is recommended to use one [script](https://www.synapse.org/#!Synapse:syn18695094) for each run. The [Reference Genomes](https://www.synapse.org/#!Synapse:syn18694208) and [Adapter](https://www.synapse.org/#!Synapse:syn18694218) files are also needed.

The outputs of the Alignment step are [BAM Files](https://www.synapse.org/#!Synapse:syn18694439), the consensus sequence and a detailed alignment log for each sample are also saved in the same folder.
An overview log file - [AlignUtils.log](https://www.synapse.org/#!Synapse:syn18695095) - is also created. 

## Codon Frequency Inference


# Contact
If you have problems running DISSEQT, please open an issue in the Issue Tracker or contact rasmus.henningsson@med.lu.se.
