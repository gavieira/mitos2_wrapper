# mitos2_wrapper

Wrapper for MITOS2 annotation tool - generates a .gbk file from plain fasta sequences

This script still needs MITOS to be installed in your machine. In order to install MITOS, first install the [Anaconda](https://www.anaconda.com/distribution/) distribution, then configure the [Bioconda](https://bioconda.github.io/user/install.html#set-up-channels) channel. Then install MITOS in a new environment (which we'll name 'mitos_local', but could have any other name):

```
conda install mitos=2.0.3 -c bioconda -m -n mitos_local
```

Every time you want to run this wrapper, you will have to activate this new environment:

```
conda activate mitos_local
```

This wrapper does not use the default parameters for MITOS2. It uses several custom options to help clean the annotation output files (no origin of replication is annotated, no plots are created and only the copy of each feature will be annotated). Also, hmmer methods and ncbi's start/stop codons are used in the annotation.

OBS: DO NOT copy only the 'mitos2_wrapper.py' file from this repository. The program won't work without reference data.

```
usage: mitos2_wrapper.py [-h] [-k] [-c GENETIC CODE] [-r REFERENCE DATA]
                         [FASTA [FASTA ...]]

A wrapper around MITOS (help on how to install and run here:
https://gitlab.com/Bernt/MITOS) to annotate metazoan mitochondrial contigs
using hmm (--alarab option)

positional arguments:
  FASTA                 Fasta file(s) to be annotated

optional arguments:
  -h, --help            show this help message and exit
  -k, --keep            Keeps MITOS annotation files. Default: False
  -c GENETIC CODE, --gencode GENETIC CODE
                        NCBI's codon table. Default: 2 (Vertebrate
                        Mitochondrial)
  -r REFERENCE DATA, --refdir REFERENCE DATA
                        Custom path to RefSeq63m directory (downloadable from
                        https://zenodo.org/record/2672835). Default:
                        script path
```
