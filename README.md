# GenoFan - Genome functional annotation pipeline

This is a Snakemake workflow for functional annotation of a genome assembly using [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper/wiki), [dbCAN2](https://bcb.unl.edu/dbCAN2/), [KOFams](https://www.genome.jp/tools/kofamkoala/), [AMRFinder](https://github.com/ncbi/amr), [gutSMASH](https://github.com/victoriapascal/gutsmash) and [antiSMASH](https://github.com/antismash).

## Installation

1. Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

2. Clone the [gutSMASH repo](https://github.com/victoriapascal/gutsmash) and install its [databases](https://gutsmash.bioinformatics.nl/download.html). Note: the instructions do not make it clear that the ClusterBlast database needs to be stored in the folder `gutsmash/antismash/databases/clusterblast`.

3. Clone this repository
```
git clone https://github.com/alexmsalmeida/genofan.git
```

4. Download main database containing eggNOG, CAZy and KOFam annotations
```
wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/protein_dbs.tar.gz
tar -xzvf protein_dbs.tar.gz
```

## How to run

1. Edit the configuration file [`config/config.yml`](config/config.yml).
    - `input_file`: Text file with each line indicating the path to the input genome FASTA files.
    - `db_dir`: Location of the database directory.
    - `gutsmash_dir`: Location where the gutSMASH repo was cloned.
    - `ncores`: Number of cores to use for the analyses.

2. (option 1) Run the pipeline locally (adjust `-j` based on the number of available cores)
```
snakemake --use-conda -k -j 4
```
2. (option 2) Run the pipeline on a cluster (e.g., SLURM)
```
snakemake --use-conda -k -j 100 --profile config/slurm --latency-wait 120
```

3. Results will be stored in the directory of each fasta file under the suffix `*_annotations`.
