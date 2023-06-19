# The Placement Workbench
A highly automated, configurable pair of pipelines, for creation and engineering of reference trees, and for their subsequent use with phylogenetic placement

## Setup
### 1) install conda/mamba
The pipelines are built on top of `snakemake` and `conda`. Thankfully, we can get `snakemake` from `conda`, which means we only have one real dependency! In principle everything should work with just the normal `conda` installation, however I've had much better results using `mamba`.
You can find the instructions for installing `mamba` [here](https://github.com/mamba-org/mamba), though to summarize, the current reccomendation is:
1) install [miniconda](https://docs.conda.io/en/latest/miniconda.html) (by downloading the installer and running it)
2) use `conda` to install mamba into the base environment:
```
conda install mamba -n base -c conda-forge
```

### 2) create the environment
Next, we create our special "place-workbench" environment (in the  repository root directory):
```
mamba env create -f environment.yml
```
### 3) enable the environment
Then, we activate that environment.
```
conda env plwb
```

This step has to be repeated every time you re-open your terminal or re-connect to a server. By default, mamba will add information about which is the active environment, to your command line:
```
(base) user@computer$ ... # by default, we are in the base environment
# now we enable our environment:
conda env plwb
(plwb) user@computer$ ... # tada!
```

### 4) thats it! ready to go
No, seriously, thats it.

## Basic Usage

### Examples
Both pipelines come with a convenient python script that should cover the majority of use-cases. 

By default, these scripts will output to a directory called `run-` followed by a timestamp. You can change the output directory by adding the `--out-dir ...` option. Use the `--help` command for a full list of options.

Some examples:
#### reftree-pipe
Infer phylogenetic trees from a set of unaligned sequences (`data/sequences.fasta`). It's DNA/Nucleotide data (`nt`), and we want to just use 4 threads:
```
./search.py --fasta-paths data/sequences.fasta  --align --datatype nt --threads 4
```

Same as the previous, except now we get some sequences by specifying their genbank accession labels in a file (`data/accessions.csv`), which are then automatically downloaded.
```
./search.py --csv-paths data/accessions.csv  --align --datatype nt --threads 4
```

You can also combine both commands, add multiple files each, and even add all fasta/csv files by the directory they're in:
```
./search.py --fasta-paths a.fasta b.fa data/ --csv-paths c.csv data/ ...
```

Finally, the pipeline also offers access to the Phylogenetic Automatic Reference Tree (PhAT) algorithm. In this mode, the input must be a database of aligned reference sequences (optionally gzipped, ending in `.gz`), and a fitting reference taxonomy must be provided (config: `data/taxonomy`, CLI: `--taxonomy-file`).
An example of the data can be found under `reftree-pipe/data/phat_test/`, and a call to the command line may look like this:
```
./search.py --fasta-paths data/phat_test/cyano.afa.gz --phat --taxonomy-file data/phat_test/cyano.tsv --datatype nt --threads 4
```
Note that, as the sequences already have to be aligned, it does not make sense to supply the `--align` flag.

#### place-pipe
```
./place.py --fasta-paths data/query.fa --reference-tree data/reference/tree.newick --reference-msa data/reference/seqs.fa --model-file data/reference/model --threads 4
```

If you also have a tab-separated file containing a mapping from the reference labels to their taxonomic paths ([example](place-pipe/data/reference/taxopaths.tsv)):
```
AY919771_clone_LG25_05_Alveolata	Eukaryota;ALVEOLATA;Alveolata_X;Alveolata_XX;Alveolata_XXX;Alveolata_XXXX;Alveolata_XXXXX;Alveolata_XXXXX+sp.;Uncultured;Uncultured+freshwater
HM245049_Alveolata_sp_CCMP3155	Eukaryota;ALVEOLATA;Alveolata_X;Alveolata_XX;Alveolata_XXX;Alveolata_XXXX;Chromerida;Chromerida+sp.;Chromerida;Chromerida+sp.
...
```

then you can add this file to the input:
```
./place.py [...] --taxonomy-file data/reference/taxopaths.tsv
```

...and taxonomic assignment ([explained here](https://github.com/lczech/gappa/wiki/Subcommand:-assign#final-output)) of the placed sequences will be added to the output.


Finally, the pipeline also supports option clustering of input query sequences:
```
./place.py [...] --sequence-clustering swarm
```

## Advanced Usage and Configuration
Both `search.py` and `place.py` are just convenience functions to manipulate and run the Snakemake configuration file. The pipelines can also be run directly from their respective folders. If you choose to do so, please consult the [Snakemake documentation](https://snakemake.readthedocs.io/en/v5.4.0/executable.html).

You can also find default configuration files for [reftree-pipe](reftree-pipe/config.yaml) and [place-pipe](place-pipe/config.yaml) that themselves are extensively documented, explaining the most commonly used options and parameters.
