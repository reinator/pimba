# PIMBA
<p align="center">
<img src="Figures/PIMBA.png" alt="PIMBA Logo" width="50%">
</p>
PIMBA, a PIpeline for MetaBarcoding Analysis which allows the use of customized databases, as well as other reference databases.

## PIMBA 3.0, the new version of PIMBA with Snakemake
PIMBA was originally implemented in Bash, which posed limitations in structure and processing speed. To enhance usability, reproducibility, and scalability, we have developed a Snakemake-based version of PIMBA. Snakemake provides structured workflow management, automated parallelization, and seamless integration with containerization technologies, making it significantly faster and more efficient than traditional Bash scripting. This new Snakemake-based pipeline (PIMBA 3.0) optimizes metabarcoding analyses, offering a powerful tool for biodiversity research, ecological monitoring, and health sciences applications. Our benchmarking analysis (comparison between version 3.0 and 2.0) demonstrated a substantial reduction in execution time, particularly in the PIMBA Prepare mode, which is crucial for processing large-scale metabarcoding projects with numerous samples. PIMBA 3.0 can be accessed here: https://github.com/itvgenomics/pimba_smk

## How to cite?
The peer-reviewed version of the paper can be found at https://doi.org/10.1007/978-3-030-91814-9_10
~~~
OLIVEIRA, R. R. M. et al. PIMBA: A PIpeline for MetaBarcoding Analysis. Advances in Bioinformatics and Computational Biology. 1ed. Switzerland: Springer, 2021, v. 13063, p. 106–116, 2021.
~~~

~~~
OLIVEIRA, RENATO RENISON MOREIRA; SILVA, R. L. ; NUNES, GISELE LOPES ; OLIVEIRA, GUILHERME .PIMBA: a PIpeline for MetaBarcoding Analysis. In: Stadler P.F., Walter M.E.M.T., Hernandez-Rosales M., Brigido M.M.. (Org.).Advances in Bioinformatics and Computational Biology. 1ed. Switzerland: Springer, 2021, v. 13063, p. 106-116
~~~

## How to install?
To run PIMBA, you just need to have docker (see https://docs.docker.com) installed in your operational system. 
~~~
sudo apt-get install docker.io
~~~
And that's all! Now you can run PIMBA on your data!

## Prepare your data (pimba_prepare.sh)
The first step to run PIMBA is to prepare your data. PIMBA can be used with paired-end or single-end reads (the latter being single-index or dual-index).
The output will be a fasta file that can be used in the next step.
### paired-end reads:
Please, place all your forward and reverse reads in one directory and make sure that forward reads contain "_R1" and reverse reads contain "_R2" in the file's name.
~~~
./pimba_prepare.sh illumina <rawdata_dir> <output_reads> <num_threads> <adapters.txt> <min_length> <min_phred>
~~~
<rawdata_dir> = path with all the R1 and R2 reads file;\
<output_reads> = name for the output file;\
<num_threads> = number of threads;\
<adapters.txt> = tab separated 2-column file with all adapters and primers used for sequencing;\
<min_lenght> = The minimum length of the read after quality treatment;\
<min_phred> = Minimum PHRED score of a read after quality treatment.

Example:
~~~
./pimba_prepare.sh ilumina rawdata/ AllSamples 24 adapters.txt 100 20
~~~

### single-end reads with dual-index:
In case your single-end reads have been multiplexed with dual-index, use the following command:
~~~
./pimba_prepare.sh iontorrent-dualindex  <rawdata.fastq> <barcodes.txt> <barcodes_reverse.txt> <barcodes.fasta> <barcodes_for_dir> <Primer_forward> <Primer_reverse> <num_threads> <output_name> <min_length> <min_phred>
~~~
<rawdata.fastq> = single file with all the reads to demultiplex;\
<barcodes.txt> = barcodes used as index in the 3' of the fragment;\
<barcodes_reverse.txt> = reverse complement of <barcodes.txt>;\
<barcodes.fasta> = fasta file for the <barcodes.txt>;\
<barcodes_for_dir> = path to all the barcodes.fasta and barcodes.txt used in the 5'  of the fragment. Each 3' barcode must have a fasta and txt file with all the associated 5' barcodes;\
<Primer_forward> = sequence of the forward primer;\
<Primer_reverse> = sequence of the reverse primer;\
<num_threads> = number of threads;\
<output_name> = name for the output fastq file;\
<min_length> = The minimum length of the read after quality treatment;\
<min_phred> = Minimum PHRED score of a read after quality treatment.

Example:
~~~
./pimba_prepare.sh iontorrent-dualindex rawdata_chip-3-4.fastq barcodes.txt barcodes_reverse.txt barcodes.fasta barcode_for/ TCCACTAATCACAAAGANATNGGNAC AGAAAATCATAATNAANGCNTGNGC 24 AllSamplesCOI.fastq 100 20
~~~

### single-end reads with single-index:
In case your single-end reads have been multiplexed with single-index, use the following command:
~~~
./pimba_prepare.sh iontorrent-singleindex <rawdata.fastq> <prefix> <barcodes.txt> <barcodes.fasta> <primer> <num_threads> <output_name> <min_length> <min_phred>
~~~
<rawdata.fastq> = single file with all the reads to demultiplex;\
<prefix> = name that will precede the barcodes names;\
<barcodes.txt> = barcodes used as index in the 5' of the fragment;\
<barcodes.fasta> = fasta file for the <barcodes.txt>;\
<primer> = primer sequence;\
<num_threads> = number of threads;\
<output_name> = name for the output fastq file;\
<min_lenght> = The minimum length of the read after quality treatment;\
<min_phred> = Minimum PHRED score of a read after quality treatment.

Example:
~~~
./pimba_prepare.sh iontorrent-singleindex SN1-45.fastq SN1-45-ITS barcodes.txt barcodes.fasta ATGCGATACTTGGTGTGAAT 24 AllSamples
~~~

## Run your metabarcoding analysis (pimba_run.sh)
The output generated by pimba_prepare.sh is a fasta file that will be used by pimba_run.sh in the following command:
~~~
./pimba_run.sh -i <input_reads> -o <output_dir> -w <approach> -s <otu_similarity> -a <assign_similarity> -c <coverage> -l <otu_length> -h <hits_per_subject> -g <marker_gene> -t <num_threads> -e <E-value> -d <databases.txt> -x <run_lulu>
~~~
-i <input_reads> = FASTA file with reads output generated by pimba_prepare.sh;\
-o <output_dir> = Directory where the results will be stored;\
-w <approach> = Analysis strategy to be used. It can be 'otu' or 'asv'. If 'otu', pimba uses vsearch. If 'asv', pimba uses swarm. Default: 'otu';\
-s <otu_similarity> = Percentage of similarity used in the otu clustering. The default is 0.97;\
-a <assign_similarity> = Percentage of similarity used in the taxonomy assignment. The default is 0.9;\
-c <coverage> = minimum coverage for the alignment. The default is 0.9;\
-l <otu_length> = Length to trim the reads. If 0, then no reads are trimmed;\
-h <hits_per_subject> = if 1, choose the best hit. If > 1, choose by majority. Default is 1;\
-g <marker_gene> = Marker gene and Database of the analysis. It can be: (16S-SILVA, 16S-GREENGENES, 16S-RDP, 16S-NCBI, ITS-FUNGI-NCBI, ITS-FUNGI-UNITE, ITS-PLANTS-NCBI, COI-NCBI). The path for each database must be configured in the <databases_file.txt>;\
-t <num_threads> = Number of threads to use in the blast step. Default is 1;\
-e <E-value> = Expected value used by Blast. The default is 0.00001;\
-d <databases_file.txt> = File with the path of the databases.\
-x <lulu> = if set as 'lulu', PIMBA will discard erroneous OTUs or ASVs with LULU. The default is not to use LULU.\
  
 The <databases_file.txt> must be properly configured. The fixed databases in the <databases_file.txt> are the following:
 
 ~~~
 #!/bin/bash

SILVA_DB_16S=/your/path/to/Silva_132_release/SILVA_132_QIIME_release/
GG_DB_16S=/your/path/to/gg_13_5_otus/
RDP_DB_16S=/your/path/to/RDP/
TAXDUMP=/your/path/to/taxdump/
ITS_UNITE_DB=/your/path/to/sh_refs_qiime_ver8.2/
~~~

You can download the databases above in the following links:
SILVA for 16S: https://www.arb-silva.de/no_cache/download/archive/qiime/ \
Greengenes for 16S: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz \
RDP for 16S: https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/RDPClassifier_16S_trainsetNo18_QiimeFormat.zip/download \
UNITE for fungal ITS: https://plutof.ut.ee/#/doi/10.15156/BIO/786385 \

For NCBI databases, make sure you download the following files to /your/path/to/taxdump/, using your terminal:
~~~
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
~~~
Then, uncompress:
~~~
tar -xzvf new_taxdump.tar.gz
~~~

### Configure your personalized database
Suppose you want to use a personalized database. In that case, you will only need a fasta file with the reference sequences and their identification, and a two-column tax.txt file with the sequence ID and the full taxonomy written for every reference sequence in the fasta file. Put them in the same directory, e.g.: /path/to/your/database/.
Example of FASTA file:

![](https://github.com/reinator/pimba/blob/main/Figures/fasta_example.png?raw=true)

Example of Taxonomy file:

![](https://github.com/reinator/pimba/blob/main/Figures/tax_example.png?raw=true)

Then, install blastn on your computer and run makeblastdb in your fasta file:
~~~
sudo apt-get install ncbi-blast+
makeblastdb -in <your_fasta.fasta> -dbtype nucl -parse_seqids
~~~
After that, all you need is to set the /path/to/your/database/ in the -g parameter when running pimba_run.sh. Example:
~~~
./pimba_run.sh -i AllSamplesCOI_chip1234_good.fasta -o AllSamplesCOI_98clust90assign -w otu -s 0.98 -a 0.9 -c 0.9 -l 130 -h 1 -g /path/to/your/database/ -t 24 -e 0.1 -d databases.txt
~~~

## Use files from the previous run to search in other reference databases (pimba_tax.sh)
Suppose you have already run pimba_run in your dataset from 16S genemarker, using the SILVA database. 
Now you can easily run the analysis by choosing a different database (e.g, 16S-GREENGENES) and reusing the previously generated FASTA file with the OTUs/ASVs and the text file with the Taxon table;

~~~
./pimba_tax.sh -i otus.fasta -u otu_table.txt -o AllSamplesCOI_98clust90assign -a 0.9 -c 0.9 -h 1 -g 16S-GREENGENES -t 24 -e 0.1 -d databases.txt
~~~

-i <otus_fasta> = FASTA file with otus;\
-u <otu_table> = otu/asv table;\
-o <output_dir> = Directory where the results will be stored;\
-a <assign_similarity> = Percentage of similarity used in the taxonomy assignment. The default is 0.9;\
-c <coverage> = minimum coverage for the alignment. The default is 0.9;\
-h <hits_per_subject> = if 1, choose the best hit. If > 1, choose by majority. Default is 1;\
-g <marker_gene> = Marker gene and Database of the analysis. It can be: (16S-SILVA, 16S-GREENGENES, 16S-RDP, 16S-NCBI, ITS-FUNGI-NCBI, ITS-FUNGI-UNITE, ITS-PLANTS-NCBI, COI-NCBI, COI-BOLD);\
-t <num_threads> = Number of threads to use in the blast step. Default is 1;\
-e <E-value> = Expected value used by Blast. The default is 0.00001;\
-d <databases_file.txt> = File with the path of the databases.\


## Plot your results (pimba_plot.sh)
When finished with pimba_run.sh, you will be able to generate some basic plots for your results, such as PCoA, rarefaction curves, alpha and beta diversity plots.
All you will need are two files that pimba_run.sh will generate and one metadata file that you will have to provide.

Example of OTU table (otu_table.txt) generated by pimba_run.sh:

![](https://github.com/reinator/pimba/blob/main/Figures/otutable_example.png?raw=true)

Example of OTU Tax assignment (tax_assignment.txt):

![](https://github.com/reinator/pimba/blob/main/Figures/taxresult_example.png?raw=true)

Example of Metadata file (metadata.csv):

![](https://github.com/reinator/pimba/blob/main/Figures/metadata_example.png?raw=true)
  
 The first column must always be “SampleID” and the second column must always be “SampleName”

Then, you can run the following command:
~~~
./pimba_plot.sh -t <otu_table> -a <tax_assignment> -m <metadata> -g <group_by>
~~~
-t <otu_table> = OTU table generated by pimba_run;\
-a <tax_assignment> = Tax assignment file generated by pimba_run;\
-m <metadata> = CSV file with columns "SampleID" and "SampleName", and other attributes related to each sample;\
-g <groupby> = A column from the metadata that will group the results. E.g., Description. If one does not want to group the results, do not specify it.\
  
 Example:
 ~~~
 ./pimba_plot.sh -t unionTriplicatas_otu_table.txt -a unionTriplicatas_otus_tax_assignments.txt -m mapping_file.csv -g Description
 ~~~
 
The list of plots that pimba_plot.sh will generate:

alpha_diversity_dotplot.svg\
rarefaction_curve2.svg\
cluster_dendogram.svg\
phylum_barplots.svg\
class_barplots.svg\
order_barplots.svg\
family_barplots.svg\
genus_barplots.svg\

## Place your taxa with Phylogenetic Placement (pimba_plwb.sh)
Once you have your OTUs/ASVs generated by pimba_run.sh, you can use a Constraint Tree and its alignments files to generate a Jplace file.

Example:
~~~
./pimba_plwb.sh -i data/AllSamples_otus_plants.fasta -c data/rbcL_ITS2_alignment.newick  -a data/rbcL_ITS2_alignment.fasta -x data/taxonpath.tsv -t 5 -d nt -o output
~~~

#-i = Inputa FASTA file with query sequences to be placed; \
#-c = Constraint tree in Newick format; \
#-a = Aligned FASTA file of sequences in constraint tree; \
#-x = File separated by tab of sequences in constraint tree and their associated taxon; \
#-t = Number of threads to be used; \
#-d = Type of sequences to be analyzed 'nt' or 'aa'; \
#-o = Output directory to store the results; 

The Jplace file will be saved to the path output/no_clustering/placed/<input_name>.jplace\

