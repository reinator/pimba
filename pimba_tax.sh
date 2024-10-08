#Authors: Renato Oliveira, Tiago Leão, Gisele Nunes, Raíssa Oliveira
#version: 2.0.13
#Date: 30-09-2024

###    Copyright (C) 2021  Renato Oliveira
###
###    This program is free software: you can redistribute it and/or modify
###    it under the terms of the GNU General Public License as published by
###    the Free Software Foundation, either version 3 of the License, or
###    any later version.
###
###    This program is distributed in the hope that it will be useful,
###    but WITHOUT ANY WARRANTY; without even the implied warranty of
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###    GNU General Public License for more details.
###
###    You should have received a copy of the GNU General Public License
###    along with this program.  If not, see <http://www.gnu.org/licenses/>.

###Contacts:
###    guilherme.oliveira@itv.org
###    renato.renison@gmail.com

#!/bin/bash
##usage: ./pimba_tax.sh -i <otus_fasta> -u <otu/asv_table> -o <output_dir> -s <otu_similarity> -a <assign_similarity> -c <coverage> -h <hits_per_subject> -g <marker_gene> -t <num_threads> -e <E-value> -d <databases.txt> -r <yes or no> -b <NCBI database prefix>
#-i <otus_fasta> = FASTA file with otus.
#-u <otu_table> = otu/asv table.
#-o <output_dir> = Directory where the results will be stored.
#-a <assign_similarity> = Percentage of similarity that will be used in the taxonomy assignment. Default is 0.9
#-c <coverage> = minimum converage for the alignment. Default is 0.9.
#-h <hits_per_subject> = if 1, choose the best hit. If > 1, choose by majority. Default is 1
#-g <marker_gene> = Marker gene and Database of the analisys. It can be: (16S-SILVA, 16S-GREENGENES, 16S-RDP, 16S-NCBI, ITS-FUNGI-NCBI, ITS-FUNGI-UNITE, ITS-PLANTS-NCBI, COI-NCBI, COI-BOLD)
#-t <num_threads> = Number or threads to use in the blast step. Default is 1.
#-e <E-value> = Expected value used by blast. Dafault is 0.00001.
#-d <databases_file.txt> = File with the databases path. Default is /bio/pimba_metabarcoding/databases.txt
#-r <yes or no> = set 'yes' to use the remote NCBI database. Default is 'no'
#-b < NCBI database prefix> = It can be 'nt', 'nt_euk', 'core_nt','nr' depending on the database you have downloaded.
#USAGE: ./pimba_tax.sh -i otus.fasta -u otu_table.txt -o AllSamplesCOI_98clust90assign -a 0.9 -c 0.9 -h 1 -g COI-ALL -t 24 -e 0.1 -d databases.txt -r no

#source activate qiime1

HITS_PER_SUBJECT=1
SIMILARITY=0.97
SIMILARITY_ASSIGN=0.9
SIMILARITY_INT=$(bc -l <<<"${SIMILARITY}*100")
SIMILARITY_INT=${SIMILARITY_INT%.*}
SIMILARITY_INT_ASG=$(bc -l <<<"${SIMILARITY_ASSIGN}*100")
SIMILARITY_INT_ASG=${SIMILARITY_INT_ASG%.*}
COVERAGE=0.9
COVERAGE_INT=$(bc -l <<<"${COVERAGE}*100")
COVERAGE_INT=${COVERAGE_INT%.*}
THREADS=1
EVALUE=0.00001
TIMESTAMP=$(date +"%Y%m%d%H%M%S")
REMOTE=no
DB_PREFIX='nt'
#DB_FILE=/bio/pimba_metabarcoding/databases.txt

while getopts "i:u:o:s:a:c:h:g:t:e:d:r:b:" opt; do
	case $opt in
		i) RAWDATA="$OPTARG"
		;;
		u) TAXTABLE="$OPTARG"
		;;
		o) OUTPUT="$OPTARG"
		;;
		s) SIMILARITY="$OPTARG"
		;;
		a) SIMILARITY_ASSIGN="$OPTARG"
		;;
		c) COVERAGE="$OPTARG"
		;;
		h) HITS_PER_SUBJECT="$OPTARG"
		;;
		g) GENE="$OPTARG"
		;;
		t) THREADS="$OPTARG"
        ;;
        e) EVALUE="$OPTARG"
        ;;
        d) DB_FILE="$OPTARG"
        ;;
        r) REMOTE="$OPTARG"
		;;
		b) DB_PREFIX="$OPTARG"
		;;
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

SIMILARITY_INT=$(bc -l <<<"${SIMILARITY}*100")
SIMILARITY_INT=${SIMILARITY_INT%.*}

COVERAGE_INT=$(bc -l <<<"${COVERAGE}*100")
COVERAGE_INT=${COVERAGE_INT%.*}

SIMILARITY_INT_ASG=$(bc -l <<<"${SIMILARITY_ASSIGN}*100")
SIMILARITY_INT_ASG=${SIMILARITY_INT_ASG%.*}

source $DB_FILE

# SILVA_DB_16S=/bio/share_bio/databases/Silva_132_release/SILVA_132_QIIME_release
# GG_DB_16S=/bio/share_bio/databases/gg_13_5_otus
# RDP_DB_16S=/bio/share_bio/databases/RDP
# NCBI_DB=/bio/share_bio/databases/NCBI/blast_nt/nt
# NCBI_DB_EXP=/bio/share_bio/databases/NCBI/blast_nt
# ITS_UNITE_DB=/bio/share_bio/databases/its_fungi/sh_refs_qiime_ver8.2
# ITS_BIOBASE_DB=/bio/share_bio/databases/metabarcoding_db/ITS_plants/ITS_plants
# COI_BIOBASE_DB=/bio/share_bio/databases/metabarcoding_db/COI_ITV/coi_itv
# COI_BOLD_DB=/bio/share_bio/databases/metabarcoding_db/COI_BOLD/coi_bold
# COI_ALL_DB=/bio/share_bio/databases/metabarcoding_db/COI_all/coi_all

CURRENT_PATH=$(pwd)

# DIR_NAME_RAW=$(dirname $RAWDATA)
# cd $DIR_NAME_RAW
# FULL_PATH_RAW=$(pwd)
# cd $CURRENT_PATH

# pathlist=$(echo $FULL_PATH_RAW; echo $CURRENT_PATH)
#COMMON_PATH=$(i=2; while [ $i -lt 500 ]; do   path=`echo "$pathlist" | cut -f1-$i -d/ | uniq -d`;   if [ -z "$path" ];   then      echo $prev_path;      break;   else      prev_path=$path;   fi;   i=`expr $i + 1`; done);


#COMMON_PATH=$({ echo $FULL_PATH_RAW; echo $CURRENT_PATH;} | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D')

#echo Common Path: $COMMON_PATH
echo Current Path: $CURRENT_PATH

mkdir $OUTPUT; chmod -R 755 $OUTPUT;
echo "RAWDATA=$RAWDATA"
cp $RAWDATA $TAXTABLE $OUTPUT
cd $OUTPUT

RAWDATA=$(basename $RAWDATA)
TAXTABLE=$(basename $TAXTABLE)
echo "RAWDATA=$RAWDATA"
newfile="$(basename $RAWDATA .fasta)"

#Dereplication <<<USING USEARCH 7>>>
echo "Creating a VSEARCH Container: "
docker run -id -v $CURRENT_PATH:/output/ --name vsearch_run_$TIMESTAMP itvdsbioinfo/pimba_vsearch:v2.15.2

#Convert UC to otu-table.txt <<< BMP SCRIPT>>>
echo "Creating a QiimePipe Container: "
docker run -id -v $CURRENT_PATH:/output/ --name qiimepipe_run_$TIMESTAMP itvdsbioinfo/pimba_qiimepipe:v2

if [ $GENE = "ITS-FUNGI-NCBI" ] || [ $GENE = "16S-NCBI" ] || [ $GENE = "COI-NCBI" ] || [ $GENE = "ITS-PLANTS-NCBI" ] || [ $GENE = "ALL-NCBI" ];
then
	PARAM_BLAST="/blastdb/$DB_PREFIX"
	EXPORTBLAST='export BLASTDB=/blastdb/;'
	PARAM_THREADS="-num_threads $THREADS"

	if [ $REMOTE = "yes" ];
	then
		echo "BLAST running in remote mode"
		PARAM_BLAST="nt -remote"
		PARAM_THREADS=''
		
	 	echo "Creating a BLAST Container: "
		docker run -id -v $CURRENT_PATH:/output/ --name blast_run_$TIMESTAMP itvdsbioinfo/pimba_blast:latest

		EXPORTBLAST=''
	else
		echo "BLAST running in local mode"
		export BLASTDB=$NCBI_DB_EXP
		echo "BLASTDB=$BLASTDB"
		echo "PARAMBLAST=$PARAM_BLAST"

		echo "Creating a BLAST Container: "
		docker run -id -v $CURRENT_PATH:/output/ -v $BLASTDB:/blastdb/ --name blast_run_$TIMESTAMP itvdsbioinfo/pimba_blast:latest

	fi
fi


if [ $GENE = "16S-SILVA" ];
then
	echo "Creating a Qiime Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $SILVA_DB_16S:/database/ --name qiime_run_$TIMESTAMP itvdsbioinfo/pimba_qiime:latest

	#Assign taxonomy to OTUS using blast method on QIIME
	echo "Running the Qiime Container - assign_taxonomy.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	assign_taxonomy.py -i '$RAWDATA' -o output -t /database/taxonomy/16S_only/'${SIMILARITY_INT}'/taxonomy_7_levels.txt \
	-r /database/rep_set/rep_set_16S_only/'${SIMILARITY_INT}'/silva_132_'${SIMILARITY_INT}'_16S.fna --similarity='$SIMILARITY_ASSIGN'; \
	chmod -R 777 output'
	#assign_taxonomy.py -i ${newfile}_otus.fasta -o output -t ${SILVA_DB_16S}/taxonomy/16S_only/${SIMILARITY_INT}/taxonomy_7_levels.txt -r ${SILVA_DB_16S}/rep_set/rep_set_16S_only/${SIMILARITY_INT}/silva_132_${SIMILARITY_INT}_16S.fna --similarity=$SIMILARITY_ASSIGN
	
	#Align sequences on QIIME, using greengenes reference sequences (use the file “otus.fa” from UPARSE as input file)
	echo "Running the Qiime Container - align_seqs.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	align_seqs.py -i '$RAWDATA' -o rep_set_align -t /database/rep_set_aligned/'${SIMILARITY_INT}'/'${SIMILARITY_INT}'_alignment.fna; \
	chmod -R 777 rep_set_align'
	#align_seqs.py -i ${newfile}_otus.fasta -o rep_set_align -t ${SILVA_DB_16S}/rep_set_aligned/${SIMILARITY_INT}/${SIMILARITY_INT}_alignment.fna

	#Filter alignments on QIIME
	echo "Running the Qiime Container - filter_alignment.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	filter_alignment.py -i rep_set_align/'${newfile}'_otus_aligned.fasta -o filtered_alignment; \
	chmod -R 777 filtered_alignment'
	#filter_alignment.py -i rep_set_align/${newfile}_otus_aligned.fasta -o filtered_alignment

	#Make the reference tree on QIIME
	echo "Running the Qiime Container - make_phylogeny.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	make_phylogeny.py -i filtered_alignment/'${newfile}'_otus_aligned_pfiltered.fasta -o rep_set.tre; \
	chmod -R 777 rep_set.tre'
	#make_phylogeny.py -i filtered_alignment/${newfile}_otus_aligned_pfiltered.fasta -o rep_set.tre

	mkdir diversity_by_sample
	cd diversity_by_sample

	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiime

	echo "Running the QiimePipe Container - createAbundanceFile.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample; \
	python3.6 /qiimepipe/createAbundanceFile.py ../output/'${newfile}'_tax_assignments.txt ../'${TAXTABLE}';\
	chmod -R 777 ../diversity_by_sample'
	#python ${SCRIPT_PATH}/createAbundanceFile.py ../output/${newfile}_tax_assignments.txt ../${newfile}_otu_table.txt
	cd ..

	docker stop qiime_run_$TIMESTAMP
	docker rm qiime_run_$TIMESTAMP


elif [ $GENE = "16S-GREENGENES" ];
then

	echo "Creating a Qiime Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $GG_DB_16S:/database/ --name qiime_run_$TIMESTAMP itvdsbioinfo/pimba_qiime:latest

	#Assign taxonomy to OTUS using blast method on QIIME
	echo "Running the Qiime Container - assign_taxonomy.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	assign_taxonomy.py -i '$RAWDATA' -o output -t /database/taxonomy/'${SIMILARITY_INT}'_otu_taxonomy.txt \
	-r /database/rep_set/'${SIMILARITY_INT}'_otus.fasta --similarity='$SIMILARITY_ASSIGN'; \
	chmod -R 777 output'
	#assign_taxonomy.py -i ${newfile}_otus.fasta -o output -t ${GG_DB_16S}/taxonomy/${SIMILARITY_INT}_otu_taxonomy.txt -r ${$GG_DB_16S}/rep_set/${SIMILARITY_INT}_otus.fasta --similarity=$SIMILARITY_ASSIGN

	#Align sequences on QIIME, using greengenes reference sequences (use the file “otus.fa” from UPARSE as input file)
	echo "Running the Qiime Container - align_seqs.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	align_seqs.py -i '$RAWDATA' -o rep_set_align -t /database/rep_set_aligned/'${SIMILARITY_INT}_otus.fasta'; \
	chmod -R 777 rep_set_align'
	#align_seqs.py -i ${newfile}_otus.fasta -o rep_set_align -t ${GG_DB_16S}/rep_set_aligned/${SIMILARITY_INT}_otus.fasta

	#Filter alignments on QIIME
	echo "Running the Qiime Container - filter_alignment.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	filter_alignment.py -i rep_set_align/'${newfile}'_otus_aligned.fasta -o filtered_alignment; \
	chmod -R 777 filtered_alignment'
	#filter_alignment.py -i rep_set_align/${newfile}_otus_aligned.fasta -o filtered_alignment

	#Make the reference tree on QIIME
	echo "Running the Qiime Container - make_phylogeny.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	make_phylogeny.py -i filtered_alignment/'${newfile}'_otus_aligned_pfiltered.fasta -o rep_set.tre; \
	chmod -R 777 rep_set.tre'
	#make_phylogeny.py -i filtered_alignment/${newfile}_otus_aligned_pfiltered.fasta -o rep_set.tre

	mkdir diversity_by_sample
	cd diversity_by_sample

	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiime
	echo "Running the QiimePipe Container - createAbundanceFile.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample; \
	python3.6 /qiimepipe/createAbundanceFile.py ../output/'${newfile}'_tax_assignments.txt ../'${TAXTABLE}';\
	chmod -R 777 ../diversity_by_sample'
	#python ${SCRIPT_PATH}/createAbundanceFile.py ../output/${newfile}_tax_assignments.txt ../${newfile}_otu_table.txt
	cd ..

	docker stop qiime_run_$TIMESTAMP
	docker rm qiime_run_$TIMESTAMP

elif [ $GENE = "16S-RDP" ];
then

	echo "Creating a Qiime Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $RDP_DB_16S:/database/ --name qiime_run_$TIMESTAMP itvdsbioinfo/pimba_qiime:latest


	#Assign taxonomy to OTUS using blast method on QIIME
	echo "Running the Qiime Container - assign_taxonomy.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	assign_taxonomy.py -i '$RAWDATA' -o output -t /database/*.t* \
	-r /database/rep_set/*.fa* --similarity='$SIMILARITY_ASSIGN'; \
	chmod -R 777 output'
	#assign_taxonomy.py -i ${newfile}_otus.fasta -o output -t ${RDP_DB_16S}/trainset16_022016.rdp.tax -r ${RDP_DB_16S}/trainset16_022016.rdp.fasta --similarity=$SIMILARITY_ASSIGN

	#Align sequences on QIIME, using RDP reference sequences (use the file “otus.fa” from UPARSE as input file)
	echo "Running the Qiime Container - align_seqs.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	align_seqs.py -i '$RAWDATA' -o rep_set_align -t /database/*align*; \
	chmod -R 777 rep_set_align'
	#align_seqs.py -i ${newfile}_otus.fasta -o rep_set_align -t ${RDP_DB_16S}/trainset16_022016.rdp.align.fasta

	#Filter alignments on QIIME
	echo "Running the Qiime Container - filter_alignment.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	filter_alignment.py -i rep_set_align/'${newfile}'_otus_aligned.fasta -o filtered_alignment; \
	chmod -R 777 filtered_alignment'
	#filter_alignment.py -i rep_set_align/${newfile}_otus_aligned.fasta -o filtered_alignment

	#Make the reference tree on QIIME
	echo "Running the Qiime Container - make_phylogeny.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	make_phylogeny.py -i filtered_alignment/'${newfile}'_otus_aligned_pfiltered.fasta -o rep_set.tre; \
	chmod -R 777 rep_set.tre'
	#make_phylogeny.py -i filtered_alignment/${newfile}_otus_aligned_pfiltered.fasta -o rep_set.tre

	mkdir diversity_by_sample
	cd diversity_by_sample

	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiimeecho "Running the QiimePipe Container - createAbundanceFile.py: "
	echo "Running the QiimePipe Container - createAbundanceFile.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample; \
	python3.6 /qiimepipe/createAbundanceFile.py ../output/'${newfile}'_tax_assignments.txt ../'${TAXTABLE}';\
	chmod -R 777 ../diversity_by_sample'
	#python ${SCRIPT_PATH}/createAbundanceFile.py ../output/${newfile}_tax_assignments.txt ../${newfile}_otu_table.txt
	cd ..

	docker stop qiime_run_$TIMESTAMP
	docker rm qiime_run_$TIMESTAMP

elif [ $GENE = "16S-NCBI" ];
then
	docker stop qiimepipe_run_$TIMESTAMP
	docker rm qiimepipe_run_$TIMESTAMP

	#Assign taxonomy to OTUS using blast. The blast database is needed.

	# echo "Creating a BLAST Container: "
	# docker run -id -v $CURRENT_PATH:/output/ -v $BLASTDB:/blastdb/ --name blast_run_$TIMESTAMP itvdsbioinfo/pimba_blast:latest
	echo "EXPORTBLAST=$EXPORTBLAST"
	echo "Running the BLAST Container - blastn: "
	echo "RAWDATA=$RAWDATA"
	#-u $(id -u)
	docker exec -u $(id -u) -i blast_run_$TIMESTAMP /bin/bash -c "cd /output/$OUTPUT; $EXPORTBLAST\
		blastn -query $RAWDATA -task megablast -db $PARAM_BLAST -perc_identity \
		$SIMILARITY_INT_ASG -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT \
		-max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE \
		-outfmt \"6 qseqid sscinames sseqid staxids stitle pident qcovs evalue\" $PARAM_THREADS > \
		'${newfile}'_blast.log; chmod 777 '${newfile}'_blast.log"
	#blastn -query ${newfile}_otus.fasta -task megablast -db $NCBI_DB -perc_identity $SIMILARITY_INT -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT -max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE -parse_deflines -num_threads $THREADS -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > ${newfile}_blast.log

	mkdir diversity_by_sample
	cd diversity_by_sample

	#Convert UC to otu-table.txt <<< BMP SCRIPT>>>
	echo "Creating a QiimePipe Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $TAXDUMP:/taxdump/ --name qiimepipe_run_$TIMESTAMP itvdsbioinfo/pimba_qiimepipe:v2
	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiime
	echo "Running the QiimePipe Container - createTaxonTable_singleFile.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample; \
	python3.6 /qiimepipe/createTaxonTable_singleFile.py ../'${newfile}'_blast.log \
	../'${TAXTABLE}' /taxdump/rankedlineage.dmp; \
	chmod -R 777 ../diversity_by_sample; chmod 777 ../'${newfile}'_tax_assignments.txt'
	#python ${SCRIPT_PATH}/createTaxonTable_singleFile.py ../${newfile}_blast.log ../${newfile}_otu_table.txt
	cd ../

	mkdir output
	chmod -R 777 output
	mv ${newfile}_tax_assignments.txt output/
	

	docker stop blast_run_$TIMESTAMP
	docker rm blast_run_$TIMESTAMP

elif [ $GENE = "ITS-FUNGI-NCBI" ];
then

	docker stop qiimepipe_run_$TIMESTAMP
	docker rm qiimepipe_run_$TIMESTAMP

	#Filter out contaminants basing on ncbi/nt blast

	echo "Running the BLAST Container - blastn: "
	docker exec -u $(id -u) -i blast_run_$TIMESTAMP /bin/bash -c "cd /output/$OUTPUT; $EXPORTBLAST\
		blastn -query $RAWDATA -task megablast -db $PARAM_BLAST -perc_identity \
		$SIMILARITY_INT_ASG -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT \
		-max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE \
		-outfmt \"6 qseqid sscinames sseqid staxids stitle pident qcovs evalue\" $PARAM_THREADS > \
		'${newfile}'_blast_ncbi.log; chmod 777 '${newfile}'_blast_ncbi.log"
	#blastn -query ${newfile}_otus.fasta -task megablast -db  $NCBI_DB -perc_identity $SIMILARITY_INT -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT -max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE -parse_deflines -num_threads $THREADS -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > ${newfile}_blast_ncbi.log

	mkdir diversity_by_sample_ncbi
	cd diversity_by_sample_ncbi

	echo "Creating a QiimePipe Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $TAXDUMP:/taxdump/ --name qiimepipe_run_$TIMESTAMP itvdsbioinfo/pimba_qiimepipe:v2

	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiime
	echo "Running the QiimePipe Container - create_otuTaxAssignment.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample_ncbi; \
	python3.6 /qiimepipe/create_otuTaxAssignment.py ../'${newfile}'_blast_ncbi.log \
	../'${TAXTABLE}' ncbi_otus_tax_assignments.txt; chmod -R 777 ../diversity_by_sample_ncbi'
	#python ${SCRIPT_PATH}/create_otuTaxAssignment.py ../${newfile}_blast_ncbi.log ../${newfile}_otu_table.txt ncbi_otus_tax_assignments.txt

	cd ..

	#Obtaining only the OTUs classified as Fungi
	echo "Running the QiimePipe Container - filterOTUs.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	python3.6 /qiimepipe/filterOTUs.py '$RAWDATA' diversity_by_sample_ncbi/ncbi_otus_tax_assignments.txt; \
	chmod 777 k__Fungi.* *filtered.fasta *contigIDs.txt'
	#python ${SCRIPT_PATH}/filterOTUs.py ${newfile}_otus.fasta diversity_by_sample_ncbi/ncbi_otus_tax_assignments.txt

	#Assign taxonomy to OTUS using blast. The blast database is needed.
	echo "Running the BLAST Container - blastn: "
	docker exec -i blast_run_$TIMESTAMP /bin/bash -c "cd /output/$OUTPUT; $EXPORTBLAST\
		blastn -query k__Fungi.fasta -task megablast -db $PARAM_BLAST -perc_identity \
		$SIMILARITY_INT_ASG -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT \
		-max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE \
		-outfmt \"6 qseqid sscinames sseqid staxids stitle pident qcovs evalue\" $PARAM_THREADS > \
		'${newfile}'_blast_fungi.log; chmod 777 '${newfile}'_blast_fungi.log"
	#blastn -query k__Fungi.fasta -task megablast -db $NCBI_DB -perc_identity $SIMILARITY_INT -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT -max_target_seqs $HITS_PER_SUBJECT -parse_deflines -evalue $EVALUE -num_threads $THREADS -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > ${newfile}_blast_fungi.log

	#Map reads back to OTU database <<<VSEARCH script>>>
	echo "Running the VSEARCH Container - --usearch_global: "
	docker exec -u $(id -u) -i vsearch_run_$TIMESTAMP  /bin/bash -c 'cd /output/'$OUTPUT'; \
	vsearch --usearch_global ../'${RAWDATA}' --db k__Fungi.fasta --strand both \
	--id '$SIMILARITY' --uc '${newfile}'_map_fungi.uc; \
	chmod 777 '${newfile}'_map_fungi.uc;'
	#vsearch  --usearch_global ../${RAWDATA} --db k__Fungi.fasta --strand both --id $SIMILARITY --uc ${newfile}_map_fungi.uc

	#Convert UC to otu-table.txt <<< BMP SCRIPT>>>
	echo "Running the QiimePipe Container - uc2otutab: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	python3.6 /qiimepipe/uc2otutab.py '${newfile}'_map_fungi.uc > '${newfile}'_otu_table_fungi.txt;\
	chmod 777 '${newfile}'_otu_table_fungi.txt'
	#python ${BMP_PATH}/uc2otutab.py ${newfile}_map_fungi.uc > ${newfile}_otu_table_fungi.txt

	mkdir diversity_by_sample
	cd diversity_by_sample
	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiime

	echo "Running the QiimePipe Container - createTaxonTable_singleFile.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample; \
	python3.6 /qiimepipe/createTaxonTable_singleFile.py ../'${newfile}'_blast_fungi.log \
	../'${newfile}'_otu_table_fungi.txt /taxdump/rankedlineage.dmp; \
	chmod -R 777 ../diversity_by_sample; chmod 777 ../'${newfile}'_tax_assignments.txt'
	#python ${SCRIPT_PATH}/createTaxonTable_singleFile.py ../${newfile}_blast_fungi.log ../${newfile}_otu_table_fungi.txt
	cd ../
	mkdir output
	chmod -R 777 output
	mv ${newfile}_tax_assignments.txt output/

	docker stop blast_run_$TIMESTAMP
	docker rm blast_run_$TIMESTAMP


elif [ $GENE = "ITS-FUNGI-UNITE" ];
then

	echo "Creating a Qiime Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $ITS_UNITE_DB:/database/ --name qiime_run_$TIMESTAMP itvdsbioinfo/pimba_qiime:latest

	#Assign taxonomy to OTUS using blast method on QIIME. Use the file .otus.fa. from UPARSE as input file and UNITE as reference database (Download UNITE database HERE)
	echo "Running the Qiime Container - assign_taxonomy.py: "
	docker exec -u $(id -u) -i qiime_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	assign_taxonomy.py -i '$RAWDATA' -o output -t /database/*'${SIMILARITY_INT}'*.txt \
	-r /database/*'${SIMILARITY_INT}'*.fasta --similarity='$SIMILARITY_ASSIGN'; \
	chmod -R 777 output'
	#assign_taxonomy.py -i ${newfile}_otus.fasta -o output -r ${ITS_UNITE_DB}/sh_refs_qiime_ver8_${SIMILARITY_INT}_s_04.02.2020.fasta -t  ${ITS_UNITE_DB}/sh_taxonomy_qiime_ver8_${SIMILARITY_INT}_s_04.02.2020.txt --similarity=$SIMILARITY_ASSIGN
	mkdir diversity_by_sample
	cd diversity_by_sample

	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiimeecho "Running the QiimePipe Container - createAbundanceFile.py: "
	echo "Running the QiimePipe Container - createAbundanceFile.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample; \
	python3.6 /qiimepipe/createAbundanceFile.py ../output/'${newfile}'_tax_assignments.txt ../'${TAXTABLE}';\
	chmod -R 777 ../diversity_by_sample'
	#python ${SCRIPT_PATH}/createAbundanceFile.py ../output/${newfile}_tax_assignments.txt ../${newfile}_otu_table.txt
	cd ..

	docker stop qiime_run_$TIMESTAMP
	docker rm qiime_run_$TIMESTAMP

elif [ $GENE = "ITS-PLANTS-NCBI" ];
then

	docker stop qiimepipe_run_$TIMESTAMP
	docker rm qiimepipe_run_$TIMESTAMP
	

	#Filter out contaminants basing on ncbi/nt blast

	echo "Running the BLAST Container - blastn: "
	docker exec -u $(id -u) -i blast_run_$TIMESTAMP /bin/bash -c "cd /output/$OUTPUT; $EXPORTBLAST\
		blastn -query $RAWDATA -task megablast -db $PARAM_BLAST -perc_identity \
		$SIMILARITY_INT_ASG -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT \
		-max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE \
		-outfmt \"6 qseqid sscinames sseqid staxids stitle pident qcovs evalue\" $PARAM_THREADS > \
		'${newfile}'_blast_ncbi.log; chmod 777 '${newfile}'_blast_ncbi.log"
	#blastn -query ${newfile}_otus.fasta -task megablast -db  $NCBI_DB -perc_identity $SIMILARITY_INT -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT -max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE -parse_deflines -num_threads $THREADS -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > ${newfile}_blast_ncbi.log

	mkdir diversity_by_sample_ncbi
	cd diversity_by_sample_ncbi

	echo "Creating a QiimePipe Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $TAXDUMP:/taxdump/ --name qiimepipe_run_$TIMESTAMP itvdsbioinfo/pimba_qiimepipe:v2

	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiime
	echo "Running the QiimePipe Container - create_otuTaxAssignment.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample_ncbi; \
	python3.6 /qiimepipe/create_otuTaxAssignment.py ../'${newfile}'_blast_ncbi.log \
	../'${TAXTABLE}' ncbi_otus_tax_assignments.txt; chmod -R 777 ../diversity_by_sample_ncbi'
	#python ${SCRIPT_PATH}/create_otuTaxAssignment.py ../${newfile}_blast_ncbi.log ../${newfile}_otu_table.txt ncbi_otus_tax_assignments.txt

	cd ..

	#Filterin out the OTUs classified as Fungi
	echo "Running the QiimePipe Container - filterOTUs.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	python3.6 /qiimepipe/filterOTUs.py '$RAWDATA' diversity_by_sample_ncbi/ncbi_otus_tax_assignments.txt; \
	chmod 777 k__Fungi.* '${newfile}'_otus_filtered.fasta'
	#python ${SCRIPT_PATH}/filterOTUs.py ${newfile}_otus.fasta diversity_by_sample_ncbi/ncbi_otus_tax_assignments.txt

	#Assign taxonomy to OTUS using blast. The blast database is needed.
	#blastn -query ${newfile}_otus_filtered.fasta -task megablast -db $DATABASE -max_target_seqs 1 -parse_deflines -num_threads 24 -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > ${newfile}_blast.log
	#blastn -query ${newfile}_otus.fasta -task megablast -db $DATABASE -max_target_seqs 1 -parse_deflines -num_threads 24 -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > ${newfile}_blast.log
	echo "Running the BLAST Container - blastn: "
	docker exec -u $(id -u) -i blast_run_$TIMESTAMP /bin/bash -c "cd /output/$OUTPUT; $EXPORTBLAST\
		blastn -query '${newfile}'_otus_filtered.fasta -task megablast -db $PARAM_BLAST -perc_identity \
		$SIMILARITY_INT_ASG -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT \
		-max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE \
		-outfmt \"6 qseqid sscinames sseqid staxids stitle pident qcovs evalue\" $PARAM_THREADS > \
		'${newfile}'_blast_plants.log; chmod 777 '${newfile}'_blast_plants.log"
	#blastn -query ${newfile}_otus.fasta -task megablast -db $NCBI_DB -perc_identity $SIMILARITY_INT -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT -max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE -parse_deflines -num_threads $THREADS -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > ${newfile}_blast.log

	#Map reads back to OTU database <<<VSEARCH script>>>
	echo "Running the VSEARCH Container - --usearch_global: "
	docker exec -u $(id -u) -i vsearch_run_$TIMESTAMP  /bin/bash -c 'cd /output/'$OUTPUT'; \
	vsearch --usearch_global ../'${RAWDATA}' --db '${newfile}'_otus_filtered.fasta --strand both \
	--id '$SIMILARITY' --uc '${newfile}'_map_plants.uc; \
	chmod 777 '${newfile}'_map_plants.uc;'
	#vsearch  --usearch_global ../${RAWDATA} --db k__Fungi.fasta --strand both --id $SIMILARITY --uc ${newfile}_map_fungi.uc

	#Convert UC to otu-table.txt <<< BMP SCRIPT>>>
	echo "Running the QiimePipe Container - uc2otutab: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
	python3.6 /qiimepipe/uc2otutab.py '${newfile}'_map_plants.uc > '${newfile}'_otu_table_plants.txt;\
	chmod 777 '${newfile}'_otu_table_plants.txt'
	#python ${BMP_PATH}/uc2otutab.py ${newfile}_map_fungi.uc > ${newfile}_otu_table_fungi.txt


	mkdir diversity_by_sample
	cd diversity_by_sample

	echo "Running the QiimePipe Container - createTaxonTable_singleFile.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample; \
	python3.6 /qiimepipe/createTaxonTable_singleFile.py ../'${newfile}'_blast_plants.log \
	../'${newfile}'_otu_table_plants.txt /taxdump/rankedlineage.dmp; \
	chmod -R 777 ../diversity_by_sample; chmod 777 ../'${newfile}'_tax_assignments.txt'
	#python ${SCRIPT_PATH}/createTaxonTable_singleFile.py ../${newfile}_blast.log ../${newfile}_otu_table.txt
	cd ../
	mkdir output
	
	mv *_otus_tax_assignments.txt ${newfile}_tax_assignments.txt
	mv ${newfile}_tax_assignments.txt output/

	chmod -R 777 output

	cp ${newfile}_otus_filtered.fasta ${newfile}_otus_plants.fasta

	docker stop blast_run_$TIMESTAMP
	docker rm blast_run_$TIMESTAMP
	

elif [ $GENE = "COI-NCBI" ];
then

	docker stop qiimepipe_run_$TIMESTAMP
	docker rm qiimepipe_run_$TIMESTAMP

	#Assign taxonomy to OTUS using blast. The blast database is needed.

	echo "Running the BLAST Container - blastn: "
	docker exec -u $(id -u) -i blast_run_$TIMESTAMP /bin/bash -c "cd /output/$OUTPUT; $EXPORTBLAST\
		blastn -query $RAWDATA -task megablast -db $PARAM_BLAST -perc_identity \
		$SIMILARITY_INT_ASG -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT \
		-max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE \
		-outfmt \"6 qseqid sscinames sseqid staxids stitle pident qcovs evalue\" $PARAM_THREADS > \
		'${newfile}'_blast.log; chmod 777 '${newfile}'_blast.log"
	#blastn -query ${newfile}_otus.fasta -task megablast -db $NCBI_DB -perc_identity $SIMILARITY_INT -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT -max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE -parse_deflines -num_threads $THREADS -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > ${newfile}_blast.log
	
	mkdir diversity_by_sample
	cd diversity_by_sample

	echo "Creating a QiimePipe Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $TAXDUMP:/taxdump/ --name qiimepipe_run_$TIMESTAMP itvdsbioinfo/pimba_qiimepipe:v2

	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiime
	echo "Running the QiimePipe Container - createTaxonTable_singleFile.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample; \
	python3.6 /qiimepipe/createTaxonTable_singleFile.py ../'${newfile}'_blast.log \
	../'${TAXTABLE}' /taxdump/rankedlineage.dmp; \
	chmod -R 777 ../diversity_by_sample; chmod 777 ../'${newfile}'_tax_assignments.txt'
	#python ${SCRIPT_PATH}/createTaxonTable_singleFile.py ../${newfile}_blast.log ../${newfile}_otu_table.txt
	cd ../
	mkdir output
	chmod -R 777 output
	mv ${newfile}_tax_assignments.txt output/

	docker stop blast_run_$TIMESTAMP
	docker rm blast_run_$TIMESTAMP
	

elif [ $GENE = "COI-BOLD" ];
then

	docker stop qiimepipe_run_$TIMESTAMP
	docker rm qiimepipe_run_$TIMESTAMP
	


	echo "Creating a BLAST Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $COI_BOLD_DB:/coibold/ --name blast_run_$TIMESTAMP itvdsbioinfo/pimba_blast:latest


	#Assign taxonomy to OTUS using blast. The blast database is needed.
	#assign_taxonomy.py -i ${newfile}_otus.fasta -o output -t ${COI_BOLD_DB}_tax.txt -r ${COI_BOLD_DB}.fasta --similarity=$SIMILARITY_ASSIGN
	echo "Running the BLAST Container - blastn: "
	docker exec -u $(id -u) -i blast_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
		blastn -query '$RAWDATA' -task megablast -db /coibold/*.fasta -perc_identity \
		'$SIMILARITY_INT_ASG' -qcov_hsp_perc '$COVERAGE_INT' -max_hsps '$HITS_PER_SUBJECT' \
		-max_target_seqs '$HITS_PER_SUBJECT' -evalue '$EVALUE' -parse_deflines -num_threads \
		'$THREADS' -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > \
		'${newfile}'_blast.log; chmod 777 '${newfile}'_blast.log'
	#blastn -query ${newfile}_otus.fasta -task megablast -db ${COI_BOLD_DB}.fasta -perc_identity $SIMILARITY_INT -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT -max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE -parse_deflines -num_threads $THREADS -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > ${newfile}_blast.log

	mkdir diversity_by_sample
	cd diversity_by_sample

	echo "Creating a QiimePipe Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $COI_BOLD_DB:/coibold/ --name qiimepipe_run_$TIMESTAMP itvdsbioinfo/pimba_qiimepipe:v2

	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiime
	echo "Running the QiimePipe Container - createTaxonTable_singleFile_flex.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample; \
	python3.6 /qiimepipe/createTaxonTable_singleFile_flex.py ../'${newfile}'_blast.log \
	../'${TAXTABLE}' /coibold/*_tax.txt; \
	chmod -R 777 ../diversity_by_sample; chmod 777 ../'${newfile}'_tax_assignments.txt ../'${newfile}'_taxon_red_flagged.txt'
	#python ${SCRIPT_PATH}/createTaxonTable_singleFile_flex.py ../${newfile}_blast.log ../${newfile}_otu_table.txt ${COI_BOLD_DB}_tax.txt
	cd ..
	mkdir output
	chmod -R 777 output
	mv ${newfile}_tax_assignments.txt ${newfile}_taxon_red_flagged.txt output/
	
	
	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiime
	#python ${SCRIPT_PATH}/createTaxonTable_singleFile_flex.py ${newfile}_blast.log ${newfile}_otu_table.txt
	#mkdir output
	#mv ${newfile}_tax_assignments.txt output/	

	docker stop blast_run_$TIMESTAMP
	docker rm blast_run_$TIMESTAMP
	
elif [ $GENE = "ALL-NCBI" ];
then
	docker stop qiimepipe_run_$TIMESTAMP
	docker rm qiimepipe_run_$TIMESTAMP

	#Assign taxonomy to OTUS using blast. The blast database is needed.

	echo "Running the BLAST Container - blastn: "
	docker exec -u $(id -u) -i blast_run_$TIMESTAMP /bin/bash -c "cd /output/$OUTPUT; $EXPORTBLAST\
		blastn -query $RAWDATA -task megablast -db $PARAM_BLAST -perc_identity \
		$SIMILARITY_INT_ASG -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT \
		-max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE \
		-outfmt \"6 qseqid sscinames sseqid staxids stitle pident qcovs evalue\" $PARAM_THREADS > \
		'${newfile}'_blast.log; chmod 777 '${newfile}'_blast.log"
	#blastn -query ${newfile}_otus.fasta -task megablast -db $NCBI_DB -perc_identity $SIMILARITY_INT -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT -max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE -parse_deflines -num_threads $THREADS -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > ${newfile}_blast.log

	mkdir diversity_by_sample
	cd diversity_by_sample

	#Convert UC to otu-table.txt <<< BMP SCRIPT>>>
	echo "Creating a QiimePipe Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $TAXDUMP:/taxdump/ --name qiimepipe_run_$TIMESTAMP itvdsbioinfo/pimba_qiimepipe:v2
	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiime
	echo "Running the QiimePipe Container - createTaxonTable_singleFile.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample; \
	python3.6 /qiimepipe/createTaxonTable_singleFile.py ../'${newfile}'_blast.log \
	../'${TAXTABLE}' /taxdump/rankedlineage.dmp; \
	chmod -R 777 ../diversity_by_sample; chmod 777 ../'${newfile}'_tax_assignments.txt'
	#python ${SCRIPT_PATH}/createTaxonTable_singleFile.py ../${newfile}_blast.log ../${newfile}_otu_table.txt
	cd ../

	mkdir output
	chmod -R 777 output
	mv ${newfile}_tax_assignments.txt output/
	

	docker stop blast_run_$TIMESTAMP
	docker rm blast_run_$TIMESTAMP
else
	docker stop qiimepipe_run_$TIMESTAMP
	docker rm qiimepipe_run_$TIMESTAMP
	

	# echo "Creating a BLAST Container: "
	# docker run -id -v $CURRENT_PATH:/output/ -v $GENE:/gene/ --name blast_run_$TIMESTAMP itvdsbioinfo/pimba_blast:latest

	echo "Running the BLAST Container - blastn: "
	docker exec -u $(id -u) -i blast_run_$TIMESTAMP /bin/bash -c "cd /output/$OUTPUT; $EXPORTBLAST\
		blastn -query $RAWDATA -task megablast -db $PARAM_BLAST -perc_identity \
		$SIMILARITY_INT_ASG -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT \
		-max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE \
		-outfmt \"6 qseqid sscinames sseqid staxids stitle pident qcovs evalue\" $PARAM_THREADS > \
		'${newfile}'_blast.log; chmod 777 '${newfile}'_blast.log"
	#blastn -query ${newfile}_otus.fasta -task megablast -db $GENE -perc_identity $SIMILARITY_INT -qcov_hsp_perc $COVERAGE_INT -max_hsps $HITS_PER_SUBJECT -max_target_seqs $HITS_PER_SUBJECT -evalue $EVALUE -parse_deflines -num_threads $THREADS -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > ${newfile}_blast.log

	mkdir diversity_by_sample
	cd diversity_by_sample
	#Generate individual diversity information for each sample in the data and convert the blast file to otu_tax_assignment file from Qiime
	echo "Creating a QiimePipe Container: "
	docker run -id -v $CURRENT_PATH:/output/ -v $GENE:/gene/ --name qiimepipe_run_$TIMESTAMP itvdsbioinfo/pimba_qiimepipe:v2

	echo "Running the QiimePipe Container - createTaxonTable_singleFile_flex.py: "
	docker exec -u $(id -u) -i qiimepipe_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'/diversity_by_sample; \
	python3.6 /qiimepipe/createTaxonTable_singleFile_flex.py ../'${newfile}'_blast.log \
	../'${TAXTABLE}' /gene/*_tax.txt; \
	chmod -R 777 ../diversity_by_sample; chmod 777 ../'${newfile}'_tax_assignments.txt ../'${newfile}'_taxon_red_flagged.txt'

	#python ${SCRIPT_PATH}/createTaxonTable_singleFile_flex.py ../${newfile}_blast.log ../${newfile}_otu_table.txt <colocar tax>
	cd ../
	mkdir output
	chmod -R 777 output
	mv ${newfile}_tax_assignments.txt ${newfile}_taxon_red_flagged.txt output/

	docker stop blast_run_$TIMESTAMP
	docker rm blast_run_$TIMESTAMP
fi

#Convert otu_table.txt to otu-table.biom, used by QIIME <<< BIOM SCRIPT>>>
echo "Creating a Biom-Format Container: "
docker run -id -v $CURRENT_PATH:/output/ --name biomformat_run_$TIMESTAMP itvdsbioinfo/pimba_biom:v2.1.10

echo "Running the Biom-Format Container - convert: "
docker exec -u $(id -u) -i biomformat_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
biom convert -i '${TAXTABLE}' -o '${newfile}'_table.biom \
--table-type="OTU table" --to-json; chmod 777 '${newfile}'_table.biom;'
#biom convert -i ${newfile}_otu_table.txt -o ${newfile}_table.biom --table-type="OTU table" --to-json

#Add metadata (taxonomy) to OTU table
echo "Running the Biom-Format Container - add-metadata: "
docker exec -u $(id -u) -i biomformat_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
biom add-metadata -i '${newfile}'_table.biom -o '${newfile}'_table_tax.biom \
--observation-metadata-fp output/'${newfile}'_tax_assignments.txt --observation-header OTUID,taxonomy,confidence \
--sc-separated taxonomy --float-fields confidence; chmod 777 '${newfile}'_table_tax.biom;'
#biom add-metadata -i ${newfile}_table.biom -o ${newfile}_table_tax.biom --observation-metadata-fp output/${newfile}_tax_assignments.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy --float-fields confidence

# Check OTU Table  on QIIME.
echo "Running the Biom-Format Container - summarize-table: "
docker exec -u $(id -u) -i biomformat_run_$TIMESTAMP /bin/bash -c 'cd /output/'$OUTPUT'; \
biom summarize-table -i '${newfile}'_table_tax.biom -o '${newfile}'_biom_table; \
chmod 777 '${newfile}'_biom_table;'
#biom summarize-table -i ${newfile}_table_tax.biom -o ${newfile}_biom_table

#conda deactivate

echo "Stopping Containeres: "
docker stop vsearch_run_$TIMESTAMP
docker stop perl_run_$TIMESTAMP
docker stop qiimepipe_run_$TIMESTAMP
docker stop biomformat_run_$TIMESTAMP


echo "Removing Containeres: "
docker rm vsearch_run_$TIMESTAMP
docker rm perl_run_$TIMESTAMP
docker rm qiimepipe_run_$TIMESTAMP
docker rm biomformat_run_$TIMESTAMP

#Get the OTU amount from the sample with the minimum value
#min="$(grep "Min: " ${newfile}_biom_table)"
#arrmin=(${min//;/ })
#min="$(echo "${arrmin[1]}")"
#min=${min%.*}

#Run diversity analyses on QIIME by applying non-phylogenetic metrics
#core_diversity_analyses.py -i ${newfile}_table_tax.biom -m ../mapping_file.txt -e $min -o core_output --nonphylogenetic_diversity