#!/bin/bash
#Authors: Renato Oliveira. Gisele Nunes, Raíssa Oliveira
#version: 1.8
#Date: 07-02-2023

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


#usage: ./pimba_prepare.sh illumina <rawdata_dir> <output_reads> <num_threads> <adapters.txt> <min_length> <min_phred>
#<rawdata_dir> = path with all the R1 and R2 reads file;
#<output_reads> = name for the output file;
#<num_threads> = number of threads;
#<adapters.txt> = txt list file with all adapters and primers used for sequencing;
#<min_lenght> = The minimum lenght of the read after quality treatment;
#<min_phred> = Minimum PHRED score of a read after quality treatment;


#####################################################################################################################################
#or usage: ./pimba_prepare.sh iontorrent-dualindex  <rawdata.fastq> <barcodes.txt> <barcodes_reverse.txt> <barcodes.fasta> <barcodes_for_dir> <Primer_forward> <Primer_reverse> <num_threads> <output_name> <min_length> <min_phred>
#<rawdata.fastq> = single file with all the reads to demultiplex;
#<barcodes.txt> = barcodes used as index in the 3' of the fragment;
#<barcodes_reverse.txt> = reverse complement of <barcodes.txt>;
#<barcodes.fasta> = fasta file for the <barcodes.txt>;
#<barcodes_for_dir> = path to all the barcodes.fasta and barcodes.txt used in the 5'  of the fragment. Each 3' barcode must have a fasta and txt file with all the associated 5' barcodes;
#<Primer_forward> = sequence of the forward primer;
#<Primer_reverse> = sequence of the reverse primer;
#<num_threads> = number of threads;
#<output_name> = name for the output fastq file;
#<min_lenght> = The minimum lenght of the read after quality treatment;
#<min_phred> = Minimum PHRED score of a read after quality treatment;

#####################################################################################################################################
#or usage: ./pimba_prepare.sh iontorrent-singleindex <rawdata.fastq> <prefix> <barcodes.txt> <barcodes.fasta> <primer> <num_threads> <output_name> <min_length> <min_phred>
#<rawdata.fastq> = single file with all the reads to demultiplex;
#<prefix> = name that will precede the barcodes names;
#<barcodes.txt> = barcodes used as index in the 5' of the fragment;
#<barcodes.fasta> = fasta file for the <barcodes.txt>;
#<primer> = primer sequence;
#<num_threads> = number of threads;
#<output_name> = name for the output fastq file;
#<min_lenght> = The minimum lenght of the read after quality treatment;
#<min_phred> = Minimum PHRED score of a read after quality treatment;

CURRENT_PATH=$(pwd)
TIMESTAMP=$(date +"%Y%m%d%H%M%S")

if [ $1 = "illumina" ];
then
	SEQUENCER=$1
	RAWDATADIR=$2
	OUTPUTNAME=$3
	NUM_THREADS=$4
	ADAPTERS=$5
	MINLENGTH=$6
	MINPHRED=$7

	#DIR_NAME_RAW=$(dirname ${RAWDATADIR}/*)
	cd $RAWDATADIR
	FULL_PATH_RAW=$(pwd)
	cd $CURRENT_PATH

	REALPATH_ADAP=$(realpath $ADAPTERS)
	DIR_NAME_ADAP=$(dirname $REALPATH_ADAP)
	FILENAME_ADAP=$(basename $REALPATH_ADAP)
	cd $CURRENT_PATH

	#pathlist=$(echo $FULL_PATH_RAW; echo $FULL_PATH_ADAP)
	#COMMON_PATH=$(i=2; while [ $i -lt 500 ]; do   path=`echo "$pathlist" | cut -f1-$i -d/ | uniq -d`;   if [ -z "$path" ];   then      echo $prev_path;      break;   else      prev_path=$path;   fi;   i=`expr $i + 1`; done);

	#COMMON_PATH=$({ echo $FULL_PATH_RAW; echo $FULL_PATH_ADAP;} | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D')

	#echo Common Path: $COMMON_PATH
	echo Current Path: $CURRENT_PATH


	mkdir prepare_output
	cd prepare_output
	mkdir R1; mkdir R2;

	cd R1
	cp ${FULL_PATH_RAW}/*_R1* .
	cd ../R2
	cp ${FULL_PATH_RAW}/*_R2* .

	cd ../../
	chmod -R 777 prepare_output

	echo DIR_NAME_ADAP: $DIR_NAME_ADAP

	#Removing adapters and filtering sequences by quality
	echo "Creating an AdapterRemoval Container: "
	docker run -id -v $DIR_NAME_ADAP:/adapter/ -v $CURRENT_PATH:/output/ --name adapter_removal_prepare_$TIMESTAMP itvdsbioinfo/pimba_adapterremoval:v2.2.3

	echo "Running the AdapterRemoval Container: "
	docker exec -u $(id -u) -i adapter_removal_prepare_$TIMESTAMP  /bin/bash -c 'ADAPTERS='$FILENAME_ADAP'; cd /output/prepare_output/;\
	 for i in R1/*.fastq*; do newfile=${i%%_*}; newfile=${newfile##*/}; echo $newfile;\
	 AdapterRemoval --file1 $i --file2 R2/${newfile}* --threads '$NUM_THREADS' --mate-separator " " --adapter-list /adapter/${ADAPTERS} \
	  --trimwindows 10 --minquality '$MINPHRED' --minlength '$MINLENGTH' --qualitymax 64 --basename ${newfile}_good --mm 5; done;\
	  rm -r R1/ R2/; chmod -R 777 /output/prepare_output/;'

	mkdir -p assemblies/pear
	cd assemblies/pear

	echo "Creating a PEAR Container: "
	docker run -id -v $CURRENT_PATH:/output/ --name pear_prepare_$TIMESTAMP itvdsbioinfo/pimba_pear:v0.9.10

	echo "Running the PEAR Container: "
	docker exec -u $(id -u) -i pear_prepare_$TIMESTAMP  /bin/bash -c 'cd /output/assemblies/pear/;\
	 for i in ../../prepare_output/*good.pair1.truncated; do newfile=$(basename $i _good.pair1.truncated);echo $newfile; \
	 mkdir $newfile; pear -j '$NUM_THREADS' -f $i -r ../../prepare_output/${newfile}_good.pair2.truncated -o ${newfile}/${newfile}; done;\
	 for i in */*assembled.fastq; do newfile=$(basename $i .assembled.fastq); echo $newfile;  sed -i "s/ /_/g" $i; done;\
	 chmod -R 777 /output/assemblies/'

    echo "Creating and running a Prinseq Container: "
	for i in */*assembled.fastq; do newfile="$(basename $i .assembled.fastq)"; echo $newfile; docker run -u $(id -u) -i -v $CURRENT_PATH:/output/ itvdsbioinfo/pimba_prinseq:v0.20.4 -fastq /output/assemblies/pear/${i} -out_format 1 -seq_id Seq -out_good /output/assemblies/pear/${newfile}.assembled; done;

	echo "Creating a QiimePipe Container: "
	docker run -id -v $CURRENT_PATH:/output/ --name qiimepipe_prepare_$TIMESTAMP itvdsbioinfo/pimba_qiimepipe:v2

	echo "Running the QiimePipe Container: "
	docker exec -u $(id -u) -i qiimepipe_prepare_$TIMESTAMP /bin/bash -c 'cd /output/assemblies/pear/;\
	for i in *assembled.fasta; do newfile=$(basename $i .assembled.fasta); echo $newfile; python3.6 /qiimepipe/relabelReads-v2.py $i .; done;\
	chmod -R 777 /output/assemblies/'

	cat *relabel.fasta > ${OUTPUTNAME}.fasta

	mv ${OUTPUTNAME}.fasta ../../

	chmod 777 ../../${OUTPUTNAME}.fasta

	for i in *relabel.fasta; do newfile=$(basename $i .relabel.fasta); mv $i ${newfile}_relabel_notSingleton.fasta; done;


	#NOW GENERATING OUTPUT WITH SINGLETON READS
	echo "Creating and running a Prinseq Container: "
	for i in */*assembled.fastq; do newfile="$(basename $i .assembled.fastq)"; echo $newfile; cat ${newfile}/*assembled.fastq ${newfile}/*unassembled.* ../../prepare_output/${newfile}_good.singleton.truncated > ${newfile}/${newfile}_withSingleton.fastq ; sed -i 's/ /_/g' ${newfile}/${newfile}_withSingleton.fastq; docker run -u $(id -u) -i -v $CURRENT_PATH:/output/ itvdsbioinfo/pimba_prinseq:v0.20.4 -fastq /output/assemblies/pear/${newfile}/${newfile}_withSingleton.fastq -out_format 1 -seq_id Seq -out_good /output/assemblies/pear/${newfile}.assembled.withSingleton; done;

	echo "Running the QiimePipe Container: "
	docker exec -u $(id -u) -i qiimepipe_prepare_$TIMESTAMP /bin/bash -c 'cd /output/assemblies/pear/;\
	for i in *.assembled.withSingleton.fasta; do newfile=$(basename $i .assembled.withSingleton.fasta); echo $newfile; python3.6 /qiimepipe/relabelReads-v2.py $i .; done;\
	chmod -R 777 /output/assemblies/'

	cat *relabel.fasta > ${OUTPUTNAME}_withSingleton.fasta

	mv ${OUTPUTNAME}_withSingleton.fasta ../../

	chmod 777 ../../${OUTPUTNAME}_withSingleton.fasta

	for i in *relabel.fasta; do newfile=$(basename $i .relabel.fasta); mv $i ${newfile}_relabel_withSingleton.fasta; done;

	echo "Stopping Containeres: "
	docker stop adapter_removal_prepare_$TIMESTAMP
	docker stop pear_prepare_$TIMESTAMP
	docker stop qiimepipe_prepare_$TIMESTAMP

	echo "Removing Containeres: "
	docker rm adapter_removal_prepare_$TIMESTAMP
	docker rm pear_prepare_$TIMESTAMP
	docker rm qiimepipe_prepare_$TIMESTAMP

	echo "Done!"

elif [ $1 = "iontorrent-dualindex" ]
then
	SEQUENCER=$1
	RAWDATA=$2
	BARCODES_3END_TXT=$3
	BARCODES_3END_REV=$4
	BARCODES_3END_FASTA=$5
	BARCODES_5END_DIR=$6
	FORWARD_ADAPTER=$7
	REVERSE_ADAPTER=$8
	NUM_THREADS=$9
	OUTPUT_NAME=${10}


	DIR_NAME_RAW=$(dirname $RAWDATA)
	cd $DIR_NAME_RAW
	FULL_PATH_RAW=$(pwd)
	cd $CURRENT_PATH

	DIR_NAME_BAR3END=$(dirname $BARCODES_3END_TXT)
	cd $DIR_NAME_BAR3END
	FULL_PATH_BAR3END=$(pwd)
	cd $CURRENT_PATH

	DIR_NAME_BAR3ENDREV=$(dirname $BARCODES_3END_REV)
	cd $DIR_NAME_BAR3ENDREV
	FULL_PATH_BAR3ENDREV=$(pwd)
	cd $CURRENT_PATH

	DIR_NAME_BAR3ENDFA=$(dirname $BARCODES_3END_FASTA)
	cd $DIR_NAME_BAR3ENDFA
	FULL_PATH_BAR3ENDFA=$(pwd)
	cd $CURRENT_PATH

	DIR_NAME_BAR5ENDDIR=$(dirname $BARCODES_5END_DIR)
	cd $DIR_NAME_BAR5ENDDIR
	FULL_PATH_BAR5ENDDIR=$(pwd)
	cd $CURRENT_PATH

	COMMON_PATH=$({ echo $FULL_PATH_RAW; echo $FULL_PATH_BAR3END; echo $FULL_PATH_BAR3ENDREV; echo $FULL_PATH_BAR3ENDFA; echo $FULL_PATH_BAR5ENDDIR;} | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D')

	echo $FULL_PATH_RAW
	echo $FULL_PATH_BAR3END
	echo $FULL_PATH_BAR3ENDREV
	echo $FULL_PATH_BAR3ENDFA
	echo $FULL_PATH_BAR5ENDDIR

	echo $COMMON_PATH

	mkdir prepare_output
	cd prepare_output

	chmod -R 777 ../prepare_output


	echo "Creating and running a fasxttoolkit Container: "
	cat ../${RAWDATA} | docker run -u $(id -u) -i -v $CURRENT_PATH:/output/ itvdsbioinfo/pimba_fastxtoolkit:v0.0.14 fastx_barcode_splitter.pl --bcfile /output/${BARCODES_3END_TXT} --prefix /output/prepare_output/indexforbol_ --suffix .fastq --bol --exact
	
	cat indexforbol_unmatched.fastq | docker run -u $(id -u) -i -v $CURRENT_PATH:/output/ itvdsbioinfo/pimba_fastxtoolkit:v0.0.14 fastx_barcode_splitter.pl --bcfile /output/${BARCODES_3END_REV} --prefix /output/prepare_output/indexreveol_ --suffix .fastq --eol --exact
	
	for i in indexreveol_R*; do newfile="${i#"indexreveol_"}"; docker run -u $(id -u) -i -v $CURRENT_PATH:/output/ itvdsbioinfo/pimba_fastxtoolkit:v0.0.14 fastx_reverse_complement -i /output/prepare_output/${i} -o /output/prepare_output/indexcreveol_${newfile}; done;

	rm indexreveol_R*

	for i in indexforbol_R*; do newfile="${i#"indexforbol_"}"; cat *_${newfile} > index_${newfile}; done;

	chmod -R 777 ../prepare_output

	echo "Creating a QiimePipe Container: "
	docker run -id -v $CURRENT_PATH:/output/ --name qiimepipe_prepare_$TIMESTAMP itvdsbioinfo/pimba_qiimepipe:v2

	echo "Running the QiimePipe Container: "
	docker exec -u $(id -u) -i qiimepipe_prepare_$TIMESTAMP /bin/bash -c 'cd /output/prepare_output/; \
	for i in index_*; do newfile="${i#"index_"}"; newfile2=$(basename $newfile .fastq); \ 
	python3.6 /qiimepipe/fastq_strip_barcode_relabel2.py $i '$REVERSE_ADAPTER' '/output/${BARCODES_3END_FASTA}' Ex > ${newfile2}_clipped.fastq; done;\
	chmod -R 777 /output/prepare_output/'

	echo "Creating and running a Prinseq Container: "
	for i in *_clipped.fastq; do newfile="$(basename $i _clipped.fastq)"; docker run -u $(id -u) -i -v $CURRENT_PATH:/output/ itvdsbioinfo/pimba_prinseq:v0.20.4 -fastq /output/prepare_output/${i} -min_len 50 -out_good /output/prepare_output/${newfile}_min50; done;


	echo "Creating and running a fasxttoolkit Container: "
	for i in *_min50.fastq; do newfile="$(basename $i _min50.fastq)"; docker run -u $(id -u) -i -v $CURRENT_PATH:/output/ itvdsbioinfo/pimba_fastxtoolkit:v0.0.14 fastx_reverse_complement -i /output/prepare_output/${i} -o /output/prepare_output/${newfile}_revclipped.fastq; done;

	echo "Running the QiimePipe Container: "
	docker exec -u $(id -u) -i qiimepipe_prepare_$TIMESTAMP /bin/bash -c 'cd /output/prepare_output/; \
	for i in *_revclipped.fastq; do newfile=$(basename $i _revclipped.fastq); \ 
	python3.6 /qiimepipe/fastq_strip_barcode_relabel2.py $i '$FORWARD_ADAPTER' '/output/${BARCODES_5END_DIR}'/barcodes_${newfile}.fasta Ex > samples_${newfile}.fastq; done;\
	chmod -R 777 /output/prepare_output/'

	cat samples_R* > ${OUTPUT_NAME}.fastq
	chmod 777 ${OUTPUT_NAME}.fastq

	echo "Creating an AdapterRemoval Container: "
	docker run -id -v $CURRENT_PATH:/output/ --name adapter_removal_prepare_$TIMESTAMP itvdsbioinfo/pimba_adapterremoval:v2.2.3

	echo "Running the AdapterRemoval Container: "
	docker exec -u $(id -u) -i adapter_removal_prepare_$TIMESTAMP  /bin/bash -c 'cd /output/prepare_output/;\
	 AdapterRemoval --file1 '${OUTPUT_NAME}'.fastq --threads '$NUM_THREADS' --trimwindows 10 \
	 --minquality '$MINPHRED' --minlength '$MINLENGTH' --qualitymax 64 --basename '${OUTPUT_NAME}'_good; \
	  chmod -R 777 /output/prepare_output/;'


	echo "Creating and running a Prinseq Container: "
	docker run -i -v $CURRENT_PATH:/output/ itvdsbioinfo/pimba_prinseq:v0.20.4 -fastq /output/prepare_output/${OUTPUT_NAME}_good.truncated -out_format 1 -out_good /output/prepare_output/${OUTPUT_NAME}_good

    rm indexcreveol* indexforbol* index_* *clipped* *min50* 

    docker exec -u $(id -u) -i qiimepipe_prepare_$TIMESTAMP  /bin/bash -c 'chmod -R 777 /output/prepare_output/'
    mv ${OUTPUT_NAME}_good.fasta ../

    echo "Stopping Containeres: "
	docker stop adapter_removal_prepare_$TIMESTAMP
	docker stop qiimepipe_prepare_$TIMESTAMP

	echo "Removing Containeres: "
	docker rm adapter_removal_prepare_$TIMESTAMP
	docker rm qiimepipe_prepare_$TIMESTAMP

	echo "Done!"

elif [ $1 = "iontorrent-singleindex" ]
then

	SEQUENCER=$1
    RAWDATA=$2
	PREFIX=$3
    BARCODES_5END_TXT=$4
    BARCODES_5END_FASTA=$5
    ADAPTER=$6
    NUM_THREADS=$7
    OUTPUT_NAME=$8


    DIR_NAME_RAW=$(dirname $RAWDATA)
	cd $DIR_NAME_RAW
	FULL_PATH_RAW=$(pwd)
	cd $CURRENT_PATH

	DIR_NAME_BAR5END=$(dirname $BARCODES_5END_TXT)
	cd $DIR_NAME_BAR5END
	FULL_PATH_BAR5END=$(pwd)
	cd $CURRENT_PATH

	DIR_NAME_BAR5ENDFA=$(dirname $BARCODES_5END_FASTA)
	cd $DIR_NAME_BAR5ENDFA
	FULL_PATH_BAR5ENDFA=$(pwd)
	cd $CURRENT_PATH

	COMMON_PATH=$({ echo $FULL_PATH_RAW; echo $FULL_PATH_BAR5END; echo $FULL_PATH_BAR5ENDFA;} | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D')

	mkdir prepare_output
	cd prepare_output

	chmod -R 777 ../prepare_output

	echo "Creating and running a fasxttoolkit Container: "
	cat ../${RAWDATA} | docker run -u $(id -u) -i -v $CURRENT_PATH:/output/ itvdsbioinfo/pimba_fastxtoolkit:v0.0.14 fastx_barcode_splitter.pl --bcfile /output/${BARCODES_5END_TXT} --prefix /output/prepare_output/${PREFIX}_ --suffix .fastq --bol --exact > stats.txt

	# cat $RAWDATA | fastx_barcode_splitter.pl --bcfile $BARCODES_5END_TXT --prefix ${PREFIX}_ --suffix .fastq --bol --exact > stats.txt

	mv ${PREFIX}_unmatched.fastq unmatched.fastq

	cat ${PREFIX}_* > ${OUTPUT_NAME}.fastq

	chmod -R 777 ${OUTPUT_NAME}.fastq

	echo "Creating an AdapterRemoval Container: "
	docker run -id -v $CURRENT_PATH:/output/ --name adapter_removal_prepare_$TIMESTAMP itvdsbioinfo/pimba_adapterremoval:v2.2.3

	echo "Running the AdapterRemoval Container: "
	docker exec -u $(id -u) -i adapter_removal_prepare_$TIMESTAMP  /bin/bash -c 'cd /output/prepare_output/;\
	 AdapterRemoval --file1 '${OUTPUT_NAME}.fastq' --threads '$NUM_THREADS' --minlength '$MINLENGTH' \
	 --basename '${OUTPUT_NAME}_min'; chmod -R 777 /output/prepare_output/;'

	# for i in ${PREFIX}_*fastq; do newfile="$(basename $i .fastq)"; echo working with $i; AdapterRemoval --file1 $i --threads 1 --minlength 50 --basename ${newfile}_min; 

	echo "Creating a QiimePipe Container: "
	docker run -id -v $CURRENT_PATH:/output/ --name qiimepipe_prepare_$TIMESTAMP itvdsbioinfo/pimba_qiimepipe:v2

	echo "Running the QiimePipe Container: "
	docker exec -u $(id -u) -i qiimepipe_prepare_$TIMESTAMP /bin/bash -c 'cd /output/prepare_output/; \
	python3.6 /qiimepipe/fastq_strip_barcode_relabel2.py '${OUTPUT_NAME}_min.truncated' '$ADAPTER' '/output/${BARCODES_5END_FASTA}' Seq > '${OUTPUT_NAME}_clipped.fastq';\
	chmod -R 777 /output/prepare_output/'

	# /bio/share_bio/softwares/BMP_UPARSE_Scripts/fastq_strip_barcode_relabel2.py ${newfile}_min.truncated $ADAPTER $BARCODES_5END_FASTA Seq >> ${OUTPUT_NAME}.fastq; done;

	# echo "Creating an AdapterRemoval Container: "
	# docker run -id -v $COMMON_PATH:/common/ -v $CURRENT_PATH:/output/ --name adapter_removal_prepare_$TIMESTAMP itvdsbioinfo/pimba_adapterremoval:v2.2.2

	echo "Running the AdapterRemoval Container: "
	docker exec -u $(id -u) -i adapter_removal_prepare_$TIMESTAMP  /bin/bash -c 'cd /output/prepare_output/;\
	 AdapterRemoval --file1 '${OUTPUT_NAME}_clipped.fastq' --threads '$NUM_THREADS' \
	 --trimwindows 10 --minquality '$MINPHRED' --minlength '$MINLENGTH' --qualitymax 64 \
	 --basename '${OUTPUT_NAME}_good'; chmod -R 777 /output/prepare_output/;'


	# AdapterRemoval --file1 ${OUTPUT_NAME}.fastq --threads $NUM_THREADS --trimwindows 10 --minquality 20 --minlength 50 --qualitymax 64 --basename ${OUTPUT_NAME}_good

	echo "Creating and running a Prinseq Container: "
	docker run -i -v $CURRENT_PATH:/output/ itvdsbioinfo/pimba_prinseq:v0.20.4 -fastq /output/prepare_output/${OUTPUT_NAME}_good.truncated -out_format 1 -out_good /output/prepare_output/${OUTPUT_NAME}_good

    docker exec -u $(id -u) -i qiimepipe_prepare_$TIMESTAMP  /bin/bash -c 'chmod -R 777 /output/prepare_output/'
    mv ${OUTPUT_NAME}_good.fasta ../

	# prinseq-lite.pl -fastq ${OUTPUT_NAME}_good.truncated -out_format 1 -out_good ${OUTPUT_NAME}_good
	#mv ${OUTPUT_NAME}_good.truncated ${OUTPUT_NAME}_good.fastq
	echo "Stopping Containeres: "
	docker stop adapter_removal_prepare_$TIMESTAMP
	docker stop qiimepipe_prepare_$TIMESTAMP

	echo "Removing Containeres: "
	docker rm adapter_removal_prepare_$TIMESTAMP
	docker rm qiimepipe_prepare_$TIMESTAMP

	echo "Done!"
else
	echo "Invalid parameter $1"
fi
