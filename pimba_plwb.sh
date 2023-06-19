#!/bin/bash
#Authors: Renato Oliveira
#version: 2.0
#Date: 07/06/2023

###    Copyright (C) 2023  Renato Oliveira
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
###    renato.renison@gmail.com


##usage: ./pimba_plwb.sh -i <fasta_query> -c <constraint_tree> -a <constraint_alignment> -x <taxonomy_file> -t <num_threads> -d <nt_or_aa> -o <outdir>
#-i = Inputa FASTA file with query sequences to be placed;
#-c = Constraint tree in newick format;
#-a = Aligned FASTA file of sequences in contraint tree;
#-x = File separated by tab of sequences in contraint tree and their associated taxon;
#-t = Number of threads to be used;
#-d = Type of sequences to be analyzed 'nt' or 'aa';
#-o = Output directory to store the results;

#SCRIPT_PATH=/bio/share_bio/utils/renato/QiimePipe
#SCRIPT_PATH=/home/pbd001/ITV/PIMBA/lulu_dada_phyloseq/Renato/Renato/Phyloseq_data_script

get_common_path() {
    IFS="/"
    read -ra ADDR1 <<< "$1"
    read -ra ADDR2 <<< "$2"
    common_path=""

    for i in "${!ADDR1[@]}"; do
        if [[ "${ADDR1[$i]}" == "${ADDR2[$i]}" ]]; then
            common_path="${common_path}/${ADDR1[$i]}"
        else
            break
        fi
    done

    echo $common_path
}

while getopts "i:c:a:x:t:d:o:" opt; do
	case $opt in
		i) INPUT="$OPTARG"
		;;
		c) TREE="$OPTARG"
		;;
		a) ALIGNMENT="$OPTARG"
		;;
		x) TAXON="$OPTARG"
		;;
		t) THREADS="$OPTARG"
		;;
		d) TYPE="$OPTARG"
		;;		
		o) OUTDIR="$OPTARG"
		;;
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

CURRENT_PATH=$(pwd)

DIR_NAME_INPUT=$(dirname $INPUT)
cd $DIR_NAME_INPUT
FULL_PATH_INPUT=$(pwd)
cd $CURRENT_PATH

DIR_NAME_TREE=$(dirname $TREE)
cd $DIR_NAME_TREE
FULL_PATH_TREE=$(pwd)
cd $CURRENT_PATH

DIR_NAME_ALG=$(dirname $ALIGNMENT)
cd $DIR_NAME_ALG
FULL_PATH_ALG=$(pwd)
cd $CURRENT_PATH

DIR_NAME_TAX=$(dirname $TAXON)
cd $DIR_NAME_TAX
FULL_PATH_TAX=$(pwd)
cd $CURRENT_PATH

pathlist=$(echo $FULL_PATH_INPUT; echo $FULL_PATH_TREE; echo $FULL_PATH_ALG; echo $FULL_PATH_TAX; echo $CURRENT_PATH)
uniqlist=$(echo $pathlist | tr ' ' '\n' | sort | uniq | tr '\n' ' ' | sed -e 's/[[:space:]]*$//')

COMMON_PATH=${uniqlist[0]}
for path in "${uniqlist[@]:1}"; do
    COMMON_PATH=$(get_common_path $COMMON_PATH $path)
done

COMMON_PATH=$(echo $COMMON_PATH | tr ' ' '/')
COMMON_PATH=$(echo $COMMON_PATH | sed 's/^/\//')


#COMMON_PATH=$(i=2; while [ $i -lt 500 ]; do   path=`echo "$uniqlist" | cut -f1-$i -d/ | uniq -d`;   if [ -z "$path" ];   then      echo $prev_path;      break;   else      prev_path=$path;   fi;   i=`expr $i + 1`; done);

#COMMON_PATH=$({ echo $FULL_PATH_OTU; echo $FULL_PATH_TAX; echo $FULL_PATH_META;} | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D')
echo $pathlist
echo $uniqlist
echo Common Path: $COMMON_PATH
echo Current Path: $CURRENT_PATH

INPUT=$(echo ${FULL_PATH_INPUT#"$COMMON_PATH"})/$(basename $INPUT)
TREE=$(echo ${FULL_PATH_TREE#"$COMMON_PATH"})/$(basename $TREE)
ALIGNMENT=$(echo ${FULL_PATH_ALG#"$COMMON_PATH"})/$(basename $ALIGNMENT)
TAXON=$(echo ${FULL_PATH_TAX#"$COMMON_PATH"})/$(basename $TAXON)


echo "Creating a plwb Container: "
docker run -i --rm -v $COMMON_PATH:/common/ -v $CURRENT_PATH:/current/ --name pimba_plwb itvdsbioinfo/pimba_plwb:v1.0.0 /bin/bash -c 'cd /plwb/place-pipe/; \
	./place.py --fasta-paths /common/'${INPUT}' --reference-tree /common/'${TREE}' \
	--reference-msa /common/'${ALIGNMENT}' --model-file /common/'${TREE}' --taxonomy-file /common/'${TAXON}' --threads '$THREADS' -d '$TYPE' --out-dir /current/'${OUTDIR}';'
