#Author: Renato Oliveira
#Creation Date: 15-03-18
#Update: 22-09-2019

#Usage: python relabelReads-v2.py <FASTA file> <Separator in filename>
import sys
import os

arq_name = sys.argv[1]
arq = open(arq_name, "r")

fim = arq_name.find(sys.argv[2])
sample = arq_name[:fim]

#os.system("usearch9.0.2132_i86linux32 -fastq_filter "+arq_name+" -fastq_maxee 0.5 -fastaout "+sample+"_filtered.fasta")

#arq2 = open(sample+"_prinseq.fasta", "r")
new_arq = open(sample+"_relabel.fasta", "w")

for line in arq:
	if(line[0]==">"):
		new_arq.write(line[:-1]+";barcodelabel="+sample+";\n")
	else:
		new_arq.write(line)

arq.close()
new_arq.close()
