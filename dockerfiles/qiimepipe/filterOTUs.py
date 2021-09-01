import sys

def filter(otu_tax_file,filterText):
	lines = otu_tax_file.readlines()
	#filter_list=[]
	filter_list = dict()
	for line in lines:
		line=line.split("\t")
		otu=line[0]
		taxon = line[1].split(";")
		kingdom=taxon[0]
		if(kingdom==filterText):
			#filter_list.append(otu)
			filter_list[otu]=line[1]
	
	return(filter_list)

def writeFasta(otuFile, filter_list,filterText):
	filteredOutTxt=""
	filteredKeepTxt=""
	filteredOutTxt2=""
	filteredOutTxt3=""
	lines=otuFile.readlines()
	num_lines = len(lines)
	i=0
	while(i < num_lines):
		if(lines[i][0]==">" and lines[i][1:-1] in filter_list.keys()):
			filteredOutTxt+=lines[i]
			filteredOutTxt2+=lines[i][:-1]+"\t"+filter_list[lines[i][1:-1]]+"\n"
			filteredOutTxt3+=lines[i]
			i+=1
			while(lines[i][0]!=">"):
				filteredOutTxt+=lines[i]
				if(i+1>=num_lines): 
					i+=1 
					break
				else: 
					i+=1
		elif(lines[i][0]==">" and lines[i][1:-1] not in filter_list.keys()):
			filteredKeepTxt+=lines[i]
			i+=1
			while(lines[i][0]!=">"):
				filteredKeepTxt+=lines[i]
				if(i+1>=num_lines): 
					i+=1
					break
				else: 
					i+=1
	filteredKeepFile = open(sys.argv[1].strip(".fasta")+"s_filtered.fasta","w")
	filteredKeepFile.write(filteredKeepTxt)
	
	filteredOutFile = open(filterText+".fasta","w")
	filteredOutFile.write(filteredOutTxt)

	filteredOutFile2 = open(filterText+".txt","w")
	filteredOutFile2.write(filteredOutTxt2)

	filteredOutFile3 = open(filterText+"contigIDs.txt","w")
	filteredOutFile3.write(filteredOutTxt3)
	
	filteredKeepFile.close()
	filteredOutFile.close()
	filteredOutFile2.close()
	filteredOutFile3.close()


otuFile = open(sys.argv[1],"r")
otu_tax_file = open(sys.argv[2],"r")
filterText= "k__Fungi"
filter_list = filter(otu_tax_file, filterText)
writeFasta(otuFile,filter_list,filterText)


