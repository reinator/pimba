import sys
import os
#USAGE: python createAbundanceFile.py taxTable.txt otuTable.txt

def createDicTaxAssign(taxTable):
	dicTaxAssign = dict()
	dicTaxon = dict()
	lines = taxTable.readlines()
	for line in lines:
		line = line.split("\t")
		#line[3]=line[3].split(";")
		#dicTaxAssign[line[0]]=line[2][0],line[3]
		#print(line)
		#dic[] = taxonomia,pident
		#dicTaxAssign[line[0]]=line[1],line[5],line[3][0]
		dicTaxAssign[line[0]]=line[1],line[2]
		key=line[0]
		line = line[1].split(";")
		text =""
		for level in range(0,len(line)):
			text+=line[level]+";"
		text=text[:-1]
		#dicTaxon[key]=line[0]+";"+line[1]+";"+line[2]+";"+line[3]+";"+line[4]+";"+line[5]+";"+line[6]
		dicTaxon[key] = text
		
	#dicTaxAssign["*"] = "unassigned",0,0
	return dicTaxAssign, dicTaxon


def createDicOtu(otuTable):
	dicOtu = dict() #Dictionary of dictionary. The first key is the Sample ID, whose value will be a dictionary with the OTU as keys and the abundance as values
	dicOtuIndex = dict() #The first key is the Sample ID and the value is the position in the header
	header = otuTable.readline()
	header = header.strip("\n")
	header = header.split("\t")
	num_samples = len(header)-1
	for i in range(num_samples):
		dicOtu[header[i+1]] = dict()
		dicOtuIndex[i+1] = header[i+1]

	lines = otuTable.readlines()
	for line in lines:
		line = line.strip("\n")
		line = line.split("\t")

		for i in range(num_samples):
			dic = dicOtu[header[i+1]]			
			dic[line[0]]=line[i+1]
			dicOtu[header[i+1]]=dic
			
		
	return dicOtu

def createFinalOutput(tableOutput):
	dicTaxon = dict()
	dicSimil = dict()
	tableOutput.readline()
	lines = tableOutput.readlines()
	dicI = dict()
	dicTaxId = dict()
	
	for line in lines:
		line = line.split("\t")
		#print(line)
		if(line[2] not in dicTaxon):
			i=1
			dicTaxon[line[2]]=int(line[1])
			if(line[3]=="NA"):
				dicSimil[line[2]]=float(0)
			else:
				dicSimil[line[2]]=float(line[3])
			dicI[line[2]]=i
			dicTaxId[line[2]] = line[4]
		else:
			dicTaxon[line[2]]+=int(line[1])
			if(not line[3]=="NA"):
				dicSimil[line[2]]+=float(line[3])

			dicI[line[2]]+=i

	return dicTaxon,dicSimil,dicI,dicTaxId

def normalize(dicTaxon):
	dicTaxonKeys = dicTaxon.keys()
	total_sum = 0
	total_sum_simil = 0
	for k in dicTaxonKeys:
		total_sum+=dicTaxon[k]
		total_sum_simil+=dicSimil[k]
	dicTaxonNormal = dict()
	for k in dicTaxonKeys:
		dicTaxonNormal[k] = float(dicTaxon[k])/float(total_sum)


	return dicTaxonNormal

#taxdump = open("/bio/share_bio/utils/renato/taxdump/rankedlineage.dmp","r")
#dicTaxdump = createDicTaxdump(taxdump)

taxTable = open(sys.argv[1], "r")	
dicTaxAssign, dicTaxon = createDicTaxAssign(taxTable)

otuTable = open(sys.argv[2], "r")
dicOtu = createDicOtu(otuTable)

fim = sys.argv[2].find("_otu_table.txt")
filename = sys.argv[2][:fim]

unassigned_otu = []

sample_ids = dicOtu.keys()

for s in sample_ids:

	tableOutput = open(s+".txt", "w")
	outtxt = ""
	dic = dicOtu[s]

	dicOtuKeys = dic.keys()
	print(s)
	outtxt = "OTUID\tread_count\ttaxon\tsimilarity%\ttaxonomy\n"
	for k in dicOtuKeys:

		if(k in dic and  k in dicTaxAssign):
			if(dic[k]!="0"):
				taxon = dicTaxAssign[k][0].split(";")
				while (taxon[-1] == ""):
					taxon.pop(-1)
				outtxt+=k+"\t"+dic[k]+"\t"+taxon[-1]+"\t"+str(dicTaxAssign[k][1])+"\t"+dicTaxon[k]+"\n"

		elif(k in dic and k not in dicTaxAssign):
			if(dic[k]!="0"):
				outtxt+=k+"\t"+dic[k]+"\tUNASSIGNED\tNA\tUNASSIGNED\n"
			if(k not in unassigned_otu):
				unassigned_otu.append(k)


	tableOutput.write(outtxt)
	tableOutput.close()

taxTable.close()
otuTable.close()

#convertBlastToQiimeTax(dicTaxAssign,filename,unassigned_otu)

for s in sample_ids:

	finalOutput = open(s+"_final.txt", "w")
	tableOutput = open(s+".txt", "r")
	dicTaxon,dicSimil,dicI,dicTaxId = createFinalOutput(tableOutput)
	dicTaxonNormal = normalize(dicTaxon)
	dicTaxonKeys = dicTaxon.keys()
	finalouttxt = "taxon\tread_count\tabundance%\tsimilarity%\ttaxonomy\n"

	for k in dicTaxonKeys:
		if(dicTaxonNormal[k]*100 >= 1.0):
			finalouttxt+=k+"\t"+str(dicTaxon[k])+"\t"+str(dicTaxonNormal[k]*100)+"\t"+str(dicSimil[k]/dicI[k])+"\t"+str(dicTaxId[k])

	finalouttxt+="\n-------------------Below 1%------------------------\n"
	for k in dicTaxonKeys:
		if(dicTaxonNormal[k]*100 < 1.0):
			finalouttxt+=k+"\t"+str(dicTaxon[k])+"\t"+str(dicTaxonNormal[k]*100)+"\t"+str(dicSimil[k]/dicI[k])+"\t"+str(dicTaxId[k])

	finalOutput.write(finalouttxt)

	tableOutput.close()
	os.remove(s+".txt")
	finalOutput.close()
