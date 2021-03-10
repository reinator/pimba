import sys
import os
#USAGE: python createTaxonTable_singleFile_flex.py blastTable.txt otuTable.txt tax_assignment.txt

def createDicTaxAssignId(taxAssignId):
	dicTaxAssignId = dict()
	for line in taxAssignId:
		line = line.split("\t")
		dicTaxAssignId[line[0]] = line[1]

	return dicTaxAssignId

def is_uncultered(organism):
	if(organism != "uncultured bacterium" and organism!= "uncultured organism" and organism != "unidentified microorganism" and organism != "bacterium" and organism != "organism" and organism.count("uncultured")==0 and organism.count("unidentified")==0):
		return False

	return True

def check_majority(organism_list, identity_list, taxid_list):
	max_vote = 0
	max_i = 0
	for i in range(len(taxid_list)):
		taxid = taxid_list[i]
		vote = taxid_list.count(taxid)
		if(vote>max_vote and not is_uncultered(organism_list[i]) ):
			max_vote = vote
			max_i = i

	return organism_list[max_i], identity_list[max_i], taxid_list[max_i]

def check_uncultered(organism_list, identity_list, taxid_list):
	max_vote = 0
	max_i = 0
	for i in range(len(taxid_list)):
		if(not is_uncultered(organism_list[i])):
			max_i = i
			break

	return organism_list[max_i], identity_list[max_i], taxid_list[max_i]

def createDicBlast(blastTable):
	dicBlast = dict()
	lines = blastTable.readlines()
	i=0

	while(i < len(lines)):
		line = lines[i]
		line = line.split("\t")
		otuid = line[0]
		if(otuid.find("gb|") != -1):
			otuid = otuid[3:-1]
		organism_list = [line[4]]
		identity_list = [line[5]]
		taxid_list = [line[2]]
		if(i+1!=len(lines)):
			while(lines[i+1].split("\t")[0] == otuid):
				i+=1
				line = lines[i]
				line = line.split("\t")
				organism_list.append(line[1])
				identity_list.append(line[5])
				taxid_list.append(line[2])
				if(i+1==len(lines)): break
			
		#dicBlast[otuid] = check_majority(organism_list, identity_list, taxid_list)
		dicBlast[otuid] = check_uncultered(organism_list, identity_list, taxid_list)
		#dicBlast[line[0]]=line[2][0],line[3]
		#dicBlast[line[0]]=line[1],line[5],line[3][0] # OTUId = organism_list, identity_list, TaxID
		i+=1
		if(i+1==len(lines)): break
	#dicBlast["*"] = "unassigned",0,0
	return dicBlast

def convertBlastToQiimeTax(dicBlast, dicTaxAssignId, filename, unassigned_otu):
	qiimeTaxFile = open(filename+"_otus_tax_assignments.txt","w")
	qiimeTaxTxt = ""
	otus = dicBlast.keys()

	for otu in otus:
		taxid = dicBlast[otu][2]
		taxon = dicTaxAssignId[str(taxid)]
		taxon = taxon.split(";")
		species = taxon[6][taxon[6].find(" ")+1:-1]
		
		qiimeTaxTxt+=otu+"\t"+"k__"+taxon[0]+"; p__"+taxon[1]+"; c__"+taxon[2]+"; o__"+taxon[3]+"; f__"+taxon[4]+"; g__"+taxon[5]+"; s__"+species+"\t"+str(float(dicBlast[otu][1])/100)+"\n"
		#qiimeTaxTxt+=otu+"\t"+"k__"+taxon[0]+"; p__"+taxon[1]+"; c__"+taxon[2]+"; o__"+taxon[3]+"; f__"+taxon[4]+"; g__"+taxon[5]+"; s__"+species.replace(" ","")+"\t"+str(float(dicBlast[otu][1])/100)+"\n"

	for unass_otu in unassigned_otu:
		qiimeTaxTxt+=unass_otu+"\tUnassigned\t1.00\n"


	qiimeTaxFile.write(qiimeTaxTxt)

	qiimeTaxFile.close()

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
		if(line[2] not in dicTaxon.keys()):
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

blastTable = open(sys.argv[1], "r")	
dicBlast = createDicBlast(blastTable)

otuTable = open(sys.argv[2], "r")
dicOtu = createDicOtu(otuTable)

taxAssignId = open(sys.argv[3], "r", encoding="utf-8")
dicTaxAssignId = createDicTaxAssignId(taxAssignId)

fim = sys.argv[2].find("_otu_table.txt")
filename = sys.argv[2][:fim]

unassigned_otu = []

sample_ids = dicOtu.keys()

for s in sample_ids:

	tableOutput = open(s+".txt", "w")
	outtxt = ""
	dic = dicOtu[s]

	dic_keys = dic.keys()
	dic_blast_keys = dicBlast.keys()

	dicOtuKeys = dic.keys()

	outtxt = "OTUID\tread_count\ttaxon\tsimilarity%\ttaxonomy\n"
	for k in dicOtuKeys:

		if(k in dic_keys and k in dic_blast_keys):
			if(dic[k]!="0"):
				taxid_blast = dicBlast[k][2]
				taxon_levels = dicTaxAssignId[taxid_blast].strip("\n").split(";")
				level = len(taxon_levels)-1
				while taxon_levels[level] == "":
					level-=1

				outtxt+=k+"\t"+dic[k]+"\t"+taxon_levels[level]+"\t"+str(dicBlast[k][1])+"\t"+dicTaxAssignId[taxid_blast].strip("\n")+"\n"

		elif(k in dic_keys and k not in dic_blast_keys):
			if(dic[k]!="0"):
				outtxt+=k+"\t"+dic[k]+"\tUNASSIGNED\tNA\tUNASSIGNED\n"
			if(k not in unassigned_otu):
				unassigned_otu.append(k)


	tableOutput.write(outtxt)
	tableOutput.close()

blastTable.close()
otuTable.close()

convertBlastToQiimeTax(dicBlast,dicTaxAssignId, filename,unassigned_otu)

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
