args<-commandArgs(TRUE)
#Rscript file.R <otu_table> <match_list> <output_otu_table>
library(lulu)

otu_file = args[1]
matchlist_file = args[2]
output_otutable_file = args[3]

### onde está -- diretorio/otu_table.txt -- colocar o diretório e nome do arquivo com a tabela de OTUs
otutab <- read.csv(otu_file ,sep='\t',header=TRUE,as.is=TRUE, row.names = 1)


### onde está -- match_list.txt -- colocar o diretório e nome do arquivo gerado no blast da etapa anterior
matchlist <- read.table(matchlist_file , header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)

curated_result <- lulu(otutab, matchlist)
curated_result
table_otu_new <- curated_result$curated_table

write.table(table_otu_new, output_otutable_file, sep="\t")
