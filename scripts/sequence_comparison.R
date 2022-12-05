# https://rdrr.io/bioc/biomaRt/man/getSequence.html
#https://www.uniprot.org/uniprotkb/P0C7U0/entry
library(biomaRt)
chromosome = 7
start = 1665745
end = 1747946 
id = "ENST00000424383.5" # ENSEMBL
getSequence(chromosome, start, end, id, type, seqType, 
            upstream, downstream, mart, verbose = FALSE)

# OR:
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

hsapiens_seq = getSequence(id = "ELFN1", 
                  type = "hgnc_symbol", 
                  seqType = "3utr",# "gene_exon", 
                  mart = mart)
show(hsapiens_seq)

dataset="mmusculus_gene_ensembl"
dataset = "rnorvegicus_gene_ensembl"
mart <- useMart("ensembl", dataset=dataset)
#mart = useMart(host="useast.ensembl.org", 
#               biomart="ENSEMBL_MART_ENSEMBL", 
#               dataset="mmusculus_gene_ensembl")
chromosomse = 5
start = 139893698
end = 139960477
getSequence(chromosome, start, end, mart=mart, type="peptide")
mmusculus_seq = getSequence(id = "ELFN1", 
                           type = "hgnc_symbol", 
                           seqType = "3utr", 
                           mart = mart)
listFilters(mart)
