library(taxonomizr)
library(taxize)

## downloaded from ENA (EBI) via "advanced search / makers (limited
## to taxon Eukaryota and daughter taxa)" on March 24 2018
ENAseq28S <- Biostrings::readDNAStringSet("/SAN/db/RDP/28S_ena.fasta")
ENAseq18S <- Biostrings::readDNAStringSet("/SAN/db/RDP/18S_ena.fasta")

## remove sequences below 300 bases
ENAseq28S <- ENAseq28S[width(ENAseq28S)>300]
ENAseq18S <- ENAseq18S[width(ENAseq18S)>300]

accession28S <- gsub("ENA\\|(.*?):.*", "\\1", names(ENAseq28S))
accession18S <- gsub("ENA\\|(.*?):.*", "\\1", names(ENAseq18S))

## FIX ME
## spc <- gsub("ENA\\|.*? (\\w+ \\w+\\.?) .*", "\\1", names(ENAseq28S))

## only to update the database
## getAccession2taxid("/SAN/db/taxonomy/")

taxdb <- "/SAN/db/taxonomy/accessionTaxa.sql"

## only to build the database 
## read.accession2taxid(list.files('/SAN/db/taxonomy','accession2taxid.gz$', full.names=TRUE),
##                      taxdb)

## wow I mean this is fast... amazing
taxaId28S <- accessionToTaxa(accession28S, taxdb)
taxaId18S <- accessionToTaxa(accession18S, taxdb)

## and has such a great coverage. Only two/111 taxa unassigned!
table(is.na(taxaId28S))
table(is.na(taxaId18S))

## only to update the database
## getNamesAndNodes("/SAN/db/taxonomy/")

taxaNodes <- read.nodes("/SAN/db/taxonomy/nodes.dmp")
taxaNames <- read.names("/SAN/db/taxonomy/names.dmp")

taxonomy28S <- getTaxonomy(taxaId28S,taxaNodes,taxaNames)
taxonomy18S <- getTaxonomy(taxaId18S,taxaNodes,taxaNames)

## how many sequences are duplicated
table(duplicated(ENAseq18S))
table(duplicated(ENAseq28S))

## how many taxonomic paths are duplicated
table(duplicated(taxonomy18S))
table(duplicated(taxonomy28S))

## rename sequences
names(ENAseq18S) <- apply(taxonomy18S[, 1:6], 1, paste, collapse=";")
names(ENAseq28S) <- apply(taxonomy28S[, 1:6], 1, paste, collapse=";")

## remove where both sequence and taxonomy path are duplicated this is
## not perfect: could be duplicates of different first occurences...
ENAseq18S <- ENAseq18S[!(duplicated(ENAseq18S)&duplicated(names(ENAseq18S)))]
ENAseq28S <- ENAseq28S[!(duplicated(ENAseq28S)&duplicated(names(ENAseq28S)))]
                       
## write the remaining to a file
Biostrings::writeXStringSet(ENAseq18S, "/SAN/db/RDP/annotated_18S_ena.fasta")
Biostrings::writeXStringSet(ENAseq28S, "/SAN/db/RDP/annotated_28S_ena.fasta")
