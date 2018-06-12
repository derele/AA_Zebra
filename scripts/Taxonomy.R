## Folling this:
## 
## Trying to use PR2 and Silva plus added own blast

## The home made part of the Taxonomy... 

## now need to run blast and the blast2alltax script

FUZZ.tax <- read.csv("/SAN/Zebra/all_asv_vs_nt.taxtable",  as.is=TRUE)

## first see that we have them more then once
FUZZ.tax <- FUZZ.tax[duplicated(FUZZ.tax$subject), ]
## then exclude duplicates
FUZZ.tax <- FUZZ.tax[!duplicated(FUZZ.tax$subject), ]
## then exclude undef at genus level
FUZZ.tax <- FUZZ.tax[!grepl("undef", FUZZ.tax$genus), ]

write.table(FUZZ.tax$subject,
            "/SAN/Zebra/selected_dada_nr_hits.acc",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

## now run blastcmd to get the fasta file

## blastdbcmd -db /SAN/db/blastdb/nt/nt -entry_batch selected_dada_nr_hits.acc -outfmt '%a|%s' > dada_nr_hits.fasta

## plus a little bit of formatting magic
## tr '|' '\n' < dada_nr_hits.fasta > tmp; mv tmp dada_nr_hits.fasta
## awk '{if(NR%2){print ">"$0} else {print}}' dada_nr_hits.fasta > tmp; mv tmp dada_nr_hits.fasta

FUZZ <- Biostrings::readDNAStringSet("/SAN/Zebra//dada_nr_hits.fasta")

names(FUZZ) <- gsub("\\.\\d+$", "", names(FUZZ))

newname <- sapply(names(FUZZ), function(acc){
    focus <- FUZZ.tax[FUZZ.tax$subject%in%acc,]
    ## names corresponding to 
    ## outlevels <- c("domain","phylum","class","order","family","genus")
    paste(focus[,c("superkingdom", "phylum", "class", 
                   "order", "family", "genus")], collapse=";")
})

names(FUZZ) <- newname

## now get rid of whole genomes and other really long sequences for
## efficiency of RDP also get rid of fragment sequences
FUZZ <- FUZZ[(width(FUZZ)<5000 & width(FUZZ)>1400)]

## properly named
FUZZ <- FUZZ[(nchar(names(FUZZ))>10)]
FUZZ <- FUZZ[!grepl("character\\(0\\)", names(FUZZ))]

allFUZZ <- DNAStringSet(FUZZ)

Biostrings::writeXStringSet(allFUZZ, "/SAN/db/RDP/NR_ALLMERGE.fasta",
                            format="fasta")
