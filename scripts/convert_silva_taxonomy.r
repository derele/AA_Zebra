o## Folling this:
## http://blog.mothur.org/2018/01/10/SILVA-v132-reference-files/

## And this: https://zenodo.org/record/1172783#.Wp7yuD0o9hE

## download the SILVA taxa mapping
## wget
## https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_132_SSURef_tax_silva_trunc.fasta.gz

## wget
## https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_132_LSURef_tax_silva_trunc.fasta.gz

map.in <- read.table("/SAN/db/RDP/Silva_132/tax_slv_ssu_132.txt", header=F,
                     sep="\t",stringsAsFactors=F)
map.in <- map.in[,c(1,3)]
colnames(map.in) <- c("taxlabel","taxlevel")

tax.ranks <- c("domain", "kingdom", "phylum", "class", "order", "family", 
               "genus")

map <- map.in[map.in$taxlevel%in%tax.ranks, ]

map$deepest.ann <- apply(map, 1, function (x) {
    all.path <- strsplit(x["taxlabel"], ";")[[1]]
    all.path[length(all.path)]
})

library(Biostrings)
FAS <- Biostrings::readRNAStringSet("/SAN/db/RDP/Silva_132/SILVA_132_SSURef_tax_silva_trunc.fasta.gz")

Nstring <- strsplit(names(FAS), " ")
nnames <- unlist(lapply(Nstring, "[", 1))
nnames <- strsplit(nnames, "\\.")
nnames <- unlist(lapply(nnames, "[", 1))

namePath <- unlist(lapply(Nstring, "[", 2))

name.map <- data.frame(acc=nnames, path=namePath)
name.map$path <- as.character(name.map$path)

library(parallel)

tax.list <- mclapply(1:nrow(name.map), function (row) {
    r.list <- as.vector(rep(NA, times=length(tax.ranks)))
    names(r.list) <- tax.ranks
    for(taxon in strsplit(name.map[row, "path"], ";")[[1]]) {
        submap <- map[map$deepest.ann%in%taxon, c("taxlevel", "deepest.ann")]
        r.list[submap[[1]]] <- submap[[2]]
    }
    r.list
}, mc.cores=20)

tax.mat <- do.call(rbind, tax.list)
name.map <- cbind(name.map, tax.mat)

new.name <- apply(name.map, 1, function(x) paste(x[c(3, 5:9)], collapse=";"))

names(FAS) <- new.name

## oops something wrong here... these are Fungi!!!
head(grep("Mitochondria", new.name, value=TRUE), n=100)
## is it enough to exclude these?

## why is nematoda not ther
grep("Nematoda", new.name, value=TRUE)

## it should be here
head(grep("Nematoda", map.in[map.in$taxlevel%in%"phylum", "taxlabel"], value=TRUE))

FAS.taxed <- FAS[rowSums(is.na(name.map[, c(3, 5:9)]))<4]

FAS.taxed <- DNAStringSet(FAS.taxed)

## This would be a very first quite "incomplete" database.
## Biostrings::writeXStringSet(FAS.taxed,
## "/SAN/db/RDP/Silva_123/SILVA_123_dada2.fasta", format="fasta")

## now need to run blast and the blast2alltax scirpt

FUZZ.tax <- read.csv("/SAN/Metabarcoding/AA_combi/all_dada_vs_nt.taxtable",
                     as.is=TRUE)
## first see that we have them more then once
FUZZ.tax <- FUZZ.tax[duplicated(FUZZ.tax$subject), ]
## then exclude duplicates
FUZZ.tax <- FUZZ.tax[!duplicated(FUZZ.tax$subject), ]
## then exclude undef at genus level
FUZZ.tax <- FUZZ.tax[!grepl("undef", FUZZ.tax$genus), ]

write.table(FUZZ.tax$subject,
            "/SAN/Metabarcoding/AA_combi/selected_dada_nr_hits.acc",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

## now run blastcmd to get the fasta file

## blastdbcmd -db /SAN/db/blastdb/nt/nt -entry_batch selected_dada_nr_hits.acc -outfmt '%a|%s' > dada_nr_hits.fasta

## plus a little bit of formatting magic
## tr '|' '\n' < dada_nr_hits.fasta > tmp; mv tmp dada_nr_hits.fasta
## awk '{if(NR%2){print ">"$0} else {print}}' dada_nr_hits.fasta > tmp; mv tmp dada_nr_hits.fasta

FUZZ <- Biostrings::readDNAStringSet("/SAN/Metabarcoding/AA_combi/dada_nr_hits.fasta")

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
## efficiency of RDP
FUZZ <- FUZZ[(width(FUZZ)<5000 & width(FUZZ)>1400)]

## properly named
FUZZ <- FUZZ[(nchar(names(FUZZ))>10)]
FUZZ <- FUZZ[!grepl("character\\(0\\)", names(FUZZ))]

allFUZZ <- DNAStringSet(c(FAS.taxed, FUZZ))

## Biostrings::writeXStringSet(FUZZ,
##                            "/SAN/Metabarcoding/AA_combi/dada2_WHexp.fasta",
##                            format="fasta")


## Biostrings::writeXStringSet(allFUZZ,
##                             "/SAN/db/RDP/Silva_123/SILVA_123_dada2_WHexp.fasta",
##                             format="fasta")
