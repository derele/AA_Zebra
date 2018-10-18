library(ggplot2)
library(MultiAmplicon)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)



## Filter # only run when new filtered data is needed
FILTER <- FALSE
newMA <- FALSE
newDeDa <- FALSE
newTax <- TRUE

## first pool file paths
files <- list.files(path="/SAN/Zebra/raw_all_fastqs",
                    pattern=".fastq.gz", full.names=TRUE)

## first pool file names
Ffq.file <- files[grepl("R1", files)]
Rfq.file <- files[grepl("R2", files)]

## first pool sample names
samples <- gsub("/SAN/Zebra/raw_all_fastqs/(.*?)_S\\d{,3}_L001_R1_001.fastq\\.gz", "\\1", Ffq.file)

filt_path <- "/SAN/Zebra/DaDaFilt"
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

if (FILTER){
   filter.track <- lapply(seq_along(Ffq.file),  function (i) {
       filterAndTrim(Ffq.file[i], filtFs[i], Rfq.file[i], filtRs[i],
                     truncLen=c(240,235), minLen=c(240,235), 
                     maxN=0, maxEE=4, truncQ=2, 
                     compress=TRUE, verbose=TRUE)
   })
   saveRDS(filter.track, file="/SAN/Zebra/filter.Rds")
} else {
    filter.track <- readRDS(file="/SAN/Zebra/filter.Rds")
}

filter <- do.call(rbind, filter.track)
colSums(filter)[2]/colSums(filter)[1]
## only 69% of data overall cept (decrease length?)

colSums(filter[1:(nrow(filter)/2),])
colSums(filter[1:(nrow(filter)/2),])[2]/colSums(filter[1:(nrow(filter)/2),])[1]
## 70% in first run

colSums(filter[(nrow(filter)/2):nrow(filter),])
colSums(filter[(nrow(filter)/2):nrow(filter),])[2]/colSums(filter[(nrow(filter)/2):nrow(filter),])[1]
## 56% in second


file.ptable <- "/home/ele/Documents/Peter_Zebras/Zebra_Primer_List_Final.csv"

ptable <- read.csv(file.ptable, sep=",", header=TRUE)

## just have a look but still keep the processing general
ptable[ptable$Target.gene%in%c("28S", "18S"), ]

primerF <- gsub(" ", "", ptable[,3])
primerR <- gsub(" ", "", ptable[,4])

names(primerF) <- as.character(ptable[,1])
names(primerR) <- as.character(paste(ptable[,2], ptable[, 5]), sep=";")

files <- PairedReadFileSet(filtFs, filtRs)

primers <- PrimerPairsSet(primerF, primerR)

rownames(ptable) <- names(primers)

MA <- MultiAmplicon(primers, files)

if(newMA) {
    if(dir.exists("./stratified_files/")){
        unlink("./stratified_files/",
               recursive=TRUE)
        }
    MA1 <- sortAmplicons(MA, starting.at=1, max.mismatch=4)
    saveRDS(MA1, file="/SAN/Zebra/MA1.Rds")
} else {
    if(!newMA){
        MA1 <- readRDS(file="/SAN/Zebra/MA1.Rds")
    } else {stop("Whant new sorting or not? Set newMA to TRUE")}
} 


pdf("figures/primers_MA_sorted.pdf", width=46)
plotAmpliconNumbers(MA1)
dev.off()

if(newDeDa){
    MA2 <- derepMulti(MA1, mc.cores=20)
    MA3.1 <- dadaMulti(MA2[, which(grepl("^P1", colnames(MA2)))],
                       Ferr=NULL, Rerr=NULL, selfConsist=TRUE,
                   multithread=FALSE, mc.cores=20, verbose=0, MAX_CONSIST=20)
    MA3.2 <- dadaMulti(MA2[, which(grepl("^P2", colnames(MA2)))],
                       Ferr=NULL, Rerr=NULL, selfConsist=TRUE,
                       multithread=FALSE, mc.cores=20, verbose=0, MAX_CONSIST=20)
    MA3 <- concatenateDadaMulti(MA3.1, MA3.2)    
    MA4.merged <- mergeMulti(MA3, justConcatenate=FALSE, mc.cores=20)
    prop.merged <- calcPropMerged(MA4.merged)
    prop.merged[is.na(prop.merged)] <- 0
    MA4 <- mergeMulti(MA3, justConcatenate=prop.merged<0.8,
                      mc.cores=20)
    ## two amplicons give error further down, they are empty anyways so
    ## lets remove them
    MA4 <- MA4[which(!rownames(MA4)%in%c("Hadz_1200CR.wang1624CR6L 18S",
                                         "LCO1490.HCO2198_5Mod COI")), ]
    ## a more general procedure (unless I kill this bug in the
    ## package) would be:
    ## non.empty <- unlist(lapply(getMergers(MA4), length))>0
    ## MA4 <- MA4[which(non.empty), ]
    MA5 <- makeSequenceTableMulti(MA4, mc.cores=20, orderBy="nsamples")
    MA6 <- removeChimeraMulti(MA5, mc.cores=20)
    saveRDS(MA6, file="/SAN/Zebra/MA6_mix.Rds")
} else {
    MA6 <- readRDS(file="/SAN/Zebra/MA6_mix.Rds")
}

##### TAXONOMY assignment

## Write out the sequences to blast them...

STnoC <- getSequenceTableNoChime(MA6)


if(newTax){
    sequences <- unlist(lapply(STnoC, colnames))
    names(sequences) <- paste0("asv_", 1:length(sequences))
    Biostrings::writeXStringSet(DNAStringSet(unlist(sequences)),
                                "/SAN/Zebra/all_seq_final.fasta")

    ## We blast this file against NR with a gi-list excluding all
    ## uncultured sequences

    ## create the gi-list as a download from an NCBI Entrez Nucleotide
    ## search '"environmental samples"[organism] OR metagenomes[orgn]'

    ## or via command line: esearch -db nucleotide -query '"environmental
    ## samples"[organism] OR metagenomes[orgn]' | efetch -format gi -mode
    ## text > /SAN/db/blastdb/uncultured_gilist.txt

    ## we use this to limit nr blast to usable entries blastn
    ##  -negative_gilist /SAN/db/blastdb/uncultured.gi -query
    ##  all_seq_final.fasta -db /SAN/db/blastdb/nt/nt -outfmt 11 -evalue
    ##  1e-5 -num_threads 10 -out all_seq_final_vs_nt.asn

    ## to be more flexible with the output of that time consuming blast I
    ## generated an archive that contains all possible information, we
    ## generate a tabular format including taxonomy ids using:
    ## blast_formatter -archive all_seq_final_vs_nt.asn -outfmt "10
    ## qaccver saccver pident length mismatch gapopen qstart qend sstart
    ## send evalue bitscore staxid" > all_seq_final_vs_nt.blttax

    ## we read that ouput into R blast <-
    blast <- read.csv("/SAN/Zebra/all_seq_final_vs_nt.blttax", header=FALSE)

    names(blast) <- c("query", "subject", "pident", "length", "mismatch",
                      "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                      "bitscore", "staxid")

    taxaNodes <- read.nodes("/SAN/db/taxonomy/nodes.dmp")
    taxaNames <- read.names("/SAN/db/taxonomy/names.dmp")

    blast.tax <- getTaxonomy(unique(blast$staxid),taxaNodes,taxaNames)
    blast.tax <- as.data.frame(blast.tax)
    blast.tax$staxid <- unique(blast$staxid)
    rownames(blast.tax) <- NULL

    saveRDS(blast.tax, file="/SAN/Zebra/blast_tax.Rds")

    blastT <- merge(blast, blast.tax, by="staxid")

    blt <- as.data.table(blastT)

    blt <- blt[,.(bitsum=sum(bitscore),
                  superkingdom, phylum, class, order, family, genus, species),
               by=c("query", "subject")]

    blt <- unique(blt)

    blt <- blt[,.(bitdiff= bitsum - max(bitsum),
                  superkingdom, phylum, class, order, family, genus, species),
               by=c("query")]

    get.opt.gen <- function(para, what){
        get.unique.or.na <- function (x){
            ## unique taxa at that level excluding potential NA's 
            ux <- unique(as.character(x[!is.na(x)]))
            ## but return NA if they are not unique
            if(length(ux)==1){return(ux)} else {as.character(NA)}
        }

        df <- blt[bitdiff > para, .(get.unique.or.na(eval(as.name((what))))),
                     by=query]
        setnames(df, "V1", what)
        length(df[, eval(as.name(what))][is.na(df[, eval(as.name(what))])])
    }

    get.opt.gen(-2.8, "genus")
    
    optimize(get.opt.gen, upper = -1, lower=-50, what="genus")


    genus <- blt[bitdiff > -2.8, .(get.unique.or.na(genus)),
                 by=query]
    
    family <- blt[bitdiff>-20, .(family=get.unique.or.na(family)),
                  by=query]

    order <- blt[bitdiff>-30, .(order=get.unique.or.na(order)),
                 by=query]

    class <- blt[bitdiff>-40, .(class=get.unique.or.na(class)),
                 by=query]

    phylum <- blt[bitdiff>-50, .(phylum=get.unique.or.na(phylum)),
                  by=query]

    superkingdom <- blt[bitdiff>-100, .(superkingdom=get.unique.or.na(superkingdom)),
                        by=query]

    annot <- cbind(superkingdom[,c("query", "superkingdom")],
                   phylum[,"phylum"],
                   class[,"class"],
                   order[,"order"],
                   family[,"family"],
                   genus[,"genus"])

    ## now we have to break this up into an annotation list for each
    ## amplicon
    seqnametab <- as.data.table(cbind(query=names(sequences), sequences))
    seqnametab <- merge(seqnametab, annot)

    dupseq <- seqnametab$sequences[duplicated(seqnametab$sequences)]
    ## seqnametab[sequences%in%dupseq,]

    seqnametab <- seqnametab[!duplicated(seqnametab$sequences),]

    annot.list <- lapply(STnoC, function (x) {
        setkey(seqnametab, sequences)
        seqnametab[colnames(x),
                   c("superkingdom", "phylum", "class", "order", "family", "genus")]
    })
    saveRDS(annot.list, file="/SAN/Zebra/annot_list.Rds")
} else {
    annot.list <- readRDS(file="/SAN/Zebra/annot_list.Rds")
}

## whatch out for this creating bugs
nrow(annot.list[["Mach1.Nem_0425_4 18S"]])
nrow(STnoC[["Mach1.Nem_0425_4 18S"]])

## drop all amplicons ones with no annotation
keep <- unlist(lapply(annot.list, nrow))>0

annot.list <- annot.list[keep]
STnoC <- STnoC[keep]

## now they are in sync
cbind(cumsum(unlist(lapply(annot.list, nrow))), cumsum(unlist(lapply(STnoC, ncol))))

## tabulate Phyla for each amplicon
lapply(annot.list, function (x) table(x[, "phylum"]))

tabulate.genera <- function(tax, subset){
    t <- tax[phylum%in%subset, ]
    if(!is.null(ncol(t))){
        table(t[,genus])
    } else {NULL}
}

## tabulate generaa for some phyla
lapply(annot.list, function (x) tabulate.genera(x,  "Nematoda"))
lapply(annot.list, function (x) tabulate.genera(x,  "Apicomplexa"))
lapply(annot.list, function (x) tabulate.genera(x,  "Platyhelminthes"))
lapply(annot.list, function (x) tabulate.genera(x,  "Streptophyta"))

## now get the real samples
sample.data <- read.csv("./sample_list_egg_counts_fgm.csv",
                        dec=",", stringsAsFactors=FALSE)

rownames(sample.data) <- make.unique(sample.data$Animal.No)

sample.data$FAA_index_I <- paste0("P1_", sample.data$FAA_index_I)

sample.data$FAA_index_II <- paste0("P2_", sample.data$FAA_index_II)

num.vars <- c("StrongyEPG", "CryptoEPG", "AscaridEPG", "AnoploEPG",
              "FGM", "rawash", "ndf", "adf", "adl", "cellul", "hemic",
              "energy", "N", "Prot", "Ca", "Cu", "Fe", "K", "Mg",
              "Mn", "Mo", "Na", "P", "S", "Zn")

fac.vars <- c("Date", "Sex", "Age", "Repro", "Season", "dens", "hab",
              "veg")

sample.data[, num.vars] <- apply(sample.data[, num.vars], 2, as.numeric)
sample.data[, fac.vars] <- apply(sample.data[, fac.vars], 2, as.factor)

samples.long <- reshape(sample.data,
                        varying=list(c("FAA_index_I", "FAA_index_II")),
                        direction="long", timevar = "pool",
                        idvar = "running_idx",
                        times=c("P_I", "P_II")) 

samples.long <-  samples.long[!samples.long$FAA_index_I%in%"P2_",]
rownames(samples.long) <- samples.long$FAA_index_I

pdf("figures/primers_MA_sorted_POS.pdf", 
    width=45, height=15, onefile=FALSE)
clust <- plotAmpliconNumbers(MA6[, which(colnames(MA6)%in%
                                         samples.long$FAA_index_I)])
dev.off()

two.clusters.row <- cutree(clust$tree_row, k=2)
two.clusters.col <- cutree(clust$tree_col, k=2)

keep.prime <- names(two.clusters.row)[two.clusters.row==1]
keep.sample <- names(two.clusters.col)[two.clusters.col==1]

MA <- MA6[which(rownames(MA6)%in%keep.prime),
          which(colnames(MA6)%in%
                keep.sample &
                colnames(MA6)%in% samples.long$FAA_index_I)]

pdf("figures/primers_MA_sorted_VAL.pdf", width=45, height=15, onefile=FALSE)
plotAmpliconNumbers(MA)
dev.off()

fill <- fillSampleTables(MA)
MA@sequenceTableFilled <- fill@sequenceTableFilled

## Analyse all at once for now
ALL <- Reduce(cbind, fill@sequenceTableFilled)

## same for tax
all.tax <- as.data.frame(Reduce(rbind, annot.list[rownames(MA)]))
all.tax <- as.matrix(all.tax)

## bring in same order and remove tax from bad amplicons (sorted out
## above) ## not necessary
## all.tax <- all.tax[colnames(ALL), ]
rownames(all.tax) <- colnames(ALL)

PS <- phyloseq(otu_table(ALL, taxa_are_rows=FALSE),
               sample_data(samples.long[rownames(ALL), ]),
               tax_table(all.tax))

Zeb.tab <- as.data.table(cbind(otu_table(PS),
                               Animal=sample_data(PS)$Animal.No))

Ccols <- colnames(Zeb.tab)[!colnames(Zeb.tab)%in%"Animal"]

Zeb <- Zeb.tab[, lapply(.SD, function(x) sum(as.numeric(x))),
               by=Animal, .SDcols=Ccols]

Zebra <- as.data.frame(Zeb[, ..Ccols])

rownames(Zebra) <- Zeb[, Animal]

PM <- phyloseq(otu_table(Zebra, taxa_are_rows=FALSE),
               sample_data(sample.data[rownames(Zebra), ]),
               tax_table(tax_table(all.tax)))
