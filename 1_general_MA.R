library(ggplot2)
library(MultiAmplicon)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)


## Filter # only run when new filtered data is needed
FILTER <- FALSE
newMA <- FALSE
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
                     truncLen=c(250,250), minLen=c(250,250), 
                     maxN=0, maxEE=4, truncQ=2, 
                     compress=TRUE, verbose=TRUE)
   })
   saveRDS(filter.track, file="/SAN/Zebra/filter.Rds")
} else {
    filter.track <- readRDS(file="/SAN/Zebra/filter.Rds")
}

filter <- do.call(rbind, filter.track)
colSums(filter)[2]/colSums(filter)[1]
## 60% of data overall cept (no need to decrease length for quality?!)

colSums(filter[1:(nrow(filter)/2),])
colSums(filter[1:(nrow(filter)/2),])[2]/colSums(filter[1:(nrow(filter)/2),])[1]
## 60 % in first run

colSums(filter[(nrow(filter)/2):nrow(filter),])
colSums(filter[(nrow(filter)/2):nrow(filter),])[2]/colSums(filter[(nrow(filter)/2):nrow(filter),])[1]
## 48% in second


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
    pdf("figures/primers_MA_sorted.pdf", width=46)
    plotAmpliconNumbers(MA1)
    dev.off()
    MA2 <- derepMulti(MA1, mc.cores=20)
    MA3.1 <- dadaMulti(MA2[, which(grepl("^P1", colnames(MA2)))],
                       Ferr=NULL, Rerr=NULL, selfConsist=TRUE,
                   multithread=FALSE, mc.cores=20, verbose=0, MAX_CONSIST=20)
    MA3.2 <- dadaMulti(MA2[, which(grepl("^P2", colnames(MA2)))],
                       Ferr=NULL, Rerr=NULL, selfConsist=TRUE,
                       multithread=FALSE, mc.cores=20, verbose=0, MAX_CONSIST=20)
    MA3 <- concatenateMultiAmplicon(MA3.1, MA3.2)    
    MA4.merged <- mergeMulti(MA3, justConcatenate=FALSE, minOverlap=8, maxMismatch=1,
                             mc.cores=20)
    prop.merged <- calcPropMerged(MA4.merged)
    prop.merged[is.na(prop.merged)] <- 0
    MA4 <- mergeMulti(MA3, justConcatenate=prop.merged<0.5,
                      minOverlap=8, maxMismatch=1,
                      mc.cores=20)
    ## removed a bug here, now working with empty amplicons
    MA5 <- makeSequenceTableMulti(MA4, mc.cores=20, orderBy="nsamples")
    MA6 <- removeChimeraMulti(MA5, mc.cores=20)
    saveRDS(MA6, file="/SAN/Zebra/MA6.Rds")
} else {
    MA6 <- readRDS(file="/SAN/Zebra/MA6.Rds")
}

## When loading an old MA object that lacks sample data, simply:
MA6 <- addSampleData(MA6)

if (newTax){
    MA7 <- blastTaxAnnot(MA6,
                         infasta="/SAN/Zebra/all_in.fasta",
                         outblast="/SAN/Zebra/blast_out.fasta", 
                         taxonSQL="/SAN/db/taxonomy/taxonomizr.sql")
    saveRDS(MA7, file="/SAN/Zebra/MA7.Rds")
} else {
    MA7 <- readRDS(file="/SAN/Zebra/MA7.Rds")
}


## ## Better way to do this needed!!!!!! Within package?
## ## tabulate Phyla for each amplicon
lapply(getTaxonTable(MA7), function (x) table(as.vector(x[, "phylum"])))

## lapply(getTaxonTable(MA7), function (x){
##     table(as.vector(x[x[, "phylum"]%in%"Nematoda", "genus"]))
## })


## lapply(getTaxonTable(MA7), function (x){
##     dim(x[as.vector(x[, "phylum"])%in%"Nematoda", "genus"])
## })


## lapply(getTaxonTable(MA7), function (x){
##     table(x[x[, "phylum"]%in%"Apicomplexa", "genus"])
## })

## lapply(getTaxonTable(MA7), function (x){
##     table(x[x[, "phylum"]%in%"Streptophyta", "genus"])
## })

## lapply(getTaxonTable(MA7), function (x){
##     table(x[x[, "phylum"]%in%"Streptophyta", "family"])
## })


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

## now add more content in  sample data
MAsample <- addSampleData(MA7, samples.long)


pdf("figures/primers_MA_sorted_POS.pdf", 
    width=45, height=15, onefile=FALSE)
clust <- plotAmpliconNumbers(MAsample[, which(colnames(MAsample)%in%
                                         samples.long$FAA_index_I)])
dev.off()

two.clusters.row <- cutree(clust$tree_row, k=2)
three.clusters.col <- cutree(clust$tree_col, k=3)

keep.prime <- names(two.clusters.row)[two.clusters.row==1]
keep.sample <- names(three.clusters.col)[three.clusters.col==1]

MA8 <- MAsample[which(rownames(MAsample)%in%keep.prime),
          which(colnames(MAsample)%in%
                keep.sample &
                colnames(MAsample)%in% samples.long$FAA_index_I)]

pdf("figures/primers_MA_sorted_VAL.pdf", width=45, height=15, onefile=FALSE)
plotAmpliconNumbers(MA8)
dev.off()

PS <- toPhyloseq(MA8, samples=colnames(MA8))

## This should acutally work with phyloseq's merge_samples function
## but doesn't as this messes up sample_data.
## CANDIDATE FOR INCLUSION IN PACKAGE...!
sumTecRep <- function (PS, by.sample){
    otab <- setDT(apply(otu_table(PS), 2, as.list))
    ## the columns giving numbers for sequences
    numcols <- colnames(otab)[nchar(colnames(otab))>10]
    sdat <- sample_data(PS, errorIfNULL = FALSE)
    otab[, (numcols):=lapply(.SD, as.numeric), .SDcols=numcols]
    otab[, sfac := as.factor(sdat[[by.sample]])]
    setkey(otab, sfac)
    otabN <- otab[, lapply(.SD, sum), by=sfac]
    setkey(otabN, sfac)
    OTN <- as.matrix(otabN, rownames=TRUE)
    ## now select the entries from colums that have the same values in
    ## the sample table...
    sdatN <- by(sdat, sfac, function(x){
        sapply(x, function (y){
            uy <- unique(y)
            if(length(uy)==1) uy else paste(uy, collapse=";")
        })
    })
    sdatN <- as.data.frame(do.call(rbind, sdatN))
    phyloseq(otu_table(OTN, taxa_are_rows=FALSE),
             sample_data(sdatN),
             tax_table(PS))
}



PM <- sumTecRep(PS, by.sample="Animal.No")


### Will see whether we need this and then make it possible WITHIN THE
### PACKAGE to generate them....

### the individual amplicons to phyloseq
## PS.l <- lapply(seq_along(rownames(MA)), function(i){
##     samples <- rownames(getSequenceTableNoChime(MA[i,]))
##     toPhyloseq(MA[i, ], samples)
## })

## names(PS.l) <- rownames(MA)

## PM.l <- lapply(PS.l, function(ps){
##     sumTecRep(ps, by.sample="Animal.No")
## })


############### Agglomerating by genus / OR NOT  ##########################

##### Raw non-merged technical replicates -> objects named "PSG": "PSG"
##### for merged amplicons "PSG.l" for the raw list per amplicosn
PSG <- tax_glom(PS, "genus")

## PSG.l <- mclapply(PS.l, tax_glom, "genus", mc.cores=12)

##### Merged technical replicates -> objects named "PMG": "PMG" for
##### merged amplicons "PMG.l" for the raw list per amplicosn

PGM <- tax_glom(PM, "genus")

## PGM.l <- mclapply(PM.l, tax_glom, "genus", mc.cores=12)

## sumSeqByTax <- function (Phy, tax) {
##     counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
##     counts$asvCount <- as.numeric(as.character(counts$asvCount))
##     tapply(counts$asvCount, counts[, tax], sum)
## }


## readNumByPhylum <- lapply(PS.l, sumSeqByTax, "phylum")
## names(readNumByPhylum) <- names(STnoC)


## readNumByGenus <- lapply(PS.l, sumSeqByTax, "genus") ## Change "text" in order to get counts per a different taxonomic level
## names(readNumByGenus) <- names(STnoC)


## ## in how many amplicons are genera observed
## all.genera <- unique(unlist(lapply(readNumByGenus, names)))

## foo <- lapply(all.genera, function (x) {
##     lapply(readNumByGenus, function(y){
##         x%in%names(y)
##     })
## })

## num.prim <- unlist(lapply(foo, function(x) sum(unlist(x))))
## names(num.prim) <- all.genera

## tail(num.prim[order(num.prim)], n=100)



## ## in how many amplicons are genera observed at a high 

## bar <- lapply(all.genera, function (x) {
##     lapply(readNumByGenus, function(y){
##         total <- sum(y)
##         prop <- y/total
##         x%in%(names(prop)[prop>0.001])
##     })
## })


## num.H.prim <- unlist(lapply(bar, function(x) sum(unlist(x))))

## names(num.H.prim) <- all.genera

## tail(num.H.prim[order(num.H.prim)], n=100)

## ## Checking out correlations
## cor.l <- lapply(PS.l, function(x) cor(otu_table(x)))


## pdf("figures/corHist.pdf")
## lapply(cor.l, function(x){
##     value <- x[upper.tri(x)]
##     ggplot(as.data.frame(value), aes(value)) +
##         geom_density() +
##         scale_y_log10()
## })
## dev.off()

