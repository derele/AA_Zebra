library(ggplot2)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)
library(dada2)
library(magrittr)

## was using the local devel version for some time
## devtools::load_all("../MultiAmplicon")

## having devel merged into main now:
## devtools::install_github("derele/MultiAmplicon")

library(MultiAmplicon)

## Filter # only run when new filtered data is needed
FILTER <- FALSE
newMAsort <- FALSE
newMApipe <- FALSE
newTax <- FALSE

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

if(newMAsort){
    fd <- "devStratFiles"
    if(dir.exists(fd)){
        unlink(fd, recursive=TRUE)
    }
    MA1 <- sortAmplicons(MA, filedir=fd, starting.at=1, max.mismatch=4)
    saveRDS(MA1, file="/SAN/Zebra/MA1.Rds")
} else {
    MA1 <- readRDS(file="/SAN/Zebra/MA1.Rds")
}

pdf("figures/primers_MA_sorted.pdf", width=46)
plotAmpliconNumbers(MA1)
dev.off()

if(newMApipe) {
    ## dereplication
    MAderep <- derepMulti(MA1, mc.cores=20)
    errF <-  learnErrors(unlist(getStratifiedFilesF(MAderep)), 
                         verbose=0, multithread = 64)
    errR <- learnErrors(unlist(getStratifiedFilesR(MAderep)), 
                        verbose=0, multithread = 64)
    
    MA3 <- dadaMulti(MAderep, Ferr=errF, Rerr=errR, pool=FALSE,
                     verbose=0, mc.cores=64)
    
    MA4.merged <- mergeMulti(MA3, justConcatenate=FALSE, mc.cores=20)
    
    prop.merged <- calcPropMerged(MA4.merged)
    prop.merged[is.na(prop.merged)] <- 0

    MA4 <- mergeMulti(MA3, justConcatenate=prop.merged<0.5,
                      mc.cores=64)
    
    MA5 <- makeSequenceTableMulti(MA4, mc.cores=64, orderBy="nsamples")
    MA6 <- removeChimeraMulti(MA5, mc.cores=64)
    saveRDS(MA6, file="/SAN/Zebra/MA6.Rds")
} else {
    MA6 <- readRDS(file="/SAN/Zebra/MA6.Rds")
}

if (newTax){
    MA7 <- blastTaxAnnot(MA6,
                         infasta="tmp_all_in.fasta",
                         outblast="tmp_blast_out.fasta",
                         negative_gilist="/SAN/db/blastdb/uncultured_gilist.txt", 
                         taxonSQL="/SAN/db/taxonomy/taxonomizr2.sql",
                         num_threads=64)
    saveRDS(MA7, file="/SAN/Zebra/MA7.Rds")
} else {
    MA7 <- readRDS(file="/SAN/Zebra/MA7.Rds")
}

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

## now add more content in  sample data, 
addSampleData(MA7, samples.long)

## In addSampleData(MA7, samples.long) : 242 samples are missing from
##   your sampleData but seem to have sequence data reported. They
##   will be omitted if you continue

## ->But we'd lose 242 samples which are not in the metadata table.

missingSamples <- rownames(MultiAmplicon::getSampleData(MA7))[
    !rownames(MultiAmplicon::getSampleData(MA7))%in%
    rownames(samples.long)]

table(rownames(MultiAmplicon::getSampleData(MA7))%in%
      rownames(samples.long))

## Good news is that these are basically like negative controls. They
## were no amplified, but the barcodes still missidentified.
summary(sapply(getSequenceTableNoChime(MA7[,which(colnames(MA7)%in%missingSamples)]), sum)/
sapply(getSequenceTableNoChime(MA7), sum) , na.rm=TRUE)
## -> about 5% average accross amplicons bad... 

sum(sapply(getSequenceTableNoChime(MA7[,which(colnames(MA7)%in%missingSamples)]), sum))/
sum(sapply(getSequenceTableNoChime(MA7), sum))
## -> about 6% of overall reads missassigend... 

## We can add them as basically negative controls.

samples.long[nrow(samples.long) + seq_along(missingSamples), "FAA_index_I"] <-
    missingSamples
samples.long[is.na(samples.long$Animal.No), "Animal.No"] <- "NEG_seq"

samples.long[, "pool"] <- ifelse(grepl("^P1_", samples.long$FAA_index_I), "P_I", "P_II")

rownames(samples.long) <- samples.long$FAA_index_I

MAsample <- addSampleData(MA7, samples.long)

## probably something to fix in the package??
table(sapply(getSequenceTableNoChime(MAsample), is.numeric))

dim(getSequenceTableNoChime(MAsample)[["rbcL-a_f_5Mod.rbcL-a_r_5Mod rbcL"]])
## [1] 0 0  -> an empty amplicon that should be removed
table(sapply(getSequenceTableNoChime(MAsample), nrow) > 0)


## we now have one amplicon without data?

## for annotation_col identifying negative samples
annoCols <- getSampleData(MAsample)[, c("pool", "Animal.No")]
annoCols$Animal.No <- ifelse(grepl("NEG_seq", annoCols$Animal.No), "NEG_seq",
                             ifelse(grepl("neg", annoCols$Animal.No), "NEG_pcr", "POS"))

attr(annoCols, "class") <- "data.frame"

pdf("figures/primers_MA_sorted_POS.pdf", 
    width=45, height=15, onefile=FALSE)
clust <- plotAmpliconNumbers(MAsample,
                             annotation_col = annoCols)
dev.off()


three.clusters.row <- cutree(clust$tree_row, k=3)
three.clusters.col <- cutree(clust$tree_col, k=3)

## keeping all primers for now
keep.prime <- names(three.clusters.row)[three.clusters.row < 3]
keep.sample <- names(three.clusters.col)[three.clusters.col==2]

MA8 <- MAsample[which(rownames(MAsample)%in%keep.prime),
          which(colnames(MAsample)%in%
                keep.sample &
                colnames(MAsample)%in% samples.long$FAA_index_I)]

## just sorting out primers whithout any taxannot
MA8 <- MA8[which( !unlist(lapply(MA8@taxonTable, is.null))), ]

pdf("figures/primers_MA_sorted_VAL.pdf", width=45, height=15, onefile=FALSE)
plotAmpliconNumbers(MA8)
dev.off()

PS <- toPhyloseq(MA8, samples=colnames(MA8))
## remove taxa with zero reports
PS <- prune_taxa(taxa_sums(PS) > 0, PS)

