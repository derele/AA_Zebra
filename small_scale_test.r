library(ggplot2)
library(MultiAmplicon)

## Filter # only run when new filtered data is needed
FILTER <- FALSE
newMA <- FALSE
newAlign <- FALSE
newTree <- FALSE
newTax <- FALSE

a.files <- list.files(path="/SAN/Zebra/raw_fastqs_Prelim/",
                      pattern=".fastq.gz", full.names=TRUE)

Ffq.file <- a.files[grepl("R1", a.files)]
Rfq.file <- a.files[grepl("R2", a.files)]

samples <- gsub("/SAN/Zebra/raw_fastqs_Prelim//(.*?)\\_L001_R1_001.fastq\\.gz", "\\1", Ffq.file)

filt_path <- "/SAN/Zebra/DaDaFiltN_Prelim/"
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples

filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

if (FILTER){
    for(i in seq_along(Ffq.file)) {
        fastqPairedFilter(c(Ffq.file[i], Rfq.file[i]), c(filtFs[i], filtRs[i]),
                          truncLen=c(170,170), 
                          maxN=0, maxEE=2, truncQ=2, 
                          compress=TRUE, verbose=TRUE)
    }
}

file.ptable <- "/home/ele/Documents/Peter_Zebras/Zebra_Primer_List_Final.csv"

ptable <- read.csv(file.ptable, sep=",", header=TRUE)

primerF <- gsub(" ", "", ptable[,3])
primerR <- gsub(" ", "", ptable[,4])

names(primerF) <- as.character(ptable[,1])
names(primerR) <- as.character(paste(ptable[,2], ptable[, 5]), sep=";")

files <- PairedReadFileSet(filtFs, filtRs)

primers <- PrimerPairsSet(primerF, primerR)

if (newMA){
    MA <- MultiAmplicon(primers, files)
    MA1 <- sortAmplicons(MA, starting.at=1, max.mismatch=4)    
    MA2 <- derepMulti(MA1, mc.cores=20)
    MA3 <- dadaMulti(MA2, err=NULL, selfConsist=TRUE,
                     multithread=TRUE)
    MA4 <- mergeMulti(MA3, justConcatenate=TRUE)
    MA5 <- sequenceTableMulti(MA4)
    MA6 <- noChimeMulti(MA5, mc.cores=20)
    pdf("figures/primers_MA_sorted.pdf", 
        width=45, height=15, onefile=FALSE)
    cluster <- plot_Amplicon_numbers(rawCounts(MA1))
    dev.off()
    names(MA6@sequenceTableNoChime) <- rownames(MA6)
    STnoC <- MA6@sequenceTableNoChime
    save(STnoC, file="/SAN/Zebra/table.Rdata")
} else {load(file="/SAN/Zebra/table.Rdata")} ## -> STnoC

primer_seq_numbers <- unlist(lapply(STnoC, sum))
primer_seq_numbers[order(primer_seq_numbers)]

primer_asv_numbers <- unlist(lapply(STnoC, ncol))
primer_asv_numbers[order(primer_asv_numbers)]

primer.df <- data.frame(primer_asv_numbers,
                        primer_seq_numbers)

primer.df$target <- gsub(".*? +(.*$)", "\\1", rownames(primer.df))

primer.df$target_cat <- primer.df$target
primer.df$target_cat <- ifelse(grepl("\\w{3}(L|H)", primer.df$target_cat),
                               "chloro", primer.df$target_cat)

pdf("figures/seq_vs_asv.pdf")
ggplot(primer.df, aes(primer_asv_numbers, primer_seq_numbers,
                      color=target_cat, shape=target_cat)) +
    geom_point()
dev.off()

pdf("figures/seq_vs_asv_log.pdf")
ggplot(primer.df, aes(primer_asv_numbers, primer_seq_numbers,
                      color=target_cat, shape=target_cat)) +
    geom_point() +
    scale_y_log10()
dev.off()

PlantSeq <- DNAStringSet(colnames(STnoC[["trnL(UAA)c_5Mod.trnL(UAA)d_3Mod trnL(UAA)"]]))

writeFasta(PlantSeq, "Plants_trnL.fasta")


COIseq <- DNAStringSet(colnames(STnoC[["COI_10F_3Mod.COI_500R COI"]]))
writeFasta(COIseq, "COI.fasta")
