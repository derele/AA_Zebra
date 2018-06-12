library(ggplot2)
library(MultiAmplicon)
library(parallel)
library(reshape)
library(vegan)

## Filter # only run when new filtered data is needed
FILTER <- FALSE
newMA <- FALSE
newAlign <- FALSE
newTree <- FALSE
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
                     truncLen=c(225,225), minLen=c(225,225), 
                     maxN=0, maxEE=2, truncQ=2, 
                     compress=TRUE, verbose=TRUE)
   })
   saveRDS(filter.track, file="/SAN/Zebra/filter.Rdata")
} else {
    filter.track <- readRDS(file="/SAN/Zebra/filter.Rdata")
}

filter <- do.call(rbind, filter.track)
colSums(filter)[2]/colSums(filter)[1]
## only 65% of data cept (decrease length?)

colSums(filter[1:(nrow(filter)/2),])
colSums(filter[1:(nrow(filter)/2),])[2]/colSums(filter[1:(nrow(filter)/2),])[1]

colSums(filter[(nrow(filter)/2):nrow(filter),])
colSums(filter[(nrow(filter)/2):nrow(filter),])[2]/colSums(filter[(nrow(filter)/2):nrow(filter),])[1]

file.ptable <- "/home/ele/Documents/Peter_Zebras/Zebra_Primer_List_Final.csv"

ptable <- read.csv(file.ptable, sep=",", header=TRUE)

primerF <- gsub(" ", "", ptable[,3])
primerR <- gsub(" ", "", ptable[,4])

names(primerF) <- as.character(ptable[,1])
names(primerR) <- as.character(paste(ptable[,2], ptable[, 5]), sep=";")

files <- PairedReadFileSet(filtFs, filtRs)

primers <- PrimerPairsSet(primerF, primerR)

rownames(ptable) <- names(primers)

getN <- function(x) sum(getUniques(x))

if (newMA){
    MA <- MultiAmplicon(primers, files)
    MA1 <- sortAmplicons(MA, starting.at=1, max.mismatch=4)
    MA2 <- derepMulti(MA1, mc.cores=20)

    MA3 <- dadaMulti(MA2, err=NULL, selfConsist=TRUE,
                     multithread=TRUE) ##  pool=TRUE)

    MA4.merged <- mergeMulti(MA3, justConcatenate=FALSE)

    before.merge <- unlist(lapply(MA4.merged@dada, function (x)
        sum(unlist(sapply(slot(x, "dadaF"), getN)))))
    after.merge <- unlist(lapply(MA4.merged@mergers, function (x)
        sum(unlist(sapply(x, getN)))))
    MA4 <- mergeMulti(MA3, justConcatenate=after.merge/before.merge<0.5)
    
    MA5 <- sequenceTableMulti(MA4)
    MA6 <- noChimeMulti(MA5, mc.cores=20)
    
    STnoC <- MA6@sequenceTableNoChime
    saveRDS(STnoC, file="/SAN/Zebra/table_STnoC_Npooled_mix.Rdata")
    saveRDS(MA6, file="/SAN/Zebra/MA6_Npooled_mix.Rdata")
} else {
    STnoC <- readRDS(file="/SAN/Zebra/table_STnoC_Npooled_mix.Rdata")
    MA6 <- readRDS(file="/SAN/Zebra/MA6_Npooled_mix.Rdata")
} 

pdf("figures/primers_MA_sorted.pdf", 
        width=45, height=15, onefile=FALSE)
    cluster <- plot_Amplicon_numbers(rawCounts(MA1))
dev.off()


## write sequences for taxonomy
if (FALSE){
    all.seq <- lapply(STnoC, colnames)
    all.seq <- unlist(all.seq)
    writeFasta(DNAStringSet(all.seq), "/SAN/Zebra/merged_seqs.fasta")
}

getN <- function(x) sum(getUniques(x))

track.l <- lapply(seq_along(MA6@dada), function (i) {
    track <- cbind(
        ## something strange like this is be needed to get the
        ## numbers for samples actually processed for that amplicon
        MA6@rawCounts[i, ][MA6@rawCounts[i, ]>0],
        sapply(slot(MA6@dada[[i]], "dadaF"), getN),
        sapply(MA6@mergers[[i]], getN, simplify=FALSE),
        rowSums(MA6@sequenceTable[[i]]),
        rowSums(MA6@sequenceTableNoChime[[i]]))
    if(ncol(track)==5){
    colnames(track) <- c("sorted",
                         "denoised", "merged",
                         "tabled", "nonchim")
    track <- apply(track, 2, as.numeric)
    rownames(track) <- rownames(MA6@sequenceTable[[i]])
    return(track)
    } else {return(data.frame(sorted=0,
                              denoised=0,
                              merged=0,
                              tabled=0,
                              nochim=0))
    }
})

track.drop <- as.data.frame(
    do.call("rbind", lapply(track.l, colSums, na.rm=TRUE)))

rownames(track.drop) <- rownames(MA6)

track.drop <- transform(track.drop,
                        perc.merge=merged/denoised)[order(track.drop$sorted),]

primer.df <- merge(track.drop, ptable, by=0)

primer_asv_numbers <- unlist(lapply(STnoC, ncol))

primer.df <- merge(primer.df, as.data.frame(primer_asv_numbers), by.x="Row.names", by.y=0)

primer.df$target_cat <- as.character(primer.df$Target.gene)
primer.df$target_cat <- ifelse(grepl("\\w{3}(L|H)", primer.df$target_cat),
                               "chloro", primer.df$target_cat)

if(FALSE){
    pdf("figures/seq_vs_asv.pdf")
    ggplot(primer.df, aes(primer_asv_numbers, sorted,
                          color=target_cat, shape=target_cat)) +
        geom_point()
    dev.off()

    pdf("figures/seq_vs_asv_log.pdf")
    ggplot(primer.df, aes(primer_asv_numbers, sorted,
                          color=target_cat, shape=target_cat)) +
        geom_point() +
        scale_y_log10()
    dev.off()

    pdf("figures/length_vs_perc.merge.pdf")
    ggplot(primer.df, aes(as.character(Amplicon_size_from), perc.merge, 
                          size=sorted, color=target_cat)) +
        geom_jitter()
    dev.off()

    pdf("figures/length_vs_numreads.pdf")
    ggplot(primer.df, aes(as.character(Amplicon_size_from), sorted, 
                          size=sorted, color=target_cat)) +
        geom_jitter()
    dev.off()
}


## quick look at what those COI sequences are
if(FALSE){
    writeFasta(colnames(MA6@sequenceTableNoChime[["COI_10F_3Mod.COI_500R COI"]]),
               "/SAN/Zebra/COI_10F_3Mod.COI_500R.fasta")
}
## they are adapters and rubbish


##### TAXONOMY assignment
assign.full.tax <- function(seqtab, what){
    seqs <- getSequences(seqtab)
    taxa <- list()
    if(what%in%c("16S", "18S")){
        cat("running silva assignment for 16S or 18S\n")
        taxa[["silva"]] <- assignTaxonomy(seqs, "/SAN/db/RDP/Silva_132/silva_nr_v132_train_set.fa.gz")
    }
    if(what%in%"18S"){
        cat("running own assignment for 18S\n")
        taxa[["own"]] <- assignTaxonomy(seqs, "/SAN/db/RDP/NR_ALLMERGE.fasta")
        cat("running pr2 assignment for 18S\n")
        taxa[["pr2"]] <- assignTaxonomy(seqs, "/SAN/db/RDP/PR2/pr2_version_4.9.2_dada2.fasta.gz",
                                      taxLevels = c("Kingdom","Supergroup","Division","Class",
                                                    "Order","Family","Genus","Species"))
    }
    if(what%in%"28S"){
        cat("running own assignment for 28S\n")
        taxa[["own"]] <- assignTaxonomy(seqs, "/SAN/db/RDP/NR_ALLMERGE.fasta")
        cat("running RDP fungi assignment for 28S\n")
        taxa[["RDPfun"]] <- assignTaxonomy(seqs, "/SAN/db/RDP/RDP_LSU_fixed_train_set_v2.fa.gz")
    }
    return(taxa)
}


assign.full.tax <- function(seqtab, what){
    seqs <- getSequences(seqtab)
    taxa <- list()
    if(what%in%c("16S")){
        cat("running silva assignment for 16S\n")
        taxa <- assignTaxonomy(seqs, "/SAN/db/RDP/Silva_132/silva_nr_v132_train_set.fa.gz")
    }
    if(what%in%"18S"){
        cat("running own assignment for 18S\n")
        taxa <- assignTaxonomy(seqs, "/SAN/db/RDP/annotated_18S_ena.fasta")
    }
    if(what%in%"28S"){
        cat("running own assignment for 28S\n")
        taxa <- assignTaxonomy(seqs, "/SAN/db/RDP/annotated_28S_ena.fasta")
    } else(taxa <- NULL)
    return(taxa)
}

## running it....
if(newTax){
    tax.l.l <- mclapply(seq_along(STnoC), function (i){
        if(length(dim(STnoC[[i]]))==2){
            assign.full.tax(STnoC[[i]], gsub(".*? (\\d\\dS$)", "\\1" , names(STnoC)[[i]]))
        } else {NULL}
    }, mc.cores=20)
    saveRDS(tax.l.l, "/SAN/Zebra/tax.own.Rds")
} else{tax.l.l <- readRDS("/SAN/Zebra/tax.own.Rds")}


## why does this fail?
assignTaxonomy(getSequences(STnoC[[3]]),
               "/SAN/db/RDP/Silva_132/silva_nr_v132_train_set.fa.gz",
               multithread = 20)


names(tax.l.l) <- names(STnoC)

 <- unlist(lapply(tax.l.l, class)%in%"matrix")


tabulate.genera <- function(tax, subset){
    t <- tax[tax[, 2]%in%subset, ]
    table(t[,6])
}

tabulate.genera(tax.l.l[["NL1_5Mod.NL4_5Mod 28S"]], "Nematoda")

taxan.18s <- unlist(lapply(tax.l.l[grepl("18S$", names(tax.l.l))], length))==3
taxan.18s <- names(taxan.18s)[taxan.18s]

STnoC.ann18s <- STnoC[taxan.18s]
tax.18s.l <- tax.l.l[taxan.18s]

## Warning message:
## In mclapply(seq_along(STnoC), function(i) { :
##   scheduled cores 15, 17, 3, 18 encountered errors in user code, all values of the jobs will be affected


## have a brief look at what we got....

cbind(
     names(STnoC),
     unlist(lapply(tax.l.l, class)),
     unlist(lapply(tax.l.l, length)),
     unlist(lapply(STnoC, ncol)),
     unlist(lapply(tax.l.l, function (x) paste(names(x), collapse="_")))
 )

foo.tax <- as.data.frame(Reduce(cbind, tax.l.l[[5]]))

colnames(foo.tax) <-
    paste(colnames(foo.tax),
          rep(names(tax.l.l[[5]]), times=unlist(lapply(tax.l.l[[5]], ncol))),
          sep="_")



foo.tax <- cbind(foo.tax, Count=unname(colSums(STnoC[[5]])))

## seq.Apicomplexa <- rownames(foo.tax[foo.tax$Division_pr2%in%"Apicomplexa", c("Family_own", "Genus_own", "Genus_pr2", "Species_pr2", "Count")])
## writeFasta(DNAStringSet(seq.Apicomplexa), "/home/ele/Apicomplexa.fasta")

rownames(foo.tax) <- NULL

## interesting
foo.tax[foo.tax$Phylum_silva%in%"Nematoda", c("Family_own", "Genus_own", "Genus_pr2", "Species_pr2", "Count")]

## interesting
foo.tax[foo.tax$Phylum_silva%in%"Nematoda", c("Family_own", "Genus_own", "Genus_pr2", "Species_pr2", "Count")]
foo.tax[foo.tax$Phylum_own%in%"Nematoda", c("Family_own", "Genus_own", "Genus_pr2", "Species_pr2", "Count")]


## interesting
foo.tax[foo.tax$Phylum_own%in%"Apicomplexa", c("Family_own", "Genus_own", "Genus_pr2", "Species_pr2", "Count")]
foo.tax[foo.tax$Division_pr2%in%"Apicomplexa", c("Family_own", "Genus_own", "Genus_pr2", "Species_pr2", "Count")]

## for 28S

foo.tax <- as.data.frame(Reduce(cbind, tax.l.l[[29]]))

colnames(foo.tax) <-
    paste(colnames(foo.tax),
          rep(names(tax.l.l[[29]]), times=unlist(lapply(tax.l.l[[29]], ncol))),
          sep="_")

foo.tax <- cbind(foo.tax, Count=unname(colSums(STnoC[[29]])))

rownames(foo.tax) <- NULL

table(foo.tax$Phylum_RDPfun)
table(foo.tax$Phylum_own)

own.genus.counts <- tapply(foo.tax$Count, foo.tax$Genus_own, sum)

own.genus.counts[order(own.genus.counts)]

table(foo.tax$Phylum_own, foo.tax$Phylum_RDPfun)

nrow(foo.tax[foo.tax$Phylum_own%in%"Nematoda", c("Family_own", "Genus_own", "Family_RDPfun", "Genus_RDPfun", "Count")])

foo.tax[foo.tax$Phylum_own%in%"Nematoda", c("Family_own", "Genus_own", "Count")]

foo.tax[foo.tax$Phylum_own%in%"Streptophyta", c("Family_own", "Genus_own", "Count")]

merged.tax.l <- lapply(seq_along(tax.l.l), function (i) {
    all.tax <- as.data.frame(Reduce(cbind, tax.l.l[[i]]))
    colnames(all.tax) <-
        paste(colnames(all.tax),
              rep(names(tax.l.l[[i]]), times=unlist(lapply(tax.l.l[[i]], ncol))),
              sep="_")
    return(all.tax)
})

names(merged.tax.l) <- names(STnoC)

X28S.tax.l <- merged.tax.l[grep("28S$", names(merged.tax.l))]
X28S.tax <- Reduce(rbind, X28S.tax.l)


X18S.tax.l <- merged.tax.l[taxan.18s]
X18S.tax <- Reduce(rbind, X18S.tax.l)
X18S.tax <- as.matrix(X18S.tax)

## Same for 18S is harder as some failed

### Here we get all tables filled for all samples 
all.reps <- unique(unlist(lapply(STnoC, rownames)))

STnoC.filled <- lapply(STnoC, function (ampST){
    missing.samples <- all.reps[!all.reps%in%rownames(ampST)]
    if(length(missing.samples)>0){
        bar <- matrix(0, nrow=length(missing.samples), ncol=ncol(ampST))
        rownames(bar) <- missing.samples
        foobar <- rbind(ampST, bar)
    } else {foobar <- ampST}
    foobar[all.reps, ]
})

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

## how many of our samples from processing files are in the sample
## table
table(samples%in%sample.data$FAA_index_I)
table(samples%in%sample.data$FAA_index_II)

samples.long <- reshape(sample.data,
                        varying=list(c("FAA_index_I", "FAA_index_II")),
                        direction="long", timevar = "pool",
                        idvar = "running_idx",
                        times=c("P_I", "P_II"))

dupes <- samples.long$FAA_index_I[duplicated(samples.long$FAA_index_I)]

### now only one "empty" sample in P2 got removed
samples.long <- samples.long[!samples.long$FAA_index_I%in%dupes, ]
rownames(samples.long) <- paste0(samples.long$FAA_index_I, "_F_filt.fastq.gz")

samples.long$is.negC <- samples.long$Animal.No%in%"neg"

## now only get samples for non-negative controls...
processed.samples <- rownames(samples.long)
non.neg.samples <- rownames(samples.long)[!samples.long$Animal.No%in%"neg"]

## intersect with everything that was sequenced
pos.samples <- intersect(non.neg.samples, all.reps)

samples.final <- samples.long

samples.final <- samples.final[pos.samples,]

RC <- rawCounts(MA6)

pdf("figures/primers_MA_sorted_POS.pdf", 
    width=45, height=15, onefile=FALSE)
plot_Amplicon_numbers(RC[, pos.samples])
dev.off()


### see which samples need to be excluded  18S
pdf("figures/primers_MA_sorted_Processed_18S.pdf", 
    width=45, height=10, onefile=FALSE)
cluster18S <- plot_Amplicon_numbers(RC[grepl("18S$", rownames(RC)), processed.samples])
dev.off()

split.2.18S <- cutree(cluster18S$tree_col, k=2)
table(split.2.18S)

X18S.pos <- names(split.2.18S)[split.2.18S==1]

### see which samples need to be excluded  16S
pdf("figures/primers_MA_sorted_Processed_16S.pdf", 
    width=45, height=4, onefile=FALSE)
cluster16S <- plot_Amplicon_numbers(RC[grepl("16S$", rownames(RC)), processed.samples])
dev.off()

split.2.16S <- cutree(cluster16S$tree_col, k=2)
table(split.2.16S)

X16S.pos <- names(split.2.16S)[split.2.16S==1]
X16S.pos <- X16S.pos[!X16S.pos%in%c("P2_FLD0265_F_filt.fastq.gz", "P2_FLD0241_F_filt.fastq.gz")]

### see which samples need to be excluded  28S
pdf("figures/primers_MA_sorted_Processed_28S.pdf", 
    width=45, height=4, onefile=FALSE)
cluster28S <- plot_Amplicon_numbers(RC[grepl("28S$", rownames(RC)), processed.samples])
dev.off()

split.2.28S <- cutree(cluster28S$tree_col, k=2)
table(split.2.28S)

X28S.pos <- names(split.2.28S)[split.2.28S==1]
X28S.pos <- X28S.pos[!X28S.pos%in%c("P2_FLD0265_F_filt.fastq.gz", "P2_FLD0241_F_filt.fastq.gz")]

## simply adding up technical replicates
library(dplyr)

STnoC.final <- lapply(STnoC.filled, function(x){
    seqs <- colnames(x)
    if(length(seqs) > 1){
        me <- merge(samples.final, x, by=0)
        summed <- summarise_at(group_by(me, Animal.No), seqs, sum)
        dat <- as.data.frame(summed)
        row.names(dat) <- dat$Animal.No
        dat$Animal.No <- NULL
        return(dat)
    } else{NULL}
})

##

## Makes not a lot of sense to analyse all at once for now
ALL <- Reduce(cbind, STnoC.final[!unlist(lapply(STnoC.final, is.null))])
dupe.otus <- colnames(ALL)[duplicated(colnames(ALL))]
## this is interesting! When merging ASVs properly, we get duplicate
## ASVs in different amplicons ... need to account for that!!!
## first a quick fix!
## loosing 13k counts with this
ALL <- ALL[, !colnames(ALL)%in%dupe.otus]
scale <- max(rowSums(ALL))/rowSums(ALL)

X18S <- Reduce(cbind, STnoC.final[taxan.18s])
scale.18S <- max(rowSums(X18S))/rowSums(X18S)
X18S.scaled <- X18S*scale.18S
X18S.scaled <- X18S.scaled[rowSums(X18S)>3000, ]
dupe.otus.18S <- colnames(X18S.scaled)[duplicated(colnames(X18S.scaled))]
X18S.scaled <- X18S.scaled[, !colnames(X18S.scaled)%in%dupe.otus.18S]

X16S <- Reduce(cbind, STnoC.final[grep("16S$", names(STnoC.final))])
## X16S <- X16S[X16S.pos,]
## not for 16S
## dupe.otus <- colnames(X16S)[duplicated(colnames(X16S))]
## X16S <- X16S[, !colnames(X16S)%in%dupe.otus]

X28S <- Reduce(cbind, STnoC.final[grep("28S$", names(STnoC.final))])
X28S <- X28S[rowSums(X28S)>3000,]
scale.28S <- max(rowSums(X28S))/rowSums(X28S)
X28S <- X28S*scale.28S

## Look here for the strainge bug...
## dupe.otus
## [1] "Animal.No" "Animal.No" "Animal.No" "Animal.No" "Animal.No"
dupe.otus.28S <- colnames(X28S)[duplicated(colnames(X28S))]
X28S <- X28S[, !colnames(X28S)%in%dupe.otus.28S]

library(phyloseq)
################################## 18S ##########################
ps18S <- phyloseq(otu_table(X18S.scaled, taxa_are_rows=FALSE),
                  sample_data(sample.data[rownames(X18S.scaled), ]),
                  tax_table(X18S.tax[colnames(X18S.scaled), ]))

ps18S.genus.own <- tax_glom(ps18S, "Genus_own")

Gen.tab <- otu_table(subset_taxa(ps18S.genus.own))

genus.names <- tax_table(ps18S.genus.own)[, 12]

gen.c.l <- by(t(Gen.tab), genus.names, colSums) 

Gen.tab <- t(do.call(rbind, gen.c.l))

round(colSums(Gen.tab)[order(colSums(Gen.tab))])

Gen.tab <- as.data.frame(merge(Gen.tab, sample.data, by=0))

get.best.cors <- function(column, cutoff) {
    foo <- cor(Gen.tab[, unique(genus.names)], Gen.tab[, column],
               use="na.or.complete", method="spearman")
    rows <- apply(foo, 1, function (x) any(abs(x)>cutoff))
    foo[which(rows), ]
}

get.best.cors("StrongyEPG", 0.1)[order(get.best.cors("StrongyEPG", 0.1))]

get.best.cors("AscaridEPG", 0.2)[order(get.best.cors("AscaridEPG", 0.2))]

get.best.cors("AnoploEPG", 0.2)[order(get.best.cors("AnoploEPG", 0.2))]

get.best.cors("CryptoEPG", 0.15)[order(get.best.cors("CryptoEPG", 0.15))]

## Nothing really...!!!! Found the bug?

devtools::source_gist("524eade46135f6348140",
                      filename = "ggplot_smooth_func.R")


pdf("figures/StrongyCor.pdf")
ggplot(Gen.tab, aes(Strongylus, StrongyEPG)) +
    geom_point() +
    scale_y_log10() + scale_x_log10() +
    geom_smooth(method="lm", se=FALSE) +
    stat_smooth_func(geom="text", method="lm", hjust=0, parse=TRUE) +
    theme_bw()
dev.off()


pdf("figures/StrongyOsterCor.pdf")
ggplot(Gen.tab, aes(Strongylus, Ostertagia)) +
    geom_point() +
    scale_y_log10() + scale_x_log10() +
    geom_smooth(method="lm", se=FALSE) +
    stat_smooth_func(geom="text", method="lm", hjust=0, parse=TRUE) +
    theme_bw()
dev.off()



pdf("figures/AscaridCor.pdf")
ggplot(Nem.tab, aes(AscaridEPG, other.nem)) +
    geom_point() +
    scale_y_log10() + scale_x_log10() +
    geom_smooth(method="lm")
dev.off()

pdf("figures/CestCor.pdf")
ggplot(Cest.tab, aes(AnoploEPG, Anoplocephala + Schistosoma + Trichobilharzia))+
    geom_point() +
    scale_y_log10() + scale_x_log10() +
    geom_smooth(method="lm")
dev.off()

pdf("figures/ApiCor.pdf")
ggplot(Api.tab, aes(CryptoEPG, Gregarina + Amoebogregarina + Leidyana))+
    geom_point() +
    scale_y_log10() + scale_x_log10() +
    geom_smooth(method="lm")
dev.off()



###  Overall analysis of variable influencing composition
##https://stats.stackexchange.com/questions/312302/adonis-in-vegan-order-of-variables-non-nested-with-one-degree-of-freedom-for

## making a data-frame from variables
ps18sdat <- sample_data(ps18S)
class(ps18sdat) <- "data.frame"


table(tax_table(subset_taxa(ps18S, Phylum_own%in%"Nematoda"))[, "Genus_pr2"])


## data without NA
ps18S.clean <- subset_samples(ps18S, !is.na(sample_data(ps18S)$cellul)&
                                     !is.na(sample_data(ps18S)$FGM))

ps18sdat.clean <- ps18sdat[!is.na(ps18sdat$cellul)&
                           !is.na(sample_data(ps18S)$FGM), ]

## now everything is gone...? But above?

## dbrda because it can do maringals  
ado.18S.all <- dbrda(otu_table(ps18S.clean)~cellul+Season+veg+
                         Sex+Age+Repro+dens+FGM,
                     data=ps18sdat.clean) ##, distance="bray")

marg.ado.18S.all <- anova(ado.18S.all, by = 'margin')

## "NMDS", "bray" does not work ... PCoA works
X18S.ord <- ordinate(ps18S, "PCoA", "bray")

X18S.ccpna <- ordinate(ps18S.clean, "CCA",
                       formula = ps18S.clean ~ Age + cellul + dens)

## insepcting this it looks like age is an artefact influenced by  outliers. 
### autoplot(X18S.ccpna)

## vegetation and density are worth plotting

pdf("figures/ordi_18S_raw_brayPCoA_densVeg.pdf")
plot_ordination(ps18S, X18S.ord, color="veg", shape="dens") +
    geom_point(size=5) 
dev.off()

## 18S: no comositional difference between sexes!?

pdf("figures/ordi_18S_raw_brayPCoA_vegcol.pdf")
plot_ordination(ps18S, X18S.ord, color="veg", shape="Season") +
        geom_point(size=5) 
dev.off()

## 18S: high comositional difference between vegetation!?



pdf("figures/ordi_18S_raw_brayPCoA_seascol.pdf")
plot_ordination(ps18S, X18S.ord, color="Season", shape="veg") +
    geom_point(size=5) 
dev.off()

## 18S: high compositional difference between seasons, maybe less than
## 16S!?

pdf("figures/ordi_18S_raw_brayPCoA_cellulcol.pdf")
plot_ordination(ps18S, X18S.ord, color="cellul", shape="Season") +
    geom_point(size=5) 
dev.off()

## 18S: composition varying with cellulose content in diet!?


ps18sdat.cellul <- ps18sdat[!is.na(ps18sdat$cellul), ]


############################################# 16 S ###############
ps16S <- phyloseq(otu_table(X16S, taxa_are_rows=FALSE),
                  sample_data(samples.final))


X16S.ord <- ordinate(ps16S, "PCoA", "bray")

pdf("figures/ordi_16S_raw_brayPCoA_sexcol.pdf")
plot_ordination(ps16S, X16S.ord, color="Sex", shape="Season")+
geom_point(size=5)
dev.off()

## 16S: no comositional difference between sexes!?

pdf("figures/ordi_16S_raw_brayPCoA_vegcol.pdf")
plot_ordination(ps16S, X16S.ord, color="veg", shape="Season") +
    geom_point(size=5) 
dev.off()

## 16S: high comositional difference between vegetation!?

ano.16S.veg <- anosim(otu_table(ps16S), sample_data(ps16S)$veg)

ps16sdat <- sample_data(ps16S)
class(ps16sdat) <- "data.frame"

pdf("figures/ordi_16S_raw_brayPCoA_seascol.pdf")
plot_ordination(ps16S, X16S.ord, color="Season", shape="veg") +
    geom_point(size=5) 
dev.off()

## 16S: high compositional difference between seasons!?

pdf("figures/ordi_16S_raw_brayPCoA_cellulcol.pdf")
plot_ordination(ps16S, X16S.ord, color="cellul", shape="Season") +
    geom_point(size=5) 
dev.off()

## 16S: composition varying with cellulose content in diet!?

ps16S.cellul <- subset_samples(ps16S, !is.na(sample_data(ps16S)$cellul))

ps16sdat.cellul <- ps16sdat[!is.na(ps16sdat$cellul), ]

###  Overall
##https://stats.stackexchange.com/questions/312302/adonis-in-vegan-order-of-variables-non-nested-with-one-degree-of-freedom-for

## dbrda because it can do maringals  
ado.16S.all <- dbrda(otu_table(ps16S.cellul)~cellul+Season+veg+
                     Sex+Age+Repro+dens,
                     data=ps16sdat.cellul)

marg.ado.16S.all <- anova(ado.16S.all, by = 'margin')

########################### 28S ##############################
ps28S <- phyloseq(otu_table(X28S, taxa_are_rows=FALSE),
                  sample_data(sample.data[rownames(X28S), ]),
                  tax_table(as.matrix(X28S.tax)))

ps28S.genus.own <- tax_glom(ps28S, "Genus_own")

Gen.tab.28S <- otu_table(subset_taxa(ps28S.genus.own))

genus.names28 <- tax_table(ps28S.genus.own)[, "Genus_own"]

gen.c.l <- by(t(Gen.tab.28S), genus.names28, colSums) 

Gen.tab.28S <- t(do.call(rbind, gen.c.l))

round(colSums(Gen.tab.28S)[order(colSums(Gen.tab.28S))])

Gen.tab.28S <- as.data.frame(merge(Gen.tab.28S, sample.data, by=0))

get.best.cors <- function(column, cutoff) {
    foo <- cor(Gen.tab.28S[, unique(genus.names28)], Gen.tab.28S[, column],
               use="na.or.complete", method="spearman")
    rows <- apply(foo, 1, function (x) any(abs(x)>cutoff))
    foo[which(rows), ]
}


get.best.cors("StrongyEPG", 0.1)[order(get.best.cors("StrongyEPG", 0.1))]

pdf("figures/Strongy28SCor.pdf")
ggplot(Gen.tab.28S, aes(Petrovinema + Cylicocyclus, StrongyEPG)) +
    geom_point() +
    scale_y_log10() + scale_x_log10() +
    geom_smooth(method="lm", se=FALSE) +
    stat_smooth_func(geom="text", method="lm", hjust=0, parse=TRUE) +
    theme_bw()
dev.off()


X28S.ord <- ordinate(ps28S, "PCoA", "bray")

pdf("figures/ordi_28S_raw_brayPCoA_sexcol.pdf")
plot_ordination(ps28S, X28S.ord, color="Sex", shape="Season")+
geom_point(size=5)
dev.off()

## 28S: no comositional difference between sexes!?

pdf("figures/ordi_28S_raw_brayPCoA_vegcol.pdf")
plot_ordination(ps28S, X28S.ord, color="veg", shape="Season") +
    geom_point(size=5) 
dev.off()

## 28S: comositional difference between vegetation (but less than
## 16S)!?

ano.28S.veg <- anosim(otu_table(ps28S), sample_data(ps28S)$veg)

ps28sdat <- sample_data(ps28S)
class(ps28sdat) <- "data.frame"

pdf("figures/ordi_28S_raw_brayPCoA_seascol.pdf")
plot_ordination(ps28S, X28S.ord, color="Season", shape="veg") +
    geom_point(size=5) 
dev.off()

## 28S: not a strong comositional difference between seasons?!

pdf("figures/ordi_28S_raw_brayPCoA_cellulcol.pdf")
plot_ordination(ps28S, X28S.ord, color="cellul", shape="Season") +
    geom_point(size=5) 
dev.off()

ps28S.cellul <- subset_samples(ps28S, !is.na(sample_data(ps28S)$cellul))

ano.28S.cell <- anosim(otu_table(ps28S.cellul),
                       sample_data(ps28S.cellul)$cellul)


ps28sdat.cellul <- ps28sdat[!is.na(ps28sdat$cellul), ]

###  Overall
##https://stats.stackexchange.com/questions/312302/adonis-in-vegan-order-of-variables-non-nested-with-one-degree-of-freedom-for

## dbrda because it can do maringals  
ado.28S.all <- dbrda(otu_table(ps28S.cellul)~cellul+Season+veg+
                     Sex+Age+Repro+dens, ###+FGM,
                     data=ps28sdat.cellul)

anova(ado.28S.all, by = 'margin')

## 28S: composition varying with cellulose content in diet!?


psALL <- phyloseq(otu_table(ALL, taxa_are_rows=FALSE),
                  sample_data(samples.final))

