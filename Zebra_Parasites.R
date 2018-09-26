library(ggplot2)
library(MultiAmplicon)
library(parallel)
library(reshape)
library(vegan)

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
    ## these two give error further down, they are empty anyways so
    ## lets remove them
    MA4 <- MA4[which(!rownames(MA4)%in%c("Hadz_1200CR.wang1624CR6L 18S",
                                         "LCO1490.HCO2198_5Mod COI")), ]
    MA5 <- makeSequenceTableMulti(MA4, mc.cores=20, orderBy="nsamples")
    MA6 <- removeChimeraMulti(MA5, mc.cores=20)
    saveRDS(MA6, file="/SAN/Zebra/MA6_mix.Rds")
} else {
    MA6 <- readRDS(file="/SAN/Zebra/MA6_mix.Rds")
}

##### TAXONOMY assignment
assign.full.tax <- function(seqtab, what){
    seqs <- getSequences(seqtab)
    if(what%in%"18S"){
        taxa <- assignTaxonomy(seqs,
                               "/SAN/db/RDP/annotated_18S_ena.fasta")
    }
    if(what%in%"28S"){
        taxa <- assignTaxonomy(seqs, "/SAN/db/RDP/annotated_28S_ena.fasta")
    }
    return(taxa)
}

STnoC <- getSequenceTableNoChime(MA6)

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

## we use this to limit nr blast to usable entries
##  blastn -negative_gilist /SAN/db/blastdb/uncultured.gi -query all_seq_final.fasta -db /SAN/db/blastdb/nt/nt -outfmt 11 -evalue 1e-5 -num_threads 10 -max_target_seqs 5 -out all_seq_final_vs_nt.asn

## to be more flexible wit the output of that time consuming blast I
## generated an archive that contains all possible information, we
## generate a tabular format including taxonomy ids using:
## blast_formatter -archive all_seq_final_vs_nt.asn -outfmt "10
## qaccver saccver pident length mismatch gapopen qstart qend sstart
## send evalue bitscore staxid" > | all_seq_final_vs_nt.blttax

## we read that ouput into R blast <-
blast <- read.csv("/SAN/Zebra/all_seq_final_vs_nt.blttax", header=FALSE)

names(blast) <- c("query", "subject", "pident", "length", "mismatch",
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                  "bitscore", "staxid")

## now we can "taxonomerize" the output (add the full paths)
library(taxonomizr)
library(taxize)

taxaNodes <- read.nodes("/SAN/db/taxonomy/nodes.dmp")
taxaNames <- read.names("/SAN/db/taxonomy/names.dmp")

blast.tax <- getTaxonomy(blast$staxid,taxaNodes,taxaNames)
blast.tax <- as.data.frame(blast.tax)
rownames(blast.tax) <- NULL

## saveRDS(blast.tax, file="/SAN/Zebra/blast_tax.Rds")
## blast.tax <- readRDS(file="/SAN/Zebra/blast_tax.Rds")

blast <- cbind(blast, blast.tax)

library(data.table)

blt <- as.data.table(blast)

blt <- blt[,.(bitsum=sum(bitscore),
                  superkingdom, phylum, class, order, family, genus, species),
               by=c("query", "subject")]

blt <- unique(blt)


blt <- blt[,.(bitdiff= bitsum - max(bitsum),
              superkingdom, phylum, class, order, family, genus, species),
           by=c("query")]


get.unique.or.na <- function (x){
    ## unique taxa at that level excluding potential NA's 
    ux <- unique(as.character(x[!is.na(x)]))
    ## but return NA if they are not unique
    if(length(ux)==1){return(ux)} else {as.character(NA)}
}


genus <- blt[bitdiff>-10, .(genus=get.unique.or.na(genus)),
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
seqnametab[sequences%in%dupseq,]

seqnametab <- seqnametab[!duplicated(seqnametab$sequences),]

annot.list <- lapply(STnoC, function (x) {
    setkey(seqnametab, sequences)
    seqnametab[colnames(x),
               c("superkingdom", "phylum", "class", "order", "family", "genus")]
})

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

library(phyloseq)
PS <- phyloseq(otu_table(ALL, taxa_are_rows=FALSE),
               sample_data(samples.long[rownames(ALL), ]),
               tax_table(all.tax))


################################## 18S and 28S  ##########################

library(DESeq2)

PS.para <- subset_taxa(PS, Phylum%in%c("Nematoda",
                                       "Apicomplexa",
                                       "Platyhelminthes"))

PG.para <- tax_glom(PS.para, "Genus", NArm = TRUE)

table(tax_table(PG.para)[,6])


Zeb <- by(otu_table(PG.para), sample_data(PG.para)$Animal.No, colSums)
Zebra <- do.call(cbind, Zeb)
Zebra <- Zebra[, !colnames(Zebra)%in%"neg"]


PSM <- phyloseq(otu_table(Zebra, taxa_are_rows=TRUE),
                sample_data(sample.data[colnames(Zebra), ]),
                tax_table(tax_table(PG.para)))


library(pheatmap)

mat <- as.matrix(otu_table(PSM.para.Ages))
rownames(mat) <- make.names(tax_table(PSM.para.Ages)[, 6])

colanot <- data.frame(lapply(sample_data(PSM.para.Ages)[, c("veg", "Sex", "Age")],
                             factor))

ann_colors <- list(
    veg = c("orange", "green", "darkgreen"),
    Sex = c("white", "red", "blue"),
    Age = c("blue", "red", "green", "orange")
)


pdf("figures/parasite_genera_heat.pdf")
pheatmap(log10(mat+1),  labels_col=rep("", times=ncol(mat)),
         annotation_col = colanot)
dev.off()


PSM.para.Ages <- subset_samples(PSM,
                                Age%in%c("fl", "juv", "mat", "sa") &
                                Season%in%c("d", "w") &
                                hab%in%c("lg", "sgp")
                                )


diagdds <- phyloseq_to_deseq2(PSM.para.Ages, ~ Age*Season)
diagdds <- estimateSizeFactors(diagdds, type="poscounts") # type="iterate"
diagdds.interact <- DESeq(diagdds, test="LRT",
                          reduced = ~ Age+Season, fitType="parametric")
diagdds.Season <- DESeq(diagdds, test="LRT", reduced = ~ Age, fitType="parametric")
diagdds.Age <- DESeq(diagdds, test="LRT", reduced = ~ Season, fitType="parametric")
diagdds.null <- DESeq(diagdds, test="LRT", reduced = ~ 1, fitType="parametric")


res.interact <- results(diagdds.interact, cooksCutoff = FALSE)
res.Season <- results(diagdds.Season, cooksCutoff = FALSE)
res.Age <- results(diagdds.Age, cooksCutoff = FALSE)
res.null <- results(diagdds.null, cooksCutoff = FALSE)

alpha <- 0.5

sigtab.interact <- res.interact[which(res.interact$padj < alpha), ]
sigtab.interact <- cbind(as(sigtab.interact, "data.frame"),
                as(tax_table(PG.para.Ages)[rownames(sigtab.interact), ], "matrix"))
rownames(sigtab.interact) <- NULL


sigtab.Age <- res.Age[which(res.Age$padj < alpha), ]
sigtab.Age <- cbind(as(sigtab.Age, "data.frame"),
                as(tax_table(PG.para.Ages)[rownames(sigtab.Age), ], "matrix"))
rownames(sigtab.Age) <- NULL

sigtab.Season <- res.Season[which(res.Season$padj < alpha), ]
sigtab.Season <- cbind(as(sigtab.Season, "data.frame"),
                as(tax_table(PG.para.Ages)[rownames(sigtab.Season), ], "matrix"))
rownames(sigtab.Season) <- NULL


sigtab.null <- res.null[which(res.null$padj < alpha), ]
sigtab.null <- cbind(as(sigtab.null, "data.frame"),
                as(tax_table(PG.para.Ages)[rownames(sigtab.null), ], "matrix"))
rownames(sigtab.null) <- NULL

### An interaction effect of Season and Age on Ostertagia and Cylicoccus
foo <- cbind(sample_data(PSM),
             t(unname(otu_table(subset_taxa(PSM, Genus%in%"Cylicocyclus")))))

sort(tapply(foo$sp1, as.factor(foo$Season):as.factor(foo$Age), median))

PSM.para.AdultF <- subset_samples(PSM,
                                  Age%in%"mat" &
                                  Repro%in%c("b", "h", "lact", "preg", "u") &
                                  Sex%in%c("f") &
                                  Season%in%c("d", "w") &
                                  hab%in%c("lg", "sgp")
                                  )


diagdds <- phyloseq_to_deseq2(PSM.para.AdultF, ~ Repro*Season)
diagdds <- estimateSizeFactors(diagdds, type="poscounts") # type="iterate"
diagdds.interact <- DESeq(diagdds, test="LRT",
                          reduced = ~ Repro+Season, fitType="parametric")

diagdds.Season <- DESeq(diagdds, test="LRT", reduced = ~ Repro, fitType="parametric")
diagdds.Repro <- DESeq(diagdds, test="LRT", reduced = ~ Season, fitType="parametric")
diagdds.null <- DESeq(diagdds, test="LRT", reduced = ~ 1, fitType="parametric")


res.interact <- results(diagdds.interact, cooksCutoff = FALSE)
res.Season <- results(diagdds.Season, cooksCutoff = FALSE)
res.Repro <- results(diagdds.Repro, cooksCutoff = FALSE)
res.null <- results(diagdds.null, cooksCutoff = FALSE)

alpha <- 0.5

sigtab.interact <- res.interact[which(res.interact$padj < alpha), ]
sigtab.Repro <- res.Repro[which(res.Repro$padj < alpha), ]
sigtab.Season <- res.Season[which(res.Season$padj < alpha), ]
sigtab.null <- res.null[which(res.null$padj < alpha), ]
## empty

PSM.para.AdultM <- subset_samples(PSM,
                                 Age%in%"mat" &
                                 Repro%in%c("b", "h", "lact", "preg", "u") &
                                 Sex%in%c("m") &
                                 Season%in%c("d", "w") &
                                 hab%in%c("lg", "sgp")
                                 )


diagdds <- phyloseq_to_deseq2(PSM.para.AdultM, ~ Repro*Season)
diagdds <- estimateSizeFactors(diagdds, type="poscounts") # type="iterate"
diagdds.interact <- DESeq(diagdds, test="LRT",
                          reduced = ~ Repro+Season, fitType="parametric")

diagdds.Season <- DESeq(diagdds, test="LRT", reduced = ~ Repro, fitType="parametric")
diagdds.Repro <- DESeq(diagdds, test="LRT", reduced = ~ Season, fitType="parametric")
diagdds.null <- DESeq(diagdds, test="LRT", reduced = ~ 1, fitType="parametric")


res.interact <- results(diagdds.interact, cooksCutoff = FALSE)
res.Season <- results(diagdds.Season, cooksCutoff = FALSE)
res.Repro <- results(diagdds.Repro, cooksCutoff = FALSE)
res.null <- results(diagdds.null, cooksCutoff = FALSE)

alpha <- 0.5

sigtab.interact <- res.interact[which(res.interact$padj < alpha), ]
sigtab.Repro <- res.Repro[which(res.Repro$padj < alpha), ]
sigtab.Season <- res.Season[which(res.Season$padj < alpha), ]
sigtab.null <- res.null[which(res.null$padj < alpha), ]

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
## Phylum order
x <- tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x <- sort(x, TRUE)
sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels=names(x))

## Genus order
x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab$Genus <- factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
    geom_point(size=6) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



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

